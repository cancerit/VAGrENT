package Sanger::CGP::VagrentSV::VagrentSV;

##########LICENCE##########
# Copyright (c) 2014 Genome Research Ltd.
# 
# Author: Cancer Genome Project cgpit@sanger.ac.uk
# 
# This file is part of VAGrENT.
# 
# VAGrENT is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
##########LICENCE##########


use strict;

use Log::Log4perl;

use Data::Dumper;

use Const::Fast qw(const);
use Try::Tiny qw(try catch);

use Sanger::CGP::Vagrent qw($VERSION);
use Sanger::CGP::VagrentSV::Base;

use FindBin qw($Bin);
Log::Log4perl->init("$Bin/../config/log4perl.vagrentsv.conf");
my $log = Log::Log4perl->get_logger(__PACKAGE__);


# reference data storage
use Sanger::CGP::VagrentSV::Data::GenomeSeq;
# bookmarkers
use Sanger::CGP::Vagrent::Bookmarkers::RepresentativeTranscriptBookmarker;
use Sanger::CGP::Vagrent::Bookmarkers::MostDeleteriousBookmarker;

# annotators
use Sanger::CGP::Vagrent::Annotators::AnnotatorCollection;
use Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator;
use Sanger::CGP::Vagrent::Annotators::InsertionAnnotator;
use Sanger::CGP::Vagrent::Annotators::DeletionAnnotator;
use Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator;

use Sanger::CGP::Vagrent::Data::Substitution;
use Sanger::CGP::Vagrent::Data::Transcript;




# Transcripts 
use Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource;



# SV modules
use Sanger::CGP::VagrentSV::Data::BreakPoint;
use Sanger::CGP::VagrentSV::Data::StructuralVariation;

# SV Annotators

use Sanger::CGP::VagrentSV::SVConstants;
use Sanger::CGP::VagrentSV::Annotators::SVAnnotator;
use Sanger::CGP::VagrentSV::Annotators::BreakPointTranscripts;
use Sanger::CGP::VagrentSV::Annotators::FusionGeneAnnotator;
use Sanger::CGP::VagrentSV::Results::TranscriptResults;
use Sanger::CGP::VagrentSV::Data::FileBasedAnnotationSource;


const my $REPRE_BM => Sanger::CGP::Vagrent::Bookmarkers::RepresentativeTranscriptBookmarker->new();
const my $WORST_BM => Sanger::CGP::Vagrent::Bookmarkers::MostDeleteriousBookmarker->new();

my @header=qw/chr lStart lEnd chr rStart rEnd name score lStrand rStrand microHomoLen cancer_type patientId sampleId
 lExon\/rExon lTranscript\/rTranscript lGeneStrand\/rGeneStrand 
 lAnnotation\/rAnnotation lGene\/rGene SpannedGenes
 lPromoterL\/lPromoterR rPromoterL\/rPromoterR 
 lEnhancerL\/lEnhancerR rEnhancerL\/rEnhancerR/; 


sub new {
  my ($class,$options) = @_;
  my $self = { };
  bless $self, $class;
  $self->_read_sv_annotation($options);
  return $self;
}

=head2 get_sv_data
Get breakpoint information for given SV 
Split breakpoints into two separate data sets with same keys
Inputs
=over 2
=item input_file 

chr[0]  bp1_start[1]    bp1_end[2]*     chr[3]  bp2[4]          bp2[5]*         Event ID:breakpoint ID:Event classification[6] 
14      55115029        55115030        14      64823314        64823315        0009b464-b376-4fbc_E000001:0009b464-b376-4fbc_B000433:del

score[7]    strand_1[8] strand2[9] MicroHLen[10] CancerType[11]    Patient_ID[12]  Sample/FileID[13]
	    1       +           -           -4              OV-AU           AOCS-117        0009b464-b376-4fbc
 
* actual break points to use for analysis
11. Length of microhomology (if negative) or non-templated sequence insertion (if positive) at the breakpoint junction.

=back
=cut

sub _read_sv_annotation {
	my ($self,$opts)=@_;
	my ($genome)= Sanger::CGP::VagrentSV::Data::GenomeSeq->new( 'fasta' => $opts->{'genome'}, 'index' => $opts->{'genome'}.'.fai' );
	my($outfile)=Sanger::CGP::VagrentSV::Base->open_to_write($opts->{'input'}.'.outttt');
	my($rfh)=Sanger::CGP::VagrentSV::Base->open_to_read($opts->{'input'});
	#my ($fai)=$self->_get_genome_object($opts->{'genome'});
	#my ($chr_lengths)=$self->_get_chromosome_length($opts->{'genome'}.'.fai');
	my ($annotator)=get_annotator($opts);	
	# get chache
	my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => $opts->{'cache'}, 'search_pad' => $Sanger::CGP::VagrentSV::SVConstants::PADDING_BUFFER_TR );
	my $bpa = Sanger::CGP::VagrentSV::Annotators::SVAnnotator->new(transcriptSource => $ts);
	my $bptr = Sanger::CGP::VagrentSV::Annotators::BreakPointTranscripts->new(transcriptSource => $ts);
	my $fga = Sanger::CGP::VagrentSV::Annotators::FusionGeneAnnotator->new(transcriptSource => $ts);
	my $promoter = Sanger::CGP::VagrentSV::Data::FileBasedAnnotationSource->new(annotationSource => '/nfs/users/nfs_s/sb43/scratch_tmp_storage_not_backed_up/vagrentSV_data/refdata/genome_features/gc19.flat.prom_per_gene.pc.bed');
	my $enhancer = Sanger::CGP::VagrentSV::Data::FileBasedAnnotationSource->new(annotationSource => '/nfs/users/nfs_s/sb43/scratch_tmp_storage_not_backed_up/vagrentSV_data/refdata/genome_features/roadmap_stringent_enhancers_with_Ensembl_ID_and_name.bed.txt');
	
	
	my $counter;
	
	print $outfile join("\t", @header)."\n";
	
	while(<$rfh>) {
		my ($chr_flag,$lparams,$lparams,$offset);
		$counter ++;
		chomp;
		my($bp1_chr,$bp1_pos,$bp2_chr,$bp2_pos,$name,$bp1_strand,$bp2_strand,$length)=(split "\t", $_ )[0,2,3,5,6,8,9,10];
		#check if SV is on same chromosome
		if(	$bp1_chr eq $bp2_chr) { $chr_flag=1; }
		if ($counter < 10 && $bp1_chr && $bp1_pos && $bp2_chr && $bp2_pos && $name && $bp1_strand && $bp2_strand ) {
				$log->debug("$counter : Format is BEDPE");
		}
		elsif ($counter < 10 ) {
			$log->logcroak('File format is not BEDPE, please convert your file to BEDPE before running this script');
		} 
		my ($svtype)=(split ":",$name)[2]; 

		next if $svtype ne 'del';
		
		# left hand break point
	
		my $lhb = Sanger::CGP::VagrentSV::Data::BreakPoint->new(	'species'				=> $opts->{'species'},
																															'genomeVersion' => $opts->{'assembly'},
																															'chr'	          => $bp1_chr,
																															'minpos'        => $bp1_pos, # include reference base 
																															'maxpos'        => $bp1_pos,
																															'strand'				=> $bp1_strand,
																																);
		#right hand break point
		
		my $rhb = Sanger::CGP::VagrentSV::Data::BreakPoint->new(	'species'				=> $opts->{'species'},
																															'genomeVersion' => $opts->{'assembly'},
																															'chr'	          => $bp2_chr,
																															'minpos'        => $bp2_pos,
																															'maxpos'        => $bp2_pos,
																															'strand'				=> $bp2_strand,
																															);
	  
	 
	  
	 my $sv = Sanger::CGP::VagrentSV::Data::StructuralVariation->new('lhb' => $lhb, 'rhb' => $rhb,  'info' => { ('name'  	=> $name,
																																																							'length'  => $length,
																																																							'svtype'  => $svtype,
																																																							'locflag' => $chr_flag,
																																																							 )} );
																																																							 
		my($lhbPromoterGenes,$rhbPromoterGenes)=$promoter->getRegulatoryAnnotations($sv);																																																					 
		my($lhbEnhancerGenes,$rhbEnhancerGenes)=$enhancer->getRegulatoryAnnotations($sv);																																																						 
		my ($lhb_genes,$spanned_genes,$rhb_genes)=$bptr->getOverlappingGenes($sv);
		#my $bp_genes = Sanger::CGP::VagrentSV::Results::TranscriptResults->new(%$bp_genes);
		
		my ($ltr,$rtr) = $bpa->getBreakpointTranscripts($sv,$bptr);	
		my ($fused_transcripts)	= $fga->getFusedTranscripts($ltr,$rtr,$sv);
					
		foreach my $ft (@$fused_transcripts) {
			if($ft){
				print $outfile $_."\t".$ft."\t".$lhb_genes.'/'.$rhb_genes."\t".$spanned_genes."\t".
				$lhbPromoterGenes->{'lgenes'}.'/'.$lhbPromoterGenes->{'rgenes'}."\t".$rhbPromoterGenes->{'lgenes'}.'/'.$rhbPromoterGenes->{'rgenes'}."\t".
				$lhbEnhancerGenes->{'lgenes'}.'/'.$lhbEnhancerGenes->{'rgenes'}."\t".$rhbEnhancerGenes->{'lgenes'}.'/'.$rhbEnhancerGenes->{'rgenes'}."\n";
				
				# output in defuse format as required by pegasus
				
			
			}
	  }
				
		 print "Annotated brekpoint: $counter\n"; 
		#$self->_writeBedpe($_,$sv_data);
	}
  
  close($rfh,$outfile);
}


sub input_file {
	shift->{'input_file'};
}

sub _get_overlapping_features {
	my ($self,$annotation_object, $chr,$start,$stop,$check_exon)=@_;
	print "****** $chr,$start,$stop:\n";
		$chr=~s/chr//g;
		my $genes;
		my $tmp_hash;
		foreach my $ann_obj_name (keys %{$annotation_object->annotation_obj}) {
 	 		my $tabix=$annotation_object->annotation_obj->{$ann_obj_name};
 	 		my $res = $tabix->query("chr$chr",$start,$stop);
 	 		if(defined $res->get) {
				while(my $record = $tabix->read($res)){	
					my ($chr,$exon,$start,$stop,$line) = (split '\t', $record) [0,2,3,4,8];
					#next if ($exon ne 'exon' && $check_exon);
					my ($tmp_gene)= (split '\"', $line) [1];
					print $record;
					if(!exists $tmp_hash->{$tmp_gene} ) {
						push(@$genes,$tmp_gene);
					}
					$tmp_hash->{$tmp_gene}++;
				}
 	 		}
 		}
 		#if (keys %$tmp_hash > 0) {
 			return $genes;
 		#}
 		#else {
 			#return "NA";
 		#}
}


sub _get_chromosome_length {
	my($self,$index_file)=@_;
	my $chr_length;
	my ($rfh_index)=Sanger::CGP::VagrentSV::Base->open_to_read($index_file);
	while(<$rfh_index>) {
		my($chr,$length)=(split "\t", $_)[0,1];
		$chr_length->{$chr}=$length;
	}
	close($rfh_index);
	return $chr_length;
}


=head2  _get_genome_object
create genome object from fasta file
Inputs
=item reference fasta file
=over 2
=back
=cut

sub _get_genome_object {
	my ($self,$genome)=@_;
	my $fai = Bio::DB::Sam::Fai->load($genome);
	return $fai;
}


sub get_annotator {
	my $options = shift;

  # creating an EnsemblTranscriptSource using the Ensembl registry
	my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => $options->{'cache'},'search_pad' => $Sanger::CGP::VagrentSV::SVConstants::PADDING_BUFFER_TR);

	# creating an AnnotatorCollection
	my $annotator = Sanger::CGP::Vagrent::Annotators::AnnotatorCollection->new( debug => 1
							#bookmarker => [$REPRE_BM,$WORST_BM],
							#only_bookmarked => 1
							);

	# creating and adding a SimpleSubstitutionAnnotator to the AnnotatorCollection
	$annotator->addAnnotator(Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts, ontologySymmary => $annotator->getSummaryCache));
	
	return $annotator;
}


sub get_annotator_complete {
	my $options = shift;

  # creating an EnsemblTranscriptSource using the Ensembl registry
	my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => $options->{'cache'});

	# creating an AnnotatorCollection
	my $annotator = Sanger::CGP::Vagrent::Annotators::AnnotatorCollection->new(
							bookmarker => [$REPRE_BM,$WORST_BM],
							only_bookmarked => 1);

	# creating and adding a SimpleSubstitutionAnnotator to the AnnotatorCollection
	$annotator->addAnnotator(Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts, ontologySymmary => $annotator->getSummaryCache));
	$annotator->addAnnotator(Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts, ontologySymmary => $annotator->getSummaryCache));
	$annotator->addAnnotator(Sanger::CGP::Vagrent::Annotators::DeletionAnnotator->new(transcriptSource => $ts, ontologySymmary => $annotator->getSummaryCache));
	$annotator->addAnnotator(Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts, ontologySymmary => $annotator->getSummaryCache));

	return $annotator;
}


sub _writeBedpe {
	my($self,$line,$sv_data)=@_;
	chomp $line;
	my $fused_genes=$self->_getString($sv_data->getFusedGenes);
	my $inb_genes=$self->_getString($sv_data->getInbGenes);
	print $line."\t".$fused_genes."\t".$inb_genes."\n";
}

sub _getString {
	my($self,$hash)=@_;
	my $array;
	foreach my $key (keys %$hash) {
		push(@$array,$key);
	}
	if($array && @$array > 0) {
  	return join(",",@$array);
  }else{
  	return "NA";
  }
}


1;


