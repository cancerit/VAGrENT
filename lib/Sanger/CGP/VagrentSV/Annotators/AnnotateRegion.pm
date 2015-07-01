package Sanger::CGP::VagrentSV::Annotators::AnnotateRegion;

use strict;
use Data::Dumper;
use Const::Fast qw(const);

use Sanger::CGP::Vagrent qw($VERSION);
use Sanger::CGP::VagrentSV::VagrentSVConstants;
use Sanger::CGP::VagrentSV::Data::SVData;
use base qw(Sanger::CGP::Vagrent::Data::AbstractVariation);




1;


sub new {
	my($class, %sv)=@_;
	my $self = { };
	bless $self, $class;
	foreach my $key(keys(%sv)){
		 $self->{"_$key"} = $sv{$key};
		}
	return $self;
}

sub isValid {
	my $self = shift;
	return 0 unless(defined($self->{_lbreak}) && defined($self->{_rbreak}) && defined($self->{_name}));
	return 1;
}




sub getDisruptedGenes{
	my($self,$ts,$fai,$chr_length)=@_;
	my ($fusedtr,$fused_genes);
	my $ext_maxpos=  $self->getInb->getMaxPos + $Sanger::CGP::VagrentSV::VagrentSVConstants::PADDING_BUFFER_SEQ;
	my $ext_minpos=  $self->getInb->getMinPos - $Sanger::CGP::VagrentSV::VagrentSVConstants::PADDING_BUFFER_SEQ;
  
	#get transcripts overlapping SV regions
	my($lhb_genes,$lhbtr)=_getSvFeatures($ts->getTranscripts($self->getLhb));
	my($rhb_genes,$rhbtr)=_getSvFeatures($ts->getTranscripts($self->getRhb));
		
	my($inb_genes,$inbtr)=_getSvFeatures($ts->getTranscripts($self->getInb));
	
	my($delseq)=_getSeq($fai,$self->getInb->getChr,$self->getInb->getMinPos,$self->getInb->getMaxPos);
	
	# check if $EXTENSION_BUFFER exceeds chromosome length
	if($chr_length->{$self->getInb->getChr} < $ext_maxpos ) {
		$ext_maxpos= $self->getInb->getMaxPos;
	}
	if( $ext_minpos < 0  ) {
		$ext_minpos= $self->getInb->getMinPos;
	}
		
	my($rhb_insseq)=_getSeq($fai,$self->getInb->getChr,$self->getInb->getMaxPos,$ext_maxpos);
	my($lhb_insseq)=_getSeq($fai,$self->getInb->getChr,$ext_minpos,$self->getInb->getMinPos );

	#get transcripts/genes exclusively within breakpoint
	$inbtr=_getUniqueFeatures($rhbtr,$inbtr);
	$inbtr=_getUniqueFeatures($lhbtr,$inbtr);
	
	$inb_genes=_getUniqueFeatures($rhb_genes,$inb_genes);
	$inb_genes=_getUniqueFeatures($lhb_genes,$inb_genes);

	
	if($lhbtr && $rhbtr) {
		$fusedtr={%$lhbtr,%$rhbtr};
		$fused_genes={%$lhb_genes,%$rhb_genes};
	}
	elsif($lhbtr) {
		$fusedtr=$lhbtr;
		$fused_genes=$lhb_genes;
	}
	elsif($rhbtr) {
		$fusedtr=$rhbtr;
		$fused_genes=$rhb_genes;
	}
	
	my ($sv_data)=Sanger::CGP::VagrentSV::Data::SVData->new('fusedtr' 		=>$fusedtr, 
								'fusedgenes'	=>$fused_genes,
								'lhbgenes'		=> $lhb_genes,
								'rhbgenes'		=> $rhb_genes,
								'inbtr'				=>$inbtr,
								'inbgenes' 		=>$inb_genes,
								#'inb_delseq' 	=>$delseq,
								'lhb_insseq' 	=>$lhb_insseq,
								'rhb_insseq' 	=>$rhb_insseq,
								'lhb_trobj' 	=>$ts->getTranscripts($self->getLhb),
								'rhb_trobj' 	=>$ts->getTranscripts($self->getRhb)     
								);		
		
	return $sv_data;
}


sub _getSvFeatures {
	my (@tr)=@_;
	my ($tr_list,$gene_list);
	foreach my $line (@tr) {
		if($line->getAccession) {
			$tr_list->{$line->getAccession}++;
			$gene_list->{$line->getGeneName}++;
	  }		
	}
	return ($gene_list,$tr_list);
}

sub _getUniqueFeatures {
	my ($break,$inbreak)=@_;
	foreach my $key(keys %$break) {
		if($inbreak->{$key}) {
			delete $inbreak->{$key};
		}
	}
	return $inbreak;
}


=head2  getSeq
get dan sequence for given interval
Inputs
=item Bio::DB::Sam index object
=over 2
=back
=cut

sub _getSeq {
	my ($fai,$chr,$start,$end)=@_;
	
	my $dna = $fai->fetch("$chr:$start-$end");
	return $dna;
}
