#!/usr/bin/env perl

use strict;
use English qw(-no_match_vars);
use warnings FATAL => 'all';
use Carp;

use Getopt::Long;
use Pod::Usage;
use Try::Tiny;

use Data::Dumper;
use File::Type;
use Readonly qw(Readonly);

use Bio::DB::Sam;

use Sanger::CGP::Vagrent::Data::Transcript;
use Sanger::CGP::Vagrent::Data::Exon;

$Data::Dumper::Indent = 0;

Readonly my @TEXT_TYPES => qw(application/octet-stream);

Readonly my @ZIP_TYPES => qw(application/x-gzip);

Readonly my $PROTEIN_CODING => 'protein_coding';
Readonly my $EXON_TYPE => 'exon';
Readonly my $CDS_START_TYPE => 'start_codon';
Readonly my $CDS_STOP_TYPE => 'stop_codon';
Readonly my $CDS_TYPE => 'CDS';

try {
	my $opts = option_builder();
	my $ccds = undef;
	my $lookup = parseFaiFile($opts->{'fai'});
	$ccds = parseCCDSFile($opts->{'c'}) if defined $opts->{'c'};
	convertGtf($opts,$lookup,$ccds);
} catch {
  warn "An error occurred while building reference support files\:\n\t$_"; # not $@
};

sub convertGtf {
	my ($opts,$lookup,$ccds) = @_;
	my $type = getInputFileType($opts->{'g'});
	my $fasta = Bio::DB::Sam::Fai->load($opts->{'f'});
	my ($inFh,$outFh);
	if($type eq 'text') {
		open $inFh, '<'.$opts->{'g'} || croak("unable to open gtf input file:".$opts->{'g'});
	} elsif($type eq 'zip'){
		open $inFh, 'zcat '.$opts->{'g'}.' |' || croak("unable to open gtf input file:".$opts->{'g'});
	}
	open $outFh, ">".$opts->{'o'} || croak("unable to open output file:".$opts->{'o'});
	my $c = 0;
	my $wip;
	my $currentChr = undef;
	while(<$inFh>){
		my @data = split /\s+/;
		my $acc;
		if($data[11] =~ m/^\"(ENST\d+)\"\;$/){
			$acc = $1;
			next unless exists $lookup->{$acc};
		} else {
			next;
		}
		$currentChr = $data[0] unless defined $currentChr;
		if($currentChr ne $data[0]) {
			# changed chr, any transcript on previous chr must be finished
			foreach my $t (processCompleted($opts,$fasta,$lookup,$ccds,$wip,$currentChr)){
				writeTranscript($outFh,$t,$wip->{$t->getAccession});
				delete $wip->{$t->getAccession};
			}
			$currentChr = $data[0];
		}
		my ($ltype,$ttype,@store) = ($data[2],$data[1],$data[0],$data[3],$data[4],$data[6],unquoteValue($data[13]));
		push(@{$wip->{$acc}->{'lines'}->{$ltype}},\@store);
		unless(exists $wip->{$acc}->{'type'}){
			$c++;
			$wip->{$acc}->{'type'} = $ttype;
			$wip->{$acc}->{'acc'} = $acc;
			$wip->{$acc}->{'gene'} = unquoteValue($data[15]);
			
		}
		if($ltype eq $CDS_TYPE && !defined $wip->{$acc}->{'protacc'}){
			$wip->{$acc}->{'protacc'} = unquoteValue($data[21]);
		}
	}
	foreach my $t (processCompleted($opts,$fasta,$lookup,$ccds,$wip,$currentChr)){
		writeTranscript($outFh,$t,$wip->{$t->getAccession});
		delete $wip->{$t->getAccession};
	}
	close $outFh;
	close $inFh;
}

sub processCompleted {
	my ($opts,$fasta,$lookup,$ccds,$wip,$currentChr) = @_;
	my @trans;
	foreach my $acc (keys %$wip){
		#warn Dumper($wip->{$acc});
		if(exists $wip->{$acc}->{'lines'}->{$EXON_TYPE} && $wip->{$acc}->{'lines'}->{$EXON_TYPE}->[0]->[0] eq $currentChr){
			my $t = convertTranscript($opts,$fasta,$lookup->{$acc},$ccds->{$acc},$wip->{$acc});
			push @trans, $t;
		}
	}
	return @trans;
}

sub convertTranscript {
	my ($opts,$fasta,$tlength,$ccds,$data) = @_;
	my $type = getGeneTypeForTranscript($data->{'type'});
	my @exons;
	my $strand = undef;
	my $rnaLengthSum = 0;
	my @sortedExonLines = sort {$a->[4] <=> $b->[4]} @{$data->{'lines'}->{$EXON_TYPE}};
	foreach my $rawE (@sortedExonLines){
		my $elength = ($rawE->[2] - $rawE->[1]) + 1;
		my $rmin = $rnaLengthSum + 1;
		my $rmax = $rnaLengthSum + $elength;
		$rnaLengthSum += $elength;
		my $convE = Sanger::CGP::Vagrent::Data::Exon->new(
							species => $opts->{'s'},
							genomeVersion => $opts->{'v'},
							chr => $rawE->[0],
							minpos => $rawE->[1],
							maxpos => $rawE->[2],
							rnaminpos => $rmin,
							rnamaxpos => $rmax,);
												
		unless(defined $strand){
			if($rawE->[3] eq '+'){
				$strand = 1;
			} elsif($rawE->[3] eq '-'){
				$strand = -1;
			} else {
				croak "Expecting strand of + or -, recieved: ".$rawE->[3];
			}
		}					
		push(@exons,$convE);
	}
	
	my ($protAcc,$protAccVer,$cdsMin,$cdsMax,$cdsPhase);
	
	if($type eq Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType()){
		# protein coding
		$protAcc = $data->{'protacc'};
		$protAccVer = 1;
		my $cdsLength = 0;
		my @sortedCdsLines = sort {$a->[4] <=> $b->[4]} @{$data->{'lines'}->{$CDS_TYPE}};
		foreach my $cds (@sortedCdsLines){
			$cdsLength += ($cds->[2] - $cds->[1]) + 1;
		}
		# find CDS start
		my $fivePrimeUtrExonLength = 0;
		foreach my $e (@sortedExonLines){
			if($sortedCdsLines[0]->[4] eq $e->[4]){
				if($strand > 0){
					$cdsMin = $fivePrimeUtrExonLength + (($sortedCdsLines[0]->[1] - $e->[1]) + 1);
					$cdsMax = ($cdsMin + $cdsLength) - 1;
				} else {
					$cdsMin = $fivePrimeUtrExonLength + (($e->[2] - $sortedCdsLines[0]->[2]) + 1);
					$cdsMax = ($cdsMin + $cdsLength) - 1;
				}
				last;
			} else {
				$fivePrimeUtrExonLength += ($e->[2] - $e->[1]) + 1;
			}
		}
		$cdsMax += 3 if exists $data->{'lines'}->{$CDS_STOP_TYPE};
		$cdsPhase = 0;
	} else {
		# non-coding
		$protAcc = undef;
		$protAccVer = undef;
		$cdsMin = undef;
		$cdsMax = undef;
		$cdsPhase = -1;
	}

	my @sortedExons = sort {$a->getMinPos <=> $b->getMinPos} @exons;
	my $convT = Sanger::CGP::Vagrent::Data::Transcript->new(
																		db => 'Ensembl',
																		dbversion => $opts->{'d'},
																		acc => $data->{'acc'},
																		accversion => 1,
																		proteinacc => $protAcc,
																		proteinaccversion => $protAccVer,
																		ccds => $ccds,
																		genename => $data->{'gene'},
																		genetype => $type,
																		strand => $strand,
																		cdnaseq => $fasta->seq($data->{'acc'},1,$tlength),
																		cdsminpos => $cdsMin,
																		cdsmaxpos => $cdsMax,
																		cdsphase => $cdsPhase,
																		genomicminpos => $sortedExons[0]->getMinPos,
																		genomicmaxpos => $sortedExons[-1]->getMaxPos,
																		exons => \@exons,);	
	return $convT;
}

sub writeTranscript {
	my ($fh,$t,$rawT) = @_;
	print $fh join("\t",$rawT->{'lines'}->{$EXON_TYPE}->[0]->[0],$t->getGenomicMinPos - 1,
											$t->getGenomicMaxPos,$t->getAccession,$t->getGeneName,length $t->getcDNASeq);																			
	$t->{_cdnaseq} = undef;
	print $fh "\t",Dumper($t),"\n";								
}

sub getGeneTypeForTranscript {
	my ($t) = @_;
	if($t eq 'protein_coding'){
		return Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType();
	} elsif($t eq 'miRNA'){
		return Sanger::CGP::Vagrent::Data::Transcript::getMicroRnaType();
	} elsif($t eq 'lincRNA'){
		return Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType();
	} elsif($t eq 'snoRNA'){
		return Sanger::CGP::Vagrent::Data::Transcript::getSnoRnaType();
	} elsif($t eq 'snRNA'){
		return Sanger::CGP::Vagrent::Data::Transcript::getSnRnaType();
	} elsif($t eq 'rRNA'){
		return Sanger::CGP::Vagrent::Data::Transcript::getRRnaType();
	} else {
		croak("Unhandled type: $t");
	}
}

sub unquoteValue {
	my $val = shift;
	$val =~ s/[\";]//g;
	return $val;
}

sub parseCCDSFile {
	my $file = shift;
	my $out = undef;
	open my $fh, "<$file" || croak ("unable to open $file for reading");
	while(<$fh>){
		my @data = split /\s+/, $_;
		next unless $data[4] =~ m/^ENST/;
		next unless $data[2] == 1;
		$out->{$data[4]} = $data[0];
	}
	close $fh;
	return $out;
}

sub parseFaiFile {
	my $fai = shift;
	open my $fh, "<$fai" || croak("unable to open fai file: $fai");
	my $out;
	while (<$fh>) {
		my ($enst,$length) = split /\s+/;
		$out->{$enst} = $length;
	}
	return $out;
}

sub getInputFileType {
	my $infile = shift;
	my $fileType = File::Type->new->checktype_filename($infile);
	foreach my $t (@TEXT_TYPES){
		return 'text' if($t eq $fileType);
	}
	foreach my $t (@ZIP_TYPES) {
		return 'zip' if($t eq $fileType);
	}
	return $fileType;
}

sub option_builder {
	my ($factory) = @_;
	my %opts = ();
	my $result = &GetOptions (
		'h|help' => \$opts{'h'},
		'f|fasta=s' => \$opts{'f'},
		'g|gtf=s' => \$opts{'g'},
		'o|output=s' => \$opts{'o'},
		's|species=s' => \$opts{'s'},
		'v|genome=s' => \$opts{'v'},
		'd|database=s' => \$opts{'d'},
		'c|ccds=s' => \$opts{'c'},
  );

  pod2usage() if($opts{'h'});
  pod2usage('Must specify the input fasta reference file') unless defined $opts{'f'} && -e $opts{'f'};
  my $fai = $opts{'f'} .'.fai';
  pod2usage("fasta index reference file not found: $fai") unless -e $fai;
  $opts{'fai'} = $fai;
  pod2usage('Must specify the input gtf reference file') unless defined $opts{'g'} && -e $opts{'g'};
  pod2usage('Must specify the output file to use') unless defined $opts{'o'};
  pod2usage('Must specify the species') unless defined $opts{'s'};
  pod2usage('Must specify the genome version') unless defined $opts{'v'};
  pod2usage('Must specify the ensembl core database version') unless defined $opts{'d'};
  if(defined($opts{'c'})){
  	pod2usage('CCDS file unreadable') unless -e $opts{'c'} && -r $opts{'c'};
  }
  return \%opts;
}

__END__

=head1 NAME

Admin_EnsemblGtf2CacheConverter.pl - Generates the Vagrent cache file from the ensembl GTF file

=head1 SYNOPSIS

Admin_EnsemblGtf2CacheConverter.pl [-h] [-f /path/to/ensembl.fa] [-g /path/to/ensembl.gtf] [-o /path/to/output.file] [-s human] [-v GRCh37] [-d homo_sapiens_core_74_37p] [-c /path/to/CCDS2Sequence.version.txt]

  General Options:

    --help         (-h)     Brief documentation

    --fasta        (-f)     Fasta file containing the transcripts to be imported, fai file must also be present.
    
    --gtf          (-g)     Ensembl GTF file for converting
    
    --output       (-o)     Output file
    
    --species      (-s)     Species (ie human, mouse)
    
    --genome       (-v)     Genome version (ie GRCh37, GRCm38)
    
    --database     (-d)     Ensembl core database version number (ie homo_sapiens_core_74_37p)
    
    --ccds         (-c)     (Optional, but strongly advised) The CCDS2Sequence file from the relevant CCDS release, see http://www.ncbi.nlm.nih.gov/CCDS
   
=cut