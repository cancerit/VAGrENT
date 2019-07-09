#!/usr/bin/perl

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
use warnings;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use English qw( -no_match_vars );
use Carp;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Capture::Tiny qw(capture);
use Const::Fast qw(const);

use Sanger::CGP::Vagrent::Data::GenomicRegion;
use Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource;
use Sanger::CGP::Vagrent::GenomicRegionDumper::SubstitutionGenomicRegionDumper;
use Sanger::CGP::Vagrent::GenomicRegionDumper::IndelGenomicRegionDumper;
use Sanger::CGP::Vagrent::IO::GenomicRegionWriter::BedWriter;

const my @RAW_DUMP_FILES => ('raw.gene_regions.bed', 'raw.exon_regions.sub.bed','raw.codingexon_regions.sub.bed','raw.exon_regions.indel.bed','raw.codingexon_regions.indel.bed');
const my @SORTED_DUMP_FILES => ('sorted.gene_regions.bed', 'sorted.exon_regions.sub.bed','sorted.codingexon_regions.sub.bed','sorted.exon_regions.indel.bed','sorted.codingexon_regions.indel.bed');
const my @FINAL_DUMP_FILES => ('gene_regions.bed', 'exon_regions.sub.bed','codingexon_regions.sub.bed','exon_regions.indel.bed','codingexon_regions.indel.bed');

my $options = option_builder(@ARGV);

my @rawFileNames;
my @sortedFileNames;
my @finalFileNames;
my $outHandles;

my $prefix = $options->{'species'}.'.'.$options->{'assembly'}.'.';

for(my $i = 0 ; $i < 5 ; $i++){
  push @rawFileNames, $prefix.$RAW_DUMP_FILES[$i];
  push @sortedFileNames, $prefix.$SORTED_DUMP_FILES[$i];
  push @finalFileNames, $prefix.$FINAL_DUMP_FILES[$i];
}

eval {
  my @faiData = parseFai($options);
  my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => $options->{'cache'});
  my $c = 0;
  my $dumpers;
  foreach my $outFile (@rawFileNames) {
    $c++;
    open my $fh,">$outFile" or croak('Unable to open file for writing: '.$outFile);
    push(@$outHandles,$fh);
    my $writer = Sanger::CGP::Vagrent::IO::GenomicRegionWriter::BedWriter->new(fh => $fh);
    if($c < 4){
			push(@$dumpers,Sanger::CGP::Vagrent::GenomicRegionDumper::SubstitutionGenomicRegionDumper->new('transcriptSource' => $ts, 'writer' => $writer));
		} else {
			push(@$dumpers,Sanger::CGP::Vagrent::GenomicRegionDumper::IndelGenomicRegionDumper->new('transcriptSource' => $ts, 'writer' => $writer));
		}
  }

  my $dc = scalar(@$dumpers);

  foreach my $c (@faiData){
 		my @d = split(/\s+/,$c);
    my $chr = Sanger::CGP::Vagrent::Data::GenomicRegion->new(
 				'species'				=> $options->{'species'},
 				'genomeVersion' => $options->{'assembly'},
 				'chr' 					=> $d[0],
 				'minpos'				=> 1,
 				'maxpos'				=> $d[1],);

    for(my $i = 0 ; $i < $dc ; $i++){
 			if($i==0){
 				$dumpers->[$i]->dumpGeneRegions($chr);
 			} elsif($i % 2 == 0){
 				$dumpers->[$i]->dumpCodingExonRegions($chr);
 			} else {
				$dumpers->[$i]->dumpExonRegions($chr);
 			}
 		}
  }
  foreach my $fh(@$outHandles){
	  close $fh;
	}
  for(my $i = 0 ; $i < scalar(@rawFileNames) ; $i++){
		my $rawFile = $rawFileNames[$i];
		my $sortedFile = $sortedFileNames[$i];
		my $finalFile = $finalFileNames[$i];
		my $sortCmd = 'bedtools sort -i '.$rawFile.' > '.$sortedFile;
		my $mergeCmd = 'bedtools merge -i '.$sortedFile.' > '.$finalFile;

    my ($sortOUT,$sortERR,$sortEXIT) = capture{ system($sortCmd) };
    croak("bedtools sort command failed : $sortCmd\n $sortERR") unless $sortEXIT == 0;

    my ($mergeOUT,$mergeERR,$mergeEXIT) = capture{ system($mergeCmd) };
    croak("bedtools merge command failed : $mergeCmd\n $mergeERR") unless $mergeEXIT == 0;

		croak("unable to remove $rawFile") unless unlink $rawFile;
    croak("unable to remove $sortedFile") unless unlink $sortedFile;

		if(defined $options->{'index'}){
		  my $sortfile = $finalFile . ".sort";
      my $gzfile = $finalFile .".gz";

      my $sortCmd = 'sort -k 1,1 -k2,3n '.$finalFile.' > '.$sortfile;
      my $zipCmd = 'bgzip -c '.$sortfile.' > '.$gzfile;
      my $tabixCmd = 'tabix -p bed '.$gzfile;

      my ($sortOUT,$sortERR,$sortEXIT) = capture{ system($sortCmd) };
      croak("sort command failed : $sortCmd\n $sortERR") unless $sortEXIT == 0;

      my ($zipOUT,$zipERR,$zipEXIT) = capture{ system($zipCmd) };
      croak("zip command failed : $zipCmd\n $zipERR") unless $zipEXIT == 0;

      my ($tabixOUT,$tabixERR,$tabixEXIT) = capture{ system($tabixCmd) };
      croak("tabix command failed : $tabixCmd\n $tabixERR") unless $tabixEXIT == 0;

      croak("unable to remove $finalFile") unless unlink $finalFile;
      croak("unable to remove $sortfile") unless unlink $sortfile;
		}
	}

  1;
} or do {
  warn "EVAL_ERROR: $EVAL_ERROR\n" if($EVAL_ERROR);
  warn "CHILD_ERROR: $CHILD_ERROR\n" if($CHILD_ERROR);
  warn "OS_ERROR: $OS_ERROR\n" if($OS_ERROR);
  croak 'A problem occurred';
};

sub parseFai {
  my $opts = shift;
  my $faiFH;
  open($faiFH, "<".$opts->{'fai'}) or die('Unable to open fai file: '.$opts->{'fai'});
  my @faiData = <$faiFH>;
  close $faiFH;
  return @faiData;
}


sub option_builder {
	my @args = @_;

	my %opts = ();

	my $result = &GetOptions (
		'h|help' => \$opts{'help'},
    'v|version' => \$opts{'version'},
		'd|debug' => \$opts{'debug'},
    'i|build_index' => \$opts{'index'},
		'sp|species=s' => \$opts{'species'},
    'as|assembly=s' => \$opts{'assembly'},
    'c|cache=s' => \$opts{'cache'},
		'fai|fai_file=s' => \$opts{'fai'},
	);

  pod2usage() unless(scalar(@args));

	pod2usage() if($opts{'help'});

  if($opts{'version'}){
    print 'Version: '.Sanger::CGP::Vagrent->VERSION."\n";
    exit;
  }

  pod2usage(q{species must be defined}) unless(defined $opts{'species'});
  pod2usage(q{assembly must be defined}) unless(defined $opts{'assembly'});

  pod2usage(q{fai file must be defined}) unless($opts{'fai'});
  pod2usage(q{fai file must exist}) unless(-e $opts{'fai'});
  pod2usage(q{fai file must be a file}) unless(-f $opts{'fai'});
  pod2usage(q{fai file is an empty file}) unless(-s $opts{'fai'});

  pod2usage(q{cache file must be defined}) unless($opts{'cache'});
  pod2usage(q{cache file must exist}) unless(-e $opts{'cache'});
	pod2usage(q{cache file must be a file}) unless(-f $opts{'cache'});
	pod2usage(q{cache file is an empty file}) unless(-s $opts{'cache'});

	return \%opts;
}


__END__

=head1 NAME

Admin_GeneRegionBedDumper.pl - Creates bed files containing regions Vagrent would consider annotatable

=head1 SYNOPSIS

Admin_GeneRegionBedDumper.pl [-h] -sp <SPECIES> -as <ASSEMBLY> -c <VAGRENT_> -fai <GENOME_FAI_FILE>

  General Options:

    --help        (-h)      Brief documentation

    --species     (-sp)     Species

    --assembly    (-as)     Genome assembly version

    --fai_file    (-fai)    Genome fai file, for chromosome names and lengths

    --cache       (-c)      Vagrent reference data cache file

  Optional

    --version     (-v)      Output version number

    --build_index (-i)      Generate a zipped and tabix indexed version of the bed file


  Examples:

    EnsemblGeneRegionBedFileDumper.pl -s human -g GRCh37 -fai /path/to/genome.fai

=cut
