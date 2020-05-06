#!/usr/bin/env perl

use strict;

use FindBin qw($Bin);
use lib "$Bin/../lib"; # the magic development hack

use Bio::DB::HTS::Faidx;
use Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource;
use Sanger::CGP::Vagrent::Data::GenomicRegion;

if(@ARGV < 4) {
  print STDERR "USAGE: genome.fa vagrent.cache.gz bases_context loci_in.tab [pad_gene]\n\n";
  print STDERR "  genome.fa     = reference genome fasta with fai index.\n";
  print STDERR "  bases_context = Number of surrounding context bases (N +/-).\n";
  print STDERR "  loci_in.tab   = tab delimited file beginning: 'chr\\tpos...\\n'.\n";
  print STDERR "  pad_gene      = 1 indicates that coordinates for gene overlap are padded with 'bases_context'\n";
  print STDERR "\nIf first line of input prefixed # understood to be header and new headings added\n";
  print STDERR "Ouput is original line with '\\tWT\\tCONTEXT\\tGENE\\tSTRAND' appended\n";
  exit 1;
}

my ($fasta, $cache, $bp_context, $loci_in, $pad_gene) = @ARGV;

$pad_gene = $pad_gene ? 1 : 0;

my $faidx = Bio::DB::HTS::Faidx->new($fasta);
my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new(cache => $cache, sorted => 1);

open my $ifh, '<', $loci_in or die "Can't open $loci_in for reading $!";
while(my $l = <$ifh>) {
  chomp $l;
  if($l =~ m/^#/) {
    # only a header if first line
    printf "%s\n", join("\t", $l, qw(WT CONTEXT GENE STRAND)) if($. == 1);
    next;
  }
  my ($chr, $pos) = split /\t/, $l;
  my $gr;
  if($pad_gene) {
    $gr = Sanger::CGP::Vagrent::Data::GenomicRegion->new(species => '.', genomeVersion => '.', chr => $chr, minpos => $pos - $bp_context, maxpos => $pos + $bp_context);
  }
  else {
    $gr = Sanger::CGP::Vagrent::Data::GenomicRegion->new(species => '.', genomeVersion => '.', chr => $chr, minpos => $pos, maxpos => $pos);
  }
  my $context = $faidx->get_sequence_no_length(sprintf '%s:%d-%d', $chr, $pos - $bp_context, $pos + $bp_context);
  my $wt = substr($context, $bp_context, 1);
  my ($gene, $strand) = qw(No_gene No_gene);
  my @trans = $ts->getTranscripts($gr);
  if(@trans) {
    $gene = $trans[0]->getAccession;
    $strand = $trans[0]->getStrand;
  }
  printf "%s\n", join("\t", $l, $wt, $context, $gene, $strand);
}
close $ifh;
