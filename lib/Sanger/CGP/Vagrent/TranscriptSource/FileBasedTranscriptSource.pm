package Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource;

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

use Carp;
use Log::Log4perl;
use Data::Dumper;
use Const::Fast qw(const);

use Bio::DB::HTS;
use Bio::DB::HTS::Tabix;
use Set::IntervalTree;

use Sanger::CGP::Vagrent::Data::Transcript;
use Sanger::CGP::Vagrent::Data::Exon;
use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource);

my $log = Log::Log4perl->get_logger(__PACKAGE__);

const my $SEARCH_BUFFER => 10000;

const my $GENE_GRAB_TEMPLATE => 'tabix %s %s | cut -s -f5 | uniq |';
const my $TRANSCRIPT_GRAB_TEMPLATE => q{tabix %s %s | awk '$5 == "%s" { print $0;}' | cut -s -f7 |};

1;

sub getTranscripts {
  my ($self,$gp) = @_;
  unless(defined($gp) && $gp->isa('Sanger::CGP::Vagrent::Data::AbstractGenomicPosition')){
    $log->error("Did not recieve a Sanger::CGP::Vagrent::Data::AbstractGenomicPosition object");
    return undef;
  }
  my $trans = $self->_getTranscriptsFromCache($gp);
  return @$trans if defined $trans;
  return;
}

sub getTranscriptsForNextGeneInDumpRegion {
	my ($self) = @_;
  my $gr = $self->getDumpRegion;
	$log->logcroak('must define a dump region before trying to loop over the genes') unless defined $gr;
	if($gr->getMinPos == $gr->getMaxPos && $gr->getMaxPos == 0){
    # full sequence region scan
    $self->{_dumpInfo}->{_fullSeq} = 1;
  } else {
    $self->{_dumpInfo}->{_fullSeq} = 0;
  }
  $self->_populateGeneList($gr) unless defined $self->{_dumpInfo}->{_geneList};
  $self->_getDumpTranscriptStore($gr) unless defined $self->{_dumpInfo}->{_transcriptStore};

my $geneListSize =  scalar @{$self->{_dumpInfo}->{_geneList}};
  for(my $i = $self->{_dumpInfo}->{_counter} ; $i <= $geneListSize ; $i++){
		$self->{_dumpInfo}->{_counter}++;
		next unless(defined $self->{_dumpInfo}->{_geneList}->[$i]);
    my $transList = $self->_getTranascriptsForGeneName($gr,$self->{_dumpInfo}->{_geneList}->[$i]);
    next unless(defined $transList && scalar(@$transList) > 0);
    return @$transList;
	}
  return undef;
}

sub _getTranascriptsForGeneName {
  my ($self,$gr,$genename) = @_;
  my $out;
  my $cmd = sprintf $TRANSCRIPT_GRAB_TEMPLATE, $self->{_cache}, $self->_generateLocationString($gr), $genename;
  open my $fh, $cmd or $log->logcroak("unable to run transcript lookup for gene name $genename: $cmd");
  while(<$fh>){
    my $VAR1;
    eval $_;
    push(@$out,$VAR1);
  }
  close $fh;
  return $out;
}

sub _getDumpTranscriptStore {
  my ($self,$gr) = @_;


}



sub _populateGeneList {
  my ($self,$gr) = @_;
  my $cmd = sprintf $GENE_GRAB_TEMPLATE, $self->{_cache}, $self->_generateLocationString($gr);
  open my $fh, $cmd or $log->logcroak('unable to run gene name lookup for region dump');
  @{$self->{_dumpInfo}->{_geneList}} = <$fh>;
  close $fh;
  chomp @{$self->{_dumpInfo}->{_geneList}};
  $self->{_dumpInfo}->{_counter} = 0;
}

sub _generateLocationString {
  my ($self,$gr) = @_;
  return $gr->getChr.':'.$gr->getMinPos.'-'.$gr->getMaxPos;
}

sub _tabix_to_interval_tree {
  my ($self, $chr) = @_;
  return 1 if defined $self->{_cache_iTree}->{$chr};

  my $full_tree = {};
  if ($self->{_sorted}) {
    $self->{_cache_iTree} = $full_tree;
  }
  else {
    $full_tree = $self->{_cache_iTree};
  }

  my %collated;
  $self->{_cache_tbx} = Bio::DB::HTS::Tabix->new(filename => $self->{_cache}) unless defined $self->{_cache_tbx};
  my $iter = $self->{_cache_tbx}->query_full($chr);
  return 1 unless defined $iter;
    while(my $line = $iter->next) {
    my ($chr, $s, $e, $object) = (split /\t/, $line)[0,1,2,6];
    # +1 on end to convert to 1 bases, tabix module would handle this
    my $this_loci = sprintf '%s:%d:%d', $chr, $s, $e+1;
    push @{$collated{sprintf '%s:%d:%d', $chr, $s, $e+1}}, $object;
  }

  my $chr_tree = Set::IntervalTree->new();
  for my $loci(keys %collated) {
    my ($chr, $s, $e) = split ':', $loci;
    $chr_tree->insert($collated{$loci}, $s, $e);
    delete $collated{$loci};
  }
  $full_tree->{$chr} = $chr_tree;
  return 1;
}

sub _getTranscriptsFromCache {
  my ($self,$gp) = @_;
  my $chr = $gp->getChr();
  $self->_tabix_to_interval_tree($chr);
  my $min;
  my $max = $gp->getMaxPos + $SEARCH_BUFFER;
  if($gp->getMinPos() < $SEARCH_BUFFER){
    $min = 0;
  } else {
    $min = ($gp->getMinPos - $SEARCH_BUFFER);
  }
  my @data = ();
  @data = @{$self->{_cache_iTree}->{$chr}->fetch($min,$max)};
  return undef unless(@data);
  my @out;
  for my $overlap(@data){
    for my $item(@{$overlap}) {
      unless(ref $item) { # turn string into object
        my $VAR1;
        eval $item;
        $VAR1->{_cdnaseq} = $self->_getTranscriptSeq($VAR1);
        $item = $VAR1;
      }
      push @out, $item;
    }
  }
  @out = sort _sort_itree @out;
  return \@out;
}

sub _sort_itree {
  if($a->{_genomicminpos} != $b->{_genomicminpos}) {
    return $a->{_genomicminpos} <=> $b->{_genomicminpos};
  }
  if($a->{_genomicmaxpos} != $b->{_genomicmaxpos}) {
    return $a->{_genomicmaxpos} <=> $b->{_genomicmaxpos};
  }
  return 0;
}

sub _init {
	my ($self, %vars) = @_;
  $self->_setCacheFile($vars{'cache'});
  $self->{_sorted} = $vars{'sorted'};
}

sub _setCacheFile {
  my ($self,$cache) = @_;

  unless(-e $cache && -f $cache && -r $cache){
    $log->logcroak("Specified cache file is unreadable: $cache");
  }
  my $cache_index = $cache .".tbi";
  unless(-e $cache_index && -f $cache_index && -r $cache_index){
    $log->logcroak("cache index file is unreadable: $cache_index");
  }
  my $fa = $cache;
  $fa =~ s/\.cache.+$/.fa/;
  unless(-e $fa && -f $fa && -r $fa){
    $log->logcroak("cache fasta file is unreadable: $fa");
  }
  my $fai = $fa . ".fai";
  unless(-e $fai && -f $fai && -r $fai){
    $log->logcroak("cache fasta index file is unreadable: $fai");
  }
  $self->{_cache} = $cache;
  $self->{_cache_fa} = $fa;
}

sub _getTranscriptSeq {
  my ($self,$trans) = @_;
  unless(defined $self->{_fai_obj}){
    $self->{_fai_obj} = Bio::DB::HTS::Fai->load($self->{_cache_fa});
  }
  return $self->{_fai_obj}->fetch($trans->getAccession);
}
