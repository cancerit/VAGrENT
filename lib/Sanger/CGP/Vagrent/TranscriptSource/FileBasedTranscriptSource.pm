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

use Bio::DB::Sam;

# prevent the init warnings from Tabix
BEGIN {
  $SIG{__WARN__} = sub {warn $_[0] unless( $_[0] =~ m/^Subroutine Tabix.* redefined/)};
};

use Tabix;

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

sub _getTranscriptsFromCache {
  my ($self,$gp) = @_;
  $self->{_cache_tbx} = Tabix->new('-data' => $self->{_cache}) unless defined $self->{_cache_tbx};
  my $min;
  my $max = $gp->getMaxPos + $SEARCH_BUFFER;
  if($gp->getMinPos() < $SEARCH_BUFFER){
    $min = 0;
  } else {
    $min = ($gp->getMinPos - $SEARCH_BUFFER) - 1;
  }
  my $res = $self->{_cache_tbx}->query($gp->getChr(),$min,$max);
  return undef unless defined $res;
  my $out = undef;
  if(defined $res->get){
    while(my $ret = $self->{_cache_tbx}->read($res)){
      my $raw = (split("\t",$ret))[6];
      my $VAR1;
      eval $raw;
      $VAR1->{_cdnaseq} = $self->_getTranscriptSeq($VAR1);
      push @$out, $VAR1;
    }
  }
  return $out;
}

sub _init {
	my $self = shift;
  my %vars = @_;
  foreach my $k(keys(%vars)){
    if($k eq 'cache'){
      $self->_setCacheFile($vars{$k});
    }
  }
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
    $self->{_fai_obj} = Bio::DB::Sam::Fai->load($self->{_cache_fa});
  }
  return $self->{_fai_obj}->fetch($trans->getAccession);
}
