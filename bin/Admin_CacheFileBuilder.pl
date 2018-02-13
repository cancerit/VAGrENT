#!/usr/bin/perl

##########LICENCE##########
# Copyright (c) 2018 Genome Research Ltd.
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
use English qw(-no_match_vars);
use warnings FATAL => 'all';
use Carp;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Getopt::Long;
use Pod::Usage;
use Try::Tiny;
use Capture::Tiny qw(capture capture_stderr);
use List::Util qw(first);
use Scalar::Util qw(looks_like_number);

use Data::Dumper;
use File::Type;
use File::Temp qw(tempdir tempfile);
use File::Spec;
use File::Copy qw(copy);
use Const::Fast qw(const);

use Bio::DB::HTS::Faidx;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::GFF;

use Sanger::CGP::Vagrent::Data::Transcript;
use Sanger::CGP::Vagrent::Data::Exon;

$Data::Dumper::Indent = 0;

const my @TEXT_TYPES => qw(text/plain application/octet-stream);

const my @ZIP_TYPES => qw(application/x-gzip);

const my $TAG_GENE_ID => 'gene_id';
const my $TAG_TRANSCRIPT_ID => 'transcript_id';
const my $TAG_PROTEIN_ID => 'protein_id';
const my $TAG_VERSION => 'version';
const my $TAG_PARENT => 'Parent';
const my $TAG_EXON_ORDER_GFF => 'rank';
const my $TAG_EXON_ORDER_GTF => 'exon_number';
const my @TAG_EXON_ORDER => ($TAG_EXON_ORDER_GTF,$TAG_EXON_ORDER_GFF);
const my $TAG_GENE_BIOTYPE => 'gene_biotype';
const my $TAG_BIOTYPE => 'biotype';
const my @TAG_BIOTYPE => ($TAG_BIOTYPE,$TAG_GENE_BIOTYPE);

const my $TAG_GENE_NAME_GFF => 'Name';
const my $TAG_GENE_NAME_GTF => 'gene_name';
const my @TAG_GENE_NAME => ($TAG_GENE_NAME_GFF,$TAG_GENE_NAME_GTF);

const my @TAG_CCDS => qw(ccds_id ccdsid);

const my $PRIMARY_UTR5P => 'five_prime_UTR';
const my $PRIMARY_CDS => 'CDS';
const my $PRIMARY_UTR3P => 'three_prime_UTR';
const my $PRIMARY_STOP => 'stop_codon';

const my @MITO_ALIASES => qw(M MT Mito MtDNA);

const my @PRIMARY_EXON => qw(exon);
const my @PRIMARY_REGIONS => ($PRIMARY_UTR5P,$PRIMARY_CDS,$PRIMARY_STOP,$PRIMARY_UTR3P);

const my $PROTEIN_CODING => 'protein_coding';
const my $EXON_TYPE => 'exon';
const my $CDS_START_TYPE => 'start_codon';
const my $CDS_STOP_TYPE => 'stop_codon';
const my $CDS_TYPE => 'CDS';

const my $UNZIP => 'gunzip -c %s > %s';
const my $BGZIP => 'bgzip -c %s > %s';
const my $SORT_N_BGZIP => 'sort -k 1,1 -k 2n,2 -k 3n,3 %s | bgzip > %s';
const my $TABIX => 'tabix -p bed %s';

const my $SEQ_LOCATION => '%s:%s-%s';

const my $OUT_CACHE_FILE => 'vagrent.%s.%s.%s.cache.gz';
const my $OUT_SEQ_FILE => 'vagrent.%s.%s.%s.fa';

my $tmpDir = tempdir("VagrentBuildCacheXXXXX", TMPDIR => 1, CLEANUP => 1);

try {
	my $opts = option_builder();
	my $transList = parseTranscriptList($opts->{'transcripts'});
	my $seqFiles = openSeqFiles($tmpDir,@{$opts->{'f'}});
	my $chrList = parseFaiFile($opts->{'fai'});
	my $ccdsList = parseCCDSfile($opts->{'c'});
	my ($tmp_cache_file,$tmp_seq_file) = processFeatureFile($opts,$tmpDir,$chrList,$seqFiles,$transList,$ccdsList);
	finaliseOutputFiles($opts,$tmp_cache_file,$tmp_seq_file);
	
} catch {
  die "An error occurred while building reference support files\:\n\t$_"; # not $@
};

sub finaliseOutputFiles {
  my ($opts,$tmp_cache,$tmp_seq) = @_;
  my $cache_file = makeCacheFilePath($opts);
  my $seq_file = makeSeqFilePath($opts);
  
  copy($tmp_seq,$seq_file) or croak("unable to copy sequence file to desination: $!");
  capture_stderr {my $fa_index = Bio::DB::HTS::Faidx->new($seq_file)};

  
  my $snz_cmd = sprintf($SORT_N_BGZIP, $tmp_cache, $cache_file);
  my ($snz_stdout, $snz_stderr, $snz_exit) = capture {system($snz_cmd);};
  croak('Error in sorting and bgzipping the cache file: '.$snz_stderr) if $snz_exit > 0;
  my $tabix_cmd = sprintf $TABIX, $cache_file;
  my ($tbx_stdout, $tbx_stderr, $tbx_exit) = capture {system($tabix_cmd);};
  croak('Error in indexing the cache file: '.$tbx_stderr) if $tbx_exit > 0;
}

sub makeCacheFilePath {
  my $opts = shift;
  my $file = sprintf $OUT_CACHE_FILE, $opts->{'sp'},$opts->{'as'},$opts->{'d'};
  return File::Spec->catfile($opts->{'o'},$file);
}

sub makeSeqFilePath {
  my $opts = shift;
  my $file = sprintf $OUT_SEQ_FILE, $opts->{'sp'},$opts->{'as'},$opts->{'d'};
  return File::Spec->catfile($opts->{'o'},$file);
}

sub processFeatureFile {
  my ($opts,$tmpDir,$chrList,$seqFiles,$transList,$ccdsList) = @_;
  my ($cache_fh, $cache_tmp_file) = tempfile('VagrentBuildCacheXXXXX', DIR => $tmpDir, SUFFIX => '.cache');
  my ($seq_tmp_fh, $fa_tmp_file) = tempfile('VagrentBuildCacheXXXXX', DIR => $tmpDir, SUFFIX => '.fa');
  my $seq_fh = Bio::SeqIO->new(-fh => $seq_tmp_fh, -format => 'fasta');
  my $features = openFeatureFile($opts->{'features'});
  my $c = 0;
  my $store = undef;
  my $gene_names = undef;
  my $current_chr = undef;
  my $current_chr_translation = undef;
  while(my $feat = $features->next_feature()) {
    my $chr = $feat->seq_id();
    if(!defined $current_chr){
      $current_chr = $chr;
      $current_chr_translation = checkChromosome($current_chr,$chrList);
      if (!defined $current_chr_translation){
        print "Sequence $current_chr not found in fai file: no alternative found, skipping\n" if $opts->{'debug'};
      } elsif ($current_chr ne $current_chr_translation){
        print "Sequence $current_chr not found in fai file: using $current_chr_translation as alternative\n" if $opts->{'debug'};
      }
    }
    next if $current_chr eq $chr && !defined $current_chr_translation;
    my $tid = undef;
    if($feat->has_tag($TAG_GENE_ID)){
      my ($gene_id) = $feat->get_tag_values($TAG_GENE_ID);
      unless(exists $gene_names->{$gene_id}){
        foreach my $tag (@TAG_GENE_NAME){
          next unless $feat->has_tag($tag);
          my ($gn) = $feat->get_tag_values($tag);
          $gene_names->{$gene_id} = $gn;
          last;
        }
      }
    }
    if($feat->has_tag($TAG_TRANSCRIPT_ID)){
      ($tid) = $feat->get_tag_values($TAG_TRANSCRIPT_ID);
    } else {
      if($feat->has_tag($TAG_PARENT)){
        my ($parent) = $feat->get_tag_values($TAG_PARENT);
        if($parent =~ m/^transcript:(.+)$/){
          $tid = $1;
        } else {
          next;
        }
      } else {
        next;
      }
    } 
    next unless defined $tid;
    if($current_chr ne $chr) {
      # chromosome has ended, have to process transcripts
      if(scalar(keys %$store) > 0){
        foreach my $obj (@{processCompleteTranscripts($opts,$store,$current_chr_translation,$gene_names,$chrList,$ccdsList,$seqFiles)}){
          writeOutput($obj,$cache_fh,$seq_fh);
        }
        $store = undef;
      }
      $current_chr = $chr;
      $current_chr_translation = checkChromosome($current_chr,$chrList);
      if (!defined $current_chr_translation){
        print "Sequence $current_chr not found in fai file: no alternative found, skipping\n" if $opts->{'debug'};
      } elsif ($current_chr ne $current_chr_translation){
        print "Sequence $current_chr not found in fai file: using $current_chr_translation as alternative\n" if $opts->{'debug'};
      }
      next if !defined $current_chr_translation;
    }  
    
    if(exists $store->{$tid}){
      push @{$store->{$tid}}, $feat;
    } elsif(exists $transList->{$tid}) {
      push @{$store->{$tid}}, $feat;
      $c++;
    } else {
      next;
    }
  }
  if(scalar(keys %$store) > 0){
    foreach my $obj (@{processCompleteTranscripts($opts,$store,$current_chr_translation,$gene_names,$chrList,$ccdsList,$seqFiles)}){
      writeOutput($obj,$cache_fh,$seq_fh);
    }
  }
  close $cache_fh;
  close $seq_tmp_fh;
  return ($cache_tmp_file,$fa_tmp_file);
}

sub writeOutput {
  my ($trans,$cache,$fa) = @_;
  my $seq = Bio::Seq->new(-display_id => $trans->getAccession(),-seq => $trans->getcDNASeq);
  my ($e) = $trans->getExons();
  my $cache_row;
  try{
    $cache_row = join("\t",$e->getChr,$trans->getGenomicMinPos - 1,$trans->getGenomicMaxPos,$trans->getAccession,$trans->getGeneName,length $trans->getcDNASeq);
  } catch {
    croak('An error occured outputting transcript '.$trans->getAccession.', please check input data');
  };
  $trans->{_cdnaseq} = undef;
  $cache_row .= "\t".Dumper($trans)."\n";
  
  $fa->write_seq($seq);
  print $cache $cache_row;
}

sub processCompleteTranscripts {
  my ($opts,$store,$chr,$gene_names,$chrList, $ccdsList,$seqFiles) = @_;
  my @out = ();
  foreach my $trans_id (keys %$store){
    my $converted = convertTranscript($opts,$store->{$trans_id},$trans_id,$gene_names,$ccdsList,$seqFiles,$chr);
    push(@out,$converted) if defined $converted;
  }
  return \@out;
}

sub convertTranscript {
  my ($opts,$raw_transcript,$transcript_id,$gene_names,$ccdsList,$seqFiles,$chr) = @_;
  my $type = getGeneTypeForTranscript($opts,$raw_transcript);
  return undef unless defined $type;
  my $rnaLengthSum = 0;
  my $strand = undef;
  my @exons;
  
  my @raw_exons;
  my @raw_cds;
  my $gene_id = undef;
  my $transcript_version = 1;
  foreach my $e(@$raw_transcript){
    if($e->has_tag($TAG_TRANSCRIPT_ID) && $e->has_tag($TAG_VERSION)){
      ($transcript_version) = $e->get_tag_values($TAG_VERSION);
    }
    unless(defined $gene_id){
      if($e->has_tag($TAG_GENE_ID)){
        ($gene_id) = $e->get_tag_values($TAG_GENE_ID);
      } elsif ($e->has_tag($TAG_PARENT)){
        my ($parent) = $e->get_tag_values($TAG_PARENT);
        my ($pt,$pv) = split ':',$parent;
        if($pt eq 'gene'){
          $gene_id = $pv;
        }
      }
    }
    push @raw_exons, $e if first { $_ eq $e->primary_tag } @PRIMARY_EXON;
    push @raw_cds, $e if first { $_ eq $e->primary_tag } @PRIMARY_REGIONS;
  }
  return undef if scalar(@raw_exons) == 0;
  return undef if scalar(@raw_cds) == 0 && $type eq Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType();
  my $gene_name;
  if(exists $gene_names->{$gene_id} && defined $gene_names->{$gene_id}){
    $gene_name = $gene_names->{$gene_id};
  } else {
    $gene_name = $gene_id
  }
  @raw_exons = sort featureOrderSortFunction @raw_exons; # sorting exons into transcript order
  @raw_cds = sort featureOrderSortFunction @raw_cds; # sorting CDS fragments into transcript order
  foreach my $e (@raw_exons){
    my $rmin = $rnaLengthSum + 1;
    my $rmax = $rnaLengthSum + $e->length;
    $rnaLengthSum += $e->length;
    my $convE = Sanger::CGP::Vagrent::Data::Exon->new(
							  species => $opts->{'sp'},
							  genomeVersion => $opts->{'as'},
							  chr => $chr,
							  minpos => $e->start,
							  maxpos => $e->end,
							  rnaminpos => $rmin,
							  rnamaxpos => $rmax,);
		
		if(defined $strand){
		  croak('Inconsistant strand in transcript') unless $strand == $e->strand;
		} else {
		  $strand = $e->strand;
		}
		push @exons, $convE;
  }
  my ($protAcc,$protAccVer,$cdsMin,$cdsMax,$cdsPhase,$ccds);
  if($type eq Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType()){
    ($cdsMin,$cdsMax,$cdsPhase) = calculateCDS($opts,\@raw_cds,\@raw_exons, $rnaLengthSum, $strand,$transcript_id);
    return undef unless(defined $cdsMin && defined $cdsMax && defined $cdsPhase);
    $ccds = getCCDS($raw_transcript,$transcript_id,$ccdsList);
    foreach my $cds(@raw_cds){
      if($cds->has_tag($TAG_PROTEIN_ID)){
        ($protAcc) = $cds->get_tag_values($TAG_PROTEIN_ID);
        last;
      }
    }
    $protAccVer = 1;
    
  } else {
		# non-coding
		$protAcc = undef;
		$protAccVer = undef;
		$cdsMin = undef;
		$cdsMax = undef;
		$cdsPhase = -1;
	}
  
  my @sortedExons = sort {$a->getMinPos <=> $b->getMinPos} @exons;
  my $transcript = Sanger::CGP::Vagrent::Data::Transcript->new(
																		db => 'Ensembl',
																		dbversion => $opts->{'d'},
																		acc => $transcript_id,
																		accversion => $transcript_version,
																		proteinacc => $protAcc,
																		proteinaccversion => $protAccVer,
																		ccds => $ccds,
																		genename => $gene_name,
																		genetype => $type,
																		strand => $strand,
																		cdnaseq => getSequence($transcript_id,$transcript_version,$seqFiles),
																		cdsminpos => $cdsMin,
																		cdsmaxpos => $cdsMax,
																		cdsphase => $cdsPhase,
																		genomicminpos => $sortedExons[0]->getMinPos,
																		genomicmaxpos => $sortedExons[-1]->getMaxPos,
																		exons => \@exons,);
    
  return $transcript; 
}

sub getSequence {
  my ($id,$version,$seqFiles) = @_;
  my $out = undef;
  my $idv = join('.',$id,$version);
  foreach my $fai(@$seqFiles){
    my $loc;
    if($fai->has_sequence($id)){
      my $length = $fai->length($id);
      $loc = sprintf($SEQ_LOCATION,$id,1,$length);
    } elsif ($fai->has_sequence($idv)){
      my $length = $fai->length($idv);
      $loc = sprintf($SEQ_LOCATION,$idv,1,$length);
    } else {
      next;
    }
    $out = $fai->get_sequence_no_length($loc);
  }
  return $out;
}

sub getCCDS {
  my ($raw_transcript,$transcript_id,$ccdsList) = @_;
  my $ccds = undef;
  foreach my $feat(@$raw_transcript){
    foreach my $tag(@TAG_CCDS){
      if($feat->has_tag($tag)){
        ($ccds) = $feat->get_tag_values($tag);
        last;
      }
    }
    last if defined $ccds;
  }
  unless(defined $ccds){
    if(defined $ccdsList){
      if(exists $ccdsList->{$transcript_id}){
        $ccds = $ccdsList->{$transcript_id};
      }
    }
  }
  return $ccds;
}

sub calculateCDS {
  my ($opts,$feats,$exons,$rnaLength, $strand,$transcript_id) = @_;
  my $utr5Length = 0;
  my $cdsLength = 0;
  my $utr3Length = 0;
  my $stopCodonLength = 0;
  my @cds;
  foreach my $f (@$feats){
    if($f->primary_tag eq $PRIMARY_UTR5P){
      $utr5Length += $f->length;
    } elsif($f->primary_tag eq $PRIMARY_CDS){
      $cdsLength += $f->length;
      push @cds, $f;
    } elsif($f->primary_tag eq $PRIMARY_UTR3P){
      $utr3Length += $f->length;
    } elsif($f->primary_tag eq $PRIMARY_STOP){
      $stopCodonLength += $f->length;
    }
  }  
  my $utr5Length_e = 0;
  my $cdsLength_e = 0;
  my $utr3Length_e = 0;
  
  my $cds_first = $cds[0];
  my $cds_last = $cds[-1];
  
  my $break = 0;
  my $mess = '';
  
  foreach my $e(@$exons){
    croak('Recieved undefined exon: '.Dumper($exons)) unless defined($e);
    # remember start = min, end = max.  start and end DO NOT infer orientation in BIOPERL
    if($strand > 0){
      # forward strand start = min = start of feature, end = max = end of feature
      if($e->end < $cds_first->start){
        # exon ends before the CDS starts
        # 5 prime UTR exon
        $utr5Length_e += $e->length;
      } elsif ($e->end < $cds_last->start && $e->start <= $cds_first->start && $e->end >= $cds_first->start){
        # exon ends before last CDS fragment starts
        # and first CDS fragment starts within this exon.
        $utr5Length_e += $cds_first->start - $e->start;
        $cdsLength_e += $e->length - ($cds_first->start - $e->start);
      } elsif ($e->start >= $cds_first->start && $e->end <= $cds_last->end){
        # exon starts after the first CDS fragment starts
        # and exon ends before the last CDS fragment ends
        # exon is entirely within the CDS
        $cdsLength_e += $e->length; 
      } elsif ($cds_first->start == $cds_last->start && $cds_first->end == $cds_last->end && $e->start <= $cds_first->start && $e->end >= $cds_last->end){
        # single CDS fragment contained entirely within this exon  
        $utr5Length_e += $cds_first->start - $e->start;
        $cdsLength_e += $e->length - (($cds_first->start - $e->start) + ($e->end - $cds_last->end));
        $utr3Length_e += $e->end - $cds_last->end;
      } elsif ($e->start > $cds_first->end && $e->start <= $cds_last->end && $e->end >= $cds_last->end){  
        # exon starts after the first CDS fragment ends
        # and the last CDS fragment ends within this exon
        $cdsLength_e += $e->length - ($e->end - $cds_last->end);
        $utr3Length_e += $e->end - $cds_last->end;
      } elsif ($e->start > $cds_last->end){
        # exon starts after the CDS ends
        $utr3Length_e += $e->length;
      } else {
        warn "\n",join(' - ',$utr5Length,$cdsLength,$stopCodonLength,$utr3Length),"\n";
        warn join(' - ',$utr5Length_e,$cdsLength_e,0,$utr3Length_e),"\n";
        warn "\n";
        foreach (@cds){warn "CDS : ",Dumper($_),"\n"};
        warn "\n";
        foreach (@$exons){warn "EXON : ",Dumper($_),"\n"};
        warn "\n";
        foreach ($cds_first,$cds_last){warn "CDS BOUNDRIES : ",Dumper($_),"\n"};   
        warn "\n";
        warn "PROC : ",Dumper($e),"\n";
        croak('unhandled exon');
      }
      
            
    } elsif($strand < 0){
      # reverse strand start = min = end of feature, end = max = start of feature
      if($e->start > $cds_first->end){
        # exon ends (start) before the CDS starts (end)
        # 5 prime UTR exon
        $utr5Length_e += $e->length;
      } elsif ($e->start > $cds_last->end && $e->end >= $cds_first->end && $e->start <= $cds_first->end){
        # exon ends (start) before last CDS fragment starts (end)
        # and first CDS fragment starts (end) within this exon
        $utr5Length_e += $e->end - $cds_first->end;
        $cdsLength_e += $e->length - ($e->end - $cds_first->end);
      } elsif ($e->end <= $cds_first->end && $e->start >= $cds_last->start){
        # exon starts (end) after the first CDS fragment starts (end)
        # and exon ends (start) before the last CDS fragment ends (start)
        # exon is entirely within the CDS
        $cdsLength_e += $e->length;
      } elsif ($cds_first->end == $cds_last->end && $cds_first->start == $cds_last->start && $e->end >= $cds_first->end && $e->start <= $cds_last->start){
        # single CDS fragment contained entirely within this exon   
        $utr5Length_e += $e->end - $cds_first->end;
        $cdsLength_e += $e->length - (($e->end - $cds_first->end) + ($cds_last->start - $e->start));
        $utr3Length_e += $cds_last->start - $e->start;
      } elsif ($e->end < $cds_first->start && $e->end >= $cds_last->start && $e->start <= $cds_last->start){
        # exon starts (end) after the first CDS fragment ends (start)
        # and the last CDS fragment ends in this exon
        $cdsLength_e += $e->length - ($cds_last->start - $e->start);
        $utr3Length_e += $cds_last->start - $e->start;
      } elsif ($e->end < $cds_last->start){
        # exons starts (end) after CDS ends (start)
        $utr3Length_e += $e->length;
      } else {
        warn "\n",join(' - ',$utr5Length,$cdsLength,$stopCodonLength,$utr3Length),"\n";
        warn join(' - ',$utr5Length_e,$cdsLength_e,0,$utr3Length_e),"\n";
        warn "\n";
        foreach (@cds){warn "CDS : ",Dumper($_),"\n"};
        warn "\n";
        foreach (@$exons){warn "EXON : ",Dumper($_),"\n"};
        warn "\n";
        foreach ($cds_first,$cds_last){warn "CDS BOUNDRIES : ",Dumper($_),"\n"};   
        warn "\n";
        warn "PROC : ",Dumper($e),"\n";
        croak('unhandled exon');
      }
    } else {
      croak('Cannot process strand: ',$strand);
    }
  }
  
  # length calculation double check
  
  if($utr5Length + $utr3Length > 0){
    # UTR features present in input file, check all 3 values, skip if data is inconsistant
    return (undef,undef,undef) unless($utr5Length == $utr5Length_e && $cdsLength == $cdsLength_e && $utr3Length == $utr3Length_e);
  } else {
    # only compare the CDS length, skip if data is inconsistant
    return (undef,undef,undef) unless($cdsLength == $cdsLength_e);
  }
  
  my $cds_min = $utr5Length_e + 1;
  my $cds_max = $utr5Length_e + $cdsLength_e + $stopCodonLength;
  my $cds_phase = undef;
  if(looks_like_number($cds_first->frame)){
    if($cds_first->frame > 0){
      $cds_phase = 3 - $cds_first->frame;
    } else {
      $cds_phase = 0;
    }
  } else {
    return (undef,undef,undef);
  }
  return ($cds_min,$cds_max,$cds_phase);
}

sub getGeneTypeForTranscript {
  my ($opts,$t) = shift;
  
  my $type = undef;
  my $type_tag = undef;
  
  foreach my $f(@$t){ 
    foreach my $tag(@TAG_BIOTYPE){
      if($f->has_tag($tag)){
        my ($raw_type) = $f->get_tag_values($tag);
        if(!defined $type){
          $type = $raw_type;
          $type_tag = $tag;
        } elsif($type ne $raw_type){
          croak('Transcript has inconsistant biotype');
        }            
      }
    }
  }

  if($type eq 'protein_coding'){
		return Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType();
	} elsif($type eq 'miRNA'){
		return Sanger::CGP::Vagrent::Data::Transcript::getMicroRnaType();
	} elsif($type eq 'lincRNA'){
		return Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType();
	} elsif($type eq 'snoRNA'){
		return Sanger::CGP::Vagrent::Data::Transcript::getSnoRnaType();
	} elsif($type eq 'snRNA'){
		return Sanger::CGP::Vagrent::Data::Transcript::getSnRnaType();
	} elsif($type eq 'rRNA'){
		return Sanger::CGP::Vagrent::Data::Transcript::getRRnaType();
	} else {
	  if($type_tag eq $TAG_GENE_BIOTYPE){
	    # this can happen, there can be a discrepancy between gene and transcript biotype in ensembl.  
	    # If they are only reporting the gene biotype and its not one we can handle we have to skip the transcript.
	    return undef;
	  } else {
	    if($opts->{'debug'}){
        foreach (@$t){
          warn Dumper($_),"\n";
        }
      }
      croak("Unhandled biotype: $type");
	  }
	}
}

sub featureOrderSortFunction {
  my $tag = undef;
  
  if($a->strand != $b->strand){
    croak('Features strands not consistent');
  }
  
  if($a->strand == 1 && $a->strand == $b->strand){
    return $a->start <=> $b->start;
  } elsif($a->strand == -1 && $a->strand == $b->strand){
    return $b->start <=> $a->start;
  } else {
    croak('unknown strand type: ',$a->strand);
  }
}

sub checkChromosome {
  my ($chr,$chrList) = @_;
  my $use_chr = undef;
  if(first { $_ eq $chr } @$chrList){
    $use_chr = $chr;
  } else {
    # chromosome not in fai file
    # lets try chr prefix discrepancies 
    if($chr =~ m/^chr(.+)$/){
      # chromosome starts with chr, try removing and search again.
      $use_chr = $1 if first { $_ eq $1 } @$chrList;
    } elsif($chr !~ m/^chr/){
      # chromosome doesn't start with chr, try adding
      $use_chr = 'chr'.$chr if first { $_ eq 'chr'.$chr } @$chrList;
    }
    if(!defined $use_chr){
      # still not found a translation
      # check for mitochondria, its a special case
      my $mito_check;
      if($chr =~ m/^chr(.+)$/){
        $mito_check = $1;
      } else {
        $mito_check = $chr;
      }
      if(first { $_ eq $mito_check } @MITO_ALIASES){
        # this is mitochondria
        # check fai file for aliases, with and without chr prefixes.
        foreach my $m(@MITO_ALIASES){
          my $chrM = 'chr'.$m;
          if(first { $_ eq $m } @$chrList){
            $use_chr = $m;
            last;
          } elsif (first { $_ eq $chrM } @$chrList){
            $use_chr = $chrM;
            last;
          }
        }
      }
    }
  }
  return $use_chr;
}

sub openFeatureFile {
  my $file = shift;
  my $fileType = File::Type->new->checktype_filename($file);
  my $fh;
  if(first { $_ eq $fileType } @ZIP_TYPES){
    # feature file is zipped
    open $fh, "zcat $file |" or croak("unable to open feature input file: $file");
  } elsif(first { $_ eq $fileType } @TEXT_TYPES){
    # feature file is plain text
    open $fh, "<$file" or croak("unable to open feature input file: $file");
  } else {
    croak('Unsupported file type: '.$fileType);
  }
  my $gff_version;
  if($file =~ m/\.gtf(?:\.gz)?/){
    $gff_version = 2.5;
  } elsif($file =~ m/\.gff3(?:\.gz)?/){
    $gff_version = 3;
  } else {
    croak('Unable to determine format of '.$file);
  } 
  
  return Bio::Tools::GFF->new(-fh => $fh, -gff_version => $gff_version);
}

sub openSeqFiles {
  my ($tmpDir,@files) = @_;
  my $out;
  foreach my $file (@files){
    my $fileType = File::Type->new->checktype_filename($file);
    if(first { $_ eq $fileType } @ZIP_TYPES){
      # File zipped, going to have to re-zip it with bgzip
      my $unzipped = makeTmpFilePath($tmpDir,'.fa');
      my $zipped = makeTmpFilePath($tmpDir,'.fa.gz');
      my $unzip_cmd = sprintf($UNZIP, $file, $unzipped);
      my $bgzip_cmd = sprintf($BGZIP, $unzipped, $zipped);
      my ($uz_stdout, $uz_stderr, $uz_exit) = capture {system($unzip_cmd);};
      croak($uz_stderr) if $uz_exit > 0;
      my ($bgz_stdout, $bgz_stderr, $bgz_exit) = capture {system($bgzip_cmd);};
      croak($bgz_stderr) if $bgz_exit > 0;
      
      $file = $zipped;
    } elsif(first { $_ eq $fileType } @TEXT_TYPES){
      # plain text file, needs zipping
      my $zipped = makeTmpFilePath($tmpDir,'.fa.gz');
      my $bgzip_cmd = sprintf($BGZIP, $file, $zipped);
      my ($bgz_stdout, $bgz_stderr, $bgz_exit) = capture {system($bgzip_cmd);};
      croak($bgz_stderr) if $bgz_exit > 0;
      $file = $zipped;
    } else {
      croak('Unsupported file type: '.$fileType)
    }
    # HTS call prints to stderr if it has to generate an fai file, which it will here
    # feels dirty but this makes it quiet.
    capture_stderr {push @$out, Bio::DB::HTS::Faidx->new($file)};
  } 
  return $out;
}

sub parseCCDSfile {
  my $file = shift;
  return undef unless defined $file;
  open my $fh, "<$file" or croak("unable to open CCDS file: ".$file);
  my @data = <$fh>;
  close $fh;
  my $out;
  foreach my $d(@data){
    next if $d !~ m/^CCDS\d/;
    my @cols = split /\s+/,$d;
    my $ccds = $cols[0];
    my $transcript_id = $cols[4];
    $out->{$transcript_id} = $ccds;
  }
  return $out;
}

sub parseFaiFile {
  my $fai = shift;
  return undef unless defined $fai;
  open my $fh, "<$fai" or croak("unable to open fai file: ".$fai);
  my @data = <$fh>;
  close $fh;
  my @out;
  foreach my $d(@data){
    my ($name) = split /\s+/,$d;
    push @out, $name;
  }
  return \@out;
}

sub parseTranscriptList {
  my $file = shift;
  open my $fh, "<$file" or croak('Unable to open transcript list file: '.$file);
  my $out;
  while(<$fh>){
    my $line = $_;
    chomp $line;
    $out->{$line} = 1;
  }
  close $fh;
  return $out;
}

sub makeTmpFilePath {
  my ($tmpDir,$ext) = @_;
  my $file;
  (undef, $file) = tempfile('VagrentBuildCacheXXXXX', DIR => $tmpDir, OPEN => 0, SUFFIX => $ext);
  return $file;
}

sub option_builder {
	my ($factory) = @_;
	my %opts = ();
	my $result = &GetOptions (
		'h|help' => \$opts{'h'},
		'f|fasta=s@' => \$opts{'f'},
		'gf|features=s' => \$opts{'features'},
		't|transcripts=s' => \$opts{'transcripts'},
		'o|output=s' => \$opts{'o'},
		'sp|species=s' => \$opts{'sp'},
		'as|assembly=s' => \$opts{'as'},
		'd|database=s' => \$opts{'d'},
		'c|ccds=s' => \$opts{'c'},
		'fai|fai=s' => \$opts{'fai'},
		'debug|debug' => \$opts{'debug'},
  );

  pod2usage() if($opts{'h'});
  if(scalar(@{$opts{'f'}}) > 0){
    foreach my $f(@{$opts{'f'}}){
      pod2usage("Unable to read the transcript sequence file: $f") unless -e $f && -r $f;
    }
  } else {
    pod2usage('Must specify one or more fasta sequence file as input');
  }

  pod2usage('Must specify the input gtf/gff3 reference file') unless defined $opts{'features'} && -e $opts{'features'};
  pod2usage('Must specify transcript list file') unless defined $opts{'transcripts'} && -e $opts{'transcripts'};
  pod2usage('Must specify the output directory') unless defined $opts{'o'} && -e $opts{'o'} && -d $opts{'o'} ;
  pod2usage('Must specify the species') unless defined $opts{'sp'};
  pod2usage('Must specify the genome version') unless defined $opts{'as'};
  pod2usage('Must specify the data version') unless defined $opts{'d'};
  if(defined($opts{'c'})){
  	pod2usage('CCDS file unreadable') unless -e $opts{'c'} && -r $opts{'c'};
  }
  if(defined($opts{'fai'})){
  	pod2usage('Unable to read the fai file: '.$opts{'fai'}) unless -e $opts{'fai'} && -r $opts{'fai'};
  }
  return \%opts;
}



__END__

=head1 NAME

Admin_CacheFileBuilder.pl - Generates the Vagrent reference data set from the supplied GFF£/GTF and Fasta files 

=head1 SYNOPSIS

Admin_CacheFileBuilder.pl [-h] [-f /path/to/sequence.fa] [-gff /path/to/annotation.gff] [-sp human] [-as GRCh37] [-d homo_sapiens_core_74_37p] [-c /path/to/CCDS2Sequence.version.txt]

  Required Options:

    --output       (-o)     Output directory

    --fasta        (-f)     Fasta file containing Transcript sequences, can be specified multiple times.

    --features     (-gf)   GFF3 or GTF file containing genomic annotation
    
    --transcripts  (-t)     Transcript accession file, (simple list, one accession per line)

    --species      (-sp)    Species (ie human, mouse)

    --assembly     (-as)    Assembly version (ie GRCh37, GRCm38)

    --data_version (-d)     Reference version number (Ensemble example: homo_sapiens_core_74_37p)
    
    Other Options:
    
    --help         (-h)     Brief documentation   

    --ccds         (-c)     (Recommended) The CCDS2Sequence file from the relevant CCDS release, see http://www.ncbi.nlm.nih.gov/CCDS
    
    --fai          (-fai)   (Recommended) The samtools fasts index file (.fai) for your reference genome
                              This is the reference genome that your bam and vcf files will be mapped to
                              
    --debug        (-debug) This turns on the debug output, useful if you are having issues with a specific Ensembl build                  
=cut