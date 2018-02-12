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

use Getopt::Long;
use Pod::Usage;
use Try::Tiny;
use Capture::Tiny qw(capture);
use FindBin qw($Bin);
use lib "$Bin/../lib";

use File::Type;
use Const::Fast qw(const);

use Data::Dumper;

const my $FASTA_HEADER_GREP => 'zgrep "^>" %s';

try {
	my $opts = option_builder();
	my $transcripts = getTranscripts($opts);
	writeOutput($opts,$transcripts);
} catch {
  warn "An error occurred while building reference support files\:\n\t$_"; # not $@
};

sub writeOutput {
  my ($opts,$trans) = @_;
  my $fh;
  if(defined $opts->{'o'}){
    open($fh,'>'.$opts->{'o'}) or die("Unable to open output file for writing: ".$opts->{'o'});
  } else {
    $fh = *STDOUT;
  }
  foreach my $t(@$trans){
    print $fh $t,"\n";
  }
  if(defined $opts->{'o'}){
    close $fh;
  }
}

sub getTranscripts {
  my $opts = shift;
  my $good_transcripts;
  foreach my $seq_file (@{$opts->{'f'}}){
    my $cmd = sprintf $FASTA_HEADER_GREP, $seq_file;
    my ($stdout, $stderr, $exit) = capture {
      system($cmd);
    };
    if($exit > 0){
      croak($stderr);
    }
    foreach my $rec (split /\n/,$stdout){
      my ($transcript,$biotype) = parseLine($rec);
      next unless hasGoodBiotype($opts,$biotype);
      # have to store the transcript accessions twice, with and without the version
      # depending on the Ensembl version or species, the accession in the fasta file may 
      # or may not have the transcript version on the end.
      # some species gff/gtf files include the the version in the accession, some don't
      push(@$good_transcripts,$transcript);
      $transcript =~ s/\.\d+?$//;
      push(@$good_transcripts,$transcript);      
    }
  }
  return $good_transcripts;
}

sub hasGoodBiotype {
  my ($opts,$biotype) = @_;
  foreach my $type (@{$opts->{'b'}}){
    return 1 if lc($biotype) eq lc($type);
  }
  return 0;
}

sub parseLine {
  my $line = shift;
  my @cols = split /\s/, $line;
  my $full = join('|',@cols);
  my $transcript = shift @cols; 
  $transcript =~ s/^>//;
  my $biotype = undef;
  foreach my $c (@cols){
    if($c =~ m/^transcript_biotype/){
      (undef,$biotype) = split(':',$c);
    }
  }
  return ($transcript,$biotype);
}

sub option_builder {
	my ($factory) = @_;

	my %opts = ();

	my $result = &GetOptions (
		'h|help' => \$opts{'h'},
		'f|fasta=s@' => \$opts{'f'},
		'o|outfile=s' => \$opts{'o'},
		'b|biotypes=s@' => \$opts{'b'},
  );
  pod2usage() if($opts{'h'});
  pod2usage('Output file must not already exist') if defined $opts{'o'} && -e $opts{'o'};
  if(scalar(@{$opts{'f'}}) > 0){
    foreach my $f(@{$opts{'f'}}){
      pod2usage("Unable to read the transcript sequence file: $f") unless -e $f && -r $f;
    }
  } else {
    pod2usage('Must specify one or more fasta sequence file as input');
  }
  pod2usage('Must specify at least one biotype') unless scalar @{$opts{'b'}} > 0 && defined($opts{'b'}->[0]);
  return \%opts;
}

__END__

=head1 NAME

Admin_EnsemblTranscriptFilter.pl - Generates filtered list of Transcript accessions

=head1 SYNOPSIS

Admin_EnsemblTranscriptFilter.pl [-h] [-f /path/to/ensembl.fa] [-o /path/to/output.list] [-b biotype]

  Required Options: 
    
    --fasta        (-f)     Ensembl fasta file from Ensembl's FTP site, can be specified multiple times.

    --biotypes     (-b)     Ensembl transcript biotypes to filter for, can be specified multiple times.

  General Options:

    --help         (-h)     Brief documentation.

    --output       (-o)     Output file, optional will default to stdout.

=cut


