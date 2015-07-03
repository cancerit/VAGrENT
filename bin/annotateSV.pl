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
use English qw(-no_match_vars);
use warnings FATAL => 'all';
use Carp;
use Const::Fast qw(const);
use Try::Tiny;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

use List::Util qw(first);
use File::Temp qw(tempfile);
use Try::Tiny qw(try catch);

use FindBin qw($Bin);
use lib "$Bin/../lib";

#use Vcf;

use Sanger::CGP::Vagrent;

# InputObject
use Sanger::CGP::VagrentSV::VagrentSV;

#td, tandem duplication
#del, interstitial deletion
#ut, unbalanced translocation
#bt, balanced translocation
#two_jump, two polymerase jumps
#fb, foldback
#inv, inversion

const my @SV_TYPE => qw(td del ut bt two_jump fb inv);

my $header_already_parsed = 0;

eval {
	my $options = option_builder();
  Sanger::CGP::VagrentSV::VagrentSV->new($options);
   
  1;
} or do {
  warn "EVAL_ERROR: $EVAL_ERROR\n" if($EVAL_ERROR);
  warn "CHILD_ERROR: $CHILD_ERROR\n" if($CHILD_ERROR);
  warn "OS_ERROR: $OS_ERROR\n" if($OS_ERROR);
  croak 'A problem occurred';
};


sub make_process_log {
  my ($opts) = @_;
  my $params;
  foreach my $key (keys %$opts){
    $params->{$key} = $opts->{$key} if(defined $opts->{$key});
  }
  return {'key'=>'vcfProcessLog',
			InputVCF => $opts->{'i'},
			InputVCFSource => 'AnnotateVcf.pl',
			InputVCFVer => Sanger::CGP::Vagrent->VERSION,
			InputVCFParam => $params,
		};
}

sub option_builder {
  my ($factory) = @_;

  my %opts = ();

  my $result = &GetOptions (
    'h|help' => \$opts{'help'},
    'v|version' => \$opts{'version'},
    'i|input=s' => \$opts{'input'},
    'o|output=s' => \$opts{'output'},
    'c|cache=s' => \$opts{'cache'},
    't|tabix' => \$opts{'tabix'},
    'p|process=n' => \$opts{'process'},
    'sp|species=s' => \$opts{'species'},
    'as|assembly=s' => \$opts{'assembly'},
    'g|genome=s' => \$opts{'genome'},
  );

  pod2usage() if($opts{'help'});

  if($opts{'version'}){
    print 'Version: '.Sanger::CGP::Vagrent->VERSION."\n";
    exit;
  }

	pod2usage(q{'-i' must be defined}) unless($opts{'input'});
	pod2usage(q{'-i' must exist}) unless(-e $opts{'input'});
	pod2usage(q{'-i' must be a file}) unless(-f $opts{'input'});
	pod2usage(q{'-i' is an empty file}) unless(-s $opts{'input'});

	pod2usage(q{'-c' must be defined}) unless($opts{'cache'});
	pod2usage(q{'-c' must exist}) unless(-e $opts{'cache'});
	pod2usage(q{'-c' must be a file}) unless(-f $opts{'cache'});
	pod2usage(q{'-c' is an empty file}) unless(-s $opts{'cache'});

	pod2usage(q{'-o' must be defined}) unless($opts{'output'});

  return \%opts;
}

__END__

=head1 NAME

AnnotateVcf.pl - Annotate variants - Sub/Snp, Insertion, Deletion, ComplexInDel

=head1 SYNOPSIS

AnnotateVcf.pl [-h] [-t] -i <IN_FILE> -o <OUT_FILE> -c <VAGRENT_CACHE_FILE> [-sp <SPECIES> -as <GENOME_VERSON>]

  General Options:

    --help      (-h)      Brief documentation

    --input     (-i)      Input BedPE file 

    --output    (-o)      Output vcf file (plain text, add -t for zip and index)

    --cache     (-c)      Vagrent reference data cache file

  Conditional (can be specified if missing from the input VCF file)

    --species   (-sp)     Species

    --assembly  (-as)     Genome assembly version

  Optional

    --version   (-v)      Output version number

    --process   (-p)      ID_PROCESS that generated this file

    --tabix     (-t)      bgzip and tabix index the output file (will generate the .gz version of the -o option)

=cut
