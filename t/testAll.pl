#!/usr/bin/perl

use strict;

use TAP::Harness;
use Cwd qw(abs_path);
use Data::Dumper;

my @testScripts = qw(substitution.t insertion.t deletion.t complex.t);

my $prog_path = abs_path($0);
$prog_path =~ m/^(.*)testAll\.pl$/;
my $stub = $1;

my @toRun;

foreach my $t(@testScripts){
	my $ref = [$stub.$t,$t];
	push(@toRun,$ref);
}

my $harness = TAP::Harness->new({verbosity => 0,timer => 1,show_count => 1,lib=>\@INC});

$harness->runtests(@toRun);
