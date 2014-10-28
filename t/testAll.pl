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
