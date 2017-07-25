##########LICENCE##########
# Copyright (c) 2017 Genome Research Ltd.
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
use Test::More;
use Try::Tiny qw(try catch);
use Const::Fast qw(const);
use FindBin qw($Bin);
use File::Temp;
use File::Path qw(rmtree);

my $script_path = "$Bin/../bin";

const my $ANNOT => '%s %s/AnnotateVcf.pl -c %s -i %s -o %s 2>&1 /dev/null';

my $perl = $^X;

my $ft_dir = File::Temp->newdir(CLEANUP=>0);
my $tmpdir = $ft_dir->dirname;


my $tmp_a = "$tmpdir/a.vcf";
my $cmd_a = sprintf $ANNOT, $^X, $script_path,
                            "$Bin/../testData/stable.cache.gz",
                            "$Bin/../testData/stable.vcf.gz",
                            $tmp_a;
my $exit_a = system($cmd_a);
ok(!$exit_a);

# the result can be the same, or be different this ensures we get a failure
for(0..3) {
  my $tmp_b = "$tmpdir/b.vcf";
  my $cmd_b = sprintf $ANNOT, $^X, $script_path,
                              "$Bin/../testData/stable.cache.gz",
                              "$Bin/../testData/stable.vcf.gz",
                              $tmp_b;
  my $exit_b = system($cmd_b);
  ok(!$exit_b);

  my $exit_diff = system("diff -I '^#' -I '^#' $tmp_a $tmp_b");

  if($exit_diff) {
    BAIL_OUT("Diff found errors see: $tmpdir, you will need to cleanup.");
  }
}

rmtree($tmpdir);

done_testing();
