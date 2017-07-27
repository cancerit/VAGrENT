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
is($exit_a, 0, 'Constructed base comparison VCF');

# the result can be the same, or be different this ensures we get a failure
for(1..4) {
  my $tmp_b = "$tmpdir/b.vcf";
  my $cmd_b = sprintf $ANNOT, $^X, $script_path,
                              "$Bin/../testData/stable.cache.gz",
                              "$Bin/../testData/stable.vcf.gz",
                              $tmp_b;
  my $exit_b = system($cmd_b);
  is($exit_b, 0, 'Constructed comparison VCF No. '.$_);

  my $exit_diff = system("diff -I '^#' -I '^#' $tmp_a $tmp_b");

  if($exit_diff) {
    BAIL_OUT("Diff found errors see: $tmpdir, you will need to cleanup.");
  }
  ok('Completed stability iteraction '.$_);
}

# now check all instances of 'VC=' do not have ':' if they do something got messed up
my $bad_vc_message = check_vc($tmp_a);
if($bad_vc_message) {
  notok($bad_vc_message);
}
else {
  ok(q{All 'VC=' entries successfully summarised.});
}

rmtree($tmpdir);

done_testing();

sub check_vc {
  my $in = shift;
  my $message;
  open my $FH, '<', $in or die $!;
  while(my $l = <$FH>) {
    next if $l =~ m/^#/;
    chomp $l;
    my @F = split /\t/, $l;
    if($F[7] =~ m/VC=([^;]+)/) {
      my $vc = $1;
      next if($vc !~ m/:/);
      $message = 'VC not summarised: '.$vc;
      last;
    }
  }
  close $FH;
  return $message;
}
