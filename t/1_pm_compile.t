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

# this is a catch all to ensure all modules do compile

use strict;
use Data::Dumper;
use Test::More;
use List::Util qw(first);
use File::Find;
use Cwd;
use Try::Tiny qw(try finally);
use File::Spec;
use Const::Fast qw(const);

use FindBin qw($Bin);
my $lib_path = "$Bin/../lib";

# Add modules here that cannot be instantiated (should be extended and have no 'new')
# or need a set of inputs - these should be tested in own test script
const my @USE_SKIP => qw(  );
const my @NEW_SKIP => qw( Sanger::CGP::Vagrent::Data::AbstractGenomicPosition
                          Sanger::CGP::Vagrent::Data::ComplexIndel
                          Sanger::CGP::Vagrent::Data::Deletion
                          Sanger::CGP::Vagrent::Data::Exon
                          Sanger::CGP::Vagrent::Data::GenomicRegion
                          Sanger::CGP::Vagrent::Data::Insertion
                          Sanger::CGP::Vagrent::Data::Substitution
                          Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource );

#warn "\n###################### WARNING #######################\n# 1_pm_compile.pl DISABLED until vagrent core update #\n###################### WARNING #######################\n\n";
#ok(1,'Here to show willing');
#done_testing();
#exit;

my $init_cwd = getcwd;

my @modules;
try {
  chdir($lib_path);
  find({ wanted => \&build_module_set, no_chdir => 1 }, './');
} finally {
  chdir $init_cwd;
  die "The try block died with: @_\n" if(@_);
};

for my $mod(@modules) {
  next if( first {$mod eq $_} @USE_SKIP );
  use_ok($mod) or BAIL_OUT("Unable to 'use' module $mod");
}

for my $mod(@modules) {
  ok($mod->VERSION, "Check version inheritance exists ($mod)");
  if($mod->can('new')) { # only try new on things that have new defined
    new_ok($mod) unless( first {$mod eq $_} (@USE_SKIP, @NEW_SKIP) );
  }
}

done_testing();

sub build_module_set {
  if($_ =~ m/\.pm$/) {

    my ($dir_str,$file) = (File::Spec->splitpath( $_ ))[1,2];
    $file =~ s/\.pm$//;
    my @dirs = File::Spec->splitdir( $dir_str );
    shift @dirs;
    push @modules, (join '::', @dirs).$file;
  }
}
