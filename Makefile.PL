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


use ExtUtils::MakeMaker;
use File::ShareDir::Install;

install_share module => 'Sanger::CGP::Vagrent', 'share';

WriteMakefile(
  NAME          => 'Vagrent',
  LICENSE       => 'agpl_3', # http://search.cpan.org/~dagolden/CPAN-Meta-2.142690/lib/CPAN/Meta/Spec.pm#license
  VERSION_FROM  => 'lib/Sanger/CGP/Vagrent.pm',
  EXE_FILES     => [qw( bin/Admin_EnsemblGtf2CacheConverter.pl
                        bin/Admin_EnsemblReferenceFileGenerator.pl
                        bin/Admin_EnsemblTranscriptFastaFilter.pl
                        bin/Admin_GeneRegionBedDumper.pl
                        bin/AnnotateVcf.pl)],
  PREREQ_PM     => {
                    'Attribute::Abstract' => '1.02',
                    'Bio::DB::HTS' => '1.11',
                    'Bio::Root::Version' => '1.006924',# this is the only way to get Bio::Perl to behave
                                                       # and return the version correctly
                    'Capture::Tiny' => 0.24,
                    'Config::IniFiles' => '2.83',
                    'Const::Fast' => '0.014',
                    'File::ShareDir' => '1.102',
                    'File::Temp' => '0.2304',
                    'File::Type' => '0.22',
                    'IPC::System::Simple' => '1.25',
                    'List::Util' => '1.41',
                    'Log::Log4perl' => '1.44',
                    'LWP::Simple' => '6.00',
                    'Pod::Usage' => '1.64',
                    'Sort::Key' => '1.33',
                    'TAP::Harness' => '3.33',
                    'Try::Tiny' => '0.22',
                   }
);

package MY;
use File::ShareDir::Install qw( postamble );
