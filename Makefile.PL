#!/usr/bin/perl

use ExtUtils::MakeMaker;

WriteMakefile(
  NAME          => 'Vagrent',
  VERSION_FROM  => 'lib/Sanger/CGP/Vagrent.pm',
  EXE_FILES     => [qw( bin/Admin_EnsemblGtf2CacheConverter.pl
                        bin/Admin_EnsemblReferenceFileGenerator.pl
                        bin/Admin_EnsemblTranscriptFastaFilter.pl
                        bin/AnnotateVcf.pl)],
  PREREQ_PM     => {
                    'Attribute::Abstract' => '1.02',
                    'Bio::DB::Sam' => '1.39',
                    'Bio::Root::Version' => '1.006924',# this is the only way to get Bio::Perl to behave
                                                       # and return the version correctly
                    'Class::Inspector' => '1.28',
                    'Clone' => '0.37',
                    'Config::IniFiles' => '2.83',
                    'Const::Fast' => '0.014',
                    'File::ShareDir' => '1.102',
                    'HTTP::Date' => '6.02',
                    'HTTP::Message' => '6.06',
                    'List::MoreUtils' => '0.33',
#                    'List::Util' => '1.39',
                    'Log::Log4perl' => '1.44',
                    'LWP' => '6.08',
                    'overloading' => '0.02',
                    'parent' => '0.228',
                    'Pod::Escapes' => '1.06',
                    'Pod::Simple' => '3.28',
                    'Sort::Key' => '1.33',
                    'Sub::Exporter::Progressive' => '0.001011',
                    'Sub::Name' => '0.09',
                    'Try::Tiny' => '0.22',
                    'URI' => '1.64',
                   }
);
