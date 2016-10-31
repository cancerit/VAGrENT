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
use Cwd qw(getcwd abs_path);

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Getopt::Long;
use Pod::Usage;
use Try::Tiny;
use LWP::Simple;
use IPC::System::Simple qw(run);
use Net::FTP;
use File::Copy;
use File::Path qw(make_path remove_tree);
use File::Temp qw(tempdir tempfile);
use Cwd qw(abs_path);

use Const::Fast qw(const);
use Data::Dumper;

const my @ENSMBL_REF_FILE_EXTENTIONS => qw(cdna.all.fa.gz ncrna.fa.gz);
const my $FASTA_FILTER_SCRIPT => 'Admin_EnsemblTranscriptFastaFilter.pl';
const my $GTF_CONVERSION_SCRIPT => 'Admin_EnsemblGtf2CacheConverter.pl';
const my $FILTERED_FASTA_SUFFIX => 'vagrent.fa';
const my $CACHE_SUFFIX_GZ => 'vagrent.cache.gz';
const my $CACHE_SUFFIX_RAW => 'vagrent.cache.raw';
const my @TRANSCRIPT_BIOTYPES => qw(protein_coding lincRNA miRNA snoRNA rRNA snRNA);
const my $ENSEMBL_SPECIES_ASSEMBLY => qr/([^\.]+?)\.(.+?)\./;
const my $ENSEMBL_VERSION_PATTERN => qr/^ftp\:\/\/ftp\.ensembl(?:genomes)?\.org\/pub\/release\-(\d+?)\//;

try {
  my $opts = option_builder();
	my $urls = getFileUrlsForRetrival($opts);

	print "Downloading Files -------- ";
	my ($codFasta, $ncFasta, $gtf) = downloadFiles($urls);
  print "Done\n";
	# make fasta file of selected cDNA sequences
	print "Building Fasta Files ----- ";
	my $fasta = generateFilteredFasta($opts, $codFasta);
	if(defined $opts->{'n'}){
		generateFilteredFasta($opts, $ncFasta, $fasta, $opts->{'n'});
	} else {
		generateFilteredFasta($opts, $ncFasta, $fasta);
	}
	print "Done\n";
	# Create fai index file
	print "Building Fasta Index ----- ";
	my $fai = createFaiIndex($fasta);
	print "Done\n";
	# Create transcript object cache file
	print "Building Vagrent Cache --- ";
	my $cache = generateCacheFile($opts,$fasta,$gtf,$codFasta);
	print "Done\n";
} catch {
  croak "An error occurred while building reference support files\:\n\t$_"; # not $@
};

sub generateCacheFile{
	my ($opts,$fasta,$gtf,$codFasta) = @_;
	my $rawCache = createFilePathFromFasta($opts,$codFasta,$CACHE_SUFFIX_RAW);
	my $cache = createFilePathFromFasta($opts,$codFasta,$CACHE_SUFFIX_GZ);
	my $cmd = $^X.' '.getCacheFileScript();
	$cmd .= " -sp ".$opts->{'sp'};
	$cmd .= " -as ".$opts->{'as'};
	$cmd .= " -d ".$opts->{'d'};
	$cmd .= " -c ".$opts->{'c'} if defined $opts->{'c'};
	$cmd .= " -f $fasta";
	$cmd .= " -g $gtf";
	$cmd .= " -o $rawCache";

	system($cmd) == 0 || croak "unable to create cache file: $rawCache";
	system("sort -k 1,1 -k 2n,2 -k 3n,3 $rawCache | bgzip > $cache") == 0 || croak "unable to sort and zip cache file: $cache";
	unlink $rawCache;
	system("tabix -p bed $cache") == 0 || croak "unable to tabix index cache file: $cache";
	return $cache;
}

sub createFaiIndex {
	my $fa = shift;
	my $cmd = "samtools faidx $fa";
	system($cmd) == 0 || croak "unable to index $fa: $?";
	return $fa.'.fai';
}

sub generateFilteredFasta {
	my ($opts, $inFa, $outFa, $statusLookupFile) = @_;
	my $app = 0;
	if(defined $outFa){
		$app = 1;
	} else {
		$outFa = createFilePathFromFasta($opts,$inFa,$FILTERED_FASTA_SUFFIX);
	}
	my $cmd = $^X.' '.getFastaFilterScript();
	$cmd .= " -f $inFa";
	$cmd .= " -o $outFa";
	$cmd .= ' -b '.join(' -b ',@TRANSCRIPT_BIOTYPES);
	$cmd .= ' -a' if $app;
	$cmd .= " -s $statusLookupFile" if defined $statusLookupFile;
	system($cmd) == 0 || croak "unable to filter $inFa: $?";
	return $outFa;
}

sub createFilePathFromFasta {
	my ($opts, $inFa, $suffix) = @_;
	my ($vol,$dirs,$file) = File::Spec->splitpath($inFa);
	my $outFa;
  if($file =~ m/$ENSEMBL_SPECIES_ASSEMBLY/){
    $outFa = File::Spec->catfile($opts->{'o'},"$1.$2.".$opts->{'e_version'}.".$suffix");
  } else {
    croak("unable to match species and assembly from file name: $inFa");
  }
	return $outFa;
}

sub downloadFiles {
	my $urls = shift;
	my @out;
	my $tmpDir = tempdir("VagrentEnsemblRefFileGenXXXXX", TMPDIR => 1, CLEANUP => 1);
	foreach my $url(@$urls){
		my $file = $tmpDir.'/'.(split /\//, $url)[-1];
		push @out, $file;
		my $rc = getstore($url, $file);
    croak "An error occured when retrieving $url\n" if(is_error($rc));
	}
	return @out;
}

sub getFileUrlsForRetrival {
	my $opts = shift;
	my @out;
	# ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/cdna
	# ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/ncrna
	# ftp://ftp.ensembl.org/pub/release-74/gtf/homo_sapiens

	# ftp://ftp.ensembl.org/pub/release-74/fasta/danio_rerio/cdna
	# ftp://ftp.ensembl.org/pub/release-74/fasta/danio_rerio/ncrna
	# ftp://ftp.ensembl.org/pub/release-74/gtf/danio_rerio

	my $ncFaDir = $opts->{'f'};
	$ncFaDir =~ s|/cdna$|/ncrna|;
	my $gtfDir = $opts->{'f'};
	$gtfDir =~ s|/cdna$||;
	$gtfDir =~ s|/fasta/|/gtf/|;
	foreach my $dir($opts->{'f'}, $ncFaDir, $gtfDir) {
		my @comps = split /\/+/, $dir;
  	my $root = shift @comps; # remove ftp:
  	my $host = shift @comps;
  	my $path = join '/', @comps;
		my $ftp = Net::FTP->new($host, Debug => 0) or croak "Failed to connect (ftp) to $host\n\t$EVAL_ERROR";
	  $ftp->login or croak "Failed to login (ftp) to $host\n\t", $ftp->message;
  	my @files = $ftp->ls($path) or croak "Failed to list directory $path on $host (ftp)\n\t", $ftp->message;
  	$ftp->quit;
		foreach my $file (@files){
			foreach my $ext (@ENSMBL_REF_FILE_EXTENTIONS){
				if($file =~ m/$ext$/){
					push @out, $dir. '/' . (split '/', $file)[-1];
				}
			}
			push @out, $dir. '/' . (split '/', $file)[-1] if($file =~ m/[[:digit:]]\.gtf\.gz$/);
		}
	}
	return \@out;
}

sub getFastaFilterScript {
	my $progPath = abs_path($0);
	my ($vol,$dirs,$file) = File::Spec->splitpath($progPath);
  return "$dirs".$FASTA_FILTER_SCRIPT;
}

sub getCacheFileScript {
	my $progPath = abs_path($0);
	my ($vol,$dirs,$file) = File::Spec->splitpath($progPath);
	my $mods = $dirs;
	$mods =~ s|bin/$|lib|;
  return "-I $mods $dirs".$GTF_CONVERSION_SCRIPT;
}

sub option_builder {
	my ($factory) = @_;
	my %opts = ();
	my $result = &GetOptions (
		'h|help' => \$opts{'h'},
		'f|ftp=s' => \$opts{'f'},
		'o|output=s' => \$opts{'o'},
		'n|ncstatus=s' => \$opts{'n'},
		'sp|species=s' => \$opts{'sp'},
		'as|assembly=s' => \$opts{'as'},
		'd|database=s' => \$opts{'d'},
		'c|ccds=s' => \$opts{'c'},
  );
  pod2usage() if($opts{'h'});
  pod2usage('Must specify the output directory to use') unless(defined $opts{'o'});
  pod2usage("Output directory must exist and be writable: $opts{o}") unless(-e $opts{'o'} && -d $opts{'o'} && -w $opts{'o'});
  pod2usage('Must specify the cDNA path for the ensembl reference data') unless defined $opts{'f'} && $opts{'f'} =~ m|cdna/?$|;
  if(defined($opts{'n'})){
  	pod2usage('Unable to read the non-coding transcript status file') unless -e $opts{'n'} && -r $opts{'n'};
  }
  pod2usage('Must specify the species') unless defined $opts{'sp'};
  pod2usage('Must specify the genome version') unless defined $opts{'as'};
  pod2usage('Must specify the ensembl core database version') unless defined $opts{'d'};
  if(defined($opts{'c'})){
  	pod2usage('CCDS file unreadable') unless -e $opts{'c'} && -r $opts{'c'};
  }
  $opts{'f'} =~ s|/$||;
  if($opts{'f'} =~ m/$ENSEMBL_VERSION_PATTERN/){
    $opts{'e_version'} = $1;
  } else {
    pod2usage('Unable to parse Ensembl version from URL');
  }

  return \%opts;
}

__END__

=head1 NAME

Admin_EnsemblReferenceFileGenerator.pl - Generates the Vagrent reference files from the specified Ensembl URL

=head1 SYNOPSIS

Admin_EnsemblReferenceFileGenerator.pl [-h] [-s human] [-v GRCh37] [-d homo_sapiens_core_74_37p] [-f <ftp://ftp.ensembl.org/pub/release-XX/fasta/XXX_XXX/cdna/>] [-c /path/to/CCDS2Sequence.version.txt] [-o /path/to/output/directory]

  General Options:

    --help         (-h)     Brief documentation

    --ftp          (-f)     Ensembl ftp directory containing the cDNA fasta sequence files

    --output       (-o)     Output directory

    --ncstatus     (-n)     Optional, path to a lookup file for the status of non-coding transcripts

    --species      (-sp)    Species (ie human, mouse)

    --assembly     (-as)    Assembly version (ie GRCh37, GRCm38)

    --database     (-d)     Ensembl core database version number (ie homo_sapiens_core_74_37p)

    --ccds         (-c)     (Optional, but strongly advised) The CCDS2Sequence file from the relevant CCDS release, see http://www.ncbi.nlm.nih.gov/CCDS

=cut
