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
use Cwd qw(getcwd abs_path);

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Getopt::Long;
use Pod::Usage;
use Try::Tiny;
use Capture::Tiny qw(capture);
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
const my $TRANSCRIPT_FILTER_SCRIPT => 'Admin_EnsemblTranscriptFilter.pl';
const my $BUILDER_SCRIPT => 'Admin_CacheFileBuilder.pl';
const my $FILTERED_FASTA_SUFFIX => 'vagrent.fa';
const my $CACHE_SUFFIX_GZ => 'vagrent.cache.gz';
const my $CACHE_SUFFIX_RAW => 'vagrent.cache.raw';
const my @TRANSCRIPT_BIOTYPES => qw(protein_coding lincRNA miRNA snoRNA rRNA snRNA);
const my $ENSEMBL_SPECIES_ASSEMBLY => qr/([^\.]+?)\.(.+?)\./;
const my $ENSEMBL_VERSION_PATTERN => qr/^ftp\:\/\/ftp\.ensembl(?:genomes)?\.org\/pub\/release\-(\d+?)\//;

my $tmpDir = tempdir("VagrentEnsemblRefFileGenXXXXX", TMPDIR => 1, CLEANUP => 1);

try {
  my $opts = option_builder();
  my ($codFasta, $ncFasta, $features,$transList);
  
  print "Downloading Files -------- ";
  if(defined $opts->{'f'}){
    my $urls = getFileUrlsForRetrival($opts);
    ($codFasta, $ncFasta, $features) = downloadFiles($tmpDir, $urls);
    print "Done\n";
  } else {
    $codFasta = $opts->{'cdna_fa'};
    $ncFasta = $opts->{'ncrna_fa'};
    $features = $opts->{'features'};
    print "Skipped, files locally supplied\n";
  }
  
  print "Obtaining Filtered Transcript List ----- ";
  unless(defined $opts->{'trans_list'}){
    $transList = generateTranscriptListFile($tmpDir,$features, $codFasta, $ncFasta);
    print "Done\n";
  } else {
    $transList = $opts->{'trans_list'};
    print "Skipped, files locally supplied\n";
  }
  
  print "Building Cache Files ----- ";
  buildCacheFiles($tmpDir, $opts, $features, $transList, $codFasta, $ncFasta);
  print "Done\n";
    
} catch {
  croak "An error occurred while building reference support files\:\n\t$_"; # not $@
};

sub buildCacheFiles {
  my ($tmpDir, $opts, $features, $transList, @seq_files) = @_;
  my $cmd = $^X.' '.getBuilderScript();
  foreach my $s(@seq_files){
    $cmd .= " -f $s";
  }
  $cmd .= ' -o '.$opts->{'o'};
  $cmd .= " -gf $features";
  $cmd .= " -t $transList";
  $cmd .= ' -sp '.$opts->{'sp'};
  $cmd .= ' -as '.$opts->{'as'};
  $cmd .= ' -d '.$opts->{'d'};
  $cmd .= ' -fai '.$opts->{'fai'} if defined $opts->{'fai'};
  my ($stdout,$stderr,$exit) = capture {system($cmd)};
  croak('Unable to generate cache files list:- '.$stderr) if $exit;
  return;
}

sub getBuilderScript {
  my $progPath = abs_path($0);
	my ($vol,$dirs,$file) = File::Spec->splitpath($progPath);
  return "$dirs".$BUILDER_SCRIPT;
}

sub generateTranscriptListFile {
  my ($tmpDir,$features, @seq_files) = @_;
  my $cmd = $^X.' '.getFilterScript();
  my $transListFile = makeTmpFilePath($tmpDir,'.translist');
  foreach my $s(@seq_files){
    $cmd .= " -f $s";
  }
  $cmd .= " -o $transListFile";
  $cmd .= ' -b '.join(' -b ',@TRANSCRIPT_BIOTYPES);
  my ($stdout,$stderr,$exit) = capture {system($cmd)};
  croak('Unable to generate transcript list:- '.$stderr) if $exit;
  return $transListFile;
}

sub getFilterScript {
  my $progPath = abs_path($0);
	my ($vol,$dirs,$file) = File::Spec->splitpath($progPath);
  return "$dirs".$TRANSCRIPT_FILTER_SCRIPT;
}

sub getFileUrlsForRetrival {
	my $opts = shift;
	my @out;
	# ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/cdna
	# ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/ncrna
	# ftp://ftp.ensembl.org/pub/release-90/gff3/homo_sapiens
	
	# ftp://ftp.ensembl.org/pub/release-70/fasta/homo_sapiens/cdna
	# ftp://ftp.ensembl.org/pub/release-70/fasta/homo_sapiens/ncrna
	# ftp://ftp.ensembl.org/pub/release-70/gtf/homo_sapiens

	# ftp://ftp.ensembl.org/pub/release-74/fasta/danio_rerio/cdna
	# ftp://ftp.ensembl.org/pub/release-74/fasta/danio_rerio/ncrna
	# ftp://ftp.ensembl.org/pub/release-74/gtf/danio_rerio

  $opts->{'f'} =~ s|/$||; # clean training '/' off url
	my $ncFaDir = $opts->{'f'};
	$ncFaDir =~ s|/cdna$|/ncrna|;
	my $gffDir = $opts->{'f'};
	$gffDir =~ s|/cdna$||;
	$gffDir =~ s|/fasta/|/gff3/|;
	my $gtfDir = $opts->{'f'};
	$gtfDir =~ s|/cdna$||;
	$gtfDir =~ s|/fasta/|/gtf/|;
	my $ftp = undef; 
	foreach my $dir($opts->{'f'}, $ncFaDir, $gffDir, $gtfDir) {
		my @comps = split /\/+/, $dir;
  	my $root = shift @comps; # remove ftp:
  	my $host = shift @comps;
  	my $path = join '/', @comps;
  	my $have_feature_file = 0;
  	unless(defined $ftp){
  	  $ftp = Net::FTP->new($host, Debug => 0) or croak "Failed to connect (ftp) to $host\n\t$EVAL_ERROR";
	    $ftp->login or croak "Failed to login (ftp) to $host\n\t", $ftp->message;
  	}
  	my @files;
  	unless(@files = $ftp->ls($path)){
  	  unless($dir eq $gffDir){ # older ensembl releases don't have a gff directory, can be missing.
  	    croak "Failed to list directory $path on $host (ftp)\n\t", $ftp->message;
  	  }
  	}
		foreach my $file (@files){
			foreach my $ext (@ENSMBL_REF_FILE_EXTENTIONS){
				if($file =~ m/$ext$/){
					push @out, $dir. '/' . (split '/', $file)[-1];
				}
			}
			if ($dir eq $gffDir && $file =~ m/\.[[:digit:]]+?\.gff3\.gz$/ && $file !~ m/\.[[:digit:]]+?\.(?:chromosome)\..+?\.gff3\.gz$/){
			  push @out, $dir. '/' . (split '/', $file)[-1];
			  $have_feature_file = 1;
			} elsif ($dir eq $gtfDir && $file =~ m/[[:digit:]]\.gtf\.gz$/){
			  push @out, $dir. '/' . (split '/', $file)[-1]
			}
		}
		last if $have_feature_file;
	}
	$ftp->quit;
	return \@out;
}

sub makeTmpFilePath {
  my ($tmpDir,$ext) = @_;
  my $file;
  (undef, $file) = tempfile('VagrentEnsemblRefFileGenXXXXX', DIR => $tmpDir, OPEN => 0, SUFFIX => $ext);
  return $file;
}

sub downloadFiles {
	my ($tmpDir, $urls) = @_;
	my @out;
	foreach my $url(@$urls){
		my $file = $tmpDir.'/'.(split /\//, $url)[-1];
		push @out, $file;
		my $rc = getstore($url, $file);
    croak "An error occured when retrieving $url\n" if(is_error($rc));
	}
	return @out;
}

sub option_builder {
	my ($factory) = @_;
	my %opts = ();
	my $result = &GetOptions (
		'h|help' => \$opts{'h'},
		'f|ftp=s' => \$opts{'f'},
		'gf|features=s' => \$opts{'features'},
		'cf|cdna_fa=s' => \$opts{'cdna_fa'},
		'nf|ncrna_fa=s' => \$opts{'ncrna_fa'},
		'o|output=s' => \$opts{'o'},
		'tl|trans_list=s' => \$opts{'trans_list'},
		'sp|species=s' => \$opts{'sp'},
		'as|assembly=s' => \$opts{'as'},
		'd|database=s' => \$opts{'d'},
		'c|ccds=s' => \$opts{'c'},
		'fai|fai=s' => \$opts{'fai'},
  );
  pod2usage() if($opts{'h'});
  pod2usage('Must specify the output directory to use') unless(defined $opts{'o'});
  pod2usage("Output directory must exist and be writable: $opts{o}") unless(-e $opts{'o'} && -d $opts{'o'} && -w $opts{'o'});
  
  if(defined $opts{'f'}){
    if(defined $opts{'features'} || defined $opts{'cdna_fa'} || defined $opts{'ncrna_fa'}){
      pod2usage('Please only define the remote URL or the local files, not both');
    } else {
      pod2usage('Must specify the cDNA URL for the ensembl reference data') unless($opts{'f'} =~ m|cdna/?$|);
    }
  } elsif (defined $opts{'cdna_fa'} && defined $opts{'ncrna_fa'} && defined $opts{'features'}){
    if(defined $opts{'f'}){
      pod2usage('Please only define the remote URL or the local files, not both');
    } else {
      pod2usage('Please specify a valid genomic feature file (-gff)') unless -e $opts{'features'} && -r $opts{'features'};
      pod2usage('Please specify a valid coding cDNA sequence file (-cf)') unless -e $opts{'cdna_fa'} && -r $opts{'cdna_fa'};
      pod2usage('Please specify a valid non-coding cDNA sequence file (-nf)') unless -e $opts{'ncrna_fa'} && -r $opts{'ncrna_fa'};
    }
  
  }
  if(defined($opts{'fai'})){
  	pod2usage('Unable to read the fai file: '.$opts{'fai'}) unless -e $opts{'fai'} && -r $opts{'fai'};
  }
  pod2usage('Must specify the species') unless defined $opts{'sp'};
  pod2usage('Must specify the genome version') unless defined $opts{'as'};
  pod2usage('Must specify the ensembl core database version') unless defined $opts{'d'};
  if(defined($opts{'c'})){
  	pod2usage('CCDS file unreadable') unless -e $opts{'c'} && -r $opts{'c'};
  }
  if(defined($opts{'trans_list'})){
  	pod2usage('Transcript file unreadable') unless -e $opts{'trans_list'} && -r $opts{'trans_list'};
  }
  return \%opts;
}

__END__

=head1 NAME

Admin_EnsemblReferenceFileGenerator.pl - Generates the Vagrent reference files from the specified Ensembl URL

=head1 SYNOPSIS

Admin_EnsemblReferenceFileGenerator.pl [-h] [-sp Human] [-as GRCh37] [-d homo_sapiens_core_74_37p] [-f <ftp://ftp.ensembl.org/pub/release-XX/fasta/XXX_XXX/cdna/>] [-o /path/to/output/directory]

  Required Options:

    --output       (-o)     Output directory    

    --species      (-sp)    Species (ie human, mouse)

    --assembly     (-as)    Assembly version (ie GRCh37, GRCm38)

    --database     (-d)     Ensembl core database version number (ie homo_sapiens_core_74_37p)

  Dynamic Download:
  
    --ftp          (-f)     Ensembl ftp directory containing the cDNA fasta sequence files

  Or Local Files:
  
    --features     (-gf)    gff3 or gtf file containing transcript and gene information
  
    --cdna_fa      (-cf)    Fasta file containing protein coding cdna sequences
  
    --ncrna_fa     (-nf)    Fasta file containing non-coding cdna sequences 

  Optional:
  
    --help         (-h)     Brief documentation
        
    --ccds         (-c)     (Recommended) The CCDS2Sequence file from the relevant CCDS release, see http://www.ncbi.nlm.nih.gov/CCDS
    
    --fai          (-fai)   (Recommended) The samtools fasts index file (.fai) for your reference genome
                              This is the reference genome that your bam and vcf files will be mapped to
    
    --trans_list   (-tl)    List of preprepared transcript accessions, only these accesions will be included in the reference output
    
=cut
