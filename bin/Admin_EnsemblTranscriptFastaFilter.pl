#!/usr/bin/perl

use strict;
use English qw(-no_match_vars);
use warnings FATAL => 'all';
use Carp;

use Getopt::Long;
use Pod::Usage;
use Try::Tiny;

use File::Type;
use Readonly qw(Readonly);

use Data::Dumper;
use Bio::SeqIO;

Readonly my @TEXT_TYPES => qw(application/octet-stream);

Readonly my @ZIP_TYPES => qw(application/x-gzip);

try {
	my $opts = option_builder();
	my $trans = getInterestingTranscriptsFromFile($opts);
	generateFilteredFasta($opts,$trans);
} catch {
  warn "An error occurred while building reference support files\:\n\t$_"; # not $@
};

sub generateFilteredFasta {
	my ($opts,$trans) = @_;
	my $type = getInputFileType($opts->{'f'});
	my ($inFh,$outFh);

	if($type eq 'text') {
		open $inFh, '<'.$opts->{'f'} || croak("unable to open sequence input file:".$opts->{'f'});
	} elsif($type eq 'zip'){
		open $inFh, 'zcat '.$opts->{'f'}.' |' || croak("unable to open sequence input file:".$opts->{'f'});
	}

	if(defined $opts->{'a'}){
		open $outFh, '>>'.$opts->{'o'} || croak("unable to open sequence ouput file:".$opts->{'f'});
	} else {
		open $outFh, '>'.$opts->{'o'} || croak("unable to open sequence ouput file:".$opts->{'f'});
	}

	my $keep = 0;

	while (<$inFh>){
		if(m/^>/){
			my ($key) = split /\s/;
			$key =~ s/^>//;
			if(exists $trans->{$key}){
				$keep = 1;
				print $outFh ">$key\n";
			} else {
				$keep = 0;
			}
		} elsif($keep){
			my $line = $_;
			$line =~ s/\s//g;
			print $outFh "$line\n";
		}
	}
	close $outFh;
	close $inFh;
}

sub getInterestingTranscriptsFromFile {
	my $opts = shift;
	my $type = getInputFileType($opts->{'f'});
	my $out;
	my $status = undef;
	my $cmd;
	if($type eq 'text'){
		$cmd = "grep";
	} elsif($type eq 'zip'){
		$cmd = "zgrep";
	} else {
		croak "file is of an unrecognised format ($type): ".$opts->{'f'};
	}

	$cmd .= " '>' ". $opts->{'f'};
	$cmd .= ' | grep ';
	foreach my $b(@{$opts->{'b'}}){
		$cmd .= " -e transcript_biotype:$b";
	}

	if(defined $opts->{'s'}){
		open my $statfh, '<'.$opts->{'s'} || croak("unable to open status lookup file:".$opts->{'s'});
		while(<$statfh>){
			my @st = split /\s+/;
			if(uc($st[1]) eq 'KNOWN'){
				$status->{$st[0]}++;
			}
		}
	} else {
		$cmd .= ' | grep known ';
	}

	open my $fh, "$cmd |" || croak("unable to grep through fasta file with command: $cmd");
	if(defined $opts->{'s'}) {
		while (<$fh>){
			my ($key) = split /\s/;
			$key =~ s/^>//;
			$out->{$key}++ if(exists $status->{$key});
		}
	} else {
		while (<$fh>){
			my ($key) = split /\s/;
			$key =~ s/^>//;
			$out->{$key}++;
		}
	}
	close $fh;
	return $out;
}

sub getInputFileType {
	my $infile = shift;
	my $fileType = File::Type->new->checktype_filename($infile);
	foreach my $t (@TEXT_TYPES){
		return 'text' if($t eq $fileType);
	}
	foreach my $t (@ZIP_TYPES) {
		return 'zip' if($t eq $fileType);
	}
	return $fileType;
}

sub option_builder {
	my ($factory) = @_;

	my %opts = ();

	my $result = &GetOptions (
		'h|help' => \$opts{'h'},
		'f|fasta=s' => \$opts{'f'},
		'o|outfile=s' => \$opts{'o'},
		'a|append' => \$opts{'a'},
		'b|biotypes=s@' => \$opts{'b'},
		's|status=s' => \$opts{'s'},
  );

  pod2usage() if($opts{'h'});
  pod2usage('Must specify the input fasta ensembl reference file') unless defined $opts{'f'} && -e $opts{'f'};
  pod2usage('Must specify the output file to use') unless defined $opts{'o'};
  pod2usage('Must specify at least one biotype') unless scalar @{$opts{'b'}} > 0 && defined($opts{'b'}->[0]);
  if(defined($opts{'s'})){
  	pod2usage('Unable to read the transcript status file') unless -e $opts{'s'} && -r $opts{'s'};
  }
  return \%opts;
}

__END__

=head1 NAME

Admin_EnsemblTranscriptFastaFilter.pl - Generates the Vagrent reference fasta file from an Ensembl one

=head1 SYNOPSIS

Admin_EnsemblTranscriptFastaFilter.pl [-h] [-f /path/to/ensembl.fa] [-o /path/to/output.fa] [-b biotype] [-s /path/to/status/lookup.list][-a]

  General Options:

    --help         (-h)     Brief documentation

    --fasta        (-f)     Ensembl fasta file

    --output       (-o)     Output file

    --biotypes     (-b)     Ensembl transcript biotypes to filter for, supports multiple instances

    --status       (-s)     Optional, Transcript Status look up file, simple transcript id - status list, space separated, one entry per line

    --append       (-a)     Append to existing output file

=cut
