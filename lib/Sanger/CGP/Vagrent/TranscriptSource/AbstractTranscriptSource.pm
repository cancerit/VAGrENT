package Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource;

use strict;

use Log::Log4perl;
use Carp;
use Attribute::Abstract;
use Data::Dumper;

my $log = Log::Log4perl->get_logger(__PACKAGE__);

1;

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
	my $self = {};
	bless($self, $class);
	$self->_init(@_);
	return $self;
}

sub _init: Abstract;

sub getTranscripts: Abstract;

sub setDumpRegion {
	my ($self,$gr) = @_;
	unless(defined($gr) && $gr->isa('Sanger::CGP::Vagrent::Data::AbstractGenomicPosition')){
		croak("Did not recieve a Sanger::CGP::Vagrent::Data::AbstractGenomicPosition object");
	}
	$self->{_dumpInfo} = undef;
	$self->{_dumpInfo}->{_region} = $gr;
	return;
}

sub getDumpRegion {
	my $self = shift;
	if(defined $self->{_dumpInfo} && defined $self->{_dumpInfo}->{_region}){
		return $self->{_dumpInfo}->{_region};
	}
	return undef;
}

sub isDumpRegionACompleteSequence {
	my $self = shift;
	if(defined $self->{_dumpInfo} && defined $self->{_dumpInfo}->{_fullSeq}){
		return $self->{_dumpInfo}->{_fullSeq};
	}
	return undef;
}

sub getTranscriptsForNextGeneInDumpRegion: Abstract;

__END__

=head1 NAME

Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource - Abstract base class for Transcript sources

=head1 DESCRIPTION

This is a base utility class all TranscriptSources should inherit from, its little more than an interface with the vast majority of the functionality being implemented in the sub class.

Sub classes must implement a getTranscripts method and an internal _init method.

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $source = Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSourceSubClass->new();

=item Function :

Builds a new Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource inheriting object

=item Returns :

Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource object initialized with parameter values

=item Params :

Hash of parameter values, the actual keys/values required depend on the sub class.

=back

=head2 Attributes

=head3 setDumpRegion

=over

=item Usage :

 $source->setDumpRegion($genomicPosition);

=item Function :

Sets the L<GenomicRegion|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition> to seed Transcript/Gene set retrieval

=item Returns :

Nothing

=item Params :

Any L<Sanger::CGP::Vagrent::Data::AbstractGenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition> implementing object

=back

=head3 getDumpRegion

=over

=item Usage :

 my $region = $source->getDumpRegion();

=item Function :

Gets the L<GenomicRegion|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition> used to seed Transcript/Gene set retrieval

=item Returns :

A L<Sanger::CGP::Vagrent::Data::AbstractGenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition> implementing object or undef

=back

=head2 Functions

=head3 isDumpRegionACompleteSequence

=over

=item Usage :

 if($source->isDumpRegionACompleteSequence()){
 	......
 }

=item Function :

Returns boolean status for the current dump region being a complete sequence (eg. a whole chromosome/contig).  The value behind this is populated by getTranscriptsForNextGeneInDumpRegion method when processing the specified gene region

=item Returns :

A Boolean, 1 if it is a complete sequence, 0 if is isn't or undef if it doesn't know

=back

=head2 Abstract

=head3 getTranscripts

=over

=item Usage :

 my @transList = $source->getTranscripts($genomicPosition);

=item Function :

Abstract function, must be overridden by a subclass.  Retrieves a list of L<Transcripts|Sanger::CGP::Vagrent::Data::Transcript> objects overlapping the specified L<GenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition>

=item Returns :

An array of L<Sanger::CGP::Vagrent::Data::Transcript|Sanger::CGP::Vagrent::Data::Transcript> objects

=item Params :

Any L<Sanger::CGP::Vagrent::Data::AbstractGenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition> implementing object

=back

=head3 getTranscriptsForNextGeneInDumpRegion

=over

=item Usage :

 while (my @transList = $source->getTranscriptsForNextGeneInDumpRegion()){
 	....
 }

=item Function :

Abstract function, must be overridden by a subclass.  Retrieves a list of L<Transcripts|Sanger::CGP::Vagrent::Data::Transcript> objects belonging to the next gene inside the previously specified L<GenomicRegion|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition>

=item Returns :

An array of L<Sanger::CGP::Vagrent::Data::Transcript|Sanger::CGP::Vagrent::Data::Transcript> objects

=item Params :

None

=back

=head3 _init

=over

=item Usage :

 $self->_init(@params);

=item Function :

Abstract internal function, must be overridden by a subclass.  This method is used to handle constructor parameters, its called by new.

=item Returns :

None

=item Params :

The flattened array of key/value pairs passed in from new.

=back
