package Sanger::CGP::Vagrent::Annotators::AnnotatorCollection;

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
use Carp;
use Data::Dumper;

use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Sanger::CGP::Vagrent::Annotators::AbstractAnnotator);

1;

sub _getAnnotation {
	my ($self,$var) = @_;
	return undef unless( defined($var) && $var->isa('Sanger::CGP::Vagrent::Data::AbstractVariation'));
	my @ans;
	foreach my $ma(@{$self->{_collection}}){
		my @tmp = $ma->getAnnotation($var);
		foreach my $t(@tmp){
			next unless(defined($t));
			push(@ans,$t);
		}
	}
	return @ans;
}

sub addAnnotator {
	my ($self,$ma) = @_;
	if($ma->isa('Sanger::CGP::Vagrent::Annotators::AbstractAnnotator')){
		push(@{$self->{_collection}},$ma);
	} else {
		$self->logger->fatal('can only deligate to other AbstractAnnotator objects');
		croak('can only deligate to other AbstractAnnotator objects');
	}
}

sub getMessages {
	my $self = shift;
	my @msgs;
	foreach my $ma(@{$self->{_collection}}){
		push(@msgs,$ma->getMessages()) if defined($ma->getMessages());
	}
	return @msgs;
}

__END__


=head1 NAME

Sanger::CGP::Vagrent::Annotators::AnnotatorCollection - Convienience class for running annotation over a variety of input data types

=head1 DESCRIPTION

This is a aggregator class, it holds several annotators and runs them one at a time over any submitted
variation.  It does none of the annotation work itself, it simply combines the answers of the other
annotators it knows about.  It extends L<Sanger::CGP::Vagrent::Annotators::AbstractAnnotator|Sanger::CGP::Vagrent::Annotators::AbstractAnnotator>
to keep the same interface.

It uses the same constructor as the L<AbstractAnnotator|Sanger::CGP::Vagrent::Annotators::AbstractAnnotator>, however it doesn't need to be given a L<TranscriptSource|Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource> as it won't be doing any actual annotating.
However it should be used for Bookmarking. By adding L<Bookmarkers|Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker> to the AnnotatorCollection, results from each contained L<Annotator|Sanger::CGP::Vagrent::Annotators::AbstractAnnotator> will be collated together and then bookmarked as a single list.
Bookmarking separately in each child L<Annotator|Sanger::CGP::Vagrent::Annotators::AbstractAnnotator> will result in duplicate answers for the same bookmark type.


=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $source = Sanger::CGP::Vagrent::Annotators::AnnotatorCollection->new(%params);

=item Function :

Builds a new Sanger::CGP::Vagrent::Annotators::AnnotatorCollection

=item Returns :

Sanger::CGP::Vagrent::Annotators::AnnotatorCollection object, optionally initialised with L<Bookmarkers|Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker>

=item Params :

Hash of parameter values

 bookmarker       => (Optional) An array reference of, or single, Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker inheriting object
 only_bookmarked  => (Optional) Boolean, only return annotations that get bookmarked

=back

=head2 Functions

=head3 addAnnotator

=over

=item Usage :

 $AnnotatorCollection->addAnnotator($substitutionAnnotator);
 $AnnotatorCollection->addAnnotator($insertionAnnotator);
 $AnnotatorCollection->addAnnotator($deletionAnnotator);

=item Function :

Adds a new child L<Annotator|Sanger::CGP::Vagrent::Annotators::AbstractAnnotator> to the collection.  Every variant given to the collection will be sent to every child Annotator in turn.

=item Returns :

Nothing

=item Params :

A L<Sanger::CGP::Vagrent::Annotators::AbstractAnnotator|Sanger::CGP::Vagrent::Annotators::AbstractAnnotator> inheriting object

=back

=head3 getMessages

=over

=item Usage :

 my @mess = $annotator->getMessages();

=item Function :

Retrieves the messages from all child Annotators and returns them as a single message list

=item Returns :

Array of Strings

=back
