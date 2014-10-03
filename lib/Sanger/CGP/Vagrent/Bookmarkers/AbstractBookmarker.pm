package Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker;

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

use Log::Log4perl;
use Attribute::Abstract;
use Data::Dumper;

use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Sanger::CGP::Vagrent::Ontology::SequenceOntologyClassifier);

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

sub markAnnotation {
	my $self = shift;
	my $mark = $self->getAnnotation(@_);
	$mark->addBookmark($self) if defined $mark;
}

sub getAnnotation {
	my $self = shift;
	my @groups = @_;
	return undef unless scalar @groups > 0 && defined $groups[0];
	foreach my $g(@groups){
		unless(defined($g) && $g->isa('Sanger::CGP::Vagrent::Data::AnnotationGroup')){
			my $msg = 'require a Sanger::CGP::Vagrent::Data::AnnotationGroup object not a '.ref($g);
			$self->addMessage($msg);
			$log->info($msg);
			return undef;
		}
	}
	return $self->_getAnnotation(@groups);
}

sub _getAnnotation: Abstract;

sub _init {
	my $self = shift;
	my %vars = @_;
	# place holder, may need to be overridden in subclass
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker - Abstract base class for the AnnotationGroup bookmarkers

=head1 DESCRIPTION

This is an abstract base class for the AnnotationGroup bookmarkers, it provides some simple shared functionality.  All
subclasses must implement the _getAnnotation object method, and optionally the _init object method.

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $marker = Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarkerSubClass->new();

=item Function :

Builds a new Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker inheriting object

=item Returns :

Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker object initialized with parameter values

=item Params :

Hash of parameter values, the actual keys/values required depend on the sub class.

=back

=head2 Functions

=head3 getAnnotation

=over

=item Usage :

 my $markedGroup = $bookmarker->getAnnotation(@annotationGroups);

=item Function :

Returns the AnnotationGroup that matches the bookmark

=item Returns :

A L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> object

=item Params :

An array of L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> objects

=back

=head3 markAnnotation

=over

=item Usage :

 $bookmarker->markAnnotation(@annotationGroups);

=item Function :

Calls getAnnotation, but marks annotation group rather than returning it

=item Returns :

None

=item Params :

An array of L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> objects

=back

=head2 Abstract

=head3 _getAnnotation

=over

=item Usage :

 my $markedGroup = $self->_getAnnotation(@annotationGroups);

=item Function :

Abstract internal function, must be overridden by a subclass.  This method should contain the actual bookmarking logic

=item Returns :

A L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> object

=item Params :

An array of L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> objects

=back

=head3 _init

=over

=item Usage :

 $self->_init(@params);

=item Function :

Abstract internal function, optionally can be overridden by a subclass.  This method is used to handle constructor parameters, its called by new.

=item Returns :

None

=item Params :

The array of constructor paramaters to be parsed.

=back
