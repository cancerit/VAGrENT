package Sanger::CGP::Vagrent::Data::AnnotationGroup;

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
use Data::Dumper;
use List::Util qw(first);
use Carp;

use Log::Log4perl;
use Sanger::CGP::Vagrent qw($VERSION);
use Sanger::CGP::Vagrent::Data::Transcript;

use base qw(Sanger::CGP::Vagrent);

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

sub _init {
	my $self = shift;
	my %vars = @_;
	foreach my $k(keys(%vars)){
		if($k eq 'label'){
			$self->{_label} = $vars{label};
		} elsif($k eq 'accession'){
			$self->{_accession} = $vars{accession};
		} elsif($k eq 'type'){
			my $good = undef;
			foreach my $type(Sanger::CGP::Vagrent::Data::Transcript::ALL_GENE_TYPES){
				if($vars{type} eq $type){
					$good = $vars{type};
					last;
				}
			}
			if(defined($good)){
				$self->{_type} = $good;
			} else {
				croak('recieved unknown gene type: '.$vars{type});
			}
		} elsif($k eq 'ccds'){
			$self->{_ccds} = $vars{ccds};
		}
	}
}

sub addAnnotation {
	my ($self,$a) = @_;
	unless(defined($a)){
		$log->error('cannot add an undef Annotation object to an annotation group');
		return;
	}
	unless($a->isa('Sanger::CGP::Vagrent::Data::Annotation')){
		$log->error('can only add Annotation objects to an annotation group');
		return;
	}
	if(exists($self->{_anno}) && defined($self->{_anno}) && scalar(@{$self->{_anno}}) > 0){
		my $good = 1;
		foreach my $ca(@{$self->{_anno}}){
			if($ca->getContext eq $a->getContext){
				$good = 0;
				$log->error('cant add more that one annotation of the same context to the same group');
				last;
			}
		}
		if($good){
			push(@{$self->{_anno}},$a);
		}
	} else {
		push(@{$self->{_anno}},$a);
	}
}

sub getAllAnnotations {
	my $self = shift;
	if(defined($self->{_anno})){
		return $self->{_anno};
	} else {
		return undef;
	}
}

sub getAnnotationByContext {
	my ($self,$ctx) = @_;
	foreach my $a(@{$self->{_anno}}){
		return $a if($a->getContext eq $ctx);
	}
	return undef;
}

sub getLabel {
	return shift->{_label};
}

sub getCCDS {
	return shift->{_ccds};
}

sub getAccession {
	return shift->{_accession};
}

sub getType {
	return shift->{_type};
}

sub addClassification {
	my ($self,@classes) = @_;
	foreach my $c (@classes){
		if(defined($c)){
			if(defined($self->getClassifications)){
				push(@{$self->{_class}},$c) unless(first {$_ eq $c} $self->getClassifications);
			} else {
				push(@{$self->{_class}},$c);
			}
		}
	}
}

sub getClassifications {
	my ($self,$class) = @_;
	return @{$self->{_class}} if(exists $self->{_class} && defined $self->{_class});
	return undef;
}

sub hasClassification {
	my ($self,$class) = @_;
	if(exists $self->{_class} && defined($self->{_class})){
		if(first {$_ eq $class} @{$self->{_class}}){
			return 1;
		} else {
			return 0;
		}
	}
	return 0;
}

sub addBookmark {
	my ($self,$marker) = @_;
	unless($marker->isa('Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker')){
		my $msg = 'can only use Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker objects to mark an annotation group';
		$log->error($msg);
		warn $msg;
		return;
	}
	$self->{_bookmarks}->{ref($marker)} = 1;
}

sub hasBookmark {
	my ($self,$marker) = @_;
	if(defined $marker){
		unless($marker->isa('Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker')){
			my $msg = 'can only use Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker objects to look for bookmarks';
			$log->error($msg);
			warn $msg;
			return;
		}
		if(exists $self->{_bookmarks} && exists $self->{_bookmarks}->{ref($marker)} && $self->{_bookmarks}->{ref($marker)} == 1){
			return 1;
		}
	} else {
		if(exists $self->{_bookmarks} && defined $self->{_bookmarks}){
			return 1;
		}
	}
	return 0;
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::Data::AnnotationGroup - Object holding data about the effect of a single
variation on a transcript.

=head1 DESCRIPTION

This holds information about a variation within a transcript along with a collection of
L<Sanger::CGP::Vagrent::Data::Annotation|Sanger::CGP::Vagrent::Data::Annotation> objects detailing the descriptions of that variation in
different sequence contexts.

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $annoGrp = Sanger::CGP::Vagrent::Data::AnnotationGroup->new(%params);

=item Function :

Builds a new Sanger::CGP::Vagrent::Data::AnnotationGroup object

=item Returns :

L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> object initialized with parameter values

=item Params :

 accession => Accession for the transcript
 type      => Gene type (defined by type constants in L<Sanger::CGP::Vagrent::Data::Transcript|Sanger::CGP::Vagrent::Data::Transcript>)
 label     => text label for the transcript, typically the gene name
 ccds      => CCDS reference number if the transcript has one

=back

=head2 Attributes

=head3 addAnnotation

=over

=item Usage :

 $annoGrp->addAnnotation($anno);

=item Function :

Validates and adds a L<Sanger::CGP::Vagrent::Data::Annotation|Sanger::CGP::Vagrent::Data::Annotation> object to the group.
Only one annotation of each context can be added to a group

=item Params :

A L<Sanger::CGP::Vagrent::Data::Annotation|Sanger::CGP::Vagrent::Data::Annotation> object

=item Returns :

None

=back

=head3 getAllAnnotations

=over

=item Usage :

 my $annoListRef = $annoGrp->getAllAnnotations;

=item Function :

Retrieves all annotations in the group

=item Returns :

Array ref of L<Sanger::CGP::Vagrent::Data::Annotation|Sanger::CGP::Vagrent::Data::Annotation> objects, or undef

=back

=head3 getAnnotationByContext

=over

=item Usage :

 my $cdsAnno = $annoGrp->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext);

=item Function :

Retrieves annotation for the specified context

=item Params :

String - A context constant from L<Sanger::CGP::Vagrent::Data::Annotation|Sanger::CGP::Vagrent::Data::Annotation>

=item Returns :

A L<Sanger::CGP::Vagrent::Data::Annotation|Sanger::CGP::Vagrent::Data::Annotation> object, or undef

=back

=head3 getLabel

=over

=item Usage :

 my $label = $grp->getLabel;

=item Function :

The label for the transcript (Gene name or similar)

=item Returns :

String or undef

=back

=head3 getCCDS

=over

=item Usage :

 my $ccds = $grp->getCCDS;

=item Function :

The CCDS reference number for the transcript

=item Returns :

String or undef

=back

=head3 getAccession

=over

=item Usage :

 my $acc = $grp->getAccession;

=item Function :

The accession number for the transcript

=item Returns :

String or undef

=back

=head3 getType

=over

=item Usage :

 my $type = $grp->getType;

=item Function :

Gene type (defined by type constants in L<Sanger::CGP::Vagrent::Data::Transcript|Sanger::CGP::Vagrent::Data::Transcript>)

=item Returns :

String - A type constant from L<Sanger::CGP::Vagrent::Data::Transcript|Sanger::CGP::Vagrent::Data::Transcript>

=back

=head3 addClassification

=over

=item Usage :

 $grp->addClassification($class1,$class2);

=item Function

Appends one or more classification terms to the annotation group

=item Params

An array of classification terms (Strings)

=item Returns

None

=back

=head3 getClassifications

=over

=item Usage :

 my @classes = $grp->getClassifications;

=item Function

Returns the stored classification terms

=item Returns

Array of Strings, or undef if none

=back

=head3 addBookmark

=over

=item Usage :

 $g->addBookmark($bookmarker);

=item Function

Adds a bookmark to the annotation group of a type defined by the supplied bookmarker

=item Params

A L<Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker|Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker> implementing object

=item Returns

None

=back

=head2 Functions

=head3 hasClassification

=over

=item Usage :

 if($grp->hasClassification($checkClass)){
  .....
 }

=item Function

Checks an annotation to see if it already contains the specified class

=item Returns

Boolean, 1 or 0

=back

=head3 hasBookmark

=over

=item Usage :

 $g->hasBookmark;
 $g->hasBookmark($bookmarker);

=item Function

Returns true if the Annotation group has a bookmark, if a Bookmarker object is
supplied it only returns true if that specific bookmark type is present

=item Params

(Optional) A L<Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker|Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker> implementing object

=item Returns

Boolean, 1 or 0

=back
