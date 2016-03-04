package Sanger::CGP::Vagrent::Data::Annotation;

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
use Const::Fast qw(const);

use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Sanger::CGP::Vagrent);

# sequence context constant values

const my $CDS_ANNOTATION_CONTEXT => 'CDS';
const my $MRNA_ANNOTATION_CONTEXT => 'mRNA';
const my $PROTEIN_ANNOTATION_CONTEXT => 'Protein';

const my @ALL_ANNOTATION_CONTEXT => ($CDS_ANNOTATION_CONTEXT,$MRNA_ANNOTATION_CONTEXT,$PROTEIN_ANNOTATION_CONTEXT);

# variation type constant values

const my $SUBSTITUTION_ANNOTATION_TYPE => 'Substitution';
const my $DELETION_ANNOTATION_TYPE => 'Deletion';
const my $INSERTION_ANNOTATION_TYPE => 'Insertion';
const my $COMPLEX_ANNOTATION_TYPE => 'Complex';
const my $FRAMESHIFT_ANNOTATION_TYPE => 'FrameShift';
const my $UNKNOWN_ANNOTATION_TYPE => 'Unknown';

const my @ALL_ANNOTATION_TYPES => ($COMPLEX_ANNOTATION_TYPE,$FRAMESHIFT_ANNOTATION_TYPE,$INSERTION_ANNOTATION_TYPE,$DELETION_ANNOTATION_TYPE,$SUBSTITUTION_ANNOTATION_TYPE,$UNKNOWN_ANNOTATION_TYPE);

# variation subtype constant values for location data description

const my $POSITION_KNOWN_SUBTYPE => 'POS-KNOWN';
const my $POSITION_OFFSET_SUBTYPE => 'POS-OFFSET';
const my $POSITION_OFF_SEQUENCE_SUBTYPE => 'POS-OFFSEQ';

const my @ALL_ANNOTATION_SUBTYPES => ($POSITION_OFF_SEQUENCE_SUBTYPE,$POSITION_OFFSET_SUBTYPE,$POSITION_KNOWN_SUBTYPE);

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
		if($k eq 'minpos'){
			$self->{_minpos} = $vars{minpos};
		} elsif($k eq 'minOffset'){
			$self->{_minOffset} = $vars{minOffset};
		} elsif($k eq 'maxpos'){
			$self->{_maxpos} = $vars{maxpos};
		} elsif($k eq 'maxOffset'){
			$self->{_maxOffset} = $vars{maxOffset};
		} elsif($k eq 'wt'){
			$self->{_wt} = $vars{wt};
		} elsif($k eq 'mt'){
			$self->{_mt} = $vars{mt};
		} elsif($k eq 'acc'){
			$self->{_acc} = $vars{acc};
		} elsif($k eq 'accversion'){
			$self->{_accversion} = $vars{accversion};
		} elsif($k eq 'seqlength'){
			$self->{_seqlength} = $vars{seqlength};
		} elsif($k eq 'db'){
			$self->{_db} = $vars{db};
		} elsif($k eq 'dbversion'){
			$self->{_dbversion} = $vars{dbversion};
		} elsif($k eq 'description'){
			$self->{_description} = $vars{description};
		} elsif($k eq 'context'){
			my $good = undef;
			foreach my $context(@ALL_ANNOTATION_CONTEXT){
				if($vars{context} eq $context){
					$good = $vars{context};
					last;
				}
			}
			if(defined($good)){
				$self->{_context} = $good;
			} else {
				croak('recieved unknown annotation context: '.$vars{context});
			}
		} elsif($k eq 'type'){
			my $good = undef;
			foreach my $type(@ALL_ANNOTATION_TYPES){
				if($vars{type} eq $type){
					$good = $vars{type};
					last;
				}
			}
			if(defined($good)){
				$self->{_type} = $good;
			} else {
				croak('recieved unknown annotation type: '.$vars{type});
			}
		} elsif($k eq 'subtype'){
			my $good = undef;
			foreach my $subtype(@ALL_ANNOTATION_SUBTYPES){
				if($vars{subtype} eq $subtype){
					$good = $vars{subtype};
					last;
				}
			}
			if(defined($good)){
				$self->{_subtype} = $good;
			} else {
				croak('recieved unknown annotation type: '.$vars{subtype});
			}
		}

	}
}

sub getCDSAnnotationContext {
	return $CDS_ANNOTATION_CONTEXT;
}

sub getmRNAAnnotationContext {
	return $MRNA_ANNOTATION_CONTEXT;
}

sub getProteinAnnotationContext {
	return $PROTEIN_ANNOTATION_CONTEXT;
}

sub getSubstitutionAnnotationType {
	return $SUBSTITUTION_ANNOTATION_TYPE;
}

sub getDeletionAnnotationType {
	return $DELETION_ANNOTATION_TYPE;
}

sub getComplexAnnotationType {
	return $COMPLEX_ANNOTATION_TYPE;
}

sub getInsertionAnnotationType {
	return $INSERTION_ANNOTATION_TYPE;
}

sub getFrameShiftAnnotationType {
	return $FRAMESHIFT_ANNOTATION_TYPE;
}

sub getUnknownAnnotationType {
	return $UNKNOWN_ANNOTATION_TYPE;
}

sub getPositionKnownSubtype {
	return $POSITION_KNOWN_SUBTYPE;
}

sub getPositionOffsetSubtype {
	return $POSITION_OFFSET_SUBTYPE;
}

sub getPositionOffSequenceSubtype {
	return $POSITION_OFF_SEQUENCE_SUBTYPE;
}

sub getMinPos {
	return shift->{_minpos};
}

sub getMaxPos {
	return shift->{_maxpos};
}

sub getMinOffset {
	return shift->{_minOffset};
}

sub getMaxOffset {
	return shift->{_maxOffset};
}

sub getWt {
	return shift->{_wt};
}

sub getMt {
	return shift->{_mt};
}

sub getDescription {
	return shift->{_description};
}

sub getContext {
	return shift->{_context};
}

sub getType {
	return shift->{_type};
}

sub getSubtype {
	return shift->{_subtype};
}

sub getAccession {
	return shift->{_acc};
}

sub getSequenceVersion {
	return shift->{_accversion};
}

sub getDatabase {
	return shift->{_db};
}

sub getDatabaseVersion {
	return shift->{_dbversion};
}

sub getSequenceLength {
	return shift->{_seqlength};
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
	if(defined($self->{_class})){
		if(first {$_ eq $class} @{$self->{_class}}){
			return 1;
		} else {
			return 0;
		}
	}
	return 0;
}


__END__

=head1 NAME

Sanger::CGP::Vagrent::Data::Annotation - Data object containing the description of an variation
within a sequence

=head1 DESCRIPTION

This holds annotation for a variation against a defined piece of sequence

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $annot = Sanger::CGP::Vagrent::Data::Annotation->new(%params);

=item Function :

Builds a new Sanger::CGP::Vagrent::Data::Annotation object

=item Returns :

L<Sanger::CGP::Vagrent::Data::Annotation|Sanger::CGP::Vagrent::Data::Annotation> object initialized with parameter values

=item Params :

 acc         => Accession for the sequence record being annotated against
 accversion  => version of the sequence record
 db          => source DB the sequence record came from
 dbversion   => version of the source DB
 seqlength   => length of sequence in the sequence record
 context     => sequence context of the annotation (defined by context constants)
 type        => variation type the annotation describes (defined by type constants)
 subtype     => variation type the annotation describes (defined by subtype constants)
 minpos      => minimum coordinate in the sequence
 maxpos      => maximum coordinate in the sequence
 minOffset   => offset distance away from the minpos (typically used for splice site mutations)
 maxOffset   => offset distance away from the maxpos (typically used for splice site mutations)
 wt          => wild type sequence of variant
 mt          => mutant sequence of variant
 description => HGVS style description string

=back

=head2 Attributes

=head3 getMinPos

=over

=item Usage :

 my $min = $a->getMinPos;

=item Function :

Returns the lowest coordinate of the annotation on the sequence

=item Returns :

Integer

=back

=head3 getMinOffset

=over

=item Usage :

 my $minOff = $a->getMinOffset;

=item Function :

Returns the offset value from the minimim coordinate

=item Returns :

Signed integer

=back

=head3 getMaxPos

=over

=item Usage :

 my $min = $a->getMaxPos;

=item Function :

Returns the highest coordinate of the annotation on the sequence

=item Returns :

Integer

=back

=head3 getMaxOffset

=over

=item Usage :

 my $maxOff = $a->getMaxOffset;

=item Function :

Returns the offset value from the maximum coordinate

=item Returns :

Signed integer

=back

=head3 getWt

=over

=item Usage :

 my $wtSeq = $a->getWt;

=item Function :

Returns the wildtype sequence string

=item Returns :

String

=back

=head3 getMt

=over

=item Usage :

 my $mtSeq = $a->getMt;

=item Function :

Returns the mutant sequence string

=item Returns :

String

=back

=head3 getDescription

=over

=item Usage :

 my $desc = $a->getDescription;

=item Function

An HGVS style description string of the annotation

=item Returns

String

=back

=head3 getContext

=over

=item Usage :

 my $ctx = $a->getContext;

=item Function

The context of the annotation, as defined by a context constant value

=item Returns

String

=back

=head3 getType

=over

=item Usage :

 my $type = $a->getType;

=item Function

The type of the annotation, as defined by a type constant value

=item Returns

String

=back

=head3 getSubtype

=over

=item Usage :

 my $subtype = $a->getSubtype;

=item Function

The subtype of the annotation, as defined by a subtype constant value

=item Returns

String

=back

=head3 getAccession

=over

=item Usage :

 my $accession = $a->getAccession;

=item Function

The accession of the sequence record the annotation is built against

=item Returns

String

=back

=head3 getSequenceVersion

=over

=item Usage :

 my $vers = $a->getSequenceVersion;

=item Function

The version of the sequence record the annotation is built against

=item Returns

String

=back

=head3 getDatabase

=over

=item Usage :

 my $db = $a->getDatabase;

=item Function

The database the sequence record originated from

=item Returns

String

=back

=head3 getDatabaseVersion

=over

=item Usage :

 my $dbVersion = $a->getDatabaseVersion;

=item Function

The version of the database the sequence record originated from

=item Returns

String

=back

=head3 getSequenceLength

=over

=item Usage :

 my $seqLength = $a->getSequenceLength;

=item Function

The full length of the sequence in the sequence record

=item Returns

Integer

=back

=head3 addClassification

=over

=item Usage :

 $a->addClassification($class1,$class2);

=item Function

Appends one or more classification terms to the annotation

=item Params

An array of classification terms (Strings)

=item Returns

None

=back

=head3 getClassifications

=over

=item Usage :

 my @classes = $a->getClassifications;

=item Function

Returns the stored classification terms

=item Returns

Array of Strings, or undef if none

=back

=head2 Functions

=head3 hasClassification

=over

=item Usage :

 if($a->hasClassification($checkClass)){
  .....
 }

=item Function

Checks an annotation to see if it already contains the specified class

=item Returns

Boolean, 1 or 0

=back

=head2 Constants

=head3 getmRNAAnnotationContext

=over

=item Usage :

 my $ctx = $anno->getmRNAAnnotationContext();
 my $ctx = Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext();

=item Function :

Constant lookup, returns the mRNA Annotation Context value

=item Returns :

String

=back

=head3 getCDSAnnotationContext

=over

=item Usage :

 my $ctx = $anno->getCDSAnnotationContext();
 my $ctx = Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext();

=item Function :

Constant lookup, returns the CDS Annotation Context value

=item Returns :

String

=back

=head3 getProteinAnnotationContext

=over

=item Usage :

 my $ctx = $anno->getProteinAnnotationContext();
 my $ctx = Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext();

=item Function :

Constant lookup, returns the Protein Annotation Context value

=item Returns :

String

=back

=head3 getSubstitutionAnnotationType

=over

=item Usage :

 my $type = $anno->getSubstitutionAnnotationType();
 my $type = Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType();

=item Function :

Constant lookup, returns the Substitution Annotation Type value

=item Returns :

String

=back

=head3 getDeletionAnnotationType

=over

=item Usage :

 my $type = $anno->getDeletionAnnotationType();
 my $type = Sanger::CGP::Vagrent::Data::Annotation::getDeletionAnnotationType();

=item Function :

Constant lookup, returns the Deletion Annotation Type value

=item Returns :

String

=back

=head3 getInsertionAnnotationType

=over

=item Usage :

 my $type = $anno->getInsertionAnnotationType();
 my $type = Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType();

=item Function :

Constant lookup, returns the Insertion Annotation Type value

=item Returns :

String

=back

=head3 getComplexAnnotationType

=over

=item Usage :

 my $type = $anno->getComplexAnnotationType();
 my $type = Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType();

=item Function :

Constant lookup, returns the Complex Annotation Type value

=item Returns :

String

=back

=head3 getFrameShiftAnnotationType

=over

=item Usage :

 my $type = $anno->getFrameShiftAnnotationType();
 my $type = Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType();

=item Function :

Constant lookup, returns the Frameshift Annotation Type value

=item Returns :

String

=back

=head3 getUnknownAnnotationType

=over

=item Usage :

 my $type = $anno->getUnknownAnnotationType();
 my $type = Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType();

=item Function :

Constant lookup, returns the Unknown Annotation Type value

=item Returns :

String

=back

=head3 getPositionKnownSubtype

=over

=item Usage :

 my $subtype = $anno->getPositionKnownSubtype();
 my $subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype();

=item Function :

Constant lookup, returns the Position Known Subtype value

=item Returns :

String

=back

=head3 getPositionOffsetSubtype

=over

=item Usage :

 my $subtype = $anno->getPositionOffsetSubtype();
 my $subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype();

=item Function :

Constant lookup, returns the Position Offset Subtype value

=item Returns :

String

=back

=head3 getPositionOffSequenceSubtype

=over

=item Usage :

 my $subtype = $anno->getPositionOffSequenceSubtype();
 my $subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype();

=item Function :

Constant lookup, returns the Position Off-Sequence Subtype value

=item Returns :

String

=back


