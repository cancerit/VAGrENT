package Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator;

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

use Bio::Seq;

use Carp qw(cluck);
use Log::Log4perl;
use POSIX qw(ceil);
use Data::Dumper;
use Sanger::CGP::Vagrent qw($VERSION);

use base qw(Sanger::CGP::Vagrent::Annotators::DeletionAnnotator);

my $log = Log::Log4perl->get_logger(__PACKAGE__);

1;

sub _getAnnotation {
	my ($self,$var) = @_;
	$self->_clearMessages();
	unless(defined($var) && $var->isa('Sanger::CGP::Vagrent::Data::ComplexIndel')){
		my $msg = 'require a Sanger::CGP::Vagrent::Data::ComplexIndel object not a '.ref($var);
		$self->addMessage($msg);
		$log->info($msg);
		return undef;
	}
	unless($var->isValid){
		my $msg = 'complex indel not valid';
		$self->addMessage($msg);
		$log->error($msg);
		return undef;
	}
	my @trans = $self->_getTranscriptSource->getTranscripts($var);
	unless(defined($trans[0])){
		my $msg = 'No transcripts returned from transcript source';
		$self->addMessage($msg);
		$log->info($msg);
		return undef;
	}
	my @groups;
	foreach my $t(@trans){
		my $g = $self->_generateAnnotatonGroup($var,$t);
		if(defined($g)){
			push(@groups,$g);
		}
	}
	unless(scalar(@groups) > 0 && defined($groups[0])){
		my $msg = 'No annotation groups generated';
		$self->addMessage($msg);
		$log->info($msg);
		return undef;
	}
	return @groups;
}

sub _generateAnnotatonGroup {
	my ($self,$var,$tran) = @_;
	my ($rAnnot,@groupClasses) = $self->_buildRNAAnnotation($var,$tran);
	unless(defined($rAnnot)){
		my $msg = 'No mRNA annotation created';
		$self->addMessage($msg);
		$log->info($msg);
		return undef;
	}
	my $group = Sanger::CGP::Vagrent::Data::AnnotationGroup->new(	accession => $tran->getAccession,
																	label => $tran->getGeneName,
																    ccds => $tran->getCCDS,
																	type => $tran->getGeneType,);

	if($tran->isProteinCoding){
 		if(	( $rAnnot->hasClassification($self->getIntronVariantClass) ||
  		   		$rAnnot->hasClassification($self->get5KBUpStreamVariantClass) ||
  		   		$rAnnot->hasClassification($self->get2KBUpStreamVariantClass) ||
  		   		$rAnnot->hasClassification($self->get5KBDownStreamVariantClass) ||
  		   		$rAnnot->hasClassification($self->get500BPDownStreamVariantClass)) &&
 			 !$rAnnot->hasClassification($self->getEssentialSpliceSiteVariantClass) &&
 			 !$rAnnot->hasClassification($self->getSpliceRegionVariantClass) &&
 			 !$rAnnot->hasClassification($self->getFrameShiftVariantClass) &&
 			 !$rAnnot->hasClassification($self->getInFrameVariantClass) &&
 			 !$rAnnot->hasClassification($self->get5PrimeUtrVariantClass) &&
 			 !$rAnnot->hasClassification($self->get3PrimeUtrVariantClass) &&
 			 !$rAnnot->hasClassification($self->getComplexChangeVariantClass)){
 			# Inronic or up/down stream mutations don't need to get any further annotation
 			$group->addAnnotation($rAnnot);
 		} else {
 			my $cAnnot = $self->_buildCDSAnnotation($var,$tran,$rAnnot);
 			unless(defined($cAnnot)){
 				my $msg = 'No CDS annotation created';
 				$self->addMessage($msg);
 				$log->info($msg);
 				return undef;
 			}
 			my $pAnnot = $self->_buildProteinAnnotation($var,$tran,$cAnnot,$rAnnot);
 			unless(defined($pAnnot)){
 				my $msg = 'No Protein annotation created';
 				$self->addMessage($msg);
 				$log->info($msg);
 				return undef;
 			}
 			$group->addAnnotation($rAnnot);
 			$group->addAnnotation($cAnnot);
 			$group->addAnnotation($pAnnot);
 		}
	} else {
		$group->addAnnotation($rAnnot);
	}
	$group->addClassification(@groupClasses);
	return $group;
}

sub _buildRNAAnnotation {
	my ($self,$var,$tran) = @_;
	my ($mrnaMin,$mrnaMinOffset,$mrnaMax,$mrnaMaxOffset) = $self->_getmRNAPositions($var,$tran);
	unless(defined($mrnaMin) && defined($mrnaMinOffset) && defined($mrnaMax) && defined($mrnaMaxOffset)){
		my $msg = 'problem generating mrna coordinates';
		$self->addMessage($msg);
		$log->error($msg);
		return undef;
	}
	unless($self->_safetyCheck($var,$tran,$mrnaMin,$mrnaMinOffset,$mrnaMax,$mrnaMaxOffset)){
		# something has gone wrong,
		my $msg = 'The affected transcript WT sequence in the deleted sequence doesnt match the affected transcript WT sequence from the transcript';
		$self->addMessage($msg);
		$log->error($msg);
		return undef;
	}

	my @classes = ($self->getComplexIndelClass);
	my @groupClasses = ($self->classifyTranscript($tran));

	my $upstream =  $self->_upstreamVariantCheck($mrnaMin,$mrnaMinOffset,$mrnaMax,$mrnaMaxOffset,\@classes);
	if(defined $upstream){
		if($upstream){
			return ($self->_buildUnknownMRNAAnnotation($var,$tran,@classes),@groupClasses);
		}
	} else {
		return undef;
	}

	my $downstream =  $self->_downstreamVariantCheck($mrnaMin,$mrnaMinOffset,$mrnaMax,$mrnaMaxOffset,\@classes);
	if(defined $downstream){
		if($downstream){
			return ($self->_buildUnknownMRNAAnnotation($var,$tran,@classes),@groupClasses);
		}
	} else {
		return undef;
	}

	# delt with up/down stream issues, can reset coordinates to the extremes of the transcript;
	if($mrnaMin == 0 && $mrnaMinOffset <= 0){
		$mrnaMin = 1;
		$mrnaMinOffset = 0;
	}
	if($mrnaMax == 0 && $mrnaMaxOffset >= 0){
		$mrnaMax = length($tran->getcDNASeq);
		$mrnaMaxOffset = 0;
	}

	if($tran->isProteinCoding){
    if(($mrnaMax > $tran->getCdsMinPos || ($mrnaMax == $tran->getCdsMinPos && $mrnaMaxOffset >= 0)) &&
       ($mrnaMin < $tran->getCdsMaxPos || ($mrnaMin == $tran->getCdsMaxPos && $mrnaMinOffset <= 0))){
#		if($mrnaMax >= $tran->getCdsMinPos && $mrnaMin <= $tran->getCdsMaxPos){
			# coding change
			push(@groupClasses,$self->getCDSClass);
		}
		if($mrnaMin < $tran->getCdsMinPos || ($mrnaMin == $tran->getCdsMinPos && $mrnaMinOffset < 0)){
		#if($mrnaMin < $tran->getCdsMinPos){
			# 5prime UTR
			push(@groupClasses,$self->get5PrimeUtrClass);
			push(@classes,$self->get5PrimeUtrVariantClass) unless($self->_arrayHasString($self->getCDSClass,@groupClasses));
		}
		if($mrnaMax > $tran->getCdsMaxPos || ($mrnaMax == $tran->getCdsMaxPos && $mrnaMaxOffset > 0)){
		#if($mrnaMax > $tran->getCdsMaxPos){
			# 3prime UTR
			push(@groupClasses,$self->get3PrimeUtrClass);
			push(@classes,$self->get3PrimeUtrVariantClass) unless($self->_arrayHasString($self->getCDSClass,@groupClasses));
		}
	}

	my $tmpGroupClassesHash = $self->_classifyDeletion($tran,$mrnaMin,$mrnaMinOffset,$mrnaMax,$mrnaMaxOffset);
	push(@groupClasses, sort keys %$tmpGroupClassesHash);
	if(scalar(keys %$tmpGroupClassesHash) == 1 && defined($tmpGroupClassesHash->{$self->getIntronClass})){
		# its intron only
		return ($self->_buildUnknownMRNAAnnotation($var,$tran,$self->getComplexIndelClass,$self->getIntronVariantClass),@groupClasses);
	}

	my $wt = $self->_getWildTypeStringForRNAAnno($var,$tran);
	$wt =~ tr/Tt/Uu/;
	my $mt = undef;
	if($tran->getStrand == 1){
		$mt = lc($var->getInsertedSequence());
	} else {
		$mt = $self->_revcompSeq(lc($var->getInsertedSequence()));
	}
	$mt =~ tr/Tt/Uu/;

	my $desc = undef;

	if($mrnaMin == 1 && $mrnaMax == length($tran->getcDNASeq)){
		$desc = 'r.0';
	} else {
		$desc = $self->_generateDNAindelDescriptionString('r.',$mrnaMin,$mrnaMax,$mrnaMinOffset,$mrnaMaxOffset,$wt,$mt);
	}

	my $subtype = undef;

	if($mrnaMin == 0 && $mrnaMax == 0){
		$subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype();
	} elsif($mrnaMinOffset == 0 && $mrnaMaxOffset == 0){
		$subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype();
	} else {
		$subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype();
	}

	my $anno = Sanger::CGP::Vagrent::Data::Annotation->new( wt => uc($wt),
															mt => uc($mt),
															minpos => $mrnaMin,
															minOffset => $mrnaMinOffset,
															maxpos => $mrnaMax,
															maxOffset => $mrnaMaxOffset,
															acc => $tran->getAccession,
															accversion => $tran->getAccessionVersion,
															db => $tran->getDatabase,
															dbversion => $tran->getDatabaseVersion,
															seqlength => length($tran->getcDNASeq),
															description => $desc,
															context => Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
															type => Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
															subtype => $subtype);

	if(defined($tmpGroupClassesHash->{$self->getExonClass}) && defined($tmpGroupClassesHash->{$self->getEssentialSpliceSiteClass})){
		# special case, variation spans over the end of the exon.  This is a complex change in the transcript
		push(@classes,$self->getComplexChangeVariantClass);
		unless($tran->isProteinCoding){
			push(@classes,$self->getNonCodingTranscriptVariantClass);
		}

	} else {
		if(defined($tmpGroupClassesHash->{$self->getSpliceRegionClass})){
			push(@classes,$self->getSpliceRegionVariantClass);
		}
		if(defined($tmpGroupClassesHash->{$self->getEssentialSpliceSiteClass})){
			push(@classes,$self->getEssentialSpliceSiteVariantClass);
		}
		if(defined($tmpGroupClassesHash->{$self->getIntronClass})){
			push(@classes,$self->getIntronVariantClass);
		}
		if(defined($tmpGroupClassesHash->{$self->getExonClass})){
			if($tran->isProteinCoding){
				if($self->_arrayHasString($self->getCDSClass,@groupClasses)){
					if($self->_arrayHasString($self->get2KBUpStreamVariantClass,@classes) || $self->_arrayHasString($self->get5PrimeUtrClass,@groupClasses) || $self->_arrayHasString($self->get3PrimeUtrClass,@groupClasses) || $self->_arrayHasString($self->get500BPDownStreamVariantClass,@classes)){
						# it overhangs the CDS-to-UTR/upstream/downstream boundry
						push(@classes,$self->getComplexChangeVariantClass);
					} else {
						if(abs(length($wt) - length($mt)) % 3 == 0){
							push(@classes,$self->getInFrameVariantClass);
						} else {
							push(@classes,$self->getFrameShiftVariantClass);
						}
					}
				}
			} else {
				push(@classes,$self->getNonCodingTranscriptVariantClass);
			}
		}
	}

	$anno->addClassification(@classes);

	return ($anno,@groupClasses);
}

sub _getMutantStringForCDSAnno {
	my ($self,$var,$tran,$rAnnot) = @_;
  	my $mt = $rAnnot->getMt;
 	$mt =~ s/u/t/ig;
 	return $mt;
}

sub _getCDSDescriptionString {
	my ($self,$tran,$mutStart,$mutEnd,$mutStartOffset,$mutEndOffset,$wt,$mt) = @_;
  	my $desc = undef;
 	if($mutStart == 1 && $mutEnd == length($tran->getCdsSeq)){
 		$desc = 'c.0';
 	} else {
 		$desc = $self->_generateDNAindelDescriptionString('c.',$mutStart,$mutEnd,$mutStartOffset,$mutEndOffset,uc($wt),uc($mt));
 	}
  	return $desc;
}

sub _getMutatedCdsSequence {
	my ($self,$wtDna,$min,$max,$mt) = @_;
 	my $codingDelLength = ($max - $min) + 1;
 	my $mtDna = $wtDna;
 	substr($mtDna,($min - 1),$codingDelLength,$mt);
	return $mtDna
}

sub _generateDNAindelDescriptionString {
	my ($self,$pre,$mutStart,$mutEnd,$mutStartOffset,$mutEndOffset,$wt,$mt) = @_;
	my $desc = 'WIBBLE';
	if($mutStart == $mutEnd && $mutStartOffset == $mutEndOffset){
		# single base change
		if($mutStartOffset > 0){
			$desc = $pre.$mutStart.'+'.$mutStartOffset.'del'.$wt;
		} elsif($mutStartOffset < 0){
			$desc = $pre.$mutStart.$mutStartOffset.'del'.$wt;
		} elsif ($mutStartOffset == 0) {
			$desc = $pre.$mutStart.'del'.$wt;
		}
		$desc .= 'ins'.$mt;

	} else {
		# multi base change
		if($mutStartOffset > 0){
			$desc = $pre.$mutStart.'+'.$mutStartOffset.'_';
		} elsif($mutStartOffset < 0){
			$desc = $pre.$mutStart.$mutStartOffset.'_';
		} elsif($mutStartOffset == 0) {
			$desc = $pre.$mutStart.'_';
		}

		if($mutEndOffset > 0){
			$desc .= $mutEnd.'+'.$mutEndOffset.'del';
		} elsif($mutEndOffset < 0){
			$desc .= $mutEnd.$mutEndOffset.'del';
		} elsif($mutEndOffset == 0){
			$desc .= $mutEnd.'del';
		}

		if(length($wt) > 15){
			$desc .= length($wt);
		} else {
			$desc .= $wt;
		}

		$desc .= 'ins'.$mt;
	}
	return $desc;
}

sub _getDefaultCDSAnnotationType {
	return Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType();
}

sub _getDefaultProteinAnnotationType {
	return Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType();
}


__END__

=head1 NAME

Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator - Annotator for complex indel variants

=head1 DESCRIPTION

This annotates complex insertion/deletion variants, it provides L<AnnotatonGroups|Sanger::CGP::Vagrent::Data::AnnotationGroup> for each transcript returned from the L<TranscriptSource|Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource>

It will only process L<Sanger::CGP::Vagrent::Data::ComplexIndel|Sanger::CGP::Vagrent::Data::ComplexIndel> objects, if any other L<Variation|Sanger::CGP::Vagrent::Data::AbstractVariation> objects are sent in it will return an empty answer.

It inherits from L<Sanger::CGP::Vagrent::Annotators::DeletionAnnotator|Sanger::CGP::Vagrent::Annotators::DeletionAnnotator>.

=head1 METHODS

See L<Sanger::CGP::Vagrent::Annotators::AbstractAnnotator|Sanger::CGP::Vagrent::Annotators::AbstractAnnotator>
