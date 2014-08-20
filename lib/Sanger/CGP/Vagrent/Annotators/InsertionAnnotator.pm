package Sanger::CGP::Vagrent::Annotators::InsertionAnnotator;

use strict;

use Bio::Seq;

use Carp qw(cluck);
use Log::Log4perl;
use POSIX qw(ceil);
use Data::Dumper;

use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Sanger::CGP::Vagrent::Annotators::AbstractAnnotator);

my $log = Log::Log4perl->get_logger(__PACKAGE__);

1;

sub _getAnnotation {
	my ($self,$var) = @_;
	$self->_clearMessages();
	unless(defined($var) && $var->isa('Sanger::CGP::Vagrent::Data::Insertion')){
		my $msg = 'require a Sanger::CGP::Vagrent::Data::Insertion object not a '.ref($var);
		$self->addMessage($msg);
		$log->info($msg);
		return undef;
	}
	unless($var->isValid){
		my $msg = 'insertion not valid';
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
		if($rAnnot->hasClassification($self->getIntronVariantClass) ||
 		   $rAnnot->hasClassification($self->get5KBUpStreamVariantClass) ||
 		   $rAnnot->hasClassification($self->get2KBUpStreamVariantClass) ||
 		   $rAnnot->hasClassification($self->get5KBDownStreamVariantClass) ||
 		   $rAnnot->hasClassification($self->get500BPDownStreamVariantClass)){
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

sub _getMutatedCdsSequence {
	my ($self,$wtDna,$min,$max,$mt) = @_;
	my $mtDna = substr($wtDna,0,$min) . $mt . substr($wtDna,$min);
	return $mtDna
}

sub _getCdsMinPosForProteinCalculation {
	my ($self,$cAnnot) = @_;
	my $min = undef;
	if($cAnnot->getMinOffset == 0){
 		$min =  $cAnnot->getMinPos;
 	} elsif($cAnnot->getMinOffset == -1){
 		$min = $cAnnot->getMinPos - 1
 	} else {
 		my $msg = 'Unexpected cds min offset, expected 0 or -1 recieved: '.$cAnnot->getMinOffset;
 		$self->addMessage($msg);
 		$log->info($msg);
 	}
 	return $min;
}

sub _getCdsMaxPosForProteinCalculation {
	my ($self,$cAnnot) = @_;
  	return $cAnnot->getMaxPos();
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
	unless(($mrnaMax - $mrnaMin == 1 && $mrnaMinOffset == 0 && $mrnaMaxOffset == 0) ||        # exon
				 ($mrnaMax == $mrnaMin && $mrnaMaxOffset - $mrnaMinOffset == 1) ||                  # intron/splice/off transcript
				 ($mrnaMin == 0 && $mrnaMax == 1 && $mrnaMinOffset == -1 && $mrnaMaxOffset == 0) || # at start of transcript
				 ($mrnaMin == length($tran->getcDNASeq) && $mrnaMax == 0 && $mrnaMinOffset == 0 && $mrnaMaxOffset == 1)){ # at end of transcript
		my $msg = "This is an insertion so the coordinates should be adjacent to each other, coordinates returned were MIN: $mrnaMin OFFSET:$mrnaMinOffset and MAX: $mrnaMax OFFSET: $mrnaMaxOffset";
		$self->addMessage($msg);
		$log->error($msg);
		return undef;
	}
	my @groupClasses = ($self->classifyTranscript($tran));
	if($mrnaMin == 0 && $mrnaMinOffset < 0){
		if($self->_isWithin2KBUpstreamOffsetDistance($mrnaMinOffset)){
			return ($self->_buildUnknownMRNAAnnotation($var,$tran,$self->getInsertionClass,$self->get2KBUpStreamVariantClass),@groupClasses);
		} elsif($self->_isWithin5KBOffsetDistance($mrnaMinOffset)){
			return ($self->_buildUnknownMRNAAnnotation($var,$tran,$self->getInsertionClass,$self->get5KBUpStreamVariantClass),@groupClasses);
		} else {
			my $msg = "variant isnt close enough to this transcript, nothing to do";
			$self->addMessage($msg);
			$log->error($msg);
			return undef;
		}
	} elsif($mrnaMax == 0 && $mrnaMaxOffset > 0){
		if($self->_isWithin500BPDownstreamOffsetDistance($mrnaMaxOffset)){
			return ($self->_buildUnknownMRNAAnnotation($var,$tran,$self->getInsertionClass,$self->get500BPDownStreamVariantClass),@groupClasses);
		} elsif($self->_isWithin5KBOffsetDistance($mrnaMaxOffset)){
			return ($self->_buildUnknownMRNAAnnotation($var,$tran,$self->getInsertionClass,$self->get5KBDownStreamVariantClass),@groupClasses);
		} else {
			my $msg = "variant isnt close enough to this transcript, nothing to do";
			$self->addMessage($msg);
			$log->error($msg);
			return undef;
		}
	}
	my $wt = '-';
	my $mt = undef;
	if($tran->getStrand == 1){
		$mt = lc($var->getInsertedSequence());
	} else {
		$mt = $self->_revcompSeq(lc($var->getInsertedSequence()));
	}
	$mt =~ s/t/u/ig;

	my $desc = undef;
	my $subtype = undef;
	my @classes = ($self->getInsertionClass);

	if($tran->isProteinCoding){
		if($mrnaMax > $tran->getCdsMinPos && $mrnaMin < $tran->getCdsMaxPos){
			# coding change
			push(@groupClasses,$self->getCDSClass);
		} elsif($mrnaMax <= $tran->getCdsMinPos){
			# 5prime UTR
			push(@groupClasses,$self->get5PrimeUtrClass);
			push(@classes,$self->get5PrimeUtrVariantClass);
		} elsif($mrnaMin >= $tran->getCdsMaxPos){
			# 3prime UTR
			push(@groupClasses,$self->get3PrimeUtrClass);
			push(@classes,$self->get3PrimeUtrVariantClass);
		} else {
			my $msg = "strange positional info, can't work out subtype : var $mrnaMin $mrnaMinOffset $mrnaMax $mrnaMaxOffset transcript ".$tran->getCdsMinPos." - ".$tran->getCdsMaxPos;
			$self->addMessage($msg);
			$log->error($msg);
			return undef;
		}
	}

	if(($mrnaMinOffset < 0 && $self->_isIntronicOffsetDistance($mrnaMinOffset)) || ($mrnaMaxOffset > 0 && $self->_isIntronicOffsetDistance($mrnaMaxOffset))){
		# its intronic
		push(@groupClasses,$self->getIntronClass);
		return ($self->_buildUnknownMRNAAnnotation($var,$tran,$self->getInsertionClass,$self->getIntronVariantClass),@groupClasses);
	} elsif(($mrnaMinOffset == 0 && $mrnaMaxOffset == 0) || ($mrnaMinOffset == -1 && $mrnaMaxOffset == 0) || ($mrnaMinOffset == 0 && $mrnaMaxOffset == 1)){
		# its in an exon
		push(@groupClasses,$self->getExonClass);
		if($mrnaMinOffset != $mrnaMaxOffset) {
			$subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype();
		} else {
			$subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype();
		}
		if($tran->isProteinCoding){
			if($self->_arrayHasString($self->getCDSClass,@groupClasses)){
				if(length($mt) % 3 == 0){
					push(@classes,$self->getInFrameVariantClass);
				} else {
					push(@classes,$self->getFrameShiftVariantClass);
				}
			}
		} else {
			push(@classes,$self->getNonCodingTranscriptVariantClass);
		}
		$desc = 'r.';
		if($mrnaMinOffset == 0){
			$desc .= $mrnaMin.'_';
		} else {
			$desc .= $mrnaMin.$mrnaMinOffset.'_';
		}
		if($mrnaMaxOffset == 0){
			$desc .= $mrnaMax;
		} else {
			$desc .= $mrnaMax.'+'.$mrnaMaxOffset;
		}
		$desc .= 'ins'.$mt;
	} else {
		# its not in an exon, and its not intronic.  Must be related to splice site.
		$subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype();
		if(($mrnaMinOffset < 0 && $mrnaMinOffset >= $self->_getConsesnsusSpliceBeforeBoundry) || ($mrnaMinOffset > 0 && $mrnaMaxOffset <= $self->_getConsesnsusSpliceAfterBoundry)){
			# Essential Splice site change
			push(@classes,$self->getEssentialSpliceSiteVariantClass);
			push(@groupClasses,$self->getEssentialSpliceSiteClass);
		} else {
			# splice region change
			push(@classes,$self->getSpliceRegionVariantClass);
			push(@groupClasses,$self->getSpliceRegionClass);
		}
		if($mrnaMinOffset < 0){
			$desc = 'r.'.$mrnaMin.$mrnaMinOffset.'_'.$mrnaMax.$mrnaMaxOffset.'ins'.$mt;
		} else {
			$desc = 'r.'.$mrnaMin.'+'.$mrnaMinOffset.'_'.$mrnaMax.'+'.$mrnaMaxOffset.'ins'.$mt;
		}
	}

	my $anno = Sanger::CGP::Vagrent::Data::Annotation->new(wt => uc($wt),
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
																												type => Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
																												subtype => $subtype);

	$anno->addClassification(@classes);

	return ($anno,@groupClasses);
}

sub _getWildTypeStringForCDSAnno {
	my ($self,$var,$tran,$rAnnot) = @_;
	my $wt = $rAnnot->getWt();
	$wt =~ s/u/t/ig;
	return $wt;
}

sub _getMutantStringForCDSAnno {
	my ($self,$var,$tran,$rAnnot) = @_;
	my $mt = $rAnnot->getMt();
	$mt =~ s/u/t/ig;
	return $mt;
}

sub _getCDSDescriptionString {
	my ($self,$tran,$mutStart,$mutEnd,$mutStartOffset,$mutEndOffset,$wt,$mt) = @_;
  	my $desc = 'c.';
 	if($mutStartOffset == 0){
 		$desc .= $mutStart.'_';
 	} elsif ($mutStartOffset > 0){
 		$desc .= $mutStart.'+'.$mutStartOffset.'_';
 	} else {
 		$desc .= $mutStart.$mutStartOffset.'_';
 	}
 	if($mutEndOffset == 0){
 		$desc .= $mutEnd;
 	} elsif($mutEndOffset > 0) {
 		$desc .= $mutEnd.'+'.$mutEndOffset;
 	} else {
 		$desc .= $mutEnd.$mutEndOffset;
 	}
 	$desc .= 'ins'.uc($mt);
 	return $desc;
}

sub _getDefaultCDSAnnotationType {
	return Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType();
}

sub _getDefaultProteinAnnotationType {
	return Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType();
}

=head1 NAME

Sanger::CGP::Vagrent::Annotators::InsertionAnnotator - Annotator for insertion variants

=head1 DESCRIPTION

This annotates insertion variants, it provides L<AnnotatonGroups|Sanger::CGP::Vagrent::Data::AnnotationGroup> for each transcript returned from the L<TranscriptSource|Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource>

It will only process L<Sanger::CGP::Vagrent::Data::Insertion|Sanger::CGP::Vagrent::Data::Insertion> objects, if any other L<Variation|Sanger::CGP::Vagrent::Data::AbstractVariation> objects are sent in it will return an empty answer.

It inherits from L<Sanger::CGP::Vagrent::Annotators::AbstractAnnotator|Sanger::CGP::Vagrent::Annotators::AbstractAnnotator>.

=head1 METHODS

See L<Sanger::CGP::Vagrent::Annotators::AbstractAnnotator|Sanger::CGP::Vagrent::Annotators::AbstractAnnotator>
