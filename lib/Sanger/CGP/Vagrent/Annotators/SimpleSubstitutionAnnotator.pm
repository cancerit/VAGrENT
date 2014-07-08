package Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator;

use strict;

use Bio::Seq;

use Log::Log4perl;
use Carp qw(cluck);
use Log::Log4perl qw(:easy);
use POSIX qw(ceil);
use Data::Dumper;

use base qw(Sanger::CGP::Vagrent::Annotators::AbstractAnnotator);

my $log = Log::Log4perl->get_logger(__PACKAGE__);

1;

sub _getAnnotation {
 	my ($self,$var) = @_;
 	$self->_clearMessages();
 	unless(defined($var) && $var->isa('Sanger::CGP::Vagrent::Data::Substitution')){
 		my $msg = 'require a Sanger::CGP::Vagrent::Data::Substitution object, not a '.ref($var);
 		$self->addMessage($msg);
 		$log->info($msg);
 		return undef;
 	}
 	unless($var->isValid){
 		my $msg = 'substitution not valid';
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
 			# mutation is only indirectly related to the transcript, only need the mRNA annotation
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
	unless($mrnaMin == $mrnaMax && $mrnaMinOffset == $mrnaMaxOffset){
		my $msg = "This should be a single base mutation, coordinates returned were MIN: $mrnaMin OFFSET:$mrnaMinOffset and MAX: $mrnaMax OFFSET: $mrnaMaxOffset";
		$self->addMessage($msg);
		$log->error($msg);
		return undef;
	}

	# This should simplify things from here on in, only need one position and one offset
	# we have just checked that they are the same
	my $pos = $mrnaMin;
	my $offset = $mrnaMinOffset;
	my $wt = undef;
	my $mt = undef;
	my $desc = undef;
	my $subtype = undef;
	my @classes = ($self->getSubstitutionClass);
	my @groupClasses = ($self->classifyTranscript($tran));

	if($pos == 0){
		# the variant is off the transcript, have to do the up/down stream check
		if($offset < 0){
			# before start of transcript
			if($self->_isWithin2KBUpstreamOffsetDistance($offset)){
				return ($self->_buildUnknownMRNAAnnotation($var,$tran,$self->getSubstitutionClass,$self->get2KBUpStreamVariantClass),@groupClasses);
			} elsif($self->_isWithin5KBOffsetDistance($offset)){
				return ($self->_buildUnknownMRNAAnnotation($var,$tran,$self->getSubstitutionClass,$self->get5KBUpStreamVariantClass),@groupClasses);
			} else {
				my $msg = "variant isnt close enough to this transcript, nothing to do";
				$self->addMessage($msg);
				$log->error($msg);
				return undef;
			}

		} elsif($offset > 0){
			# after end of transcript
			if($self->_isWithin500BPDownstreamOffsetDistance($offset)){
				return ($self->_buildUnknownMRNAAnnotation($var,$tran,$self->getSubstitutionClass,$self->get500BPDownStreamVariantClass),@groupClasses);
			} elsif($self->_isWithin5KBOffsetDistance($offset)){
				return ($self->_buildUnknownMRNAAnnotation($var,$tran,$self->getSubstitutionClass,$self->get5KBDownStreamVariantClass),@groupClasses);
			} else {
				my $msg = "variant isnt close enough to this transcript, nothing to do";
				$self->addMessage($msg);
				$log->error($msg);
				return undef;
			}
		} else {
			my $msg = "strange positional info, a variant with a position of 0 must have an offset";
			$self->addMessage($msg);
			$log->error($msg);
			return undef;
		}
	}

	if($tran->isProteinCoding){
		if($pos >= $tran->getCdsMinPos && $pos <= $tran->getCdsMaxPos){
			# coding change
			push(@groupClasses,$self->getCDSClass);
		} elsif($mrnaMax < $tran->getCdsMinPos){
			# 5prime UTR
			push(@groupClasses,$self->get5PrimeUtrClass);
		} elsif($mrnaMin > $tran->getCdsMaxPos){
			# 3prime UTR
			push(@groupClasses,$self->get3PrimeUtrClass);
		} else {
			my $msg = "strange positional info, can't work out subtype : var $mrnaMin $mrnaMinOffset $mrnaMax $mrnaMaxOffset transcript ".$tran->getCdsMinPos." - ".$tran->getCdsMaxPos;
			$self->addMessage($msg);
			$log->error($msg);
			return undef;
		}
	}

	if($self->_isIntronicOffsetDistance($offset)){
		# its intronic
		push(@groupClasses,$self->getIntronClass);
		return ($self->_buildUnknownMRNAAnnotation($var,$tran,$self->getSubstitutionClass,$self->getIntronVariantClass),@groupClasses);
	}
	if($tran->getStrand == 1){
		# transcript is on the same strand as the genome
		$wt = $var->getWt();
		$mt = $var->getMt();
	} else {
		# transcript reversed on the genome
		$wt = $self->_revcompSeq($var->getWt());
		$mt = $self->_revcompSeq($var->getMt());
	}
	if($offset == 0){
		# its in an exon
		push(@groupClasses,$self->getExonClass);
		$subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype();
		# sanity check - does the substring on the cDNA sequence at the mut location equal wt
		my $substr = substr($tran->getcDNASeq,($pos - 1),1);
		if(lc($substr) ne lc($wt)){
			my $msg = "calculated wt doesn't match substring - $wt vs $substr";
			$self->addMessage($msg);
			$log->error($msg);
			return undef;
		}
		if($tran->isProteinCoding){
			if($self->_arrayHasString($self->get5PrimeUtrClass,@groupClasses)) {
				if($self->_isStartGained($var,$tran,$pos,$pos,$wt,$mt)){
					push(@classes,$self->getPrematureStartGainedVariantClass);
				} else {
					push(@classes,$self->get5PrimeUtrVariantClass);
				}
			} elsif($self->_arrayHasString($self->get3PrimeUtrClass,@groupClasses)){
				push(@classes,$self->get3PrimeUtrVariantClass);
			} else {
				push(@classes,$self->getCodonVariantClass);
			}
		} else {
			push(@classes,$self->getNonCodingTranscriptVariantClass);
		}
	} else {
		# its not in an exon, and its not intronic.  Must be related to splice site.
		$subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype();
		if($self->_isOffsetAConsensusSpliceDistance($offset)){
			push(@classes,$self->getEssentialSpliceSiteVariantClass);
			push(@groupClasses,$self->getEssentialSpliceSiteClass);
		} else {
			push(@classes,$self->getSpliceRegionVariantClass);
			push(@groupClasses,$self->getSpliceRegionClass);
		}
		if($self->_arrayHasString($self->get5PrimeUtrClass,@groupClasses)){
			push(@classes,$self->get5PrimeUtrVariantClass);
		} elsif($self->_arrayHasString($self->get3PrimeUtrClass,@groupClasses)){
			push(@classes,$self->get3PrimeUtrVariantClass);
		}
	}

	$wt =~ s/t/u/ig;
	$mt =~ s/t/u/ig;
	if($offset == 0){
		$desc = 'r.'.$pos . lc($wt) . '>' . lc($mt);
	} elsif($offset > 0){
		$desc = 'r.'.$pos.'+'.$offset. lc($wt).'>'.lc($mt);
	} else {
		$desc = 'r.'.$pos.$offset.lc($wt).'>'.lc($mt);
	}



	my $anno = Sanger::CGP::Vagrent::Data::Annotation->new(wt => uc($wt),
																												mt => uc($mt),
																												minpos => $pos,
																									    	minOffset => $offset,
																												maxpos => $pos,
																												maxOffset => $offset,
																												acc => $tran->getAccession,
																												accversion => $tran->getAccessionVersion,
																												db => $tran->getDatabase,
																												dbversion => $tran->getDatabaseVersion,
																												seqlength => length($tran->getcDNASeq),
																												description => $desc,
																												context => Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
																												type => Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
																												subtype => $subtype);

	$anno->addClassification(@classes);

	return ($anno,@groupClasses);
}

sub _getWildTypeStringForCDSAnno {
	my ($self,$var,$tran,$rAnnot) = @_;
	my $wt = $rAnnot->getWt();
	$wt =~ s/u/t/ig;
	return uc($wt);
}

sub _getMutantStringForCDSAnno {
	my ($self,$var,$tran,$rAnnot) = @_;
	my $mt = $rAnnot->getMt();
	$mt =~ s/u/t/ig;
	return uc($mt);
}

sub _getCDSDescriptionString {
	my ($self,$tran,$mutStart,$mutEnd,$mutStartOffset,$mutEndOffset,$wt,$mt) = @_;
	my $desc;
  	if($mutStartOffset == 0){
 		$desc = 'c.'.$mutStart . $wt . '>' . $mt;
 	} elsif($mutStartOffset > 0){
 		$desc = 'c.'.$mutStart.'+'.$mutStartOffset. $wt.'>'.$mt;
 	} else {
 		$desc = 'c.'.$mutStart.$mutStartOffset.$wt.'>'.$mt;
 	}
 	return $desc;
}

sub _getMutatedCdsSequence {
	my ($self,$wtDna,$min,$max,$mt) = @_;
	my $mtDna = substr($wtDna,0,$min - 1) . $mt . substr($wtDna,$max);
	return $mtDna
}

sub _getCdsMinPosForProteinCalculation {
	my ($self,$cAnnot) = @_;
  	return $cAnnot->getMinPos();
}

sub _getCdsMaxPosForProteinCalculation {
	my ($self,$cAnnot) = @_;
  	return $cAnnot->getMaxPos();
}

sub _isStartGained {
 	my ($self,$var,$tran,$varTranStart,$varTranEnd,$wt,$mt) = @_;
 	my $varLength = ($varTranEnd - $varTranStart) + 1;
 	my $substrStart = $varTranStart - 3;
 	my $substrLength;

 	if($substrStart < 0){
 		$substrLength = (($varTranEnd - $varTranStart) + 5) + $substrStart;
 		$substrStart = 0;
 	} else {
 		$substrLength = ($varTranEnd - $varTranStart) + 5;
 	}

 	my $wtSubstr = substr($tran->getcDNASeq,$substrStart,$substrLength);
 	my $mtSubstr = $wtSubstr;
 	$mtSubstr =~ s/\w{$varLength}(\w{2})$/$mt$1/;
 	my @wtPos;
 	my @mtPos;
 	while($wtSubstr =~ m/atg/ig){
 		push(@wtPos,pos($wtSubstr));
 	}
 	while($mtSubstr =~ m/atg/ig){
 		push(@mtPos,pos($mtSubstr));
 	}
 	if(scalar(@mtPos) > scalar(@wtPos)){
 		# more in mutant, must have created a start
 		return 1;
 	} elsif(scalar(@mtPos) == scalar(@wtPos) && scalar(@wtPos) > 0){
 		# both have same number, must check positions, if any change must have created a start
 		for(my $i = 0 ;$i < scalar(@wtPos) ; $i++){
 			if($mtPos[$i] != $wtPos[$i]){
 				return 1;
 			}
 		}
 	}
 	return 0;
}

sub _getDefaultCDSAnnotationType {
	return Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType();
}

sub _getDefaultProteinAnnotationType {
	return Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType();
}


=head1 NAME

Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator - Annotator for simple substitution variants

=head1 DESCRIPTION

This annotates substitution variants, it provides L<AnnotatonGroups|Sanger::CGP::Vagrent::Data::AnnotationGroup> for each transcript returned from the L<TranscriptSource|Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource>

It will only process L<Sanger::CGP::Vagrent::Data::Substitution|Sanger::CGP::Vagrent::Data::Substitution> objects, if any other L<Variation|Sanger::CGP::Vagrent::Data::AbstractVariation> objects are sent in it will return an empty answer.

It inherits from L<Sanger::CGP::Vagrent::Annotators::AbstractAnnotator|Sanger::CGP::Vagrent::Annotators::AbstractAnnotator>.

=head1 METHODS

See L<Sanger::CGP::Vagrent::Annotators::AbstractAnnotator|Sanger::CGP::Vagrent::Annotators::AbstractAnnotator>
