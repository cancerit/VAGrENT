package Sanger::CGP::Vagrent::Annotators::DeletionAnnotator;

use strict;

use Bio::Seq;

use Carp qw(cluck);
use Log::Log4perl;
use POSIX qw(ceil);
use Data::Dumper;

use Sanger::CGP::Vagrent::Data::Deletion;
use Sanger::CGP::Vagrent qw($VERSION);

use base qw(Sanger::CGP::Vagrent::Annotators::AbstractAnnotator);

my $log = Log::Log4perl->get_logger(__PACKAGE__);

1;

sub _getAnnotation {
	my ($self,$var) = @_;
	$self->_clearMessages();
	unless(defined($var) && $var->isa('Sanger::CGP::Vagrent::Data::Deletion')){
		my $msg = 'require a Sanger::CGP::Vagrent::Data::Deletion object not a '.ref($var);
		$self->addMessage($msg);
		$log->info($msg);
		return undef;
	}
	unless($var->isValid){
		my $msg = 'deletion not valid';
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
	my @classes = ($self->getDeletionClass);
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
		if($mrnaMax >= $tran->getCdsMinPos && $mrnaMin <= $tran->getCdsMaxPos){
			# coding change
			push(@groupClasses,$self->getCDSClass);
		}
		if($mrnaMin < $tran->getCdsMinPos){
			# 5prime UTR
			push(@groupClasses,$self->get5PrimeUtrClass);
			push(@classes,$self->get5PrimeUtrVariantClass) unless($self->_arrayHasString($self->getCDSClass,@groupClasses));
		}
		if($mrnaMax > $tran->getCdsMaxPos){
			# 3prime UTR
			push(@groupClasses,$self->get3PrimeUtrClass);
			push(@classes,$self->get3PrimeUtrVariantClass) unless($self->_arrayHasString($self->getCDSClass,@groupClasses));
		}
	}

	my $tmpGroupClassesHash = $self->_classifyDeletion($tran,$mrnaMin,$mrnaMinOffset,$mrnaMax,$mrnaMaxOffset);
	push(@groupClasses,keys %$tmpGroupClassesHash);
	if(scalar(keys %$tmpGroupClassesHash) == 1 && defined($tmpGroupClassesHash->{$self->getIntronClass})){
		# its intron only
		return ($self->_buildUnknownMRNAAnnotation($var,$tran,$self->getDeletionClass,$self->getIntronVariantClass),@groupClasses);
	}

	my $wt = $self->_getWildTypeStringForRNAAnno($var,$tran);
	my $mt = '-';

	$wt =~ tr/Tt/Uu/;

	my $desc = undef;

	if($mrnaMin == 1 && $mrnaMax == length($tran->getcDNASeq)){
		$desc = 'r.0';
	} else {
		$desc = $self->_generateDNADelDescriptionString('r.',$mrnaMin,$mrnaMax,$mrnaMinOffset,$mrnaMaxOffset,$wt);
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
																													mt => $mt,
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
																													type => Sanger::CGP::Vagrent::Data::Annotation::getDeletionAnnotationType(),
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
						if(length($wt) % 3 == 0){
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

sub _upstreamVariantCheck {
	my ($self,$mrnaMin,$mrnaMinOffset,$mrnaMax,$mrnaMaxOffset,$classes) = @_;

	my $upStartScore = 0;
	my $upEndScore = 0;
	if($mrnaMin > 0 || ($mrnaMin == 0 && $mrnaMinOffset >= 0)){
		$upStartScore = 3;
	} elsif($mrnaMin == 0 && $mrnaMinOffset < 0){
		if($self->_isWithin2KBUpstreamOffsetDistance($mrnaMinOffset)){
			$upStartScore = 2;
		} elsif($self->_isWithin5KBOffsetDistance($mrnaMinOffset)){
			$upStartScore = 1;
		}
	}
	if($mrnaMax > 0 || ($mrnaMax == 0 && $mrnaMaxOffset >= 0)){
		$upEndScore = 3;
	} elsif ($mrnaMax == 0 && $mrnaMaxOffset < 0){
		if($self->_isWithin2KBUpstreamOffsetDistance($mrnaMaxOffset)){
			$upEndScore = 2;
		} elsif($self->_isWithin5KBOffsetDistance($mrnaMaxOffset)){
			$upEndScore = 1;
		}
	}
	print "up scores: $upStartScore - $upEndScore\n" if($self->_debug);
	if($upStartScore == 0 && $upEndScore == 0){
		my $msg = "variant isnt close enough to this transcript, nothing to do";
		$self->addMessage($msg);
		$log->error($msg);
		return undef;
	}
	for(my $i = $upStartScore; $i <= $upEndScore ; $i++){
		if($i == 1){
			push(@$classes,$self->get5KBUpStreamVariantClass);
		} elsif($i == 2){
			push(@$classes,$self->get2KBUpStreamVariantClass);
		}
	}
	if($upEndScore < 3){
		return 1
	}
	return 0;
}

sub _downstreamVariantCheck {
	my ($self,$mrnaMin,$mrnaMinOffset,$mrnaMax,$mrnaMaxOffset,$classes) = @_;

	my $downStartScore = 0;
	my $downEndScore = 0;
	if($mrnaMin > 0 || ($mrnaMin == 0 && $mrnaMinOffset <= 0)){
		$downStartScore = 3;
	} elsif($mrnaMin == 0 && $mrnaMinOffset > 0){
		if($self->_isWithin500BPDownstreamOffsetDistance($mrnaMinOffset)){
			$downStartScore = 2;
		} elsif($self->_isWithin5KBOffsetDistance($mrnaMinOffset)){
			$downStartScore = 1;
		}
	}
	if($mrnaMax > 0 || ($mrnaMax == 0 && $mrnaMaxOffset <= 0)){
		$downEndScore = 3;
	} elsif ($mrnaMax == 0 && $mrnaMaxOffset > 0){
		if($self->_isWithin500BPDownstreamOffsetDistance($mrnaMaxOffset)){
			$downEndScore = 2;
		} elsif($self->_isWithin5KBOffsetDistance($mrnaMaxOffset)){
			$downEndScore = 1;
		}
	}
	print "down scores: $downStartScore - $downEndScore\n" if($self->_debug);
	if($downStartScore == 0 && $downEndScore == 0){
		my $msg = "variant isnt close enough to this transcript, nothing to do";
		$self->addMessage($msg);
		$log->error($msg);
		return undef;
	}
	for(my $i = $downEndScore ; $i <= $downStartScore; $i++){
		if($i == 1){
			push(@$classes,$self->get5KBDownStreamVariantClass);
		} elsif($i == 2){
			push(@$classes,$self->get500BPDownStreamVariantClass);
		}
	}
	if($downStartScore < 3){
		return 1;
	}
	return 0;
}

sub _classifyDeletion {
	my ($self,$tran,$mrnaMin,$mrnaMinOffset,$mrnaMax,$mrnaMaxOffset) = @_;
	# big scary logic ladder time
	my $grpCls;
	my @exons = $tran->getExons;
	@exons = sort {$a->getRnaMinPos <=> $b->getRnaMinPos} @exons;
	print "$mrnaMin,$mrnaMinOffset,$mrnaMax,$mrnaMaxOffset\n" if($self->_debug);
	for(my $i = 0 ; $i < scalar(@exons) ; $i++){
		my $e = $exons[$i];
		if($mrnaMin > $e->getRnaMaxPos){
			# exon occurs before mutation in transcript, ignore
			next;
		}
		if($mrnaMax < $e->getRnaMinPos){
			# exon occurs after mutation in transcript, ignore
			next;
		}
		# if we get here the variant is somehow linked to the current exon

		# check to see if the variant overlaps the exon
		if(($mrnaMin < $e->getRnaMaxPos || ($mrnaMin == $e->getRnaMaxPos && $mrnaMinOffset == 0)) && ($mrnaMax > $e->getRnaMinPos || ($mrnaMax == $e->getRnaMinPos && $mrnaMaxOffset == 0))){
			$grpCls->{$self->getExonClass}++;
		}
		print Dumper($e) if($self->_debug);
		# check to see if the variant effects the area before the exon
		if($mrnaMin < $e->getRnaMinPos || ($mrnaMin == $e->getRnaMinPos && $mrnaMinOffset < 0)){
			# variant starts before the exon
			if($mrnaMin < $e->getRnaMinPos || $self->_isIntronicOffsetDistance($mrnaMinOffset)){
				# variant starts in the intron
				if($mrnaMax > $e->getRnaMinPos || ($mrnaMax == $e->getRnaMinPos && $mrnaMaxOffset == 0)){
					# variant ends after or on the start of the exon, so it gets all 3 terms
					$grpCls->{$self->getIntronClass}++;
					$grpCls->{$self->getSpliceRegionClass}++;
					$grpCls->{$self->getEssentialSpliceSiteClass}++;
				} elsif($self->_isIntronicOffsetDistance($mrnaMaxOffset)){
					# variant also ends in the intron
					$grpCls->{$self->getIntronClass}++;
				} elsif($mrnaMaxOffset < 0 && $mrnaMaxOffset >= $self->_getConsesnsusSpliceBeforeBoundry){
					# variant ends in the essential splice site
					$grpCls->{$self->getIntronClass}++;
					$grpCls->{$self->getSpliceRegionClass}++;
					$grpCls->{$self->getEssentialSpliceSiteClass}++;
				} else {
					# variant ends in the splice region
					$grpCls->{$self->getIntronClass}++;
					$grpCls->{$self->getSpliceRegionClass}++;
				}
			} elsif ($mrnaMin == $e->getRnaMinPos && $mrnaMinOffset >= $self->_getConsesnsusSpliceBeforeBoundry){
				# variant starts in the essential splice site
				$grpCls->{$self->getEssentialSpliceSiteClass}++;
			} else {
				# variant starts in the splice region
				if($mrnaMax > $e->getRnaMinPos || ($mrnaMax == $e->getRnaMinPos && $mrnaMaxOffset == 0)){
					# variant ends after the start of the exon, so
					$grpCls->{$self->getSpliceRegionClass}++;
					$grpCls->{$self->getEssentialSpliceSiteClass}++;
				} elsif($mrnaMaxOffset < 0 && $mrnaMaxOffset >= $self->_getConsesnsusSpliceBeforeBoundry){
					# variant ends in the essential splice site
					$grpCls->{$self->getSpliceRegionClass}++;
					$grpCls->{$self->getEssentialSpliceSiteClass}++;
				} else {
					# variant also ends in the splice region
					$grpCls->{$self->getSpliceRegionClass}++;
				}
			}
		}

		# check to see if the variant effects the area after the exon
		if($mrnaMax > $e->getRnaMaxPos || ($mrnaMax == $e->getRnaMaxPos && $mrnaMaxOffset > 0)){
			# variant ends after the exon
			if($mrnaMax > $e->getRnaMaxPos || $self->_isIntronicOffsetDistance($mrnaMaxOffset)){
				# variant ends in the intron
				if($mrnaMin < $e->getRnaMaxPos || ($mrnaMin == $e->getRnaMaxPos && $mrnaMinOffset == 0)){
					# variant starts before or on the end of the exon, so it gets all 3 terms
					$grpCls->{$self->getIntronClass}++;
					$grpCls->{$self->getSpliceRegionClass}++;
					$grpCls->{$self->getEssentialSpliceSiteClass}++;
				} elsif($self->_isIntronicOffsetDistance($mrnaMinOffset)){
					# variant also starts in the intron
					$grpCls->{$self->getIntronClass}++;
				} elsif ($mrnaMinOffset > 0 && $mrnaMinOffset <= $self->_getConsesnsusSpliceAfterBoundry){
					# variant starts in the essential splice site
					$grpCls->{$self->getIntronClass}++;
					$grpCls->{$self->getSpliceRegionClass}++;
					$grpCls->{$self->getEssentialSpliceSiteClass}++;
				} else {
					# variant must start in the splice region
					$grpCls->{$self->getIntronClass}++;
					$grpCls->{$self->getSpliceRegionClass}++;
				}
			} elsif($mrnaMax == $e->getRnaMaxPos && $mrnaMaxOffset <= $self->_getConsesnsusSpliceAfterBoundry){
				# variant ends in the essential splice site
				$grpCls->{$self->getEssentialSpliceSiteClass}++;
			} else {
				# variant ends in the splice region
				if($mrnaMin < $e->getRnaMaxPos || ($mrnaMin == $e->getRnaMaxPos && $mrnaMinOffset == 0)){
					# variant starts before or on the end of the exon, so it gets all 3 terms
					$grpCls->{$self->getSpliceRegionClass}++;
					$grpCls->{$self->getEssentialSpliceSiteClass}++;
				} elsif ($mrnaMinOffset > 0 && $mrnaMinOffset <= $self->_getConsesnsusSpliceAfterBoundry){
					# variant starts in the essential splice site
					$grpCls->{$self->getSpliceRegionClass}++;
					$grpCls->{$self->getEssentialSpliceSiteClass}++;
				} else {
					# variant must start in the splice region
					$grpCls->{$self->getSpliceRegionClass}++;
				}
			}
		}
	}
	print Dumper($grpCls) if($self->_debug);
	return $grpCls;
}

sub _getCDSDescriptionString {
	my ($self,$tran,$mutStart,$mutEnd,$mutStartOffset,$mutEndOffset,$wt,$mt) = @_;
  	my $desc = undef;
 	if($mutStart == 1 && $mutEnd == length($tran->getCdsSeq)){
 		$desc = 'c.0';
 	} else {
 		$desc = $self->_generateDNADelDescriptionString('c.',$mutStart,$mutEnd,$mutStartOffset,$mutEndOffset,$wt);
 	}
  	return $desc;
}

sub _generateDNADelDescriptionString {
	my ($self,$pre,$mutStart,$mutEnd,$mutStartOffset,$mutEndOffset,$wt) = @_;
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
	}
	return $desc;
}

sub _getWildTypeStringForRNAAnno {
	my ($self,$var,$tran) = @_;
	my $wt = '';
	my @seq = split('',$var->getDeletedSequence);
	my @exons = sort{$a->getMinPos <=> $b->getMinPos} $tran->getExons;
	for(my $i = $var->getMinPos, my $c = 0; $c < scalar(@seq); $i++, $c++){
		my $isExonic = 0;
		my $ec = 0;
		my $use = 0;
		foreach my $e(@exons){
			$ec++;
			if($i > $e->getMaxPos){ # base after exon, continue
				if($ec == scalar(@exons)){ # if the base is after the last exon we don't want to use it
					$use = 0;
				} else {
					$use = 1;
				}
				next;
			} elsif($i >= $e->getMinPos && $i <= $e->getMaxPos){ # base in exon, quit
				$isExonic = 1;
				$use = 1;
				last;
			} elsif($i < $e->getMinPos){ # base before exon, quit
				$use = 1 unless($ec==1); # ignore the base if its before the first exon
				last;
			}
		}
		next unless($use == 1); # skip the base if its not flagged for use
		if($isExonic){
			$wt .= uc($seq[$c]);
		} else {
			$wt .= lc($seq[$c]);
		}
	}
	if($tran->getStrand == -1){
		$wt = $self->_revcompSeq($wt);
	}
	return $wt;
}

sub _getWildTypeStringForCDSAnno {
	my ($self,$var,$tran,$rAnnot) = @_;
	my $wt = '';
	if($rAnnot->hasClassification($self->getInFrameVariantClass) || $rAnnot->hasClassification($self->getFrameShiftVariantClass)){
		my @seq = split('',$var->getDeletedSequence);
		my @exons = sort{$a->getMinPos <=> $b->getMinPos} $tran->getExons;
		for(my $i = $var->getMinPos, my $c = 0; $c < scalar(@seq); $i++, $c++){
			my $isCod = 0;
			foreach my $e(@exons){
				next unless($i >= $e->getMinPos && $i <= $e->getMaxPos);  # base is in this exon
				if($tran->getStrand == 1){
					my $pos = ($i - $e->getMinPos) + $e->getRnaMinPos;
					if($pos >= $tran->getCdsMinPos && $pos <= $tran->getCdsMaxPos){
						$isCod = 1;
						last;
					}
				} else {
					my $pos = $e->getRnaMaxPos - ($i - $e->getMinPos);
					if($pos >= $tran->getCdsMinPos && $pos <= $tran->getCdsMaxPos){
						$isCod = 1;
						last;
					}
				}
			}
			if($isCod){
				$wt .= uc($seq[$c]);
			}
		}
		if($tran->getStrand == -1){
			$wt = $self->_revcompSeq($wt);
		}
	} elsif ($rAnnot->hasClassification($self->getComplexChangeVariantClass)){
		if($rAnnot->getMinPos <= $tran->getCdsMinPos || $rAnnot->getMaxPos >= $tran->getCdsMaxPos){
			# variant exists outside the CDS, have to truncate it back
			my $len = ($tran->getCdsMaxPos - $tran->getCdsMinPos) + 1;
			my ($genomicCdsLow,$genomicCdsHigh) = $self->_calculateGenomicCdsPosition($tran);
			my $newDelMin = undef;
			my $newDelMax = undef;
			my $newDelSeq = $var->getDeletedSequence;
			if($var->getMinPos < $genomicCdsLow){
				$newDelMin = $genomicCdsLow;
				substr($newDelSeq,0,($newDelMin - $var->getMinPos),'');
			} else {
				$newDelMin = $var->getMinPos;
			}
			if($var->getMaxPos > $genomicCdsHigh){
				$newDelMax = $genomicCdsHigh;
				$newDelSeq = substr($newDelSeq,0,($newDelMax - $newDelMin) + 1);
			} else {
				$newDelMax = $var->getMaxPos;
			}

			my $newDel = Sanger::CGP::Vagrent::Data::Deletion->new(
										'species'				=> $var->getSpecies,
										'genomeVersion' => $var->getGenomeVersion,
										'chr' 					=> $var->getChr,
										'minpos'				=> $newDelMin,
										'maxpos'				=> $newDelMax,
										'delseq' 				=> $newDelSeq);

			$wt = $self->_getWildTypeStringForRNAAnno($newDel,$tran);
		} else {
			$wt = $self->_getWildTypeStringForRNAAnno($var,$tran);
		}
	} else {
		$wt = $self->_getWildTypeStringForRNAAnno($var,$tran);
	}
	$wt =~ s/u/t/ig;
	return $wt;
}

sub _getMutantStringForCDSAnno {
	my ($self,$var,$tran,$rAnnot) = @_;
	return '-';
}

sub _getMutatedCdsSequence {
	my ($self,$wtDna,$min,$max,$mt) = @_;
 	my $codingDelLength = ($max - $min) + 1;
 	my $mtDna = $wtDna;
 	substr($mtDna,($min - 1),$codingDelLength,'');
	return $mtDna
}

sub _getCdsMinPosForProteinCalculation {
	my ($self,$cAnnot) = @_;
	my $min = undef;
	if($cAnnot->getMinOffset > 0){
 		$min = $cAnnot->getMinPos + 1;
 	} else {
 		$min = $cAnnot->getMinPos;
 	}
 	return $min;
}

sub _getCdsMaxPosForProteinCalculation {
	my ($self,$cAnnot) = @_;
  	my $max = undef;
  	if($cAnnot->getMaxOffset < 0){
 		$max = $cAnnot->getMaxPos - 1;
 	} else {
 		$max = $cAnnot->getMaxPos;
 	}
  	return $max;
}

sub _safetyCheck {
	my ($self,$var,$tran,$mrnaMin,$mrnaMinOffset,$mrnaMax,$mrnaMaxOffset) = @_;
	my $tranString = $self->_getTranscriptSubstringForCheck($tran,$mrnaMin,$mrnaMinOffset,$mrnaMax,$mrnaMaxOffset);
	unless(defined($tranString)){
		return 0;
	}
	my $genoString = $self->_getExonicGenomicSubstringForCheck($var,$tran);
	unless(defined($genoString)){
		return 0;
	}
	if($self->_debug){
		print 'DELSEQ = '.$var->getDeletedSequence."\n";
		print 'TRANSEQ = '.$tranString."\n";
		print 'GENOSEQ = '.$genoString."\n";
	}
	if(length($tranString) == 0 && length($genoString) == 0){
		# not exonic sequence, no transcript sequence to check, just have to hope its OK
		return 1;
	} elsif(uc($tranString) eq uc($genoString)){
		return 1;
	}
	return 0;
}

sub _getTranscriptSubstringForCheck {
	my ($self,$tran,$mrnaMin,$mrnaMinOffset,$mrnaMax,$mrnaMaxOffset) = @_;
	my $substr = undef;
	my $pos = undef;
	my $len = undef;
	print "TRAN SUB STR $mrnaMin,$mrnaMinOffset,$mrnaMax,$mrnaMaxOffset\n" if($self->_debug);
	if($mrnaMin == 0 && $mrnaMinOffset > 0){
		# variant starts after transcript, nothing to do
		return '';
	} elsif($mrnaMax == 0 && $mrnaMaxOffset < 0){
		# variant ends before transcript, nothing to do
		return '';
	} elsif($mrnaMinOffset == 0 && $mrnaMaxOffset == 0){
		$pos = $mrnaMin - 1;
		$len = ($mrnaMax - $mrnaMin) + 1
	} else {
		my $leftAnchor;
		my $rightAnchor;
		if($mrnaMinOffset > 0){
			$leftAnchor = $mrnaMin + 1;
		} else {
			$leftAnchor = $mrnaMin;
		}
		if($mrnaMaxOffset < 0){
			$rightAnchor = $mrnaMax - 1;
		} else {
			$rightAnchor = $mrnaMax;
		}
		$leftAnchor = 1 if $leftAnchor < 1;
		$rightAnchor = length($tran->getcDNASeq) if($rightAnchor == 0);
		print "ANCHOR $leftAnchor $rightAnchor\n" if($self->_debug);
		if($leftAnchor > $rightAnchor){
			# its an intronic change that doesn't span any coding bases
			return '';
		}
		$pos = $leftAnchor - 1;
		$len = ($rightAnchor - $leftAnchor) + 1;

	}
	print "FINAL $pos $len\n" if($self->_debug);
	$substr = substr($tran->getcDNASeq,$pos,$len);
	return $substr;
}

sub _getExonicGenomicSubstringForCheck {
	my ($self,$var,$tran) = @_;
	my $substr = undef;
	my $c = 0;
	my $beforeC = 0;
	my $duringC = 0;
	foreach my $e ($tran->getExonsGenomicOrder){
		$c++;
		if($e->getMaxPos < $var->getMinPos){
			# exon before deletion - skip
			print "EXON $c BEFORE\n" if($self->_debug);
			$beforeC++;
			next;
		}
		if($var->getMinPos <= $e->getMaxPos && $e->getMinPos <= $var->getMaxPos){
			# overlap
			$duringC++;
			my $start = ($e->getMinPos - $var->getMinPos) + 1;
			my $end = ($e->getMaxPos - $var->getMinPos) + 1;
			$start = 1 if($start < 1);
			$end = length($var->getDeletedSequence) if($end > length($var->getDeletedSequence));

			my $pos = $start - 1;
			my $len = ($end - $start) + 1;
			my $tmp = substr($var->getDeletedSequence,$pos,$len);
			if(defined($substr)){
				$substr .= $tmp;
			} else {
				$substr = $tmp;
			}
			print "EXON $c DURING $start $end $tmp\n" if($self->_debug);
		}
		if($var->getMaxPos < $e->getMinPos){
			print "EXON $c AFTER\n" if($self->_debug);
			# exon after deletion
			if(!defined($substr)){
				# its not a coding change;
				return '';
			}
		}
	}
	if($duringC == 0){
 		# variant doesn't overlap any exons
 		return '';
 	}
	if($tran->getStrand == -1){
		$substr = $self->_revcompSeq($substr);
	}
	return $substr;
}

sub _getDefaultCDSAnnotationType {
	return Sanger::CGP::Vagrent::Data::Annotation::getDeletionAnnotationType();
}

sub _getDefaultProteinAnnotationType {
	return Sanger::CGP::Vagrent::Data::Annotation::getDeletionAnnotationType();
}


=head1 NAME

Sanger::CGP::Vagrent::Annotators::DeletionAnnotator - Annotator for deletion variants

=head1 DESCRIPTION

This annotates deletion variants, it provides L<AnnotatonGroups|Sanger::CGP::Vagrent::Data::AnnotationGroup> for each transcript returned from the L<TranscriptSource|Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource>

It will only process L<Sanger::CGP::Vagrent::Data::Deletion|Sanger::CGP::Vagrent::Data::Deletion> objects, if any other L<Variation|Sanger::CGP::Vagrent::Data::AbstractVariation> objects are sent in it will return an empty answer.

It inherits from L<Sanger::CGP::Vagrent::Annotators::AbstractAnnotator|Sanger::CGP::Vagrent::Annotators::AbstractAnnotator>.

=head1 METHODS

See L<Sanger::CGP::Vagrent::Annotators::AbstractAnnotator|Sanger::CGP::Vagrent::Annotators::AbstractAnnotator>
