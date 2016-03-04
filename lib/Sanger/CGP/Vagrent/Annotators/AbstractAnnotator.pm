package Sanger::CGP::Vagrent::Annotators::AbstractAnnotator;

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
use POSIX qw(ceil);
use Data::Dumper;
use Attribute::Abstract;
use Const::Fast qw(const);

use Sanger::CGP::Vagrent qw($VERSION);

use Sanger::CGP::Vagrent::Data::Annotation;
use Sanger::CGP::Vagrent::Data::AnnotationGroup;

use base qw(Sanger::CGP::Vagrent::Ontology::SequenceOntologyClassifier);

my $log = Log::Log4perl->get_logger(__PACKAGE__);

# constant reference values for consensus splice site values.
const my @CONSENSUS_SPLICE_OFFSETS => (-2, -1, 1, 2, 5);
const my $CONSENSUS_SPLICE_BEFORE_BOUNDRY => -2;
const my $CONSENSUS_SPLICE_AFTER_BOUNDRY => 5;

# constant value representing the cutoff for intronic calls
const my $INTRONIC_OFFSET_CUTOFF => 11;

const my $UPDOWNSTREAM_5KB_CUTOFF => 5000;
const my $UPSTREAM_2KB_CUTOFF => -2000;
const my $DOWNSTREAM_500BP_CUTOFF => 500;

1;

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
	my $self = {};
	bless($self, $class);
	$self->_init(@_);
  $self->_ontologyInit(@_);
	return $self;
}

sub getConsensusSpliceOffsets {
  return @CONSENSUS_SPLICE_OFFSETS;
}

sub getAnnotation {
	my ($self,$var) = @_;
	my @ann = $self->_getAnnotation($var);
	$self->_assignBookmarks(@ann);
	if(exists $self->{'_only_bookmarked'} && $self->{'_only_bookmarked'}){
		my @out;
		foreach my $a(@ann){
			if($a->hasBookmark()){
				push(@out,$a);
			}
		}
		return @out;
	}
	return @ann;
}

sub _getAnnotation: Abstract;

sub addMessage {
	my ($self,$msg) = @_;
	push(@{$self->{_msg}},ref($self).": ".$msg);
}

sub _debug {
	my $self = shift;
	if(exists($self->{_debug}) && defined($self->{_debug}) && $self->{_debug}){
		return 1;
	} else {
		return 0;
	}
}

sub getMessages {
	my $self = shift;
	return @{$self->{_msg}} if defined($self->{_msg});
	return undef;
}

sub _clearMessages {
	shift->{_msg} = undef;
}

sub _assignBookmarks {
	my ($self,@groups) = @_;
	if(exists $self->{'_bookmarkers'} && defined $self->{'_bookmarkers'}){
		foreach my $bm (@{$self->{'_bookmarkers'}}){
			$bm->markAnnotation(@groups);
		}
	}
}

sub _isStartGained: Abstract;

sub _init {
	my $self = shift;
	my %vars = @_;
	foreach my $k(keys(%vars)){
		if($k eq 'transcriptSource' && $vars{'transcriptSource'}->isa('Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource')){
			$self->{'_transcriptSource'} = $vars{'transcriptSource'};
		} elsif($k eq 'bookmarker'){
			if(ref($vars{'bookmarker'}) eq 'ARRAY'){
				foreach my $tmp(@{$vars{'bookmarker'}}){
					if($tmp->isa('Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker')){
						push(@{$self->{'_bookmarkers'}},$tmp);
					}
				}
			} elsif($vars{'bookmarker'}->isa('Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker')){
				push(@{$self->{'_bookmarkers'}},$vars{'bookmarker'});
			}
		} elsif($k eq 'only_bookmarked'){
			$self->{'_only_bookmarked'} = $vars{'only_bookmarked'};
		} elsif($k eq 'debug' && $vars{'debug'}){
			$self->{'_debug'} = 1;
		}
	}
}

sub _getmRNAPositions {
	my ($self,$mut,$tran) = @_;
	my $exonBeforeInTrans = undef;
	my $exonContainingStartInTrans = undef;
	my $exonsDuringInTrans = undef;
	my $exonContainingEndInTrans = undef;
	my $exonAfterInTrans = undef;
	foreach my $e ($tran->getExons){
		if($e->getChr() ne $mut->getChr()){
			# thats not right, chromosomes don't match;
			my $msg = 'Unable to process this variant, chromosome doesnt match transcript';
			$self->addMessage($msg);
			$log->error($msg);
			return undef;
		}
		if($e->getMaxPos < $mut->getMinPos){
			# exon before mutation on genome, they don't overlap
			if($tran->getStrand == 1){
				# we are the easy way up
				$exonBeforeInTrans = $e;
			} else {
				# we are the tricky way up
				unless(defined($exonAfterInTrans)){
					$exonAfterInTrans = $e;
				}
			}
			next;
		}
		if($mut->getMinPos <= $e->getMinPos && $e->getMaxPos <= $mut->getMaxPos){
			# exon is completely covered by the mutation
			push(@$exonsDuringInTrans,$e);
			next;
		}

		if($mut->getMaxPos < $e->getMinPos){
			# exon after mutation on genome, they don't overlap
			if($tran->getStrand == 1){
				# we are the easy way up
				unless(defined($exonAfterInTrans)){
					$exonAfterInTrans = $e;
				}
			} else {
				# we are the tricky way up
				$exonBeforeInTrans = $e;
			}
			next;
		}
		if($e->getMinPos <= $mut->getMinPos && $mut->getMinPos <= $e->getMaxPos ){
			# start of mutation is in this exon on the genome
			if($tran->getStrand == 1){
				# we are the easy way up
				$exonContainingStartInTrans = $e;
			} else {
				# we are the tricky way up
				$exonContainingEndInTrans = $e;
			}
		}
		if($e->getMinPos <= $mut->getMaxPos && $mut->getMaxPos <= $e->getMaxPos){
			# end of mutation is in this exon on the genome
			if($tran->getStrand == 1){
				# we are the easy way up
				$exonContainingEndInTrans = $e;
			} else {
				# we are the tricky way up
				$exonContainingStartInTrans = $e;
			}
		}
	}

	if($self->_debug){
		print 'VAR '.Dumper($mut);
		print 'BEFORE '.Dumper($exonBeforeInTrans);
		print 'START '.Dumper($exonContainingStartInTrans);
		print 'DURING '.Dumper($exonsDuringInTrans);
		print 'END '.Dumper($exonContainingEndInTrans);
		print 'AFTER '.Dumper($exonAfterInTrans);
	}

	my $mutTranStart = undef;
	my $mutStartOffset = undef;
	my $mutTranEnd = undef;
	my $mutEndOffset = undef;

	if(!defined($exonBeforeInTrans) &&
		 !defined($exonContainingStartInTrans)){
		# variant starts before the transcript, set start to special case 0
		my $tmpAfter;
		if(defined($exonsDuringInTrans) && defined($exonsDuringInTrans->[0])){
			$tmpAfter = $exonsDuringInTrans->[0];
		} elsif(defined($exonContainingEndInTrans)){
			$tmpAfter = $exonContainingEndInTrans;
		} else {
			$tmpAfter = $exonAfterInTrans;
		}
		$mutTranStart = 0;
		if($tran->getStrand == 1){
			$mutStartOffset = $mut->getMinPos - $tmpAfter->getMinPos;
		} else {
			$mutStartOffset = $tmpAfter->getMaxPos - $mut->getMaxPos;
		}
	} elsif (defined($exonBeforeInTrans) &&
		       !defined($exonContainingStartInTrans) &&
		       !defined($exonsDuringInTrans) &&
		       !defined($exonContainingEndInTrans) &&
		       !defined($exonAfterInTrans)) {
		# variant starts after the transcript, set start to special case 0
		$mutTranStart = 0;
		if($tran->getStrand == 1){
			$mutStartOffset = $mut->getMinPos - $exonBeforeInTrans->getMaxPos;
		} else {
			$mutStartOffset = $exonBeforeInTrans->getMinPos - $mut->getMaxPos;
		}
	} elsif(defined($exonContainingStartInTrans)){
		# mutations starts within an exon;
		if($tran->getStrand == 1){
			# we are the easy way up
			$mutTranStart = ($mut->getMinPos - $exonContainingStartInTrans->getMinPos) + $exonContainingStartInTrans->getRnaMinPos;
			$mutStartOffset = 0;
		} else {
			# we are the tricky way up
			$mutTranStart = ($exonContainingStartInTrans->getMaxPos - $mut->getMaxPos) + $exonContainingStartInTrans->getRnaMinPos;
			$mutStartOffset = 0;
		}
		print "START - $mutTranStart - $mutStartOffset\n" if($self->_debug);
	} elsif(defined($exonBeforeInTrans)) {
		# mutation starts between exons, this makes things harder
		my $before = $exonBeforeInTrans;
		my $after = undef;
		my $beforeBoundry = undef;
		my $afterBoundry = undef;
		if(defined($exonsDuringInTrans)){
			$after = $exonsDuringInTrans->[0];
		} elsif(defined($exonContainingEndInTrans)){
			$after = $exonContainingEndInTrans;
		} else {
			$after = $exonAfterInTrans;
		}
		if($tran->getStrand == 1){
			$beforeBoundry->{mut} = $mut->getMinPos;
			$beforeBoundry->{exon} = $before->getMaxPos;
			$beforeBoundry->{off} = $mut->getMinPos - $before->getMaxPos;
			$beforeBoundry->{pos} = $before->getRnaMaxPos;
			$afterBoundry->{mut} = $mut->getMinPos;
			$afterBoundry->{exon} = $after->getMinPos;
			$afterBoundry->{off} =  $mut->getMinPos - $after->getMinPos;
			$afterBoundry->{pos} = $after->getRnaMinPos;
		} else {
			$beforeBoundry->{mut} = $mut->getMaxPos;
			$beforeBoundry->{exon} = $before->getMinPos;
			$beforeBoundry->{off} = $before->getMinPos - $mut->getMaxPos;
			$beforeBoundry->{pos} = $before->getRnaMaxPos;
			$afterBoundry->{mut} = $mut->getMaxPos;
			$afterBoundry->{exon} = $after->getMaxPos;
			$afterBoundry->{off} = $after->getMaxPos - $mut->getMaxPos;
			$afterBoundry->{pos} = $after->getRnaMinPos;
		}
		if($self->_debug){
			print 'START BEFORE BOUNDRY '. Dumper($beforeBoundry);
			print 'START AFTER BOUNDRY '.Dumper($afterBoundry);
		}
		if($mut->isa('Sanger::CGP::Vagrent::Data::Insertion') && $beforeBoundry->{off} + $afterBoundry->{off} == 0 && $beforeBoundry->{pos} + 1 == $afterBoundry->{pos}){
			# very special case.  we have an insertion starting at the mid point of an odd length intron.  To avoid confusion
			# we'll attach the start of the insertion to the exon after the insertion as that is going to be closer to the insertion end point
			$mutTranStart = $afterBoundry->{pos};
			$mutStartOffset = $afterBoundry->{off};
		} else {
			if(abs($beforeBoundry->{off}) > abs($afterBoundry->{off})){
				$mutTranStart = $afterBoundry->{pos};
				$mutStartOffset = $afterBoundry->{off};
			} else {
				$mutTranStart = $beforeBoundry->{pos};
				$mutStartOffset = $beforeBoundry->{off};
			}
		}
	} else {
		my $msg = 'Unable to process this variant, very strange coordinate data';
		$self->addMessage($msg);
		$log->error($msg);
		return undef;
	}
	if(!defined($exonBeforeInTrans) &&
 			!defined($exonContainingStartInTrans) &&
 			!defined($exonsDuringInTrans) &&
 			!defined($exonContainingEndInTrans) &&
 			defined($exonAfterInTrans)){
		# variant ends before the transcripts starts, set end to special case 0
		$mutTranEnd = 0;
		if($tran->getStrand == 1){
			$mutEndOffset = $mut->getMaxPos - $exonAfterInTrans->getMinPos;
		} else {
			$mutEndOffset = $exonAfterInTrans->getMaxPos - $mut->getMinPos;
		}
	} elsif(!defined($exonContainingEndInTrans) &&
 			    !defined($exonAfterInTrans)){
		# variant ends after the transcripts ends, set end to special case 0
		my $tmpBefore;
		if(defined($exonsDuringInTrans) && defined($exonsDuringInTrans->[scalar(@$exonsDuringInTrans) - 1])){
			$tmpBefore = $exonsDuringInTrans->[scalar(@$exonsDuringInTrans) - 1];
		} elsif(defined($exonContainingStartInTrans)){
			$tmpBefore = $exonContainingStartInTrans;
		} else {
			$tmpBefore = $exonBeforeInTrans;
		}
		$mutTranEnd = 0;
		if($tran->getStrand == 1){
			$mutEndOffset = $mut->getMaxPos - $tmpBefore->getMaxPos;
		} else {
			$mutEndOffset = $tmpBefore->getMinPos - $mut->getMinPos;
		}
	}	elsif(defined($exonContainingEndInTrans)){
		# mutations ends within an exon;
		if($tran->getStrand == 1){
			# we are the easy way up
			$mutTranEnd = ($mut->getMaxPos - $exonContainingEndInTrans->getMinPos) + $exonContainingEndInTrans->getRnaMinPos;
			$mutEndOffset = 0;
		} else {
			# we are the tricky way up
			$mutTranEnd = ($exonContainingEndInTrans->getMaxPos - $mut->getMinPos) + $exonContainingEndInTrans->getRnaMinPos;
			$mutEndOffset = 0;
		}
		print "END - $mutTranEnd - $mutEndOffset\n" if($self->_debug);
	} elsif(defined($exonAfterInTrans)) {
		# mutation ends between exons, this makes things harder
		my $before = undef;
		my $after = $exonAfterInTrans;
		my $beforeBoundry = undef;
		my $afterBoundry = undef;
		if(defined($exonsDuringInTrans)){
			$before = $exonsDuringInTrans->[scalar(@$exonsDuringInTrans) - 1];
		} elsif(defined($exonContainingStartInTrans)){
			$before = $exonContainingStartInTrans;
		} else {
			$before = $exonBeforeInTrans;
		}
		if($tran->getStrand == 1){
			$beforeBoundry->{mut} = $mut->getMaxPos;
			$beforeBoundry->{exon} = $before->getMaxPos;
			$beforeBoundry->{off} = $mut->getMaxPos - $before->getMaxPos;
			$beforeBoundry->{pos} = $before->getRnaMaxPos;
			$afterBoundry->{mut} = $mut->getMaxPos;
			$afterBoundry->{exon} = $after->getMinPos;
			$afterBoundry->{off} =  $mut->getMaxPos - $after->getMinPos;
			$afterBoundry->{pos} = $after->getRnaMinPos;
		} else {
			$beforeBoundry->{mut} = $mut->getMinPos;
			$beforeBoundry->{exon} = $before->getMinPos;
			$beforeBoundry->{off} = $before->getMinPos - $mut->getMinPos;
			$beforeBoundry->{pos} = $before->getRnaMaxPos;
			$afterBoundry->{mut} = $mut->getMinPos;
			$afterBoundry->{exon} = $after->getMaxPos;
			$afterBoundry->{off} = $after->getMaxPos - $mut->getMinPos;
			$afterBoundry->{pos} = $after->getRnaMinPos;
		}
		if($self->_debug){
			print 'END BEFORE BOUNDRY '. Dumper($beforeBoundry);
			print 'END AFTER BOUNDRY '.Dumper($afterBoundry);
		}
		if($mut->isa('Sanger::CGP::Vagrent::Data::Insertion') && $beforeBoundry->{off} + $afterBoundry->{off} == 1 && $beforeBoundry->{pos} + 1 == $afterBoundry->{pos}){
			# very special case.  we have an insertion dead centre of an evenly sized intron.  To avoid confusion
			# we'll attach the end of the insertion to the exon before the insertion as that is the same as the start
			$mutTranEnd = $beforeBoundry->{pos};
			$mutEndOffset = $beforeBoundry->{off};
		} else {
			if(abs($beforeBoundry->{off}) > abs($afterBoundry->{off})){
				$mutTranEnd = $afterBoundry->{pos};
				$mutEndOffset = $afterBoundry->{off};
			} else {
				$mutTranEnd = $beforeBoundry->{pos};
				$mutEndOffset = $beforeBoundry->{off};
			}
		}
	} else {
		my $msg = 'Unable to process this variant, very strange coordinate data';
		$self->addMessage($msg);
		$log->error($msg);
		return undef;
	}
	print "mRNA VAR START = $mutTranStart $mutStartOffset\nmRNA VAR END = $mutTranEnd $mutEndOffset\n" if($self->_debug);
	return ($mutTranStart,$mutStartOffset,$mutTranEnd,$mutEndOffset);
}

sub _buildProteinAnnotation {
	my ($self,$var,$tran,$cAnnot,$rAnnot) = @_;
	my @classes;
 	my $wtDna = $tran->getCdsSeq;
 	my ($prePad,$postPad) = $self->_calculateCdsTranslationPadStrings($tran);
 	my $wtProt = Bio::Seq->new(-seq => $prePad . $wtDna . $postPad)->translate->seq(); # wild type protein sequence
	unless($self->_canAnnotateToProtein($tran,$cAnnot)){
		if($cAnnot->hasClassification($self->getComplexChangeVariantClass) && $cAnnot->getDescription eq 'c.0'){
			# special case, the CDS has completely gone so the protein must follow.
			return $self->_buildRemovedProteinAnnotation($var,$tran,$wtProt,($self->getDeletionClass));
		}
		my $msg = 'CDS annotation is non-coding or too complex, cant annotate this';
		$self->addMessage($msg);
		$log->info($msg);
		my $pAnnot = $self->_buildUnknownProteinAnnotation($var,$tran,$cAnnot,length($wtProt),@classes);
		my $coversStart = 0;
		my $coversEnd = 0;

		$coversStart = 1 if($self->_coversStartCodon($tran,$cAnnot));
		$coversEnd = 1 if($self->_coversStopCodon($tran,$cAnnot));
		if($coversStart){
			$pAnnot->addClassification($self->getStartLostVariantClass);
		}
		if($coversEnd){
			$pAnnot->addClassification($self->getStopLostVariantClass);
		}
		return $pAnnot;
	}
	my $cdsMinPos = $self->_getCdsMinPosForProteinCalculation($cAnnot);
 	my $cdsMaxPos = $self->_getCdsMaxPosForProteinCalculation($cAnnot);
 	my $desc = undef;
 	my $mutProtMin = undef;
  	my $mutProtMax = undef;
  	my $wt = undef;
  	my $mt = undef;
  	my $subtype = undef;
  	my $type = undef;
 	unless(defined $cdsMinPos && defined $cdsMaxPos) {
 		# something has gone wrong
 		return undef;
 	}

	my $mtDna = $self->_getMutatedCdsSequence($wtDna,$cdsMinPos,$cdsMaxPos,$cAnnot->getMt());
	my $mtProt = Bio::Seq->new(-seq => $prePad . $mtDna . $postPad)->translate->seq(); # mutant protein sequence
	my $maxMtProt = Bio::Seq->new(-seq => $prePad . $mtDna . substr($tran->getcDNASeq,$tran->getCdsMaxPos()))->translate->seq(); # maximised protein sequence, overruns the natural stop and translates to the end of the transcript
	if($wtProt eq $mtProt){
		# wt and mt protein sequences are the same, its silent
		$mutProtMin = ceil(($cAnnot->getMinPos / 3));
  	$mutProtMax = ceil(($cAnnot->getMaxPos / 3));
  	$wt = substr($wtProt,($mutProtMin - 1),(($mutProtMax - $mutProtMin) + 1));
  	$mt = substr($mtProt,($mutProtMin - 1),(($mutProtMax - $mutProtMin) + 1));
  	if(length($wt) == 1 && length($mt) == 1 && $mutProtMin == $mutProtMax){
		  $desc = 'p.'.$wt.$mutProtMin.$mt;
  	} else {
			$desc = 'p.(=)';
		}
		$type = $self->_getDefaultProteinAnnotationType();
		$subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype();
		if($wt eq '*'){
			push(@classes,$self->getStopRetainedVariantClass);
		} else {
			push(@classes,$self->getSynonymousVariantClass);
		}
	} elsif($cAnnot->hasClassification($self->getFrameShiftVariantClass)){
		# frame shift
    if($wtProt eq $maxMtProt || $self->_sequenceStartsWithSequence($maxMtProt,$wtProt)){
 			my $msg = 'very strange frameshift, there is no change in the protein sequence so calling it unknown';
 			$self->addMessage($msg);
 			$log->info($msg);
 			return $self->_buildUnknownProteinAnnotation($var,$tran,$cAnnot,length($wtProt),@classes);
 		}
 		($subtype,$wt,$mt,$mutProtMin,$mutProtMax,$desc) = $self->_calculateFrameShiftProteinChange($wtProt,$maxMtProt);
 		if($mutProtMin == 1){
 	  		# its frame shifted the start codon, no idea what this is going to cause.
 	  		push(@classes,$self->getStartLostVariantClass);
	  		return $self->_buildUnknownProteinAnnotation($var,$tran,$cAnnot,length($wtProt),@classes);
    }
 		$type = Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType();
    push(@classes,$self->getFrameShiftVariantClass);
	} else {
		$wt = $wtProt;
		$mt = $mtProt;
		$mutProtMin = 0;
		while(substr($wt,0,1) eq substr($mt,0,1)){
			substr($wt,0,1,'');
			substr($mt,0,1,'');
			$mutProtMin++;
		}
		while(substr($wt,-1,1) eq substr($mt,-1,1)){
			substr($wt,-1,1,'');
			substr($mt,-1,1,'');
		}

		#warn "|$wt| to |$mt|\n";
		if($wt ne ''){
			# wild type residue has been changed
			if($mt ne '' && length($wt) == length($mt) && length($wt) == 1){
				# its a simple sub
				$mutProtMin++;
				$mutProtMax = $mutProtMin;
				$type = Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType();
				$subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype();
				$desc = 'p.'.$wt.$mutProtMin.$mt;
				push(@classes,$self->getSubstitutionClass);
				if($wt eq $mt){
					# this could still happen, but it should have already been caught
					$desc = 'p.(=)';
					if($wt eq '*'){
  						push(@classes,$self->getStopRetainedVariantClass);
  					} else {
  						push(@classes,$self->getSynonymousVariantClass);
  					}
				} else {
					if($mt eq '*'){
						push(@classes,$self->getStopGainedVariantClass);
					} elsif($wt eq '*'){
						push(@classes,$self->getStopLostVariantClass);
					} elsif($wt eq 'M' && $mutProtMin == 1){
						push(@classes,$self->getStartLostVariantClass);
					} else {
						push(@classes,$self->getNonSynonymousVariantClass);
					}
				}
			} elsif ($mt ne '' && (length($wt) > 1 || length($mt) > 1)){
				# complex sub
				$mutProtMin++;
				$mutProtMax = ($mutProtMin + length($wt)) - 1;
				$type = Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType();
				$subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype();
				if($mt =~ s/^(.*?\*).*$/$1/ || $mt eq '*'){
					# this variant introduces a stop codon
					push(@classes,$self->getStopGainedVariantClass);
					push(@classes,$self->getComplexIndelClass);
				} else {
					push(@classes,$self->getComplexIndelClass);
				}
				if($mutProtMin == $mutProtMax){
					# replaces single WT residue
					$desc = 'p.'.substr($wtProt,($mutProtMin - 1),1).$mutProtMin.'delins'.$mt;
				} else {
					$desc = 'p.'.substr($wtProt,($mutProtMin - 1),1).$mutProtMin.'_'.substr($wtProt,($mutProtMax - 1),1).$mutProtMax.'delins'.$mt;
				}
			} elsif($mt eq ''){
				# its an inframe deletion
				$mt = '-';
				$mutProtMin++;
 				$mutProtMax = ($mutProtMin + length($wt)) - 1;
				$type = Sanger::CGP::Vagrent::Data::Annotation::getDeletionAnnotationType();
				$subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype();
				if($mutProtMin == $mutProtMax){
 					$desc = 'p.'.$wt.$mutProtMin.'del'.$wt;
 				} else {
 					$desc = 'p.'.substr($wt,0,1).$mutProtMin.'_'.substr($wt,-1,1).$mutProtMax.'del'.$wt;
 				}
				push(@classes,$self->getInFrameCodonLossVariantClass);
			}
		} elsif($wt eq '' && $mt ne '' && length($mt) > 0){
			# This is an in frame insertion
			$wt = '-';
			$mutProtMax = $mutProtMin + 1;
			$type = Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType();
			$subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype();
			if($mt =~ s/^(.*?\*).*$/$1/){
 				# this variant introduces a stop codon
 				push(@classes,$self->getStopGainedVariantClass);
 				if(length($mt) > 1){
 					push(@classes,$self->getInFrameCodonGainVariantClass);
 				}
 			} else {
 				push(@classes,$self->getInFrameCodonGainVariantClass);
 			}
			$desc = 'p.'.substr($wtProt,($mutProtMin - 1),1).$mutProtMin.'_'.substr($wtProt,($mutProtMax - 1),1).$mutProtMax.'ins'.$mt;
		} else {
			my $msg = 'Unhandled case, reporting unknown protein annotation';
 			$self->addMessage($msg);
 			$log->info($msg);
 			return $self->_buildUnknownProteinAnnotation($var,$tran,$cAnnot,length($wtProt),@classes);
		}
	}

	my $anno = Sanger::CGP::Vagrent::Data::Annotation->new( wt => $wt,
															mt => $mt,
															minpos => $mutProtMin,
															minOffset => 0,
															maxpos => $mutProtMax,
															maxOffset => 0,
															acc => $tran->getProteinAccession,
															accversion => $tran->getProteinAccessionVersion,
															db => $tran->getDatabase,
															dbversion => $tran->getDatabaseVersion,
															seqlength => length($wtProt),
															description => $desc,
															context => Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
															type => $type,
															subtype => $subtype);
  	$anno->addClassification(@classes);
  	return $anno;
}

sub _getMutatedCdsSequence: Abstract;

sub _getCdsMinPosForProteinCalculation: Abstract;

sub _getCdsMaxPosForProteinCalculation: Abstract;

sub _buildCDSAnnotation {
	my ($self,$var,$tran,$rAnnot) = @_;
	my @classes;
	unless($self->_canAnnotateToCDS($tran,$rAnnot)){
		my $msg = 'mRNA annotation doesnt affect the CDS, cant annotate this';
		$self->addMessage($msg);
		$log->info($msg);
		return $self->_buildUnknownCDSAnnotation($var,$tran,$rAnnot,@classes);
	}
	my ($cdsMin,$cdsMinOffset,$cdsMax,$cdsMaxOffset) = (undef,undef,undef,undef);
	if($rAnnot->getMinPos < $tran->getCdsMinPos){
		$cdsMin = 1;
		$cdsMinOffset = 0;
  } elsif($rAnnot->getMinPos == $tran->getCdsMinPos && $rAnnot->getMinOffset() < 0){
    $cdsMin = 1;
    if($self->_isIntronicOffsetDistance($rAnnot->getMinOffset())){
      $cdsMinOffset = 0;
    } else {
      $cdsMinOffset = $rAnnot->getMinOffset();
    }
	} else {
		$cdsMin = ($rAnnot->getMinPos - $tran->getCdsMinPos) + 1;
		$cdsMinOffset = $rAnnot->getMinOffset();
	}

  if($rAnnot->getMaxPos > $tran->getCdsMaxPos){
    $cdsMax = length($tran->getCdsSeq);
		$cdsMaxOffset = 0;
  } elsif($rAnnot->getMaxPos == $tran->getCdsMaxPos && $rAnnot->getMaxOffset() > 0){
    $cdsMax = length($tran->getCdsSeq);
    if($self->_isIntronicOffsetDistance($rAnnot->getMaxOffset())){
      $cdsMaxOffset = 0;
    } else {
      $cdsMaxOffset = $rAnnot->getMaxOffset();
    }
  } else {
    $cdsMax = ($rAnnot->getMaxPos() - $tran->getCdsMinPos) + 1;
		$cdsMaxOffset = $rAnnot->getMaxOffset();
  }

  print "CDS: $cdsMin , $cdsMinOffset - $cdsMax, $cdsMaxOffset\n" if $self->_debug();

	my $wt = $self->_getWildTypeStringForCDSAnno($var,$tran,$rAnnot);
	my $mt = $self->_getMutantStringForCDSAnno($var,$tran,$rAnnot);
	my $desc = $self->_getCDSDescriptionString($tran,$cdsMin,$cdsMax,$cdsMinOffset,$cdsMaxOffset,$wt,$mt);
	my $subtype = $rAnnot->getSubtype();
	my $anno = Sanger::CGP::Vagrent::Data::Annotation->new( wt => uc($wt),
															mt => uc($mt),
															minpos => $cdsMin,
															minOffset => $cdsMinOffset,
															maxpos => $cdsMax,
															maxOffset => $cdsMaxOffset,
															acc => $tran->getAccession,
															accversion => $tran->getAccessionVersion,
															db => $tran->getDatabase,
															dbversion => $tran->getDatabaseVersion,
															seqlength => length($tran->getCdsSeq),
															description => $desc,
															context => Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
															type => $self->_getDefaultCDSAnnotationType(),
															subtype => $subtype);

	foreach my $rc($rAnnot->getClassifications){
		next if($rc eq $self->get5KBUpStreamVariantClass ||
						$rc eq $self->get2KBUpStreamVariantClass ||
						$rc eq $self->get5KBDownStreamVariantClass ||
						$rc eq $self->get500BPDownStreamVariantClass );
		$anno->addClassification($rc);
	}
	$anno->addClassification(@classes);
	return $anno;



}

sub _getDefaultCDSAnnotationType: Abstract;

sub _getDefaultProteinAnnotationType: Abstract;

sub _getWildTypeStringForCDSAnno: Abstract;

sub _getMutantStringForCDSAnno: Abstract;

sub _getCDSDescriptionString: Abstract;

sub _buildUnknownCDSAnnotation {
	my ($self,$mut,$tran,$rnaAnno,@classes) = @_;
	my $anno = Sanger::CGP::Vagrent::Data::Annotation->new( wt => '?',
																													mt => '?',
																													minpos => 0,
																													minOffset => 0,
																													maxpos => 0,
																													maxOffset => 0,
																													acc => $tran->getAccession,
																													accversion => $tran->getAccessionVersion,
																													db => $tran->getDatabase,
																													dbversion => $tran->getDatabaseVersion,
																													seqlength => length($tran->getCdsSeq),
																													description => 'c.?',
																													context => Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
																													type => Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
																													subtype => Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype());
	$anno->addClassification(@classes);
	$anno->addClassification($self->getUnknownVariantClass);
	return $anno;
}

sub _buildUnaffectedCDSAnnotation {
	my ($self,$mut,$tran,$rnaAnno,@classes) = @_;
	my $anno = Sanger::CGP::Vagrent::Data::Annotation->new(wt => '-',
																																		mt => '-',
																																		minpos => 0,
																																		minOffset => 0,
																																		maxpos => 0,
																																		maxOffset => 0,
																																		acc => $tran->getAccession,
																																		accversion => $tran->getAccessionVersion,
																																		db => $tran->getDatabase,
																																		dbversion => $tran->getDatabaseVersion,
																																		seqlength => length($tran->getCdsSeq),
																																		description => 'c.=',
																																		context => Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
																																		type => Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
																																		subtype => Sanger::CGP::Vagrent::Data::Annotation::getNoneAnnotationSubtype());
	$anno->addClassification(@classes);
	return $anno;
}

sub _buildUnknownProteinAnnotation {
	my ($self,$mut,$tran,$cdsAnno,$protLength,@classes) = @_;
	my $anno = Sanger::CGP::Vagrent::Data::Annotation->new(wt => '?',
																																		mt => '?',
																																		minpos => 0,
																																		minOffset => 0,
																																		maxpos => 0,
																																		maxOffset => 0,
																																		acc => $tran->getProteinAccession,
																																		accversion => $tran->getProteinAccessionVersion,
																																		db => $tran->getDatabase,
																																		dbversion => $tran->getDatabaseVersion,
																																		seqlength => $protLength,
																																		description => 'p.?',
																																		context => Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
																																		type => Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
																																		subtype => Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype());
	$anno->addClassification(@classes);
	$anno->addClassification($self->getUnknownVariantClass);
	return $anno;
}

sub _buildUnaffectedProteinAnnotation {
	my ($self,$mut,$tran,$cdsAnno,$protLength,@classes) = @_;
	my $anno = Sanger::CGP::Vagrent::Data::Annotation->new(wt => '-',
																																		mt => '-',
																																		minpos => 0,
																																		minOffset => 0,
																																		maxpos => 0,
																																		maxOffset => 0,
																																		acc => $tran->getProteinAccession,
																																		accversion => $tran->getProteinAccessionVersion,
																																		db => $tran->getDatabase,
																																		dbversion => $tran->getDatabaseVersion,
																																		seqlength => $protLength,
																																		description => 'p.=',
																																		context => Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
																																		type => Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
																																		subtype => Sanger::CGP::Vagrent::Data::Annotation::getNoneAnnotationSubtype());
	$anno->addClassification(@classes);
	return $anno;
}

sub _buildRemovedProteinAnnotation {
	my ($self,$mut,$tran,$wtProt,@classes) = @_;
	my $anno = Sanger::CGP::Vagrent::Data::Annotation->new( wt => $wtProt,
															mt => '-',
															minpos => 1,
															minOffset => 0,
															maxpos => length($wtProt),
															maxOffset => 0,
															acc => $tran->getProteinAccession,
															accversion => $tran->getProteinAccessionVersion,
															db => $tran->getDatabase,
															dbversion => $tran->getDatabaseVersion,
															seqlength => length($wtProt),
															description => 'p.0',
															context => Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
															type => $self->_getDefaultProteinAnnotationType(),
															subtype => Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype());
	$anno->addClassification(@classes);
	return $anno;
}

sub _buildUnknownMRNAAnnotation {
	my ($self,$var,$tran,@classes) = @_;
	my $anno = Sanger::CGP::Vagrent::Data::Annotation->new(wt => '?',
																																		mt => '?',
																																		minpos => 0,
																																		minOffset => 0,
																																		maxpos => 0,
																																		maxOffset => 0,
																																		acc => $tran->getAccession,
																																		accversion => $tran->getAccessionVersion,
																																		db => $tran->getDatabase,
																																		dbversion => $tran->getDatabaseVersion,
																																		seqlength => length($tran->getcDNASeq),
																																		description => 'r.?',
																																		context => Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
																																		type => Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
																																		subtype => Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype());

	$anno->addClassification(@classes);
	return $anno;
}

sub _isOffsetAConsensusSpliceDistance {
	my ($self,$offset) = @_;
  foreach my $cf(@CONSENSUS_SPLICE_OFFSETS){
  	if($offset == $cf){
  		return 1;
  	}
  }
	return 0;
}

sub _getConsesnsusSpliceBeforeBoundry {
	return $CONSENSUS_SPLICE_BEFORE_BOUNDRY;
}

sub _getConsesnsusSpliceAfterBoundry {
	return $CONSENSUS_SPLICE_AFTER_BOUNDRY;
}

sub _isIntronicOffsetDistance {
	my ($self,$offset) = @_;
	if(abs($offset) >= $INTRONIC_OFFSET_CUTOFF){
		return 1;
	}
	return 0;
}

sub _isWithin5KBOffsetDistance {
	my ($self,$offset) = @_;
	if(abs($offset) <= $UPDOWNSTREAM_5KB_CUTOFF){
		return 1;
	}
	return 0;
}

sub _isWithin2KBUpstreamOffsetDistance {
	my ($self,$offset) = @_;
	if($offset < 0 && $offset >= $UPSTREAM_2KB_CUTOFF){
		return 1;
	}
	return 0;
}

sub _isWithin500BPDownstreamOffsetDistance {
	my ($self,$offset) = @_;
	if($offset > 0 && $offset <= $DOWNSTREAM_500BP_CUTOFF){
		return 1;
	}
	return 0;
}

sub _coversStartCodon {
	my ($self,$tran,$anno) = @_;
	unless($tran->isProteinCoding){
		# if the transcript isn't protein coding it can't have a start codon
		return 0;
	}

  my ($startMin,$startMax);
  if($anno->getContext eq Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()){
    $startMin = $tran->getCdsMinPos;
    $startMax = $tran->getCdsMinPos + 2;
  } elsif($anno->getContext eq Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()){
    $startMin = 1;
    $startMax = 3;
  } else {
    # don't know, assume no
		return 0;
  }

  if($anno->hasClassification($self->getInsertionClass)){
    # insertions are a special case, coordinates are outside the variant
    if($anno->getMinPos < $startMax && $anno->getMaxPos > $startMin){
      # var started before the end of the first codon, and ended after the start of the first codon
      return 1;
    }
  } else {
    if(($anno->getMinPos < $startMax || ($anno->getMinPos == $startMax && $anno->getMinOffset == 0)) &&
        ($anno->getMaxPos > $startMin || ($anno->getMaxPos == $startMin && $anno->getMaxOffset == 0))){
      # var started before the end of the first codon OR explicitly on the last base of the first codon AND
      # var ended after the start of the first codon OR explicitly on the first base of the first codon
      return 1;
    }
  }
	return 0;
}

sub _coversStopCodon {
	my ($self,$tran,$anno) = @_;
	unless($tran->isProteinCoding){
		# if the transcript isn't protein coding it can't have a stop codon
		return 0;
	}
  my ($stopMin,$stopMax);
  if($anno->getContext eq Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()){
    $stopMin = $tran->getCdsMaxPos - 2;
    $stopMax = $tran->getCdsMaxPos;
  } elsif($anno->getContext eq Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()){
    $stopMin = $tran->getCdsLength - 2;
    $stopMax = $tran->getCdsLength;
  } else {
    # don't know, assume no
		return 0;
  }

  if($anno->hasClassification($self->getInsertionClass)){
    # insertions are a special case, coordinates are outside the variant
    if($anno->getMinPos < $stopMax && $anno->getMaxPos > $stopMin){
      # var started before the end of the last codon, and ended after the start of the last codon
      return 1;
    }
  } else {
    if(($anno->getMinPos < $stopMax || ($anno->getMinPos == $stopMax && $anno->getMinOffset == 0)) &&
        ($anno->getMaxPos > $stopMin || ($anno->getMaxPos == $stopMin && $anno->getMaxOffset == 0))){
      # var started before the end of the last codon OR explicitly on the last base of the last codon AND
      # var ended after the start of the last codon OR explicitly on the first base of the last codon
      return 1;
    }
  }
	return 0;
}

sub _getTranscriptSource {
	return shift->{_transcriptSource};
}

sub _calculateCdsTranslationPadStrings {
	my ($self,$tran) = @_;
	my $prePad = '';
	my $postPad = '';
	if($tran->getCdsPhase == 0){
		# nothing to pad
		$prePad = '';
	} elsif($tran->getCdsPhase > 0 && $tran->getCdsPhase < 3){
		if($tran->getCdsPhase == 1){
			$prePad = 'N';
		} elsif ($tran->getCdsPhase == 2){
			$prePad = 'NN';
		} else {
			$self->throw('this should be impossible');
		}
	} else {
		warn Dumper($tran);
		$self->throw("Unhandled phase");
	}

	if(length($prePad . $tran->getCdsSeq) % 3 == 0){
		# nothing to pad
		$postPad = '';
	} else {
		# yuk, its not a round number of codons long, got to pad the end.
		my $rem = 3 - (length($prePad . $tran->getCdsSeq) % 3);
		for(my $i = 0 ; $i < $rem; $i++){
			$postPad .= 'N';
		}
	}
	return ($prePad,$postPad);
}

sub _canAnnotateToCDS {
	my ($self,$tran,$anno) = @_;
	unless($tran->isProteinCoding){
		# if the transcript isn't protein coding it can't be a coding change
		return 0;
	}
	if($anno->getContext eq Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()){
		# first work out if the variant overlaps with the CDS
		if($anno->hasClassification($self->getInsertionClass)){
			# insertions are a special case.
			# Coordinates are the last WT positions, and not the first variant ones like everything else

      print 'ANNO POS: '.$anno->getMinPos.' , '.$anno->getMinOffset.' - '.$anno->getMaxPos.' , '.$anno->getMaxOffset."\n" if $self->_debug();
      print 'CDS POS: '.$tran->getCdsMinPos.' , '.$tran->getCdsMaxPos."\n" if $self->_debug();

      if($anno->getMaxPos < $tran->getCdsMinPos || $anno->getMinPos > $tran->getCdsMaxPos){
        # ends before CDS or starts afterwards
        return 0;
      } elsif($anno->getMaxPos == $tran->getCdsMinPos) {
        # potential start codon issues
        if($anno->getMinPos == $anno->getMaxPos && $anno->getMinPos == $tran->getCdsMinPos && abs($anno->getMinOffset) + abs($anno->getMaxOffset) > 0){
          # probably start coordinate issues
          unless($anno->getMaxOffset <= 0 && $self->_isIntronicOffsetDistance($anno->getMaxOffset) == 0){
            # or not
            return 0;
          }
        } else {
          return 0;
        }
      }
		} else {
			if($anno->getMaxPos < $tran->getCdsMinPos || $anno->getMinPos > $tran->getCdsMaxPos){
				# its outside the CDS
				return 0;
			}
		}
		# now we look at the classifications to workout if can be annotated to the CDS
		if($anno->hasClassification($self->getSpliceRegionVariantClass) || $anno->hasClassification($self->getEssentialSpliceSiteVariantClass)){
			return 1;
		} elsif($anno->hasClassification($self->getComplexChangeVariantClass)){
			return 1;
		} elsif($anno->hasClassification($self->getInFrameVariantClass)){
			return 1;
		} elsif($anno->hasClassification($self->getFrameShiftVariantClass)){
			return 1;
		} elsif($anno->hasClassification($self->getCodonVariantClass)){
			return 1;
		} elsif($anno->hasClassification($self->getIntronVariantClass)){
			return 0;
		} elsif($anno->hasClassification($self->getUnknownVariantClass)){
			return 0;
    } elsif($anno->hasClassification($self->getInsertionClass) && $anno->hasClassification($self->get5PrimeUtrVariantClass)){
      # odd case, insertions close to the start codons can be described on the CDS even though they don't change it.
      return 1;
		} else {
			my $msg = "Unable to calculate CDS relevance - UNKNOWN CLASSIFICATION: ".join(' ',$anno->getClassifications);
			$self->addMessage($msg);
			$log->info($msg);
			return 0;
		}
	} else {
		my $msg = "Unable to calculate CDS relevance - UNKNOWN CONTEXT: ".$anno->getContext;
		$self->addMessage($msg);
		$log->info($msg);
		return 0;
	}
}

sub _canAnnotateToProtein {
	my ($self,$tran,$anno) = @_;

	unless($tran->isProteinCoding){
		# if the transcript isn't protein coding it can't be a coding change
		return 0;
	}
	if($anno->getContext eq Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()){
		if($anno->hasClassification($self->getUnknownVariantClass)){
			return 0;
		} elsif($anno->hasClassification($self->getEssentialSpliceSiteVariantClass)){
			return 0;
		} elsif($anno->hasClassification($self->getSpliceRegionVariantClass)){
			return 0;
		} elsif($anno->hasClassification($self->getComplexChangeVariantClass)){
			return 0;
		} elsif($anno->hasClassification($self->getInFrameVariantClass)){
			return 1;
		} elsif($anno->hasClassification($self->getFrameShiftVariantClass)){
			return 1;
		} elsif($anno->hasClassification($self->getCodonVariantClass)){
			return 1;
		} else {
			my $msg = "Unable to calculate protein relevance - UNKNOWN CLASSIFICATION: ".join(' ',$anno->getClassifications);
			$self->addMessage($msg);
			$log->info($msg);
			return 0;
		}
	} else {
		my $msg = "Unable to calculate protein relevance - UNKNOWN CONTEXT: ".$anno->getContext;
		$self->addMessage($msg);
		$log->info($msg);
		return 0;
	}
}

sub _revcompSeq {
	my ($self,$seqIn) = @_;
	my $seqOut;
	foreach my $b(reverse(split('',$seqIn))){
		$b =~ tr/atcgATCG/tagcTAGC/;
		$seqOut .= $b;
	}
	return $seqOut;
}

sub _sequenceStartsWithSequence {
	my ($self,$seqRef,$seqCheck) = @_;
	warn "CHECKING IF\n$seqRef\nSTARTS WITH\n$seqCheck\n" if ($self->_debug);
	if(length($seqRef) < length($seqCheck)){
		warn "FALSE\n" if ($self->_debug);
		return 0;
	}
	warn substr($seqRef,0,length($seqCheck))."\n" if ($self->_debug);
	warn $seqCheck."\n" if ($self->_debug);

	if(substr($seqRef,0,length($seqCheck)) eq $seqCheck){
		warn "TRUE\n" if ($self->_debug);
		return 1;
	}
	warn "FALSE\n" if ($self->_debug);
	return 0;
}

sub _calculateFrameShiftProteinChange {
	my ($self,$wtProt,$mtProt) = @_;
	my $subtype = Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype();
	my $wt = undef;
	my $mt = undef;
	my $mutProtMin = undef;
	my $mutProtMax = undef;
	my $desc = undef;
	warn "\n$wtProt\n$mtProt\n" if ($self->_debug);
	for(my $i = 0 ; $i < length($mtProt); $i++){
		warn "$wt, ".substr($wtProt,$i,1)." ne ".substr($mtProt,$i,1)."\n" if ($self->_debug);
		if(!defined($wt) && substr($wtProt,$i,1) ne substr($mtProt,$i,1)){
			$wt = substr($wtProt,$i,1);
			$mt = substr($mtProt,$i,1);
			$mutProtMin = $i + 1;
			$mutProtMax = $mutProtMin;
			if($mt eq '*'){
				# variation actually creates the stop codon
				$desc = 'p.'.$wt.$mutProtMin.'fs*1';
				last;
			} else {
				# we'll reserve an I don't know description, should get overridden in further iterations.
				$desc = 'p.'.$wt.$mutProtMin.'fs?';
			}
		} elsif(!defined($wt) && length($mtProt) == $i + 1 && length($wtProt) > $i + 1 && substr($wtProt,$i,1) eq substr($mtProt,$i,1)){
			# interesting case, the current position is the last position of the mutant protein
			# no differences have been found between mt and wt, but the wt sequence is longer, ie we have lost residues from the end
			$wt = substr($wtProt,($i + 1),1);
			$mt = '?';
			$mutProtMin = $i + 2;
			$mutProtMax = $mutProtMin;
			$desc = 'p.'.$wt.$mutProtMin.'fs?';
			last;
		} elsif(defined($wt) && length($mtProt) == $i + 1 && substr($mtProt,$i,1) ne '*'){
			# this is the last base of the mutant protein and we never reach a stop.
			$mt .= substr($mtProt,$i,1);
			$desc = 'p.'.$wt.$mutProtMin.'fs*>'.length($mt);
			last;
		} elsif(defined($wt)){
			$mt .= substr($mtProt,$i,1);
			if(substr($mtProt,$i,1) eq '*'){
				$desc = 'p.'.$wt.$mutProtMin.'fs*'.length($mt);
				last;
			}
		}
	}
	warn "$subtype,$wt,$mt,$mutProtMin,$mutProtMax,$desc\n" if ($self->_debug);
	return ($subtype,$wt,$mt,$mutProtMin,$mutProtMax,$desc);
}

sub _calculateGenomicCdsPosition {
	my ($self,$tran) = @_;
  my @exons = $tran->getExonsGenomicOrder;
	my $genoCdsMin = undef;
	my $genoCdsMax = undef;
	foreach my $e(@exons){
		if($tran->getStrand == 1){
			if(!defined($genoCdsMin)){
				if($e->getRnaMinPos <= $tran->getCdsMinPos && $e->getRnaMaxPos >= $tran->getCdsMinPos){
					$genoCdsMin = $e->getMinPos + ($tran->getCdsMinPos - $e->getRnaMinPos);
				}
			}
			if(!defined($genoCdsMax)){
				if($e->getRnaMinPos <= $tran->getCdsMaxPos && $e->getRnaMaxPos >= $tran->getCdsMaxPos){
					$genoCdsMax = $e->getMinPos + ($tran->getCdsMaxPos - $e->getRnaMinPos);
				}
			}
		} else {
			if(!defined($genoCdsMin)){
				if($e->getRnaMinPos <= $tran->getCdsMaxPos && $e->getRnaMaxPos >= $tran->getCdsMaxPos){
					$genoCdsMin = $e->getMaxPos - ($tran->getCdsMaxPos - $e->getRnaMinPos);
				}
			}
			if(!defined($genoCdsMax)){
				if($e->getRnaMinPos <= $tran->getCdsMinPos && $e->getRnaMaxPos >= $tran->getCdsMinPos){
					$genoCdsMax = $e->getMaxPos - ($tran->getCdsMinPos - $e->getRnaMinPos);
				}
			}
		}
	}
	return ($genoCdsMin,$genoCdsMax);
}

sub _arrayHasString {
	my ($self,$val,@arr) = @_;
	foreach my $a (@arr){
		return 1 if $a eq $val;
	}
	return 0;
}

sub _defaultTranscriptSort {
	my ($self,@trans) = @_;
	my $customSort = sub {
		my $accds = 0;
		my $bccds = 0;
		$accds = 1 if(defined $a->getCCDS && length($a->getCCDS) > 0);
		$bccds = 1 if(defined $b->getCCDS && length($b->getCCDS) > 0);
		my $chk = $bccds <=> $accds;
		if($chk == 0){
			$chk = $b->getCdsLength <=> $a->getCdsLength;
			if($chk == 0){
				$chk = length($b->getcDNASeq()) <=> length($a->getcDNASeq());
			}
		}
		return $chk;};
	return sort $customSort @trans;
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::Annotators::AbstractAnnotator - Abstract base class for the annotation generators

=head1 DESCRIPTION

This is an abstract template class for the mutation annotators, it provides a lot of shared behind the scenes functionality.  All
subclasses must implement the _getAnnotation method.

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $source = Sanger::CGP::Vagrent::Annotators::AbstractAnnotatorSubClass->new(%params);

=item Function :

Builds a new Sanger::CGP::Vagrent::Annotators::AbstractAnnotator inheriting object

=item Returns :

Sanger::CGP::Vagrent::Annotators::AbstractAnnotator object initialized with parameter values

=item Params :

Hash of parameter values

 transcriptSource => A Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource inheriting object
 bookmarker       => (Optional) An array reference of, or single, Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker inheriting object
 only_bookmarked  => (Optional) Boolean, only return annotations that get bookmarked

=back

=head2 Functions

=head3 getAnnotation

=over

=item Usage :

 my @annoGrps = $annotator->getAnnotation($variation);

=item Function :

Annotates the supplied L<Variation|Sanger::CGP::Vagrent::Data::AbstractVariation> object and returns a list of L<AnnotationGroups|Sanger::CGP::Vagrent::Data::AnnotationGroup> objects.
If L<Bookmarkers|Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker> have been set, the AnnotationGroups will have been marked before being returned.
If 'only_bookmarked' was set to true, only AnnotationGroups that match Bookmarkers will be returned.

=item Returns :

An array of L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> objects

=item Params :

A L<Sanger::CGP::Vagrent::Data::AbstractVariation|Sanger::CGP::Vagrent::Data::AbstractVariation> implementing object

=back

=head3 addMessage

=over

=item Usage :

 $annotator->addMessage("Interesting event found");

=item Function :

Adds a text message to the message list.  All messages are reset every time C<getAnnotation> is called

=item Returns :

Nothing

=item Params :

String

=back

=head3 getMessages

=over

=item Usage :

 my @mess = $annotator->getMessages();

=item Function :

Retrieves a list of message strings about the most recent annotation attempt

=item Returns :

Array of String

=back

=head2 Abstract

=head3 _getAnnotation

=over

=item Usage :

 my $type = $annotator->_getAnnotation($variation);

=item Function :

Abstract internal function, must be implemented in subclass.  Returns a list of L<AnnotationGroups|Sanger::CGP::Vagrent::Data::AnnotationGroup>

=item Returns :

An array of L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> objects

=item Params :

A L<Sanger::CGP::Vagrent::Data::AbstractVariation|Sanger::CGP::Vagrent::Data::AbstractVariation> implementing object

=back

=head3 _getDefaultCDSAnnotationType

=over

=item Usage :

 my $type = $annotator->_getDefaultCDSAnnotationType();

=item Function :

Abstract internal function, must be implemented in subclass.  Returns default CDS annotation (L<Annotation|Sanger::CGP::Vagrent::Data::Annotation> type constant) type for the Annotator

=item Returns :

String - L<Annotation|Sanger::CGP::Vagrent::Data::Annotation> type constant

=back

=head3 _getMutatedCdsSequence

=over

=item Usage :

 my $varSeq = $annotator->_getMutatedCdsSequence();

=item Function :

Abstract internal function, must be implemented in subclass.  Returns the full variant form of the CDS sequence

=item Returns :

String, DNA sequence

=item Params :

 String - Full wildtype CDS DNA sequence
 Integer - Minimum position of the variant on the wildtype CDS sequence
 Integer - Maximum position of the variant on the wildtype CDS sequence
 String - Mutant DNA sequence for the variant

=back

=head3 _getWildTypeStringForCDSAnno

=over

=item Usage :

 my $wtseq = $annotator->_getWildTypeStringForCDSAnno($var,$tran,$mrnaAnno);

=item Function :

Abstract internal function, must be implemented in subclass.  Generates the CDS wildtype string from the mRNA annotation

=item Returns :

String, DNA sequence

=item Params :

 A Sanger::CGP::Vagrent::Data::AbstractVariation implementing object
 A Sanger::CGP::Vagrent::Data::Transcript object
 A Sanger::CGP::Vagrent::Data::Annotation object

=back

=head3 _getMutantStringForCDSAnno

=over

=item Usage :

 my $varseq = $annotator->_getMutantStringForCDSAnno($var,$tran,$mrnaAnno);

=item Function :

Abstract internal function, must be implemented in subclass.  Generates the CDS variant string from the mRNA annotation

=item Returns :

String, DNA sequence

=item Params :

 A Sanger::CGP::Vagrent::Data::AbstractVariation implementing object
 A Sanger::CGP::Vagrent::Data::Transcript object
 A Sanger::CGP::Vagrent::Data::Annotation object

=back

=head3 _getCDSDescriptionString

=over

=item Usage :

 my $desc = $annotator->_getCDSDescriptionString($tran,$cdsMin,$cdsMax,$cdsMinOffset,$cdsMaxOffset,$wt,$mt);

=item Function :

Abstract internal function, must be implemented in subclass.  Takes the plotted CDS variation data and returns the HGVS syntax describing the change

=item Returns :

String, HGVS description

=item Params :

 A Sanger::CGP::Vagrent::Data::Transcript object
 Integer - CDS minimum position
 Integer - CDS maximum position
 Integer - Offset from the CDS minimum position (signed)
 Integer - Offset from the CDS maximum position (signed)
 String - Wildtype cDNA sequence of variant
 String - Variant cDNA sequence of variant

=back

=head3 _getDefaultProteinAnnotationType

=over

=item Usage :

 my $type = $annotator->_getDefaultProteinAnnotationType();

=item Function :

Abstract internal function, must be implemented in subclass.  Returns default Protein annotation (L<Annotation|Sanger::CGP::Vagrent::Data::Annotation> type constant) type for the Annotator

=item Returns :

String - L<Annotation|Sanger::CGP::Vagrent::Data::Annotation> type constant

=back

=head3 _getCdsMinPosForProteinCalculation

=over

=item Usage :

 my $max = $annotator->_getCdsMinPosForProteinCalculation($cDNAAnnot);

=item Function :

Abstract internal function, must be implemented in subclass.  Returns the minimum cDNA location of the variant that can be used for protein annotation

=item Returns :

Integer - cDNA position

=item Params :

A L<Annotation|Sanger::CGP::Vagrent::Data::Annotation> representing the cDNA annotation

=back

=head3 _getCdsMaxPosForProteinCalculation

=over

=item Usage :

 my $max = $annotator->_getCdsMaxPosForProteinCalculation($cDNAAnnot);

=item Function :

Abstract internal function, must be implemented in subclass.  Returns the maximum cDNA location of the variant that can be used for protein annotation

=item Returns :

Integer - cDNA position

=item Params :

A L<Annotation|Sanger::CGP::Vagrent::Data::Annotation> representing the cDNA annotation

=back

=head3 _isStartGained

=over

=item Usage :

 if($annotator->_isStartGained($var,$tran,$mRNAmin,$mRNAmax,$wt,$mt)){
   .......
 }

=item Function :

Abstract internal function, must be implemented in subclass.  Returns returns true if there is a start codon created or moved in the 5' UTR

=item Returns :

Boolean

=item Params :

 A Sanger::CGP::Vagrent::Data::AbstractVariation implementing object
 A Sanger::CGP::Vagrent::Data::Transcript object
 Integer - Minimum mRNA coordinate of the variant
 Integer - Maximum mRNA coordinate of the variant
 String - Wildtype variant sequence
 String - Mutant variant sequence

=back
