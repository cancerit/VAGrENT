package Sanger::CGP::Vagrent::Bookmarkers::MostDeleteriousBookmarker;

use strict;

use Log::Log4perl qw(:easy);
use Data::Dumper;
use Carp qw(croak);

use base qw(Sanger::CGP::Vagrent::Bookmarkers::RepresentativeTranscriptBookmarker);

1;
use constant DOWNSTREAM_SCORE => 4;
use constant UPSTREAM_SCORE => 7;
use constant INTRONIC_SCORE => 10;
use constant NONCODING_GENE_SCORE => 15;
use constant COMPLEX_IN_MRNA_SCORE => 50;
use constant NONCODING_GENE_SPLICE_REGION_SCORE => 95;
use constant THREEPRIME_UTR_SPLICE_REGION_SCORE => 100;
use constant FIVEPRIME_UTR_SPLICE_REGION_SCORE => 105;
use constant CODING_SPLICE_REGION_SCORE => 200;
use constant FIVEPRIME_UTR_SCORE => 300;
use constant THREEPRIME_UTR_SCORE => 305;
use constant NONCODING_GENE_ESS_SPLICE_SCORE => 395;
use constant THREEPRIME_UTR_ESS_SPLICE_SCORE => 400;
use constant FIVEPRIME_UTR_ESS_SPLICE_SCORE => 405;
use constant START_GAINED_SCORE => 450;
use constant SYNONYMOUS_SCORE => 500;
use constant COMPLETE_NONCODING_TRANSCRIPT_LOSS_SCORE => 525;
use constant COMPLEX_IN_CDS_SCORE => 550;
use constant NON_SYNONYMOUS_SCORE => 600;
use constant STOP_LOST_SCORE => 700;
use constant INITIATOR_CHANGE_SCORE => 800;
use constant INFRAME_CODON_GAIN_SCORE => 825;
use constant INFRAME_CODON_LOSS_AND_GAIN_SCORE => 840;
use constant INFRAME_CODON_LOSS_SCORE => 850;
use constant CODING_ESS_SPLICE_SCORE => 900;
use constant STOP_GAINED_SCORE => 1000;
use constant FRAMESHIFT_SCORE => 1100;
use constant COMPLETE_PROTEIN_LOSS_SCORE => 1200;

sub _getAnnotation {
	my $self = shift;
	my @groups = @_;
	my @sorted = $self->_sortAnnotations(@groups);
	return $self->_getMostDeleterious(@sorted);
}

sub _getMostDeleterious {
	my $self = shift;
	my @groups = @_;
	my $mostScore = 0;
	my $mostGroup = undef;
	foreach my $g(@groups){
		my $score = 1;
		my $mrna = $g->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext);
		my $cds = $g->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext);
		my $prot = $g->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext);
		if($g->hasClassification($self->getProteinCodingClass)){
			# protein coding transcript
			if(defined($prot) && $prot->getType() ne $mrna->getUnknownAnnotationType){
				# Protein annotation
				if($self->COMPLETE_PROTEIN_LOSS_SCORE > $score &&
					$prot->hasClassification($self->getDeletionClass) &&
					$prot->getMinPos() == 1 && $prot->getMaxPos() == $prot->getSequenceLength()){
					# if marked as a deletion, start = 1 and end = protein length the protein is gone.
						$score = $self->COMPLETE_PROTEIN_LOSS_SCORE;
				}
				if($self->FRAMESHIFT_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$prot->hasClassification($self->getFrameShiftVariantClass)){
					# Frameshift
						$score = $self->FRAMESHIFT_SCORE;
				}
				if($self->STOP_GAINED_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$prot->hasClassification($self->getStopGainedVariantClass)){
					# non-sense / stop gained
						$score = $self->STOP_GAINED_SCORE;
				}
				if($self->INFRAME_CODON_LOSS_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$prot->hasClassification($self->getInFrameCodonLossVariantClass)){
					# in frame deletion
						$score = $self->INFRAME_CODON_LOSS_SCORE;
				}
				if($self->INFRAME_CODON_LOSS_AND_GAIN_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$prot->hasClassification($self->getComplexIndelClass)){
					# in frame complex sub
						$score = $self->INFRAME_CODON_LOSS_AND_GAIN_SCORE;
				}
				if($self->INFRAME_CODON_GAIN_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$prot->hasClassification($self->getInFrameCodonGainVariantClass)){
					# in frame insertion
						$score = $self->INFRAME_CODON_GAIN_SCORE;
				}
				if($self->INITIATOR_CHANGE_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$prot->hasClassification($self->getStartLostVariantClass)){
					# start lost
						$score = $self->INITIATOR_CHANGE_SCORE;
				}
				if($self->STOP_LOST_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$prot->hasClassification($self->getStopLostVariantClass)){
					# stop lost
						$score = $self->STOP_LOST_SCORE;
				}
				if($self->NON_SYNONYMOUS_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$prot->hasClassification($self->getNonSynonymousVariantClass)){
					# mis sense
						$score = $self->NON_SYNONYMOUS_SCORE;
				}
				if($self->SYNONYMOUS_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					($prot->hasClassification($self->getSynonymousVariantClass) || $prot->hasClassification($self->getStopRetainedVariantClass))){
					# silent including terminator silent
						$score = $self->SYNONYMOUS_SCORE;
				}
			} elsif(defined($cds) && $cds->getType() ne $mrna->getUnknownAnnotationType){
				# CDS annotation
				if($self->FRAMESHIFT_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$cds->hasClassification($self->getFrameShiftVariantClass)){
					# frame shift again, incase protein translation was too complex
						$score = $self->FRAMESHIFT_SCORE;
				}
				if($self->CODING_ESS_SPLICE_SCORE > $score &&
					$g->hasClassification($self->getEssentialSpliceSiteClass) && $g->hasClassification($self->getCDSClass) &&
					$cds->hasClassification($self->getEssentialSpliceSiteVariantClass)){
					# essential splice change in CDS
						$score = $self->CODING_ESS_SPLICE_SCORE;
				}

				if($self->COMPLEX_IN_CDS_SCORE > $score &&
					$g->hasClassification($self->getCDSClass) &&
					$cds->hasClassification($self->getComplexChangeVariantClass)){
					# complex transcript consequence involving CDS
						if($cds->getMinPos() == 1 && $cds->getMaxPos() == $cds->getSequenceLength()){
							# position 1 to CDS length effected, transcript lost
							$score = $self->COMPLETE_PROTEIN_LOSS_SCORE;
						} elsif($g->hasClassification($self->getEssentialSpliceSiteClass) && $g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass)){
							# essential splice change
							$score = $self->CODING_ESS_SPLICE_SCORE;
						} elsif($g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
									($g->hasClassification($self->get5PrimeUtrClass) || $mrna->hasClassification($self->get2KBUpStreamVariantClass)) &&
									$cds->getMinPos() == 1){
							# start codon lost
							$score = $self->INITIATOR_CHANGE_SCORE;
						} elsif($g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
									($g->hasClassification($self->get3PrimeUtrClass) || $mrna->hasClassification($self->get500BPDownStreamVariantClass)) &&
									$cds->getMaxPos() == $cds->getSequenceLength()){
							# stop codon lost
							$score = $self->STOP_LOST_SCORE;
						} else {
							# if its none of the above, its just complex in CDS
							$score = $self->COMPLEX_IN_CDS_SCORE;
						}
				}
				if($self->CODING_SPLICE_REGION_SCORE > $score &&
					$g->hasClassification($self->getSpliceRegionClass) && $g->hasClassification($self->getCDSClass) &&
					$cds->hasClassification($self->getSpliceRegionVariantClass)){
					# splice region change in CDS
						$score = $self->CODING_SPLICE_REGION_SCORE;
				}

			} elsif(defined($mrna)) {
				# cDNA annotation
				if($self->START_GAINED_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->get5PrimeUtrClass) &&
					$mrna->hasClassification($self->getPrematureStartGainedVariantClass)){
						# new start codon created in the 5' UTR
						$score = $self->START_GAINED_SCORE;
				}
				if($self->FIVEPRIME_UTR_ESS_SPLICE_SCORE > $score &&
					$g->hasClassification($self->getEssentialSpliceSiteClass) && $g->hasClassification($self->get5PrimeUtrClass) &&
					$mrna->hasClassification($self->getEssentialSpliceSiteVariantClass) && $mrna->hasClassification($self->get5PrimeUtrVariantClass)){
						# essential splice change in 5' UTR
						$score = $self->FIVEPRIME_UTR_ESS_SPLICE_SCORE;
				}
				if($self->THREEPRIME_UTR_ESS_SPLICE_SCORE > $score &&
					$g->hasClassification($self->getEssentialSpliceSiteClass) && $g->hasClassification($self->get3PrimeUtrClass) &&
					$mrna->hasClassification($self->getEssentialSpliceSiteVariantClass) && $mrna->hasClassification($self->get3PrimeUtrVariantClass)){
						# essential splice change in 3' UTR
						$score = $self->THREEPRIME_UTR_ESS_SPLICE_SCORE;
				}
				if($self->THREEPRIME_UTR_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->get3PrimeUtrClass) &&
					$mrna->hasClassification($self->get3PrimeUtrVariantClass)){
						# change in 3' UTR exon
						$score = $self->THREEPRIME_UTR_SCORE;
				}
				if($self->FIVEPRIME_UTR_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->get5PrimeUtrClass) &&
					$mrna->hasClassification($self->get5PrimeUtrVariantClass)){
						# change in 5' UTR exon
						$score = $self->FIVEPRIME_UTR_SCORE;
				}
				if($self->FIVEPRIME_UTR_SPLICE_REGION_SCORE > $score &&
					$g->hasClassification($self->getSpliceRegionClass) && $g->hasClassification($self->get5PrimeUtrClass) &&
					$mrna->hasClassification($self->getSpliceRegionVariantClass) && $mrna->hasClassification($self->get5PrimeUtrVariantClass)){
						# splice region change in 5' UTR
						$score = $self->FIVEPRIME_UTR_SPLICE_REGION_SCORE;
				}
				if($self->THREEPRIME_UTR_SPLICE_REGION_SCORE > $score &&
					$g->hasClassification($self->getSpliceRegionClass) && $g->hasClassification($self->get3PrimeUtrClass) &&
					$mrna->hasClassification($self->getSpliceRegionVariantClass) && $mrna->hasClassification($self->get3PrimeUtrVariantClass)){
						# splice region change in 3' UTR
						$score = $self->THREEPRIME_UTR_SPLICE_REGION_SCORE;
				}
				if($self->COMPLEX_IN_MRNA_SCORE > $score &&
					$g->hasClassification($self->getCDSClass) &&
					$mrna->hasClassification($self->getComplexChangeVariantClass)){
					# complex transcript consequence involving only UTR (no CDS)
						if($g->hasClassification($self->getEssentialSpliceSiteClass) && $g->hasClassification($self->getExonClass) && $g->hasClassification($self->get5PrimeUtrClass)){
							# essential splice change in 5' UTR
							$score = $self->FIVEPRIME_UTR_ESS_SPLICE_SCORE;
						} elsif($g->hasClassification($self->getEssentialSpliceSiteClass) && $g->hasClassification($self->getExonClass) && $g->hasClassification($self->get3PrimeUtrClass)){
							# essential splice change in 3' UTR
							$score = $self->THREEPRIME_UTR_ESS_SPLICE_SCORE;
						} else {
							# if its none of the above, its just complex in mrna
							$score = $self->COMPLEX_IN_MRNA_SCORE;
						}
				}
				if($self->INTRONIC_SCORE > $score &&
					$g->hasClassification($self->getIntronClass) &&
					$mrna->hasClassification($self->getIntronVariantClass)){
						# intronic change
						$score = $self->INTRONIC_SCORE;
				}
				if($self->UPSTREAM_SCORE > $score &&
					($mrna->hasClassification($self->get2KBUpStreamVariantClass) || $mrna->hasClassification($self->get5KBUpStreamVariantClass))){
						# upstream of transcript
						$score = $self->UPSTREAM_SCORE;
				}
				if($self->DOWNSTREAM_SCORE > $score &&
					($mrna->hasClassification($self->get500BPDownStreamVariantClass) || $mrna->hasClassification($self->get5KBDownStreamVariantClass))){
						# downstream of transcript
						$score = $self->DOWNSTREAM_SCORE;
				}
			}

		} else {
			# non-coding transcript
			if($self->COMPLETE_NONCODING_TRANSCRIPT_LOSS_SCORE > $score &&
				$mrna->getMinPos() == 1 && $mrna->getMaxPos() == $mrna->getSequenceLength()){
				# if marked as a variant, start = 1 and end = sequence length the transcript is gone.
				$score = $self->COMPLETE_NONCODING_TRANSCRIPT_LOSS_SCORE;
			}
			if($self->NONCODING_GENE_ESS_SPLICE_SCORE > $score &&
				$g->hasClassification($self->getEssentialSpliceSiteClass) &&
				$mrna->hasClassification($self->getEssentialSpliceSiteVariantClass)){
					# essential splice change
					$score = $self->NONCODING_GENE_ESS_SPLICE_SCORE;
			}
			if($self->NONCODING_GENE_SPLICE_REGION_SCORE > $score &&
				$g->hasClassification($self->getSpliceRegionClass) &&
				$mrna->hasClassification($self->getSpliceRegionVariantClass)){
					# splice region change
					$score = $self->NONCODING_GENE_SPLICE_REGION_SCORE;
			}

			if($self->COMPLEX_IN_MRNA_SCORE > $score && $mrna->hasClassification($self->getComplexChangeVariantClass)){
				# complex transcript consequence
				if($g->hasClassification($self->getEssentialSpliceSiteClass) && $g->hasClassification($self->getExonClass)){
					# essential splice change
					$score = $self->NONCODING_GENE_ESS_SPLICE_SCORE;
				} else {
					# if its not ess splice, its just complex (probably straddles transcript boundary)
					$score = $self->COMPLEX_IN_MRNA_SCORE;
				}
			}
			if($self->NONCODING_GENE_SCORE > $score &&
				$g->hasClassification($self->getExonClass) &&
				$mrna->hasClassification($self->getNonCodingTranscriptVariantClass)){
					# exonic change
					$score = $self->NONCODING_GENE_SCORE;
			}
			if($self->INTRONIC_SCORE > $score &&
				$g->hasClassification($self->getIntronClass) &&
				$mrna->hasClassification($self->getIntronVariantClass)){
					# intronic change
					$score = $self->INTRONIC_SCORE;
			}
			if($self->UPSTREAM_SCORE > $score &&
				($mrna->hasClassification($self->get2KBUpStreamVariantClass) || $mrna->hasClassification($self->get5KBUpStreamVariantClass))){
					# upstream of transcript
					$score = $self->UPSTREAM_SCORE;
			}
			if($self->DOWNSTREAM_SCORE > $score &&
				($mrna->hasClassification($self->get500BPDownStreamVariantClass) || $mrna->hasClassification($self->get5KBDownStreamVariantClass))){
					# downstream of transcript
					$score = $self->DOWNSTREAM_SCORE;
			}
		}

		if($score == 1){
			croak("Unable to classify\n".Dumper($g));
		}
		if($score > $mostScore){
			$mostGroup = $g;
			$mostScore = $score;
		}
	}
	return $mostGroup;
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::Bookmarkers::MostDeleteriousBookmarker - Finds the L<AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> with the most severe consequence from the supplied list

=head1 DESCRIPTION

This bookmarker will mark/return the L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> with the most catastrophic consequence from the list of supplied groups.

It inherits from L<Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker|Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker> and L<Sanger::CGP::Vagrent::Bookmarkers::RepresentativeTranscriptBookmarker|Sanger::CGP::Vagrent::Bookmarkers::RepresentativeTranscriptBookmarker>

When handed an array of L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> objects it will sort them into the same order as
the L<RepresentativeTranscriptBookmarker|Sanger::CGP::Vagrent::Bookmarkers::RepresentativeTranscriptBookmarker> does and then prioritise AnnotationGroups by matching
their ontology terms to the following hierarchy.  The highest priority AnnotationGroup is considered to be the match

=over

=item 1 Complete protein loss

=item 2 Frameshift

=item 3 Stop gained/Nonsense

=item 4 Essential splice site change in CDS

=item 5 In-frame codon loss

=item 6 Simultaneous in-frame codon loss and gain

=item 7 In-frame codon gain

=item 8 Initiator/start codon change

=item 9 Terminator/stop codon change

=item 10 Non-synonymous/Missense

=item 11 Complex change in CDS

=item 12 Complete non-coding transcript loss

=item 13 Synonymous/Silent

=item 14 Premature start codon gained in 5'UTR

=item 15 Essential splice site change in 5'UTR

=item 16 Essential splice site change in 3'UTR

=item 17 Essential splice site change in non-coding transcript

=item 18 Exonic change in 3'UTR

=item 19 Exonic change in 5'UTR

=item 20 Splice region change in CDS

=item 21 Splice region change in 5'UTR

=item 22 Splice region change in 3'UTR

=item 23 Splice region change in non-coding transcript

=item 24 Complex change in mRNA

=item 25 Change in non-coding transcript

=item 26 Intronic change

=item 27 Upstream change

=item 28 Downstream change

=back

=head2 NOTE

The rank in the hierarchy is calculated for each L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> in isolation, and comparison between L<AnnotationGroups|Sanger::CGP::Vagrent::Data::AnnotationGroup> is based solely on the ranks. There is no consideration given to variants that are annotated to multiple genes and only the L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> with the most disruptive effect will be returned.

=head1 METHODS

see L<Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker|Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker>
