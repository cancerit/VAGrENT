package Sanger::CGP::Vagrent::Bookmarkers::MostDeleteriousBookmarker;

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

use Log::Log4perl qw(:easy);
use Data::Dumper;
use Carp qw(croak);
use Const::Fast qw(const);

use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Sanger::CGP::Vagrent::Bookmarkers::RepresentativeTranscriptBookmarker);

const my $DOWNSTREAM_SCORE => 4;
const my $UPSTREAM_SCORE => 7;
const my $INTRONIC_SCORE => 10;
const my $NONCODING_GENE_SCORE => 15;
const my $COMPLEX_IN_MRNA_SCORE => 50;
const my $NONCODING_GENE_SPLICE_REGION_SCORE => 95;
const my $THREEPRIME_UTR_SPLICE_REGION_SCORE => 100;
const my $FIVEPRIME_UTR_SPLICE_REGION_SCORE => 105;
const my $CODING_SPLICE_REGION_SCORE => 200;
const my $FIVEPRIME_UTR_SCORE => 300;
const my $THREEPRIME_UTR_SCORE => 305;
const my $NONCODING_GENE_ESS_SPLICE_SCORE => 395;
const my $THREEPRIME_UTR_ESS_SPLICE_SCORE => 400;
const my $FIVEPRIME_UTR_ESS_SPLICE_SCORE => 405;
const my $START_GAINED_SCORE => 450;
const my $SYNONYMOUS_SCORE => 500;
const my $COMPLETE_NONCODING_TRANSCRIPT_LOSS_SCORE => 525;
const my $COMPLEX_IN_CDS_SCORE => 550;
const my $NON_SYNONYMOUS_SCORE => 600;
const my $STOP_LOST_SCORE => 700;
const my $INITIATOR_CHANGE_SCORE => 800;
const my $INFRAME_CODON_GAIN_SCORE => 825;
const my $INFRAME_CODON_LOSS_AND_GAIN_SCORE => 840;
const my $INFRAME_CODON_LOSS_SCORE => 850;
const my $CODING_ESS_SPLICE_SCORE => 900;
const my $STOP_GAINED_SCORE => 1000;
const my $FRAMESHIFT_SCORE => 1100;
const my $COMPLETE_PROTEIN_LOSS_SCORE => 1200;

1;

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
		my $mrna = $g->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation->getmRNAAnnotationContext);
		my $cds = $g->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation->getCDSAnnotationContext);
		my $prot = $g->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation->getProteinAnnotationContext);
		if($g->hasClassification($self->getProteinCodingClass)){
			# protein coding transcript
			if(defined($prot) && $prot->getType() ne $mrna->getUnknownAnnotationType){
				# Protein annotation
				if($COMPLETE_PROTEIN_LOSS_SCORE > $score &&
					$prot->hasClassification($self->getDeletionClass) &&
					$prot->getMinPos() == 1 && $prot->getMaxPos() == $prot->getSequenceLength()){
					# if marked as a deletion, start = 1 and end = protein length the protein is gone.
						$score = $COMPLETE_PROTEIN_LOSS_SCORE;
				}
				if($FRAMESHIFT_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$prot->hasClassification($self->getFrameShiftVariantClass)){
					# Frameshift
						$score = $FRAMESHIFT_SCORE;
				}
				if($STOP_GAINED_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$prot->hasClassification($self->getStopGainedVariantClass)){
					# non-sense / stop gained
						$score = $STOP_GAINED_SCORE;
				}
				if($INFRAME_CODON_LOSS_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$prot->hasClassification($self->getInFrameCodonLossVariantClass)){
					# in frame deletion
						$score = $INFRAME_CODON_LOSS_SCORE;
				}
				if($INFRAME_CODON_LOSS_AND_GAIN_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$prot->hasClassification($self->getComplexIndelClass)){
					# in frame complex sub
						$score = $INFRAME_CODON_LOSS_AND_GAIN_SCORE;
				}
				if($INFRAME_CODON_GAIN_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$prot->hasClassification($self->getInFrameCodonGainVariantClass)){
					# in frame insertion
						$score = $INFRAME_CODON_GAIN_SCORE;
				}
				if($INITIATOR_CHANGE_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$prot->hasClassification($self->getStartLostVariantClass)){
					# start lost
						$score = $INITIATOR_CHANGE_SCORE;
				}
				if($STOP_LOST_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$prot->hasClassification($self->getStopLostVariantClass)){
					# stop lost
						$score = $STOP_LOST_SCORE;
				}
				if($NON_SYNONYMOUS_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$prot->hasClassification($self->getNonSynonymousVariantClass)){
					# mis sense
						$score = $NON_SYNONYMOUS_SCORE;
				}
				if($SYNONYMOUS_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					($prot->hasClassification($self->getSynonymousVariantClass) || $prot->hasClassification($self->getStopRetainedVariantClass))){
					# silent including terminator silent
						$score = $SYNONYMOUS_SCORE;
				}
			} elsif(defined($cds) && $cds->getType() ne $mrna->getUnknownAnnotationType){
				# CDS annotation
				if($FRAMESHIFT_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
					$cds->hasClassification($self->getFrameShiftVariantClass)){
					# frame shift again, incase protein translation was too complex
						$score = $FRAMESHIFT_SCORE;
				}
				if($CODING_ESS_SPLICE_SCORE > $score &&
					$g->hasClassification($self->getEssentialSpliceSiteClass) && $g->hasClassification($self->getCDSClass) &&
					$cds->hasClassification($self->getEssentialSpliceSiteVariantClass)){
					# essential splice change in CDS
						$score = $CODING_ESS_SPLICE_SCORE;
				}

				if($COMPLEX_IN_CDS_SCORE > $score &&
					$g->hasClassification($self->getCDSClass) &&
					$cds->hasClassification($self->getComplexChangeVariantClass)){
					# complex transcript consequence involving CDS
						if($cds->getMinPos() == 1 && $cds->getMaxPos() == $cds->getSequenceLength()){
							# position 1 to CDS length effected, transcript lost
							$score = $COMPLETE_PROTEIN_LOSS_SCORE;
						} elsif($g->hasClassification($self->getEssentialSpliceSiteClass) && $g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass)){
							# essential splice change
							$score = $CODING_ESS_SPLICE_SCORE;
						} elsif($g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
									($g->hasClassification($self->get5PrimeUtrClass) || $mrna->hasClassification($self->get2KBUpStreamVariantClass)) &&
									$cds->getMinPos() == 1){
							# start codon lost
							$score = $INITIATOR_CHANGE_SCORE;
						} elsif($g->hasClassification($self->getExonClass) && $g->hasClassification($self->getCDSClass) &&
									($g->hasClassification($self->get3PrimeUtrClass) || $mrna->hasClassification($self->get500BPDownStreamVariantClass)) &&
									$cds->getMaxPos() == $cds->getSequenceLength()){
							# stop codon lost
							$score = $STOP_LOST_SCORE;
						} else {
							# if its none of the above, its just complex in CDS
							$score = $COMPLEX_IN_CDS_SCORE;
						}
				}
				if($CODING_SPLICE_REGION_SCORE > $score &&
					$g->hasClassification($self->getSpliceRegionClass) && $g->hasClassification($self->getCDSClass) &&
					$cds->hasClassification($self->getSpliceRegionVariantClass)){
					# splice region change in CDS
						$score = $CODING_SPLICE_REGION_SCORE;
				}
        if($FIVEPRIME_UTR_ESS_SPLICE_SCORE > $score &&
					$g->hasClassification($self->getEssentialSpliceSiteClass) && $g->hasClassification($self->get5PrimeUtrClass) &&
					$mrna->hasClassification($self->getEssentialSpliceSiteVariantClass) && $mrna->hasClassification($self->get5PrimeUtrVariantClass)){
            # essential splice change in 5' UTR, the splice site is probably directly before to the start codon
            $score = $FIVEPRIME_UTR_ESS_SPLICE_SCORE;
        }
        if($THREEPRIME_UTR_ESS_SPLICE_SCORE > $score &&
					$g->hasClassification($self->getEssentialSpliceSiteClass) && $g->hasClassification($self->get3PrimeUtrClass) &&
					$mrna->hasClassification($self->getEssentialSpliceSiteVariantClass) && $mrna->hasClassification($self->get3PrimeUtrVariantClass)){
						# essential splice change in 3' UTR, the splice site is probably directly after to the stop codon
						$score = $THREEPRIME_UTR_ESS_SPLICE_SCORE;
				}
        if($FIVEPRIME_UTR_SPLICE_REGION_SCORE > $score &&
					$g->hasClassification($self->getSpliceRegionClass) && $g->hasClassification($self->get5PrimeUtrClass) &&
					$mrna->hasClassification($self->getSpliceRegionVariantClass) && $mrna->hasClassification($self->get5PrimeUtrVariantClass)){
						# splice region change in 5' UTR, the splice site is probably directly before to the start codon
						$score = $FIVEPRIME_UTR_SPLICE_REGION_SCORE;
				}
				if($THREEPRIME_UTR_SPLICE_REGION_SCORE > $score &&
					$g->hasClassification($self->getSpliceRegionClass) && $g->hasClassification($self->get3PrimeUtrClass) &&
					$mrna->hasClassification($self->getSpliceRegionVariantClass) && $mrna->hasClassification($self->get3PrimeUtrVariantClass)){
						# splice region change in 3' UTR, the splice site is probably directly after to the stop codon
						$score = $THREEPRIME_UTR_SPLICE_REGION_SCORE;
				}
			} elsif(defined($mrna)) {
				# cDNA annotation
				if($START_GAINED_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->get5PrimeUtrClass) &&
					$mrna->hasClassification($self->getPrematureStartGainedVariantClass)){
						# new start codon created in the 5' UTR
						$score = $START_GAINED_SCORE;
				}
				if($FIVEPRIME_UTR_ESS_SPLICE_SCORE > $score &&
					$g->hasClassification($self->getEssentialSpliceSiteClass) && $g->hasClassification($self->get5PrimeUtrClass) &&
					$mrna->hasClassification($self->getEssentialSpliceSiteVariantClass) && $mrna->hasClassification($self->get5PrimeUtrVariantClass)){
						# essential splice change in 5' UTR
						$score = $FIVEPRIME_UTR_ESS_SPLICE_SCORE;
				}
				if($THREEPRIME_UTR_ESS_SPLICE_SCORE > $score &&
					$g->hasClassification($self->getEssentialSpliceSiteClass) && $g->hasClassification($self->get3PrimeUtrClass) &&
					$mrna->hasClassification($self->getEssentialSpliceSiteVariantClass) && $mrna->hasClassification($self->get3PrimeUtrVariantClass)){
						# essential splice change in 3' UTR
						$score = $THREEPRIME_UTR_ESS_SPLICE_SCORE;
				}
				if($THREEPRIME_UTR_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->get3PrimeUtrClass) &&
					$mrna->hasClassification($self->get3PrimeUtrVariantClass)){
						# change in 3' UTR exon
						$score = $THREEPRIME_UTR_SCORE;
				}
				if($FIVEPRIME_UTR_SCORE > $score &&
					$g->hasClassification($self->getExonClass) && $g->hasClassification($self->get5PrimeUtrClass) &&
					$mrna->hasClassification($self->get5PrimeUtrVariantClass)){
						# change in 5' UTR exon
						$score = $FIVEPRIME_UTR_SCORE;
				}
				if($FIVEPRIME_UTR_SPLICE_REGION_SCORE > $score &&
					$g->hasClassification($self->getSpliceRegionClass) && $g->hasClassification($self->get5PrimeUtrClass) &&
					$mrna->hasClassification($self->getSpliceRegionVariantClass) && $mrna->hasClassification($self->get5PrimeUtrVariantClass)){
						# splice region change in 5' UTR
						$score = $FIVEPRIME_UTR_SPLICE_REGION_SCORE;
				}
				if($THREEPRIME_UTR_SPLICE_REGION_SCORE > $score &&
					$g->hasClassification($self->getSpliceRegionClass) && $g->hasClassification($self->get3PrimeUtrClass) &&
					$mrna->hasClassification($self->getSpliceRegionVariantClass) && $mrna->hasClassification($self->get3PrimeUtrVariantClass)){
						# splice region change in 3' UTR
						$score = $THREEPRIME_UTR_SPLICE_REGION_SCORE;
				}
				if($COMPLEX_IN_MRNA_SCORE > $score &&
					$g->hasClassification($self->getCDSClass) &&
					$mrna->hasClassification($self->getComplexChangeVariantClass)){
					# complex transcript consequence involving only UTR (no CDS)
						if($g->hasClassification($self->getEssentialSpliceSiteClass) && $g->hasClassification($self->getExonClass) && $g->hasClassification($self->get5PrimeUtrClass)){
							# essential splice change in 5' UTR
							$score = $FIVEPRIME_UTR_ESS_SPLICE_SCORE;
						} elsif($g->hasClassification($self->getEssentialSpliceSiteClass) && $g->hasClassification($self->getExonClass) && $g->hasClassification($self->get3PrimeUtrClass)){
							# essential splice change in 3' UTR
							$score = $THREEPRIME_UTR_ESS_SPLICE_SCORE;
						} else {
							# if its none of the above, its just complex in mrna
							$score = $COMPLEX_IN_MRNA_SCORE;
						}
				}
				if($INTRONIC_SCORE > $score &&
					$g->hasClassification($self->getIntronClass) &&
					$mrna->hasClassification($self->getIntronVariantClass)){
						# intronic change
						$score = $INTRONIC_SCORE;
				}
				if($UPSTREAM_SCORE > $score &&
					($mrna->hasClassification($self->get2KBUpStreamVariantClass) || $mrna->hasClassification($self->get5KBUpStreamVariantClass))){
						# upstream of transcript
						$score = $UPSTREAM_SCORE;
				}
				if($DOWNSTREAM_SCORE > $score &&
					($mrna->hasClassification($self->get500BPDownStreamVariantClass) || $mrna->hasClassification($self->get5KBDownStreamVariantClass))){
						# downstream of transcript
						$score = $DOWNSTREAM_SCORE;
				}
			}

		} else {
			# non-coding transcript
			if($COMPLETE_NONCODING_TRANSCRIPT_LOSS_SCORE > $score &&
				$mrna->getMinPos() == 1 && $mrna->getMaxPos() == $mrna->getSequenceLength()){
				# if marked as a variant, start = 1 and end = sequence length the transcript is gone.
				$score = $COMPLETE_NONCODING_TRANSCRIPT_LOSS_SCORE;
			}
			if($NONCODING_GENE_ESS_SPLICE_SCORE > $score &&
				$g->hasClassification($self->getEssentialSpliceSiteClass) &&
				$mrna->hasClassification($self->getEssentialSpliceSiteVariantClass)){
					# essential splice change
					$score = $NONCODING_GENE_ESS_SPLICE_SCORE;
			}
			if($NONCODING_GENE_SPLICE_REGION_SCORE > $score &&
				$g->hasClassification($self->getSpliceRegionClass) &&
				$mrna->hasClassification($self->getSpliceRegionVariantClass)){
					# splice region change
					$score = $NONCODING_GENE_SPLICE_REGION_SCORE;
			}

			if($COMPLEX_IN_MRNA_SCORE > $score && $mrna->hasClassification($self->getComplexChangeVariantClass)){
				# complex transcript consequence
				if($g->hasClassification($self->getEssentialSpliceSiteClass) && $g->hasClassification($self->getExonClass)){
					# essential splice change
					$score = $NONCODING_GENE_ESS_SPLICE_SCORE;
				} else {
					# if its not ess splice, its just complex (probably straddles transcript boundary)
					$score = $COMPLEX_IN_MRNA_SCORE;
				}
			}
			if($NONCODING_GENE_SCORE > $score &&
				$g->hasClassification($self->getExonClass) &&
				$mrna->hasClassification($self->getNonCodingTranscriptVariantClass)){
					# exonic change
					$score = $NONCODING_GENE_SCORE;
			}
			if($INTRONIC_SCORE > $score &&
				$g->hasClassification($self->getIntronClass) &&
				$mrna->hasClassification($self->getIntronVariantClass)){
					# intronic change
					$score = $INTRONIC_SCORE;
			}
			if($UPSTREAM_SCORE > $score &&
				($mrna->hasClassification($self->get2KBUpStreamVariantClass) || $mrna->hasClassification($self->get5KBUpStreamVariantClass))){
					# upstream of transcript
					$score = $UPSTREAM_SCORE;
			}
			if($DOWNSTREAM_SCORE > $score &&
				($mrna->hasClassification($self->get500BPDownStreamVariantClass) || $mrna->hasClassification($self->get5KBDownStreamVariantClass))){
					# downstream of transcript
					$score = $DOWNSTREAM_SCORE;
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
