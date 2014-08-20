package Sanger::CGP::Vagrent::Bookmarkers::RepresentativeTranscriptBookmarker;

use strict;

use Log::Log4perl;
use Data::Dumper;

use Sanger::CGP::Vagrent qw($VERSION);
use Sanger::CGP::Vagrent::Data::Annotation;

use base qw(Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker);

1;

sub _getAnnotation {
	my $self = shift;
	my @groups = @_;
	my @out = $self->_sortAnnotations(@groups);
	return $out[0];
}

sub _sortAnnotations {
	my ($self,@groups) = @_;

	my @out = sort{
				my $a_ccds = 0;
				my $b_ccds = 0;
				if(defined($a->getCCDS) && $a->getCCDS ne ''){
					$a_ccds = 1;
				}
				if(defined($b->getCCDS) && $b->getCCDS ne ''){
					$b_ccds = 1;
				}
				my $ccds_cmp = $b_ccds <=> $a_ccds;
				if($ccds_cmp == 0){
					my $a_cds_len = 0;
					my $b_cds_len = 0;
					my $a_cds = $a->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext);
					my $b_cds = $b->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext);
					if(defined($a_cds)){
						$a_cds_len = $a_cds->getSequenceLength;
					}
					if(defined($b_cds)){
						$b_cds_len = $b_cds->getSequenceLength;
					}
					my $cds_len_cmp = $b_cds_len <=> $a_cds_len;
					if($cds_len_cmp == 0){
						return $b->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext)->getSequenceLength <=> $a->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext)->getSequenceLength;
					} else {
						return $cds_len_cmp;
					}
				} else {
					return $ccds_cmp;
				}
					} @groups;
	return @out;
}


sub _annotationSort {
	my $a_ccds = 0;
	my $b_ccds = 0;
	if(defined($a->getCCDS) && $a->getCCDS ne ''){
		$a_ccds = 1;
	}
	if(defined($b->getCCDS) && $b->getCCDS ne ''){
		$b_ccds = 1;
	}
	my $ccds_cmp = $b_ccds <=> $a_ccds;
	if($ccds_cmp == 0){
		my $a_cds_len = 0;
		my $b_cds_len = 0;
		my $a_cds = $a->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext);
		my $b_cds = $b->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext);
		if(defined($a_cds)){
			$a_cds_len = $a_cds->getSequenceLength;
		}
		if(defined($b_cds)){
			$b_cds_len = $b_cds->getSequenceLength;
		}
		my $cds_len_cmp = $b_cds_len <=> $a_cds_len;
		if($cds_len_cmp == 0){
			return $b->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext)->getSequenceLength <=> $a->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext)->getSequenceLength;
		} else {
			return $cds_len_cmp;
		}
	} else {
		return $ccds_cmp;
	}
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::Bookmarkers::RepresentativeTranscriptBookmarker - Finds the L<AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> on the best, most representative Transcript from the supplied list.

=head1 DESCRIPTION

This bookmarker will mark/return the L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> linked to the transcript thats the most reliable and complete representation of a gene.

It inherits from L<Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker|Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker>.

When handed an array of L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> objects it will sort them using the following criteria

=over

=item * Presence of a CCDS identifier (having one is better)

=item * CDS length (longer is better)

=item * cDNA length (longer is better)

=back

The top answer is then returned as the matching L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup>.

=head2 NOTE

This is calculated using the transcript information only, there is no consideration given to variants that are annotated to multiple genes.  Only a single L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> will be returned even if the input data contained annotation against multiple genes.

=head1 METHODS

see L<Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker|Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker>
