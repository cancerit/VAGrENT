package Sanger::CGP::Vagrent::Data::GenomicRegion;

use strict;

use base qw(Sanger::CGP::Vagrent::Data::AbstractGenomicPosition);

1;

sub _init {
	my $self = shift;
	$self->throw('must specify a minpos lower than a max pos') unless($self->{_minpos} <= $self->{_maxpos});
	$self->throw('coordinates are 1 based, a start position of 0 is not allowed') if($self->{_minpos} == 0);
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::Data::GenomicRegion - Data object representing a generic genomic region

=head1 DESCRIPTION

This is a data class to hold details about an simple genomic region.

It inherits from L<Sanger::CGP::Vagrent::Data::AbstractGenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition> and only implements internal constructor data validation, no extra methods.
