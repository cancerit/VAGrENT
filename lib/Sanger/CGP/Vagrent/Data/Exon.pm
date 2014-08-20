package Sanger::CGP::Vagrent::Data::Exon;

use strict;
use Data::Dumper;
use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Sanger::CGP::Vagrent::Data::AbstractGenomicPosition);

1;

sub _init {
	my $self = shift;
	my %vars = @_;
	foreach my $k(keys(%vars)){
		if($k eq 'rnaminpos'){
			$self->{_rnaminpos} = $vars{rnaminpos};
		} elsif($k eq 'rnamaxpos'){
			$self->{_rnamaxpos} = $vars{rnamaxpos};
		}
	}
}

sub getRnaMinPos {
	return shift->{_rnaminpos};
}

sub getRnaMaxPos {
	return shift->{_rnamaxpos};
}

sub toString {
	my $self = shift;
	print $self->getChr.':'.$self->getMinPos.'..'.$self->getMaxPos;
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::Data::Exon - Data object representing an exon

=head1 DESCRIPTION

This is a data class to hold details about an exon. This allows the exon information to be
abstracted away from its original source.

It inherits from L<Sanger::CGP::Vagrent::Data::AbstractGenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition>

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $exon = Sanger::CGP::Vagrent::Data::Exon->new(%params);

=item Function :

Builds a new Sanger::CGP::Vagrent::Data::Exon object

=item Returns :

Sanger::CGP::Vagrent::Data::Exon object initialized with parameter values

=item Params :

Same as the constructor from L<Sanger::CGP::Vagrent::Data::AbstractGenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition> plus

 rnaminpos => minimum position of the exon in the mRNA
 rnamaxpos => maximum position of the exon in the mRNA

=back

=head2 Attributes

=head3 getRnaMinPos

=over

=item Usage :

 my $rnaMin = $exon->getRnaMinPos;

=item Function :

Returns the exons minimum position in the rna sequence

=item Returns :

Integer

=back

=head3 getRnaMaxPos

=over

=item Usage :

 my $rnaMax = $exon->getRnaMaxPos;

=item Function :

Returns the exons maximum position in the rna sequence

=item Returns :

Integer

=back

=head2 Functions

=head3 toString

=over

=item Usage :

 print $exon->toString;

=item Function :

Returns a simple string representation of the genomic view of the exon (chr:min..max)

=item Returns :

String

=back
