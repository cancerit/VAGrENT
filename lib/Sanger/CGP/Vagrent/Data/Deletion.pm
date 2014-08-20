package Sanger::CGP::Vagrent::Data::Deletion;

use strict;
use Data::Dumper;
use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Sanger::CGP::Vagrent::Data::AbstractGenomicPosition Sanger::CGP::Vagrent::Data::AbstractVariation);

1;

sub _init {
	my $self = shift;
	my %vars = @_;
	foreach my $k(keys(%vars)){
		if($k eq 'delseq'){
			$self->{_delseq} = $vars{delseq};
		}
	}
}

sub isValid {
	my $self = shift;

	return 0 unless(defined($self->{_minpos}) && defined($self->{_maxpos}));
	return 0 unless($self->{_minpos} <= $self->{_maxpos});
	return 0 unless($self->{_minpos} > 0);

	return 0 unless(defined($self->{_delseq}));
	return 0 unless(length($self->{_delseq}) > 0);
	return 0 unless($self->{_delseq} =~ m/^[atcgn]+$/i);

	return 0 if(length($self->{_delseq}) != ($self->{_maxpos} - $self->{_minpos}) + 1);

	return 1;
}

sub getDeletedSequence {
	return shift->{_delseq};
}

sub toString {
	my $self = shift;
	my $out = 'chr'.$self->getChr.':g.';
	if($self->getMinPos == $self->getMaxPos){
		$out .= $self->getMinPos.'del'.$self->getDeletedSequence;
	} else {
		$out .= $self->getMinPos.'_'.$self->getMaxPos.'del';
		if(length($self->getDeletedSequence) <= 10){
			$out .= $self->getDeletedSequence;
		} else {
			$out .= length($self->getDeletedSequence);
		}
	}
	return $out;
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::Data::Deletion - Data object representing a deletion

=head1 DESCRIPTION

This is a data class describing a deletion plotted to a genome.

It inherits from L<Sanger::CGP::Vagrent::Data::AbstractGenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition> and L<Sanger::CGP::Vagrent::Data::AbstractVariation|Sanger::CGP::Vagrent::Data::AbstractVariation>

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $del = Sanger::CGP::Vagrent::Data::Deletion->new(%params);

=item Function :

Builds a new Sanger::CGP::Vagrent::Data::Deletion object

=item Returns :

Sanger::CGP::Vagrent::Data::Deletion object initialized with parameter values

=item Params :

Same as the constructor from L<Sanger::CGP::Vagrent::Data::AbstractGenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition> plus

 delseq => the deleted sequence fragment

=back

=head2 Attributes

=head3 getDeletedSequence

=over

=item Usage :

 my $seq = $del->getDeletedSequence;

=item Function :

Returns the deleted sequence fragment

=item Returns :

String - DNA sequence

=back

=head2 Functions

=head3 toString

=over

=item Usage :

 print $variant->toString;

=item Function :

Returns a simple string representation of the variant in hgvs genomic syntax

=item Returns :

String

=back
