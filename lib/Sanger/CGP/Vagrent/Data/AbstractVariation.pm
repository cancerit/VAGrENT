package Sanger::CGP::Vagrent::Data::AbstractVariation;

use strict;

use base qw(Bio::Root::Root);

1;

sub isValid {
	shift->throw_not_implemented;
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::Data::AbstractVariation - Abstract Data object representing
a variation

=head1 DESCRIPTION

This is an abstract data class designed to be extended, it provides a common base class for all
variations.  Its also privides a placeholder function isValid that must me implemented in child
classes

=head1 METHODS

=head2 Functions

=head3 isValid

=over

=item Usage :

 if($var->isValid()){ ....... }

=item Function :

Abstract validation function, checks internal data, must be implemented by subclasses

=item Returns :

1 for pass, 0 for fail

=back
