package Sanger::CGP::VagrentSV::Data::AbstractVariation;

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

use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Bio::Root::Root);

1;

sub isValid {
	shift->throw_not_implemented;
}

sub getLength {
	return shift->{_length};
}

sub getSvType {
	return shift->{_svtype};
}

sub getName {
	return shift->{_name};
}

sub getLocFlag {
	return shift->{_locflag};
}

sub getLhb {
	return shift->{_lhb};
}
sub getRhb {
	return shift->{_rhb};
}

sub getInb {
	return shift->{_inb};
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
