package Sanger::CGP::Vagrent::Data::GenomicRegion;

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
use base qw(Sanger::CGP::Vagrent::Data::AbstractGenomicPosition);

1;

sub _init {
	my $self = shift;
	$self->throw('must specify a minpos lower than a max pos') unless($self->{_minpos} <= $self->{_maxpos});
	$self->throw('coordinates are 1 based, a start position of 0 is not allowed') if($self->{_minpos} == 0);
}

sub toString {
	my $self = shift;
	return $self->getChr.':'.$self->getMinPos.'-'.$self->getMaxPos;
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::Data::GenomicRegion - Data object representing a generic genomic region

=head1 DESCRIPTION

This is a data class to hold details about an simple genomic region.

It inherits from L<Sanger::CGP::Vagrent::Data::AbstractGenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition> and only implements internal constructor data validation, no extra methods.
