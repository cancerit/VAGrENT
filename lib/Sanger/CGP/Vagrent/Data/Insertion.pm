package Sanger::CGP::Vagrent::Data::Insertion;

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
use Data::Dumper;
use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Sanger::CGP::Vagrent::Data::AbstractGenomicPosition Sanger::CGP::Vagrent::Data::AbstractVariation);

1;

sub _init {
	my $self = shift;
	my %vars = @_;
	foreach my $k(keys(%vars)){
		if($k eq 'insseq'){
			$self->{_insseq} = $vars{insseq};
		}
	}
}

sub isValid {
	my $self = shift;
	return 0 if(!defined($self->{_minpos}) || !defined($self->{_maxpos}) || $self->{_minpos} + 1 != $self->{_maxpos} || $self->{_minpos} < 1);
	return 0 if(!defined($self->{_insseq}) || length($self->{_insseq}) < 1 || $self->{_insseq} !~ m/^[atcgn]+$/i);
	return 1;
}

sub getInsertedSequence {
	return shift->{_insseq};
}

sub toString {
	my $self = shift;
	my $out = 'chr'.$self->getChr.':g.'.$self->getMinPos.'_'.$self->getMaxPos.'ins'.$self->getInsertedSequence;
	return $out;
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::Data::Insertion - Data object representing an insertion event

=head1 DESCRIPTION

This is a data class describing an insertion variant plotted to a genome.

It inherits from L<Sanger::CGP::Vagrent::Data::AbstractGenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition> and L<Sanger::CGP::Vagrent::Data::AbstractVariation|Sanger::CGP::Vagrent::Data::AbstractVariation>

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $ins = Sanger::CGP::Vagrent::Data::Insertion->new(%params);

=item Function :

Builds a new Sanger::CGP::Vagrent::Data::Insertion object

=item Returns :

Sanger::CGP::Vagrent::Data::Insertion object initialized with parameter values

=item Params :

Same as the constructor from L<Sanger::CGP::Vagrent::Data::AbstractGenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition> plus

 insseq => the inserted sequence fragment

=back

=head2 Attributes

=head3 getInsertedSequence

=over

=item Usage :

 my $seq = $ins->getInsertedSequence;

=item Function :

Returns the inserted sequence fragment

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

