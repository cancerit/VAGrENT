package Sanger::CGP::Vagrent::Data::Substitution;

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
		if($k eq 'wt'){
			$self->{_wt} = $vars{wt};
		} elsif($k eq 'mt'){
			$self->{_mt} = $vars{mt};
		}
	}
}

sub isValid {
	my $self = shift;
	return 0 unless(defined($self->{_wt}) && defined($self->{_mt}));
	return 0 unless(length($self->{_wt}) == 1 && length($self->{_mt}) == 1);
	return 0 unless($self->{_wt} =~ m/[atcg]/i && $self->{_mt} =~ m/[atcg]/i);
	return 0 if($self->{_mt} eq $self->{_wt});

	return 0 unless(defined($self->{_minpos}) && defined($self->{_maxpos}));
	return 0 unless($self->{_minpos} == $self->{_maxpos});
	return 0 unless($self->{_minpos} > 0);

	return 1;
}

sub getWt {
	return shift->{_wt};
}

sub getMt {
	return shift->{_mt};
}

sub toString {
	my $self = shift;
	return 'chr'.$self->getChr.':g.'.$self->getMinPos.$self->getWt.'>'.$self->getMt;
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::Data::Substitution - Data object representing a substitution

=head1 DESCRIPTION

This is a data class describing a substitution variant plotted to a genome.

It inherits from L<Sanger::CGP::Vagrent::Data::AbstractGenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition> and L<Sanger::CGP::Vagrent::Data::AbstractVariation|Sanger::CGP::Vagrent::Data::AbstractVariation>

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $cplx = Sanger::CGP::Vagrent::Data::Substitution->new(%params);

=item Function :

Builds a new Sanger::CGP::Vagrent::Data::Substitution object

=item Returns :

Sanger::CGP::Vagrent::Data::Substitution object initialized with parameter values

=item Params :

Same as the constructor from L<Sanger::CGP::Vagrent::Data::AbstractGenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition> plus

 wt => the wildtype sequence
 mt => the mutant/variant sequence

=back

=head2 Attributes

=head3 getWt

=over

=item Usage :

 my $seq = $sub->getWt;

=item Function :

Returns the wildtype sequence string

=item Returns :

String - DNA sequence

=back

=head3 getMt

=over

=item Usage :

 my $seq = $sub->getMt;

=item Function :

Returns the mutated/variant sequence string

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
