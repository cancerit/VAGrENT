package Sanger::CGP::Vagrent::Data::AbstractGenomicPosition;

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
use Carp;
use Attribute::Abstract;

use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Sanger::CGP::Vagrent);

1;

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
	my $self = {};
	bless($self, $class);
	$self->_localInit(@_);
	$self->_init(@_);
	return $self;
}

sub _localInit {
	my $self = shift;
	my %vars = @_;
	foreach my $k(keys(%vars)){
		if($k eq 'species'){
			$self->{_species} = $vars{species};
		} elsif($k eq 'genomeVersion'){
			$self->{_genomeVersion} = $vars{genomeVersion};
		} elsif($k eq 'chr'){
			$self->{_chr} = $vars{'chr'};
		} elsif($k eq 'minpos'){
			$self->{_minpos} = $vars{minpos};
		} elsif($k eq 'maxpos'){
			$self->{_maxpos} = $vars{maxpos};
		} elsif($k eq 'id'){
			$self->{_id} = $vars{id};
		}
	}
	croak('must specify a species') unless(exists($self->{_species}) && defined($self->{_species}));
	croak('must specify a genomeVersion') unless(exists($self->{_genomeVersion}) && defined($self->{_genomeVersion}));
	croak('must specify a chr') unless(exists($self->{_chr}) && defined($self->{_chr}));
	croak('must specify a minpos') unless(exists($self->{_minpos}) && defined($self->{_minpos}));
	croak('must specify a maxpos') unless(exists($self->{_maxpos}) && defined($self->{_maxpos}));
}

sub _init: Abstract;

sub getSpecies {
	return shift->{_species};
}

sub getGenomeVersion {
	return shift->{_genomeVersion};
}

sub getChr {
	return shift->{_chr};
}

sub getMinPos {
	return shift->{_minpos};
}

sub setMinPos {
	my ($self,$pos) = @_;
	if(defined $pos){
		$self->{_minpos} = $pos;
	}
}

sub getMaxPos {
	return shift->{_maxpos};
}

sub setMaxPos {
	my ($self,$pos) = @_;
	if(defined $pos){
		$self->{_maxpos} = $pos;
	}
}

sub getId {
	return shift->{_id};
}

sub getLength {
  my $self = shift;
  return (($self->{_maxpos} - $self->{_minpos}) + 1);
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::Data::AbstractGenomicPosition - Abstract Data object representing
Genomic position

=head1 DESCRIPTION

This is an abstract data class designed to be extended, it provides basic functionality for holding
genomic position.  Child classes must implement an _init method to handle parameter values handed to
the object when the constructor is called.

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $gp = Sanger::CGP::Vagrent::Data::AbstractGenomicPositionSubClass->new(%params);

=item Function :

Builds a new Sanger::CGP::Vagrent::Data::AbstractGenomicPosition inheriting object

=item Returns :

L<Sanger::CGP::Vagrent::Data::AbstractGenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition> sub class object initialized with parameter values

=item Params :

 species => Species string (eg Human)
 genomeVersion => Genome version String (eg GRCh37)
 chr => Chromosome/contig name
 minpos => Lowest coordinate of the feature
 maxpos => Highest coordinate of the feature
 id => Identifier (Optional)

=back

=head3 _init

=over

=item Usage :

Abstract internal initialisation method, must be implemented in subclasses. All parameters are passed through from the constructor.

=back

=head2 Attributes

=head3 getId

=over

=item Usage :

 my $id = $gp->getId;

=item Function :

Returns the value of Id

=item Returns :

String

=back

=head3 getSpecies

=over

=item Usage :

 my $species = $gp->getSpecies;

=item Function :

Returns the species name

=item Returns :

String

=back

=head3 getGenomeVersion

=over

=item Usage :

 my $gVers = $gp->getGenomeVersion;

=item Function :

Returns the genome version string

=item Returns :

String

=back

=head3 getChr

=over

=item Usage :

 my $chr = $gp->getChr;

=item Function :

Returns the chromosome/contig name

=item Returns :

String

=back

=head3 getMinPos

=over

=item Usage :

 my $min = $gp->getMinPos;

=item Function :

Returns the lowest coordinate of the feature on the sequence

=item Returns :

Integer

=back

=head3 setMinPos

=over

=item Usage :

 $gp->setMinPos($newPos);

=item Function :

Sets the lowest coordinate of the feature on the sequence

=item Params

Integer, new position

=item Returns :

None

=back

=head3 getMaxPos

=over

=item Usage :

 my $max = $gp->getMaxPos;

=item Function :

Returns the highest coordinate of the feature on the sequence

=item Returns :

Integer

=back

=head3 setMaxPos

=over

=item Usage :

 $gp->setMaxPos($newPos);

=item Function :

Sets the highest coordinate of the feature on the sequence

=item Params

Integer, new position

=item Returns :

None

=back
