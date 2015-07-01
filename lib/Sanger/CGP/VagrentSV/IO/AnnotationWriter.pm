package Sanger::CGP::Vagrent::IO::AnnotationWriter;

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
use base qw(Sanger::CGP::Vagrent);

1;

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
	my $self = {};
	bless($self, $class);
	$self->_init(@_);
	return $self;
}

sub _init {
	my $self = shift;
  my %vars = @_;
	foreach my $k(keys(%vars)){
		if($k eq 'fh'){
			$self->{_fh} = $vars{'fh'};
		}
	}
}

sub write {
	shift->throw_not_implemented;
}

sub _fh {
	return shift->{_fh};
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::IO::AnnotationWriter - Abstract class for annotation writers

=head1 DESCRIPTION

Abstract class to hold basic shared functionality for annotation writers

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $source = Sanger::CGP::Vagrent::IO::AnnotationWriterSubClass->new();

=item Function :

Builds a new Sanger::CGP::Vagrent::IO::AnnotationWriter inheriting object

=item Returns :

Sanger::CGP::Vagrent::IO::AnnotationWriter object initialized with parameter values

=item Params :

 fh => File handle to write to

=back

=head2 Abstract

=head3 write

=over

=item Usage :

 $writer->write($var,$annoGroup);

=item Function :

Abstract function, must be overwritten. Writes the annotation out to the file

=item Returns :

Nothing

=item Params :

A L<Sanger::CGP::Vagrent::Data::AbstractVariation|Sanger::CGP::Vagrent::Data::AbstractVariation> inheriting object

A L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> object

=back

