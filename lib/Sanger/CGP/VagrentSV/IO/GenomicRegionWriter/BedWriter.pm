package Sanger::CGP::Vagrent::IO::GenomicRegionWriter::BedWriter;

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
use base qw(Sanger::CGP::Vagrent::IO::GenomicRegionWriter);

1;

my $seperator = "\t";

sub write {
	my ($self,$gr) = @_;
	unless($gr->isa('Sanger::CGP::Vagrent::Data::GenomicRegion')){
		warn 'expecting an Sanger::CGP::Vagrent::Data::GenomicRegion, received a' . ref($gr);
		return;
	}
	print {$self->_fh} join($seperator,$gr->getChr,($gr->getMinPos - 1),$gr->getMaxPos,),"\n";
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::IO::GenomicRegionWriter::BedWriter - class for writing genomic regions in bed format

=head1 DESCRIPTION

This class writes out the supplied L<GenomicRegion|Sanger::CGP::Vagrent::Data::GenomicRegion> objects to a bed file.

Inherits from L<Sanger::CGP::Vagrent::IO::GenomicRegionWriter|Sanger::CGP::Vagrent::IO::GenomicRegionWriter>, view that for method details
