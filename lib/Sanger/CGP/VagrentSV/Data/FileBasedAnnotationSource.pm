package Sanger::CGP::VagrentSV::Data::FileBasedAnnotationSource;

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
use Log::Log4perl;
use Data::Dumper;
use Const::Fast qw(const);
use Sanger::CGP::VagrentSV::SVConstants;

# prevent the init warnings from Tabix
BEGIN {
  $SIG{__WARN__} = sub {warn $_[0] unless( $_[0] =~ m/^Subroutine Tabix.* redefined/)};
};

use Tabix;

use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Sanger::CGP::VagrentSV::Data::AbstractAnnotationSource);

my $log = Log::Log4perl->get_logger(__PACKAGE__);



1;


sub _locatInit {
	my $self=shift;
}


sub getRegulatoryAnnotations {
	my($self,$sv)=@_;
	 my $buffer=$Sanger::CGP::VagrentSV::SVConstants::PROMOTER_SEARCH_PADDING;
	my($lhbGenes)=$self->_getIntervalOverlap($sv->getLhb->getChr,$sv->getLhb->getMaxPos,$buffer);
	my($rhbGenes)=$self->_getIntervalOverlap($sv->getRhb->getChr,$sv->getRhb->getMaxPos,$buffer);
	
	return($lhbGenes,$rhbGenes);
	
	
}


sub _getIntervalOverlap {
	my($self,$chr,$pos,$buffer)=@_;
	my $overlap;
	my $lstart=$pos - $buffer;
	my $rstop=$pos + $buffer;
	my ($lgenes)=$self->_getOverlappingGenes($chr,$lstart,$pos);
	my ($rgenes)=$self->_getOverlappingGenes($chr,$pos,$rstop);
	$overlap->{'lgenes'}=$lgenes;
	$overlap->{'rgenes'}=$rgenes;
	return $overlap;
}



sub _getOverlappingGenes {
	my ($self,$chr,$start,$stop)=@_;
	my $tabix=$self->{'tabix'};
	$chr=~s/chr//g;
	my $genes;
	#extract all the bed lines in an interval
	my $res = $tabix->query('chr'.$chr,$start,$stop);
		while(my $record = $tabix->read($res)){	
			my($tmp_gene,$tmp_strand)=(split "\t", $record)[3,5];
			next if $tmp_gene eq '.';
			push(@$genes,split(" ",$tmp_gene));
		}
		if (defined $genes){
			my %unique_genes   = map { $_, 1 } @$genes;
			return join(",", keys %unique_genes);
		}
		else{
			return 'NA';
		}		
}

