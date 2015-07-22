package Sanger::CGP::VagrentSV::Data::GenomeSeq;
																					 
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

use Log::Log4perl;
use POSIX qw(ceil);
use Data::Dumper;

use Bio::DB::Sam;
use Sanger::CGP::VagrentSV::Base;
use Sanger::CGP::Vagrent qw($VERSION);


my $log = Log::Log4perl->get_logger(__PACKAGE__);


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
	$self->{'_fai'}	= Bio::DB::Sam::Fai->load($vars{'fasta'});
	$self->{'_chrlen'} = $self->_getChrLength($vars{'index'});
}

sub getSeq {
	my ($self,$chr,$start,$end)=@_;
	my $dna = $self->_getFai->fetch("$chr:$start-$end");
	return $dna;
}

sub _getFai {
	return shift->{'_fai'};
}


sub getChrLength {
	return shift->{'_chrlen'};
}

sub _getChrLength {
	my($self,$index_file)=@_;
	my $chr_length;
	my ($rfh_index)=Sanger::CGP::VagrentSV::Base->open_to_read($index_file);
	while(<$rfh_index>) {
		my($chr,$length)=(split "\t", $_)[0,1];
		$chr_length->{$chr}=$length;
	}
	close($rfh_index);
	return $chr_length;
}









