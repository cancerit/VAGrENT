package Sanger::CGP::VagrentSV::Results::TranscriptResults;

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

use Data::Dumper;

use Const::Fast qw(const);
use Try::Tiny qw(try catch);

use Sanger::CGP::Vagrent qw($VERSION);
use Sanger::CGP::VagrentSV::Base;

use FindBin qw($Bin);
Log::Log4perl->init("$Bin/../config/log4perl.vagrentsv.conf");
my $log = Log::Log4perl->get_logger(__PACKAGE__);

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
		$self->{"_$k"} = $vars{$k};
	}
}


sub getLExonNum {
	shift->{_lexon};
}

sub getRExonNum {
	shift->{_rexon};
}

sub getAnno {
	shift->{_anno};
}

sub getRhbSeq {
	shift->{_rseq};
}

sub getLhbSeq {
	shift->{_lseq};
}

sub getTranscript {
	shift->{_tr};
}

sub getLhbGenes {
	shift->{_lhb_genes};
}
sub getRhbGenes {
	shift->{_rhb_genes};
}

sub getSpannedGenes {
	shift->{_spanned_genes};
}
























