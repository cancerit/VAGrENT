package Sanger::CGP::VagrentSV::IO::OutputFormatter;
																					 
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
use Attribute::Abstract;

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
}


sub getPegasusInputFormat {
	my($self,$ltr,$rtr,$sv,$anno)=@_;
	my $pegasus;
	my $ldata=$self->_getRepresentativeTranscript($ltr);
	my $rdata=$self->_getRepresentativeTranscript($rtr);
		
	if (defined $ldata && defined $rdata) {
		foreach my $lgene (keys %$ldata) {
		  my $ldefuse=$self->_formatDefuseLine($ldata,$lgene);
				foreach my $rgene(keys %$rdata) {		
					my $rdefuse=$self->_formatDefuseLine($rdata,$rgene);
					push (@$pegasus,$self->_getDefuseLine($ldefuse,$rdefuse,$sv));					
				}
		}		
	}
	elsif($ldata){
		my $rdefuse=$self->_formatDefuseLine(undef,undef);
		foreach my $lgene (keys %$ldata) {
			my $ldefuse=$self->_formatDefuseLine($ldata,$lgene);
		 	push (@$pegasus,$self->_getDefuseLine($ldefuse,$rdefuse,$sv));
		 }
	}
	elsif($rdata) {
		my $ldefuse=$self->_formatDefuseLine(undef,undef);
		foreach my $rgene (keys %$rdata) {
			my $rdefuse=$self->_formatDefuseLine($rdata,$rgene);
			push (@$pegasus,$self->_getDefuseLine($ldefuse,$rdefuse,$sv));
		}
	}
	else{ 
		my $ldefuse=$self->_formatDefuseLine(undef,undef);
		my $rdefuse=$self->_formatDefuseLine(undef,undef);
		#push (@$pegasus,$self->_getDefuseLine($ldefuse,$rdefuse,$sv));
	}
	return $pegasus;
}

sub _getRepresentativeTranscript {
	my ($self,$tr)=@_;
	my $results=undef;
	my $len=0;
	my $tmp_gene='NA';
	foreach my $t (@$tr) {
		if( ($len < $t->getTranscript->getCdsLength) or  ($tmp_gene ne $t->getTranscript->getGeneName) ) {
			$results->{$t->getTranscript->getGeneName}=$t;
			$len=$t->getTranscript->getCdsLength;
			$tmp_gene=$t->getTranscript->getGeneName;
		}
	}
	return $results;
}


sub getCharStrand {
	my($self,$strand)=@_;
	my $chrStrand='-';
	if($strand > 0 ) {
		$chrStrand='+';
	}
	return $chrStrand;
}



sub _getDefuseLine {
	my($self,$ldefuse,$rdefuse,$sv)=@_;
	
	my $line= $ldefuse->{'tr'}."\t".$rdefuse->{'tr'}."\t".
					$sv->getLhb->getStrand."\t".$sv->getRhb->getStrand."\t".
					$sv->getLhb->getChr."\t".$sv->getRhb->getChr."\t".
					$ldefuse->{'gene_end'}."\t".$rdefuse->{'gene_end'}."\t".
					$ldefuse->{'anno'}."\t".$rdefuse->{'anno'}."\t".
					$ldefuse->{'gene'}."\t".$rdefuse->{'gene'}."\t".
					$ldefuse->{'gene_start'}."\t".$rdefuse->{'gene_start'}."\t".
					$ldefuse->{'strand'}."\t".$rdefuse->{'strand'}."\t0\t".
					$sv->getLhb->getMaxPos."\t".$sv->getRhb->getMaxPos;
					
	return $line;
}


sub _formatDefuseLine {
	my ($self,$data,$gene)=@_;
	my $line;
	if($data){
		$line->{'tr'}=$data->{$gene}->getTranscript->getAccession;
		$line->{'gene_end'}=$data->{$gene}->getTranscript->getGenomicMaxPos;
		$line->{'anno'}=$data->{$gene}->getAnno;
		$line->{'gene_start'}=$data->{$gene}->getTranscript->getGenomicMinPos;
		$line->{'strand'}=$self->getCharStrand($data->{$gene}->getTranscript->getStrand);
		$line->{'gene'}=$gene;
	}
	else {
		$line->{'tr'}='NA';
		$line->{'gene_end'}='NA';
		$line->{'anno'}='NA';
		$line->{'gene_start'}='NA';
		$line->{'strand'}='NA';
		$line->{'gene'}='NA';
	}
 	return $line;
}



