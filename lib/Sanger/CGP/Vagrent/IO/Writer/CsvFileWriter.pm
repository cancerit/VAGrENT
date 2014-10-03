package Sanger::CGP::Vagrent::IO::Writer::CsvFileWriter;

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
use base qw(Sanger::CGP::Vagrent::IO::AnnotationWriter);

1;

sub write {
	my ($self,$m,$g) = @_;
	if($self->_writeMutation($m)){
		$self->_writeAnnotationGroup($g);
	}
}

sub _writeSub {
	my ($self,$m) = @_;
	print {$self->_fh} join(',','ANNOGROUP',$m->getSpecies,$m->getGenomeVersion,'SUB',$m->getChr,$m->getMinPos,$m->getMaxPos,$m->getWt,$m->getMt,'');
}

sub _writeDel {
	my ($self,$m) = @_;
	my $delSeq;
	if(length($m->getDeletedSequence) > 20){
		$delSeq = length($m->getDeletedSequence);
	} else {
		$delSeq = $m->getDeletedSequence;
	}

	print {$self->_fh} join(',','ANNOGROUP',$m->getSpecies,$m->getGenomeVersion,'DEL',$m->getChr,$m->getMinPos,$m->getMaxPos,$delSeq,'-','');
}

sub _writeIns {
	my ($self,$m) = @_;
	print {$self->_fh} join(',','ANNOGROUP',$m->getSpecies,$m->getGenomeVersion,'INS',$m->getChr,$m->getMinPos,$m->getMaxPos,'-',$m->getInsertedSequence,'');
}

sub _writeMutation {
	my ($self,$m) = @_;
	if($m->isa('Sanger::CGP::Vagrent::Data::Substitution')){
		$self->_writeSub($m);
		return 1;
	} elsif($m->isa('Sanger::CGP::Vagrent::Data::Deletion')){
		$self->_writeDel($m);
		return 1;
	} elsif($m->isa('Sanger::CGP::Vagrent::Data::Insertion')){
		$self->_writeIns($m);
		return 1;
	} else {
		return 0;
	}
}

sub _writeAnnotationGroup {
	my ($self,$g) = @_;
	if(defined($g->getClassifications)){
		print {$self->_fh} join(',',$g->getType,$g->getLabel,$g->getCCDS,$g->getAccession,join('|',$g->getClassifications),"\n");
	} else {
		print {$self->_fh} join(',',$g->getType,$g->getLabel,$g->getCCDS,$g->getAccession,'',"\n");
	}

	foreach my $a(@{$g->getAllAnnotations}){
		$self->_writeAnnotation($a);
	}
}

sub _writeAnnotation {
	my ($self,$a) = @_;
	if(defined($a)){
		if(defined($a->getClassifications)){
			print {$self->_fh} join(',','ANNO',$a->getContext,$a->getType,$a->getSubtype,$a->getDatabase,$a->getDatabaseVersion,$a->getAccession,$a->getSequenceVersion,$a->getMinPos,$a->getMinOffset,$a->getMaxPos,$a->getMaxOffset,$a->getWt,$a->getMt,$a->getDescription,join('|',$a->getClassifications),"\n");
		} else {
			print {$self->_fh} join(',','ANNO',$a->getContext,$a->getType,$a->getSubtype,$a->getDatabase,$a->getDatabaseVersion,$a->getAccession,$a->getSequenceVersion,$a->getMinPos,$a->getMinOffset,$a->getMaxPos,$a->getMaxOffset,$a->getWt,$a->getMt,$a->getDescription,",\n");
		}
	}
}

=head1 NAME

Sanger::CGP::Vagrent::IO::Writer::CsvFileWriter - Simple csv file writer

=head1 DESCRIPTION

This class prints a comma seperated text representation of L<AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> and L<Annotation|Sanger::CGP::Vagrent::Data::Annotation> objects over separate lines, 1 line for the Group and 1 line for each Annotation.

Inherits from L<Sanger::CGP::Vagrent::IO::AnnotationWriter|Sanger::CGP::Vagrent::IO::AnnotationWriter>, view that for method details

