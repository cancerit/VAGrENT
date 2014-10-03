package Sanger::CGP::Vagrent::IO::Writer::SingleLineCsvFileWriter;

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

my $seperator = ',';

use constant HEADER => 'MUT ID,MUT TYPE,SPECIES,GENE TYPE,GENE,CCDS,DB SOURCE,DB VERSION,GROUP SO CLASSES,mRNA ANNO TYPE,mRNA ANNO SUBTYPE,mRNA ACCESSION,mRNA VERSION,mRNA LENGTH,mRNA MIN,mRNA MINOFFSET,mRNA MAX,mRNA MAXOFFSET,mRNA WT,mRNA MT,mRNA DESC,mRNA SO CLASSES,CDS ANNO TYPE,CDS ANNO SUBTYPE,CDS ACCESSION,CDS VERSION,CDS LENGTH,CDS MIN,CDS MINOFFSET,CDS MAX,CDS MAXOFFSET,CDS WT,CDS MT,CDS DESC,CDS SO CLASSES,PROTEIN ANNO TYPE,PROTEIN ANNO SUBTYPE,PROTEIN ACCESSION,PROTEIN VERSION,PROTEIN LENGTH,PROTEIN MIN,PROTEIN MINOFFSET,PROTEIN MAX,PROTEIN MAXOFFSET,PROTEIN WT,PROTEIN MT,PROTEIN DESC,PROTEIN SO CLASSES,';

1;

sub write {
	my ($self,$m,$g) = @_;
	unless($m->isa('Sanger::CGP::Vagrent::Data::AbstractVariation')){
		warn 'expecting an Sanger::CGP::Vagrent::Data::AbstractVariation, received a' . ref($m);
		return;
	}
	unless($g->isa('Sanger::CGP::Vagrent::Data::AnnotationGroup')){
		warn 'expecting an Sanger::CGP::Vagrent::Data::AnnotationGroup, received a' . ref($g);
		return;
	}
	unless(exists($self->{_linesWritten}) && defined($self->{_linesWritten}) && $self->{_linesWritten} > 0){
		$self->_writeHeader();
	}
	$self->{_linesWritten}++;
	$self->_writeAnnotationGroup($m,$g);
}

sub _writeHeader {
	my $self = shift;
	print {$self->_fh} HEADER."\n";
}

sub _writeAnnotationGroup {
	my ($self,$m,$g) = @_;
	my $id = 'NONE';
	if(defined($m->getId())){
		$id = $m->getId();
	}

	my $mutType = $self->_getMutType($m);

	print {$self->_fh} join($seperator,$id,$mutType,$m->getSpecies(),$g->getType,$g->getLabel,$g->getCCDS,'');

	my $mRNA = $g->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext());
	my $classes = '';
	$classes = join('|',$g->getClassifications) if(defined($g->getClassifications));

	print {$self->_fh} join($seperator,$mRNA->getDatabase,$mRNA->getDatabaseVersion,$classes,'');

	$self->_writeAnnotation($mRNA);

	foreach my $ctx ($mRNA->getCDSAnnotationContext(),$mRNA->getProteinAnnotationContext()){
		my $tmpAnno = $g->getAnnotationByContext($ctx);
		if(defined($tmpAnno)){
			$self->_writeAnnotation($tmpAnno);
		}
	}
	print {$self->_fh} "\n";
}

sub _writeAnnotation {
	my ($self,$a) = @_;
	if(defined($a)){
		my $classes = '';
		$classes = join('|',$a->getClassifications) if(defined($a->getClassifications));
		print {$self->_fh} join($seperator,$a->getType,$a->getSubtype,
		$a->getAccession,$a->getSequenceVersion,$a->getSequenceLength,$a->getMinPos,
		$a->getMinOffset,$a->getMaxPos,$a->getMaxOffset,
		$a->getWt,$a->getMt,$a->getDescription,
		$classes,'');
	}
}

sub _getMutType {
	my ($self,$m) = @_;
	if($m->isa('Sanger::CGP::Vagrent::Data::Substitution')){
		return 'SUB';
	} elsif($m->isa('Sanger::CGP::Vagrent::Data::Deletion')){
		return 'DEL';
	} elsif($m->isa('Sanger::CGP::Vagrent::Data::Insertion')){
		return 'INS';
	} elsif($m->isa('Sanger::CGP::Vagrent::Data::ComplexIndel')){
		return 'INDEL';
	} else {
		$self->throw('unknown mutation class: '.ref($m));
	}
}

=head1 NAME

Sanger::CGP::Vagrent::IO::Writer::SingleLineCsvFileWriter - Another simple csv file writer

=head1 DESCRIPTION

This class prints a comma seperated text representation of L<AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> and L<Annotation|Sanger::CGP::Vagrent::Data::Annotation> objects.
It represents a single AnnotationGroup and its Annotation objects as a single line.

Inherits from L<Sanger::CGP::Vagrent::IO::AnnotationWriter|Sanger::CGP::Vagrent::IO::AnnotationWriter>, view that for method details
