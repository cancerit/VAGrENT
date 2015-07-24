package Sanger::CGP::VagrentSV::Annotators::BreakPointTranscripts;
																					 
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
use List::Util qw(first);
use Data::Dumper;
use Attribute::Abstract;

use Sanger::CGP::Vagrent qw($VERSION);
use Sanger::CGP::VagrentSV::Data::BreakPoint;
use Sanger::CGP::Vagrent::Data::Insertion;
my $log = Log::Log4perl->get_logger(__PACKAGE__);

use base qw(Sanger::CGP::VagrentSV::Annotators::AbstractSVAnnotator);



1;


sub _localInit {
	my $self = shift;
	#my %vars = @_;
}

sub getTranscriptAnnotation {
	my ($self,$t,$pos,$ins)=@_;
	my $trObj;
	my($lexon,$lcds,$rexon,$rcds)=$self->_getBreakPointSeq($t,$pos);
	$trObj->{'lexon'} = $lexon; 
	$trObj->{'rexon'} = $rexon; 
	$trObj->{'lseq'} = $lcds;
	$trObj->{'rseq'} = $rcds;
	$trObj->{'tr'} = $t;
	$trObj->{'anno'}= $self->_stringifyAnnotation($self->getSVannotator->getTranscriptAnnotation($ins,$t));
	return $trObj;
}

sub _getBreakPointSeq {
	my ($self,$tr,$pos)=@_;
	my @exons=$tr->getExons;
	my $total=0;
	my $lcount=0;
	my $rcount=0;
	my $final_exon=0;
	my $lexon=0;
	my $rexon=0;
	my $lcds=undef;
	my $rcds=undef;
	#print Dumper $tr;
	
	#next if $ann!~/ENST00000450551/;

	foreach my $exon (@exons) {
		$total++;
		if ( ($exon->getMinPos < $pos) && ($exon->getMaxPos > $pos) ) {
			$final_exon=$total;
			last;
		}
		if ( ($exon->getMinPos < $pos) && ($exon->getMaxPos < $pos) ) {
			$lcount++;
		}
		if (($exon->getMinPos > $pos) && ($exon->getMaxPos > $pos) ) {
			$rcount++;
		}
	}
	
	# if RNA gene then set cds min =1  and max  as RnaMaxPos of Last exon
	if (!$tr->getCdsMaxPos && !$tr->getCdsMinPos) {
		$tr->setCdsMinPos(1);
		$tr->setCdsMaxPos($exons[$#exons]->getRnaMaxPos);
	}
	# if breakpoint pos is in exon 
	if( $final_exon > 0 ) {
			$rexon=$final_exon;
			$lexon=$final_exon;
		my $rnapos = ( ($exons[$final_exon - 1]->getRnaMinPos) + ($pos - ($exons[$final_exon - 1]->getMinPos) ) );
		if($tr->getStrand eq '-1') {
				$lcds=$self->getSeq($tr->getcDNASeq,$rnapos,$tr->getCdsMaxPos);
				$rcds=$self->getSeq($tr->getcDNASeq, $tr->getCdsMinPos, $rnapos);
			}
		else{
			$lcds=$self->getSeq($tr->getcDNASeq,$tr->getCdsMinPos, $rnapos);
			$rcds=$self->getSeq($tr->getcDNASeq, $rnapos, $tr->getCdsMaxPos );
		}
	}
	elsif($tr->getStrand eq '-1') {
		$lexon=$rcount + 1 ;
		$lcds=$self->getSeq($tr->getcDNASeq,$exons[$lexon - 1	]->getRnaMinPos,$tr->getCdsMaxPos);
		$rexon=$rcount;
		$rcds=$self->getSeq($tr->getcDNASeq,$tr->getCdsMinPos,$exons[$rexon - 1]->getRnaMaxPos);
	}
	else{
		$lexon=$lcount;
		$lcds=$self->getSeq($tr->getcDNASeq,$tr->getCdsMinPos, $exons[$lexon - 1]->getRnaMaxPos);
		$rexon=$lcount + 1 ;
		$rcds=$self->getSeq($tr->getcDNASeq,$exons[$rexon - 1]->getRnaMinPos,$tr->getCdsMaxPos);
	}
		
	#print "\n Exon: $final_exon : pos $pos \n $cds \n STRAND ".$tr->getStrand."\n";
	#print Dumper $tr;

	return ($lexon,$lcds,$rexon,$rcds);
}

sub getOverlappingGenes {
	my ($self,$sv)=@_;
	# if breakpoint interval is on same chromosome
			my $spanned_genes;
			$log->debug("breakpoints on same chromosome");
			my($lhbGenes)=$self->_getBreakPointGenes($sv->getLhb);
			my($rhbGenes)=$self->_getBreakPointGenes($sv->getRhb);
			if ($sv->getLhb->getChr eq $sv->getRhb->getChr) {
				my $intObj=$self->_createIntervalObject($sv->getLhb->getSpecies,$sv->getLhb->getGenomeVersion,$sv->getLhb->getChr,$sv->getLhb->getMaxPos,$sv->getRhb->getMaxPos);
				$spanned_genes=$self->_getSpannedGenes($intObj);
			}
			else {
				$spanned_genes='NA';
			}
			return ($lhbGenes,$spanned_genes,$rhbGenes);
}

sub _createIntervalObject {
		my ($self,$species,$genomeversion,$chr,$minpos,$maxpos)=@_;
		if($maxpos<$minpos) {
			my $tmp=$maxpos;
			$maxpos=$minpos;
			$minpos=$tmp;
		}
    my $intObj = Sanger::CGP::VagrentSV::Data::BreakPoint->new(
							'species'				=> $species,
							'genomeVersion' => $genomeversion,
							'chr'	          => $chr,
							'minpos'        => $minpos,
							'maxpos'        => $maxpos); 
 	
	return $intObj;
}


sub _getBreakPointGenes {
	my ($self,$bp)=@_;
	my $overlapping_genes=undef;
	foreach my $t ($self->{'_transcriptSource'}->getTranscripts($bp)) {
		if (($t->getGenomicMinPos < $bp->getMaxPos) && ($t->getGenomicMaxPos > $bp->getMaxPos)) {
			$overlapping_genes->{$t->getGeneName}++;
		}
	}
	my $genes=join (',', keys %$overlapping_genes);
	$genes='NA' if(!$genes);
	return $genes;
}


sub _getSpannedGenes {
	my ($self,$bp)=@_;
	my $spanned_genes=undef;
	foreach my $t ($self->{'_transcriptSource'}->getTranscripts($bp)) {
		if (($t->getGenomicMinPos > $bp->getMinPos) && ($t->getGenomicMaxPos < $bp->getMaxPos)) {
			$spanned_genes->{$t->getGeneName}++;
		}
	}
	my $genes=join (',', keys %$spanned_genes);
	$genes='NA' if(!$genes);
	return $genes;
}



sub _stringifyAnnotation {
  my ($self,$anno)=@_;
  return undef unless defined $anno;
  my $mrna = $anno->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext);
	#my $cds = $anno->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext);
	#my $prot = $anno->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext);
my $desc;
	#$desc = $anno->getLabel.'|';
=head
	if(defined($anno->getCCDS) && $anno->getCCDS ne ''){
		$desc .= $anno->getCCDS.'|';
	} else {
		$desc .= $anno->getAccession.'|';
	}
	$desc .= $mrna->getDescription.'|';
	if(defined($cds)){
		$desc .= $cds->getDescription.'|';
	} else {
		$desc .= '-|';
	}
	if(defined($prot)){
		$desc .= $prot->getDescription.'|';
	} else {
		$desc .= '-|';
	}
=cut
 	my @classIds;
 	my @classText;

 #	foreach my $a($anno,$mrna,$cds,$prot){
 	foreach my $a($anno,$mrna){
 		next unless defined $a;
 		foreach my $term ($a->getClassifications){
 			my ($id,$text) = (split(':',$term))[1,2];
 			if((scalar @classIds) == 0){
 				push @classIds, $id;
 				push @classText, $text;
 			}
 			else {
 				next if(first {$_ eq $id} @classIds);
 				push @classIds, $id;
 				push @classText, $text;
 			}
 		}
 	}
	$desc.= join(':',@classText[1]);
	#$desc .= join(':',@classText.'|';
	#$desc .= 'SO:'.join(':SO:',@classIds);

	return $desc;
}










