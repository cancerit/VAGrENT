package Sanger::CGP::Vagrent::GenomicRegionDumper::IndelGenomicRegionDumper;

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

use Log::Log4perl qw(:easy);
use Data::Dumper;
use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Sanger::CGP::Vagrent::GenomicRegionDumper::AbstractGenomicRegionDumper);

1;

sub _convertExonToAnnotatableExonicRegions {
	my ($self,$e,$t,$procStart,$procEnd) = @_;
	my @out;
  	my $min = $e->getMinPos;
  	my $max = $e->getMaxPos;
  	if($procStart){
		if($t->getStrand > 0){
			# forward strand
			my $off = $self->_getConsesnsusSpliceBeforeBoundry();
			$min = $min + $off;
		} else {
			# reverse strand
			my $off = $self->_getConsesnsusSpliceAfterBoundry();
			$min = $min - $off;
		}
  	}
  	if($procEnd){
  		if($t->getStrand > 0){
			# forward strand
			my $off = $self->_getConsesnsusSpliceAfterBoundry();
			$max = $max + $off;
		} else {
			# reverse strand
			my $off = $self->_getConsesnsusSpliceBeforeBoundry();
			$max = $max - $off;
		}
  	}
  	push(@out,Sanger::CGP::Vagrent::Data::GenomicRegion->new(
			'species'				=> $e->getSpecies,
			'genomeVersion'         => $e->getGenomeVersion,
			'chr' 					=> $e->getChr,
			'minpos'				=> $min,
			'maxpos'				=> $max,));
	return @out;
}

sub _convertExonToAnnotatableCodingExonicRegions {
	my ($self,$e,$t,$procStart,$procEnd) = @_;
  	my @out = ();
  	my $min = undef;
  	my $max = undef;
  	if(($t->isProteinCoding && $e->getRnaMinPos < $t->getCdsMaxPos && $e->getRnaMaxPos > $t->getCdsMinPos) || !$t->isProteinCoding){
		# overhangs with CDS,
		#print join(',',"\t",$e->getChr,$e->getMinPos,$e->getMaxPos,$e->getRnaMinPos,$e->getRnaMaxPos,'CDS'),"\n";
		$min = $e->getMinPos;
  		$max = $e->getMaxPos;
  		if($t->isProteinCoding){
  			if($e->getRnaMinPos < $t->getCdsMinPos && $e->getRnaMaxPos > $t->getCdsMinPos){
  				my $mod = $t->getCdsMinPos - $e->getRnaMinPos;
				#print "START ($mod)\n";
				if($t->getStrand > 0){
					$procStart = 0;
					$min = $min + $mod;
				} else {
					$procEnd = 0;
					$max = $max - $mod;
				}
  			}
  			if($e->getRnaMinPos < $t->getCdsMaxPos && $e->getRnaMaxPos > $t->getCdsMaxPos){
  				my $mod = $e->getRnaMaxPos - $t->getCdsMaxPos;
				#print "END ($mod)\n";
				if($t->getStrand > 0){
					$max = $max - $mod;
					$procEnd = 0;
				} else {
					$procStart = 0;
					$min = $min + $mod;
				}
  			}
  		}
  	} else {
  		# UTR, SKIP;
  		return;
  	}
  	if($procStart){
		if($t->getStrand > 0){
			# forward strand
			my $off = $self->_getConsesnsusSpliceBeforeBoundry();
			$min = $min + $off;
		} else {
			# reverse strand
			my $off = $self->_getConsesnsusSpliceAfterBoundry();
			$min = $min - $off;
		}
  	}
  	if($procEnd){
  		if($t->getStrand > 0){
			# forward strand
			my $off = $self->_getConsesnsusSpliceAfterBoundry();
			$max = $max + $off;
		} else {
			# reverse strand
			my $off = $self->_getConsesnsusSpliceBeforeBoundry();
			$max = $max - $off;
		}
  	}
  	push(@out,Sanger::CGP::Vagrent::Data::GenomicRegion->new(
			'species'				=> $e->getSpecies,
			'genomeVersion'         => $e->getGenomeVersion,
			'chr' 					=> $e->getChr,
			'minpos'				=> $min,
			'maxpos'				=> $max,));
	return @out;
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::GenomicRegionDumper::IndelGenomicRegionDumper - Class for saving regions of the genome that would provide annotations for Insertion, Deletion or Complex variation types.

=head1 DESCRIPTION

This is an implementation of the L<AbstractGenomicRegionDumper|Sanger::CGP::Vagrent::GenomicRegionDumper::AbstractGenomicRegionDumper> to handle indel type variants, specifically any kind of variant that has the potential cause a length change in the genomic sequence.

Using the supplied L<TranscriptSource|Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource>, L<Transcripts|Sanger::CGP::Vagrent::Data::Transcript> that could generate annotations are selected and the relevent genomic regions are saved to the specified L<GenomicRegionWriter|Sanger::CGP::Vagrent::IO::GenomicRegionWriter>

Inherits from L<Sanger::CGP::Vagrent::GenomicRegionDumper::AbstractGenomicRegionDumper|Sanger::CGP::Vagrent::GenomicRegionDumper::AbstractGenomicRegionDumper>.

=head1 METHODS

See L<Sanger::CGP::Vagrent::GenomicRegionDumper::AbstractGenomicRegionDumper|Sanger::CGP::Vagrent::GenomicRegionDumper::AbstractGenomicRegionDumper>

