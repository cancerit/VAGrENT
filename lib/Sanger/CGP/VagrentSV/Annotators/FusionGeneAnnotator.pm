package Sanger::CGP::VagrentSV::Annotators::FusionGeneAnnotator;
																					 
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

use base qw(Sanger::CGP::VagrentSV::Annotators::AbstractSVAnnotator);

1;


sub _localInit {
	my $self = shift;
}


sub getFusedTranscripts {
	my($self,$ltr,$rtr,$sv)=@_;
	
	# transcripts present at both the breakpoints
	my $fusionData;
	if($ltr && $rtr) {
		foreach my $l (@$ltr) {
			foreach my $r(@$rtr) {
				
				if($sv->getSvType eq 'del') {
				push(@$fusionData,	$l->getLExonNum.'/'.$r->getRExonNum."\t".
														$l->getTranscript->getAccession.'/'.$r->getTranscript->getAccession."\t".
														$l->getTranscript->getStrand.'/'.$r->getTranscript->getStrand."\t".
														$l->getAnno.'/'.$r->getAnno);
				}
				
				elsif($sv->getSvType eq 'td') {
				push(@$fusionData,	$l->getLExonNum.'/'.$r->getRExonNum."\t".
														$l->getTranscript->getAccession.'/'.$r->getTranscript->getAccession."\t".
														$l->getTranscript->getStrand.'/'.$r->getTranscript->getStrand."\t".
														$l->getAnno.'/'.$r->getAnno);
				}
				
				elsif($sv->getSvType eq 'bt') {
				push(@$fusionData,	$l->getLExonNum.'/'.$r->getRExonNum."\t".
														$l->getTranscript->getAccession.'/'.$r->getTranscript->getAccession."\t".
														$l->getTranscript->getStrand.'/'.$r->getTranscript->getStrand."\t".
														$l->getAnno.'/'.$r->getAnno);
				}																			
			}
		}
	}
	elsif($ltr && ($sv->getSvType eq 'bt' || $sv->getSvType eq 'del'|| $sv->getSvType eq 'td')){
		foreach my $l(@$ltr){
			push(@$fusionData, $l->getLExonNum."\t".
														$l->getTranscript->getAccession."\t".
														$l->getTranscript->getStrand."\t".
														$l->getAnno);	
		}
	}
	elsif($rtr && ($sv->getSvType eq 'bt' || $sv->getSvType eq 'del'|| $sv->getSvType eq 'td')){
		foreach my $r(@$rtr){
			push(@$fusionData, $r->getRExonNum."\t".
														$r->getTranscript->getAccession."\t".
														$r->getTranscript->getStrand."\t".
														$r->getAnno);	
		}
	}
	else{
		push(@$fusionData, "NA\tNA\tNA\tNA");
	}

return ($fusionData);

}


sub getLminPos {
	shift->{_lmin};
}

sub getRminPos {
	shift->{_rmin};
}

sub getLmaxPos {
	shift->{_lmax};
}

sub getRmaxPos {
	shift->{_rmax};
}

sub getLChr {
	shift->{_lchr};
}

sub getRChr {
	shift->{_rchr};
}

sub getKey {
	shift->{_key};
}

sub getLAnno {
	shift->{_lann};
}

sub getRAnno {
	shift->{_rann};
}

sub getLSeq {
	shift->{_lseq};
}

sub getRSeq {
	shift->{_rseq};
}

sub getLExon {
	shift->{_lexon};
}

sub getRExon {
	shift->{_rexon};
}




