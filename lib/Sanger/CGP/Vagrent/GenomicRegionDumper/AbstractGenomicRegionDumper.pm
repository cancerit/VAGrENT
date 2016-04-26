package Sanger::CGP::Vagrent::GenomicRegionDumper::AbstractGenomicRegionDumper;

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
use Carp;
use Attribute::Abstract;
use Data::Dumper;
use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Sanger::CGP::Vagrent::Annotators::AbstractAnnotator);

my $log = Log::Log4perl->get_logger(__PACKAGE__);

1;

sub _init {
	my $self = shift;
	my %vars = @_;
	foreach my $k(keys(%vars)){
		if($k eq 'transcriptSource' && $vars{transcriptSource}->isa('Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource')){
			$self->{_transcriptSource} = $vars{transcriptSource};
		} elsif($k eq 'writer' && $vars{writer}->isa('Sanger::CGP::Vagrent::IO::GenomicRegionWriter')){
			$self->{_writer} = $vars{writer};
		} elsif($k eq 'debug' && $vars{debug}){
			$self->{_debug} = 1;
		}
	}
}

sub dumpGeneRegions {
	my ($self,$gr) = @_;
	unless(defined($gr) && $gr->isa('Sanger::CGP::Vagrent::Data::GenomicRegion')){
		$log->error("Did not recieve a Sanger::CGP::Vagrent::Data::GenomicRegion object");
		return undef;
	}
	$self->_getTranscriptSource()->setDumpRegion($gr);
  foreach my $t($self->_getTranscriptSource()->getTranscripts($gr)){

      my $min = $t->getGenomicMinPos - $self->getUpDownStream5kbCutoff;
      my $max = $t->getGenomicMaxPos + $self->getUpDownStream5kbCutoff;
      if($min < $gr->getMinPos){
        $min = $gr->getMinPos;
      }
      if($max > $gr->getMaxPos){
        $max = $gr->getMaxPos;
      }
      my $geneR = Sanger::CGP::Vagrent::Data::GenomicRegion->new(
        'species'				=> $gr->getSpecies,
        'genomeVersion'         => $gr->getGenomeVersion,
        'chr' 					=> $gr->getChr,
        'minpos'				=> $min,
        'maxpos'				=> $max,);

      $self->_getWriter->write($geneR);
		}

#	}

	return 1;
}

sub dumpExonRegions {
	my ($self,$gr) = @_;
	return $self->_dumpRegionsWithExonConverter($gr,sub {$self->_convertExonToAnnotatableExonicRegions(@_)});
}

sub dumpCodingExonRegions {
	my ($self,$gr) = @_;
	return $self->_dumpRegionsWithExonConverter($gr,sub {$self->_convertExonToAnnotatableCodingExonicRegions(@_)});
}

sub _dumpRegionsWithExonConverter {
  my ($self,$gr,$exonConverter) = @_;
	unless(defined($gr) && $gr->isa('Sanger::CGP::Vagrent::Data::GenomicRegion')){
		$log->error("Did not recieve a Sanger::CGP::Vagrent::Data::GenomicRegion object");
		return undef;
	}
  $self->_getTranscriptSource()->setDumpRegion($gr);
  foreach my $t($self->_getTranscriptSource()->getTranscripts($gr)){
    my @exons = $t->getExonsGenomicOrder;
    my $ec = 0;
    foreach my $e(@exons){
      $ec++;
      my $editStart = 1;
      my $editEnd = 1;
      $editStart = 0 if $ec == 1;
      $editEnd = 0 if $ec == scalar(@exons);
      my @eReg = &$exonConverter($e,$t,$editStart,$editEnd);
      next unless(defined $eReg[0]);
      foreach my $r(@eReg){
        $self->_getWriter->write($r);
      }
    }
  }
  return 1;
}


sub _regionComparison {
	my ($self,$ref,$test,$loud) = @_;
	if($test->getChr eq $ref->getChr && $test->getMinPos <= $ref->getMaxPos && $test->getMaxPos >= $ref->getMinPos){
		if($test->getMinPos == $ref->getMinPos && $test->getMaxPos == $ref->getMaxPos){
			# identical region already exists, do nothing
			 print join(',',"\t",$ref->getChr,$ref->getMinPos,$ref->getMaxPos,'IDENT'),"\n" if $loud;
		} elsif($test->getMinPos >= $ref->getMinPos && $test->getMaxPos <= $ref->getMaxPos) {
			# a region that fully encompasses the current one already exists, do nothing
			 print join(',',"\t",$ref->getChr,$ref->getMinPos,$ref->getMaxPos,'CONTAINED'),"\n" if $loud;
		} else {
			# a region that overlaps the current one exists, edit that region as its already in the list
			 print join(',',"\t",$ref->getChr,$ref->getMinPos,$ref->getMaxPos,'OVERHANG'),"\t" if $loud;
			if($test->getMinPos < $ref->getMinPos){
				$ref->setMinPos($test->getMinPos);
			}
			if($test->getMaxPos > $ref->getMaxPos){
				$ref->setMaxPos($test->getMaxPos);
			}
			 print join(',',"\t",$ref->getChr,$ref->getMinPos,$ref->getMaxPos),"\n" if $loud;
		}
		# region is not new, shouldn't be added to the list
		return 1;
	} else {
		 print join(',',"\t",$ref->getChr,$ref->getMinPos,$ref->getMaxPos,'NOPE'),"\n" if $loud;
		 return 0;
	}
}

sub _convertExonToAnnotatableExonicRegions: Abstract;

sub _convertExonToAnnotatableCodingExonicRegions: Abstract;

sub _getWriter {
	return shift->{_writer};
}


__END__

=head1 NAME

Sanger::CGP::Vagrent::GenomicRegionDumper::AbstractGenomicRegionDumper - Abstract template class for dumping annotatable genomic regions

=head1 DESCRIPTION

This is an abstract template class for saving the regions of the genome that could be annotated to.

Using the supplied L<TranscriptSource|Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource>, L<Transcripts|Sanger::CGP::Vagrent::Data::Transcript> that could generate annotations are selected and the relevent genomic regions are saved to the specified L<GenomicRegionWriter|Sanger::CGP::Vagrent::IO::GenomicRegionWriter>

Inherits from L<Sanger::CGP::Vagrent::Annotators::AbstractAnnotator|Sanger::CGP::Vagrent::Annotators::AbstractAnnotator>.


=head1 METHODS

=head2 Functions

=head3 dumpGeneRegions

=over

=item Usage :

 $dumper->dumpGeneRegions($genomicRegion);

=item Function :

Finds annotatable gene footprints with in the specified L<GenomicRegion|Sanger::CGP::Vagrent::Data::GenomicRegion> from the attached L<TranscriptSource|Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource> and sends them to the L<GenomicRegionWriter|Sanger::CGP::Vagrent::IO::GenomicRegionWriter>.

=item Returns :

Boolean, true if successful

=item Params :

A L<Sanger::CGP::Vagrent::Data::GenomicRegion|Sanger::CGP::Vagrent::Data::GenomicRegion> object

=back

=head3 dumpExonRegions

=over

=item Usage :

 $dumper->dumpExonRegions($genomicRegion);

=item Function :

Finds annotatable exon footprints with in the specified L<GenomicRegion|Sanger::CGP::Vagrent::Data::GenomicRegion> from the attached L<TranscriptSource|Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource> and sends them to the L<GenomicRegionWriter|Sanger::CGP::Vagrent::IO::GenomicRegionWriter>.

=item Returns :

Boolean, true if successful

=item Params :

A L<Sanger::CGP::Vagrent::Data::GenomicRegion|Sanger::CGP::Vagrent::Data::GenomicRegion> object

=back

=head3 dumpCodingExonRegions

=over

=item Usage :

 $dumper->dumpCodingExonRegions($genomicRegion);

=item Function :

Finds annotatable protein coding exon footprints with in the specified L<GenomicRegion|Sanger::CGP::Vagrent::Data::GenomicRegion> from the attached L<TranscriptSource|Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource> and sends them to the L<GenomicRegionWriter|Sanger::CGP::Vagrent::IO::GenomicRegionWriter>.

=item Returns :

Boolean, true if successful

=item Params :

A L<Sanger::CGP::Vagrent::Data::GenomicRegion|Sanger::CGP::Vagrent::Data::GenomicRegion> object

=back

=head2 Abstract

=head3 _convertExonToAnnotatableExonicRegions

=over

=item Usage :

 $dumper->_convertExonToAnnotatableExonicRegions($exon,$transcript,$process5Psplice,$process3Psplice);

=item Function :

Abstract method, needs to be overridden in subclass.  Converts an exon into annotatable genomic regions

=item Returns :

An array of L<Sanger::CGP::Vagrent::Data::GenomicRegion|Sanger::CGP::Vagrent::Data::GenomicRegion> objects

=item Params :

A L<Sanger::CGP::Vagrent::Data::Exon|Sanger::CGP::Vagrent::Data::Exon> object

A L<Sanger::CGP::Vagrent::Data::Transcript|Sanger::CGP::Vagrent::Data::Transcript> object

A boolean switch for applying splice site calculation rules to the 5' end of the exon

A boolean switch for applying splice site calculation rules to the 3' end of the exon

=back

=head3 _convertExonToAnnotatableCodingExonicRegions

=over

=item Usage :

 $dumper->_convertExonToAnnotatableCodingExonicRegions($exon,$transcript,$process5Psplice,$process3Psplice);

=item Function :

Abstract method, needs to be overridden in subclass.  Converts an exon into annotatable genomic regions corresponding to coding sequence

=item Returns :

An array of L<Sanger::CGP::Vagrent::Data::GenomicRegion|Sanger::CGP::Vagrent::Data::GenomicRegion> objects

=item Params :

A L<Sanger::CGP::Vagrent::Data::Exon|Sanger::CGP::Vagrent::Data::Exon> object

A L<Sanger::CGP::Vagrent::Data::Transcript|Sanger::CGP::Vagrent::Data::Transcript> object

A boolean switch for applying splice site calculation rules to the 5' end of the exon

A boolean switch for applying splice site calculation rules to the 3' end of the exon

=back
