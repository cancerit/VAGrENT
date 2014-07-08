package Sanger::CGP::Vagrent::GenomicRegionDumper::SubstitutionGenomicRegionDumper;

use strict;

use Log::Log4perl qw(:easy);
use Data::Dumper;

use base qw(Sanger::CGP::Vagrent::GenomicRegionDumper::AbstractGenomicRegionDumper);

1;

sub _convertExonToAnnotatableExonicRegions {
	my ($self,$e,$t,$procStart,$procEnd) = @_;
  	my @out;
  	my $other = undef;
  	my $min = $e->getMinPos;
  	my $max = $e->getMaxPos;
  	if($procStart){
		if($t->getStrand > 0){
			# forward strand
			my $off = $self->_getLowestContinuousAcceptorEssentialSpliceSiteOffset();
			$min = $min + $off;
		} else {
			# reverse strand
			my $iso = $self->_getIsolatedDonorEssentialSpliceSiteOffset();
			$other = Sanger::CGP::Vagrent::Data::GenomicRegion->new(
				'species'				=> $e->getSpecies,
				'genomeVersion'         => $e->getGenomeVersion,
				'chr' 					=> $e->getChr,
				'minpos'				=> ($min - $iso),
				'maxpos'				=> ($min - $iso),);

			my $off = $self->_getHighestContinuousDonorEssentialSpliceSiteOffset();
			$min = $min - $off;
		}
  	}
  	if($procEnd){
  		if($t->getStrand > 0){
			# forward strand
			my $iso = $self->_getIsolatedDonorEssentialSpliceSiteOffset();
			$other = Sanger::CGP::Vagrent::Data::GenomicRegion->new(
				'species'				=> $e->getSpecies,
				'genomeVersion'         => $e->getGenomeVersion,
				'chr' 					=> $e->getChr,
				'minpos'				=> ($max + $iso),
				'maxpos'				=> ($max + $iso),);

			my $off = $self->_getHighestContinuousDonorEssentialSpliceSiteOffset();
			$max = $max + $off;
		} else {
			my $off = $self->_getLowestContinuousAcceptorEssentialSpliceSiteOffset();
			$max = $max - $off;
			# reverse strand
		}
  	}
  	push(@out,Sanger::CGP::Vagrent::Data::GenomicRegion->new(
			'species'				=> $e->getSpecies,
			'genomeVersion'         => $e->getGenomeVersion,
			'chr' 					=> $e->getChr,
			'minpos'				=> $min,
			'maxpos'				=> $max,));

	if(defined $other){
		push(@out,$other);
	}

  	return @out;
}

sub _convertExonToAnnotatableCodingExonicRegions {
	my ($self,$e,$t,$procStart,$procEnd) = @_;
  	my @out = ();
  	my $other = undef;
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
  		#print join(',',"\t",$e->getChr,$e->getMinPos,$e->getMaxPos,$e->getRnaMinPos,$e->getRnaMaxPos,'UTR'),"\n";
  		return;
  	}
	if($procStart){
		if($t->getStrand > 0){
			# forward strand
			my $off = $self->_getLowestContinuousAcceptorEssentialSpliceSiteOffset();
			$min = $min + $off;
		} else {
			# reverse strand
			my $iso = $self->_getIsolatedDonorEssentialSpliceSiteOffset();
			$other = Sanger::CGP::Vagrent::Data::GenomicRegion->new(
				'species'				=> $e->getSpecies,
				'genomeVersion'         => $e->getGenomeVersion,
				'chr' 					=> $e->getChr,
				'minpos'				=> ($min - $iso),
				'maxpos'				=> ($min - $iso),);

			my $off = $self->_getHighestContinuousDonorEssentialSpliceSiteOffset();
			$min = $min - $off;
		}
  	}
  	if($procEnd){
  		if($t->getStrand > 0){
			# forward strand
			my $iso = $self->_getIsolatedDonorEssentialSpliceSiteOffset();
			$other = Sanger::CGP::Vagrent::Data::GenomicRegion->new(
				'species'				=> $e->getSpecies,
				'genomeVersion'         => $e->getGenomeVersion,
				'chr' 					=> $e->getChr,
				'minpos'				=> ($max + $iso),
				'maxpos'				=> ($max + $iso),);

			my $off = $self->_getHighestContinuousDonorEssentialSpliceSiteOffset();
			$max = $max + $off;
		} else {
			my $off = $self->_getLowestContinuousAcceptorEssentialSpliceSiteOffset();
			$max = $max - $off;
			# reverse strand
		}
  	}
  	push(@out,Sanger::CGP::Vagrent::Data::GenomicRegion->new(
			'species'				=> $e->getSpecies,
			'genomeVersion'         => $e->getGenomeVersion,
			'chr' 					=> $e->getChr,
			'minpos'				=> $min,
			'maxpos'				=> $max,));

	if(defined $other){
		push(@out,$other);
	}
	return @out;
}

sub _getLowestContinuousAcceptorEssentialSpliceSiteOffset {
	my @offsets = sort {$b <=> $a} Sanger::CGP::Vagrent::Annotators::AbstractAnnotator::CONSENSUS_SPLICE_OFFSETS;
	my $out = 0;
	foreach my $o(@offsets){
		next if $o > 0;
		if($out - 1 == $o){
			$out = $o;
		} else {
			last;
		}
	}
	return $out;
}

sub _getHighestContinuousDonorEssentialSpliceSiteOffset {
	my @offsets = sort {$a <=> $b} Sanger::CGP::Vagrent::Annotators::AbstractAnnotator::CONSENSUS_SPLICE_OFFSETS;
	my $out = 0;
	foreach my $o(@offsets){
		next if $o < 0;
		if($out + 1 == $o){
			$out = $o;
		} else {
			last;
		}
	}
	return $out;
}

sub _getIsolatedDonorEssentialSpliceSiteOffset {
	my ($out) = sort {$b <=> $a} Sanger::CGP::Vagrent::Annotators::AbstractAnnotator::CONSENSUS_SPLICE_OFFSETS;
	return $out;
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::GenomicRegionDumper::SubstitutionGenomicRegionDumper - Class for saving regions of the genome that would provide annotations for Substitution variation types.

=head1 DESCRIPTION

This is an implementation of the L<AbstractGenomicRegionDumper|Sanger::CGP::Vagrent::GenomicRegionDumper::AbstractGenomicRegionDumper> to handle substitution type variants, ie any kind of variant that cannot have an effect on sequence length and only effect a single wildtype base

Using the supplied L<TranscriptSource|Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource>, L<Transcripts|Sanger::CGP::Vagrent::Data::Transcript> that could generate annotations are selected and the relevent genomic regions are saved to the specified L<GenomicRegionWriter|Sanger::CGP::Vagrent::IO::GenomicRegionWriter>

Inherits from L<Sanger::CGP::Vagrent::GenomicRegionDumper::AbstractGenomicRegionDumper|Sanger::CGP::Vagrent::GenomicRegionDumper::AbstractGenomicRegionDumper>.

=head1 METHODS

See L<Sanger::CGP::Vagrent::GenomicRegionDumper::AbstractGenomicRegionDumper|Sanger::CGP::Vagrent::GenomicRegionDumper::AbstractGenomicRegionDumper>
