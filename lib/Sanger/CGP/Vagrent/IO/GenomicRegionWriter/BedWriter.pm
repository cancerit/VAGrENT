package Sanger::CGP::Vagrent::IO::GenomicRegionWriter::BedWriter;

use strict;
use Data::Dumper;

use base qw(Sanger::CGP::Vagrent::IO::GenomicRegionWriter);

1;

my $seperator = "\t";

sub write {
	my ($self,$gr) = @_;
	unless($gr->isa('Sanger::CGP::Vagrent::Data::GenomicRegion')){
		warn 'expecting an Sanger::CGP::Vagrent::Data::GenomicRegion, received a' . ref($gr);
		return;
	}
	print {$self->_fh} join($seperator,$gr->getChr,($gr->getMinPos - 1),$gr->getMaxPos,),"\n";
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::IO::GenomicRegionWriter::BedWriter - class for writing genomic regions in bed format

=head1 DESCRIPTION

This class writes out the supplied L<GenomicRegion|Sanger::CGP::Vagrent::Data::GenomicRegion> objects to a bed file.

Inherits from L<Sanger::CGP::Vagrent::IO::GenomicRegionWriter|Sanger::CGP::Vagrent::IO::GenomicRegionWriter>, view that for method details
