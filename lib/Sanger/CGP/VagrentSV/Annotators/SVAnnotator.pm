package Sanger::CGP::VagrentSV::Annotators::SVAnnotator;

use strict;
use Data::Dumper;

use Const::Fast qw(const);
use Try::Tiny qw(try catch);

use Bio::Seq;
use Bio::SeqUtils;

use Sanger::CGP::Vagrent qw($VERSION);
use Sanger::CGP::VagrentSV::SVConstants;
use Sanger::CGP::Vagrent::Data::Insertion;
use Sanger::CGP::VagrentSV::Results::TranscriptResults;
use base qw(Sanger::CGP::VagrentSV::Annotators::AbstractSVAnnotator);

1;


sub _localInit {
	my $self = shift;
	#$self->{'reserved'}='reserved';
}


sub getBreakpointTranscripts {
	my($self,$sv,$bptr)=@_;
		my $ltr=$self->_getTranscriptData($bptr,$sv->getLhb);
		my $rtr=$self->_getTranscriptData($bptr,$sv->getRhb);	
	return($ltr,$rtr);
}

sub _getTranscriptData {
	my($self,$bptr,$bp)=@_;
	my $results;
	my $tmp_ins=$self->_getInsertionObject($bp);
	foreach my $t ($self->{'_transcriptSource'}->getTranscripts($bp)) {
		my $bptrData = $bptr->getTranscriptAnnotation($t,$bp->getMaxPos,$tmp_ins);
		my $res = Sanger::CGP::VagrentSV::Results::TranscriptResults->new(%$bptrData);
		push(@$results,$res);
	}
	return $results;
}

sub _getInsertionObject {
		my ($self,$bp)=@_;
    my $var = Sanger::CGP::Vagrent::Data::Insertion->new(
							'species'				=> $bp->getSpecies,
							'genomeVersion' => $bp->getGenomeVersion,
							'chr'	          => $bp->getChr,
							'minpos'        => $bp->getMinPos,
							'maxpos'        => $bp->getMinPos + 1,
							'insseq'        => 'A'); # as alt always includes the reference base prior to the change
 	
	return $var;

}



















