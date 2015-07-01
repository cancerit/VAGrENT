package Sanger::CGP::VagrentSV::Data::SVData;

use strict;
use Data::Dumper;
use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Sanger::CGP::Vagrent::Data::AbstractVariation );

1;

sub new {
	my($class, %sv_data)=@_;
	my $self = { };
	bless $self, $class;
	foreach my $key(keys %sv_data){
		 $self->{"_$key"} = $sv_data{$key};
	}
	return $self;
}

sub getFusedTr {
	return shift->{_fusedtr};
}

sub getFusedGenes {
	return shift->{_fusedgenes};
}

sub getInbTr {
	return shift->{_inbtr};
}

sub getInbGenes {
	return shift->{_inbgenes};
}

sub getRhbGenes {
	return shift->{_rhbgenes};
}

sub getLhbGenes {
	return shift->{_lhbgenes};
}

sub getInbDeletedSequence {
	return shift->{_inb_delseq};
}

sub getLhbInsertedSequence {
	return shift->{_lhb_insseq};
}

sub getRhbInsertedSequence {
	return shift->{_rhb_insseq};
}

sub getLhbTranscriptObj {
	return shift->{_lhb_trobj};
}

sub getRhbTranscriptObj {
	return shift->{_rhb_trobj};
}


