package Sanger::CGP::VagrentSV::Data::BreakPoint;

use strict;
use Data::Dumper;
use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Sanger::CGP::Vagrent::Data::AbstractGenomicPosition Sanger::CGP::Vagrent::Data::AbstractVariation );

1;

sub _init {
	my $self = shift;
	my %vars = @_;
	foreach my $k(keys(%vars)){
		if($k eq 'strand'){
			$self->{_strand} = $vars{strand};
		}
		elsif($k eq 'insseq') {
			$self->{_insseq} = $vars{insseq};
		}
	}
}


sub isValid {
	my $self = shift;
	return 0 unless(defined($self->{_minpos}) && defined($self->{_maxpos}));
	return 0 unless($self->{_minpos} <= $self->{_maxpos});
	return 0 unless($self->{_minpos} > 0);

	return 1;
}





