package Sanger::CGP::VagrentSV::Data::StructuralVariation;

use strict;
use Data::Dumper;
use Const::Fast qw(const);

use Sanger::CGP::Vagrent qw($VERSION);
use Sanger::CGP::VagrentSV::SVConstants;
use Sanger::CGP::VagrentSV::Data::SVData;
use base qw(Sanger::CGP::Vagrent::Data::AbstractVariation);




1;


sub new {
	my($class, %sv)=@_;
	my $self = { };
	bless $self, $class;
	foreach my $key(keys(%sv)){
		 $self->{"_$key"} = $sv{$key};
		}
	return $self;
}

sub isValid {
	my $self = shift;
	return 0 unless(defined($self->{_lbreak}) && defined($self->{_rbreak}) && defined($self->{_name}));
	return 1;
}


sub getSvType {
	return shift->{_info}{svtype};
}

sub getSvName {
	return shift->{_info}{name};
}

sub getLocFlag {
	return shift->{_info}{locflag};
}

sub getLhb {
	return shift->{_lhb};
}
sub getRhb {
	return shift->{_rhb};
}

sub getInfo {
	return shift->{_info};
}



=head2  getSeq
get dan sequence for given interval
Inputs
=item Bio::DB::Sam index object
=over 2
=back
=cut

sub _getSeq {
	my ($fai,$chr,$start,$end)=@_;
	
	my $dna = $fai->fetch("$chr:$start-$end");
	return $dna;
}
