package Sanger::CGP::Vagrent::IO::GenomicRegionWriter;

use strict;
use Data::Dumper;
use Sanger::CGP::Vagrent qw($VERSION);
use base qw(Sanger::CGP::Vagrent);

1;

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
	my $self = {};
	bless($self, $class);
	$self->_init(@_);
	return $self;
}

sub _init {
	my $self = shift;
  my %vars = @_;
	foreach my $k(keys(%vars)){
		if($k eq 'fh'){
			$self->{_fh} = $vars{'fh'};
		}
	}
}

sub write {
	shift->throw_not_implemented;
}

sub _fh {
	return shift->{_fh};
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::IO::GenomicRegionWriter - Abstract class for genomic region writers

=head1 DESCRIPTION

Abstract class to hold basic shared functionality for genomic region writers

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $source = Sanger::CGP::Vagrent::IO::GenomicRegionWriterSubClass->new();

=item Function :

Builds a new Sanger::CGP::Vagrent::IO::GenomicRegionWriter inheriting object

=item Returns :

Sanger::CGP::Vagrent::IO::GenomicRegionWriter object initialized with parameter values

=item Params :

 fh => File handle to write to

=back

=head2 Abstract

=head3 write

=over

=item Usage :

 $writer->write($region);

=item Function :

Abstract function, must be overwritten. Writes the GenomicRegion out to the file

=item Returns :

Nothing

=item Params :

A L<Sanger::CGP::Vagrent::Data::GenomicRegion|Sanger::CGP::Vagrent::Data::GenomicRegion> object

=back
