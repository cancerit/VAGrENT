package Sanger::CGP::Vagrent::IO::AnnotationWriter;

use strict;
use Data::Dumper;

use base qw(Sanger::CGP::Vagrent::Logger);

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

Sanger::CGP::Vagrent::IO::AnnotationWriter - Abstract class for annotation writers

=head1 DESCRIPTION

Abstract class to hold basic shared functionality for annotation writers

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $source = Sanger::CGP::Vagrent::IO::AnnotationWriterSubClass->new();

=item Function :

Builds a new Sanger::CGP::Vagrent::IO::AnnotationWriter inheriting object

=item Returns :

Sanger::CGP::Vagrent::IO::AnnotationWriter object initialized with parameter values

=item Params :

 fh => File handle to write to

=back

=head2 Abstract

=head3 write

=over

=item Usage :

 $writer->write($var,$annoGroup);

=item Function :

Abstract function, must be overwritten. Writes the annotation out to the file

=item Returns :

Nothing

=item Params :

A L<Sanger::CGP::Vagrent::Data::AbstractVariation|Sanger::CGP::Vagrent::Data::AbstractVariation> inheriting object

A L<Sanger::CGP::Vagrent::Data::AnnotationGroup|Sanger::CGP::Vagrent::Data::AnnotationGroup> object

=back

