package Sanger::CGP::VagrentSV::Data::AbstractAnnotationSource;
																					 
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
use Sanger::CGP::VagrentSV::Base;

my $log = Log::Log4perl->get_logger(__PACKAGE__);


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
	$self->{'tabix'}=$self->_getTabixIndex($vars{'annotationSource'});
}


=head2 _getTabixIndex
create tabix object using tabix indexed bed file
Inputs
=over 2
=item options -user input paramaters
=back
=cut

sub _getTabixIndex {
	my ($self,$annotation_source)=@_;
	my $cmd;
	if(-e $annotation_source){
		$cmd= "bgzip $annotation_source";
		Sanger::CGP::VagrentSV::Base->_run_cmd($cmd);
	}
	if(-e "$annotation_source.gz" && (! -e "$annotation_source.gz.tbi") ){
		$cmd="tabix -p bed $annotation_source.gz";
		Sanger::CGP::VagrentSV::Base->_run_cmd($cmd);
	}
	elsif(! -e "$annotation_source.gz.tbi"){
		$log->logcroak("Unable to find input file:".$annotation_source);
	}
	my $tabix_obj = new Tabix(-data => $annotation_source.'.gz');
	$log->debug("Tabix object created successfully");
	return $tabix_obj;
}

sub _localInit: Abstract;

sub addMessage {
	my ($self,$msg) = @_;
	push(@{$self->{_msg}},ref($self).": ".$msg);
}

sub _debug {
	my $self = shift;
	if(exists($self->{_debug}) && defined($self->{_debug}) && $self->{_debug}){
		return 1;
	} else {
		return 0;
	}
}

sub getMessages {
	my $self = shift;
	return @{$self->{_msg}} if defined($self->{_msg});
	return undef;
}

sub _clearMessages {
	shift->{_msg} = undef;
}


__END__

=head1 NAME

Sanger::CGP::VagrentSV::Annotator::AbstractSVAnnotator - Abstract base class for the SV annotation generators

=head1 DESCRIPTION

This is an abstract template class for the SV annotators, it provides a lot of shared behind the scenes functionality.  All
subclasses must implement the _getAnnotation method.

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $source = Sanger::CGP::Vagrent::Annotators::AbstractAnnotatorSubClass->new(%params);

=item Function :

Builds a new Sanger::CGP::Vagrent::Annotators::AbstractAnnotator inheriting object

=item Returns :

Sanger::CGP::Vagrent::Annotators::AbstractAnnotator object initialized with parameter values

=item Params :

Hash of parameter values

 transcriptSource => A Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource inheriting object
 bookmarker       => (Optional) An array reference of, or single, Sanger::CGP::Vagrent::Bookmarkers::AbstractBookmarker inheriting object
 only_bookmarked  => (Optional) Boolean, only return annotations that get bookmarked

=back



