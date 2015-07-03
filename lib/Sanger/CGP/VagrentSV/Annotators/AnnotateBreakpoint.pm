package Sanger::CGP::VagrentSV::Annotators::AnnotateBreakpoint;

use strict;
use Data::Dumper;
use Const::Fast qw(const);
use Try::Tiny qw(try catch);

use Sanger::CGP::Vagrent qw($VERSION);
use Sanger::CGP::VagrentSV::SVConstants;
use Sanger::CGP::VagrentSV::Data::SVData;
#use base qw(Sanger::CGP::Vagrent::Data::AbstractVariation);




1;


sub new {
	my($class,$bp,$ts,$annotator)=@_;
	my $self = { };
	bless $self, $class;
	$self->_init(@_);
	#$self->getBreakpointGenes($bp,$ts,$annotator);
	return $self;
}

sub _init {
	my $self = shift;
	my %vars = @_;
	foreach my $k(keys(%vars)){
		if($k eq 'transcriptSource' && $vars{'transcriptSource'}->isa('Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource')){
			$self->{'_transcriptSource'} = $vars{'transcriptSource'};
		} elsif($k eq 'debug' && $vars{'debug'}){
			$self->{'_debug'} = 1;
		}
	}
	$self->{_subannot} = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);
}

sub getBreakpointGenes{
	my($self,$bp,$ts,$annotator)=@_;
	#get transcripts overlapping SV regions
	my @tr=$ts->getTranscripts($bp);
	
=head	
	foreach my $t(@tr){
		my $res = $self->_generateResult($bp,$t);
		push @results, $res if defined $res;
	}
	
	
=cut	
	
	my($bp_genes,$bp_tr)=_getSvFeatures($bp,\@tr,$annotator);
	print Dumper $bp_genes;
	#return $sv_data;
}

sub _getSvFeatures {
	my ($self,$bp,$tr)=@_;
	my ($tr_list,$gene_list);
	foreach my $tr (@$tr) {
		try {
   my $g=$self->{_subannot}->getTranscriptAnnotation($bp,$tr);
  } catch {
  		warn join "\n", $annotator->getMessages;
      warn "caught error: $_\n"; # not $@
  }; 				
		
		exit;
		if($tr->getAccession) {
			$tr_list->{$tr->getAccession}++;
			$gene_list->{$tr->getGeneName}++;
	  }		
	}
	return ($gene_list,$tr_list);
}
