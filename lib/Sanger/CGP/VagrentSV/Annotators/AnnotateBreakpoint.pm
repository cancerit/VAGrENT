package Sanger::CGP::VagrentSV::Annotators::AnnotateBreakpoint;

use strict;
use Data::Dumper;
use Const::Fast qw(const);

use Sanger::CGP::Vagrent qw($VERSION);
use Sanger::CGP::VagrentSV::SVConstants;
use Sanger::CGP::VagrentSV::Data::SVData;
use base qw(Sanger::CGP::Vagrent::Data::AbstractVariation);




1;


sub new {
	my($class,$bp,$ts,$annotator)=@_;
	my $self = { };
	bless $self, $class;
	$self->getBreakpointGenes($bp,$ts,$annotator);
	return $self;
}

sub getBreakpointGenes{
	my($self,$bp,$ts,$annotator)=@_;
	#get transcripts overlapping SV regions
	my @tr=$ts->getTranscripts($bp);
	my($bp_genes,$bp_tr)=_getSvFeatures($bp,\@tr,$annotator);
	print Dumper $bp_genes;
	#return $sv_data;
}

sub _getSvFeatures {
	my ($bp,$tr,$annotator)=@_;
	my ($tr_list,$gene_list);
	foreach my $tr (@$tr) {
		try {
   my $g=$annotator->getTranscriptAnnotation($bp,$tr);
  } catch {
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
