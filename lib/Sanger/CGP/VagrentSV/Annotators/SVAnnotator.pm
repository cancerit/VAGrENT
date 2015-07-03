package Sanger::CGP::VagrentSV::Annotators::SVAnnotator;

use strict;
use Data::Dumper;
use Const::Fast qw(const);
use Try::Tiny qw(try catch);

use Sanger::CGP::Vagrent qw($VERSION);
use Sanger::CGP::VagrentSV::SVConstants;
use Sanger::CGP::VagrentSV::Data::SVData;

use base qw(Sanger::CGP::VagrentSV::Annotators::AbstractSVAnnotator);

1;


sub _myinit {
	my $self = shift;
	$self->{'reserved'}='reserved';
}


sub getSVAnnotations{
	my($self,$sv)=@_;
	#get transcripts overlapping SV regions
	my @tr= $self->{'_transcriptSource'}->getTranscripts($sv->getLhb);
	foreach my $t(@tr){
		my $res = $self->_generateResults($sv->getLhb,$t);
		#push @results, $res if defined $res;
	}

	#my($bp_genes,$bp_tr)=_getSvFeatures($bp,\@tr,$annotator);
	#print Dumper $bp_genes;
	#return $sv_data;
}

=head1

sub _getSvFeatures {
	my ($self,$bp,$tr)=@_;
	my ($tr_list,$gene_list);
	foreach my $tr (@$tr) {
		try {
   my $g=$self->{_subannot}->getTranscriptAnnotation($bp,$tr);
  } catch {
  		warn join "\n", $self->{_subannot}->getMessages;
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
=cut

sub _generateResults {
	my($self,$bp,$t)=@_;
	my $g=$self->getSVannotator->getTranscriptAnnotation($bp,$t);
	
	print Dumper $g;

}































