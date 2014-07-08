=head1 NAME

Sanger::CGP::Vagrent::Ontology::SequenceOntologyClassifier - Contains a set of rules to classify
Annotations using Sequence Ontology terms

=head1 DESCRIPTION

This class is designed to be extended, it contains a collection of test functions that classify
Sanger::CGP::Vagrent::Data::Annotation objects and assign the appropriate Sequence Ontology terms.

=cut

package Sanger::CGP::Vagrent::Ontology::SequenceOntologyClassifier;

use strict;

use Const::Fast qw(const);
use Cwd qw(abs_path);
use File::Basename;
use File::Path;
use List::Util qw(first);

use Log::Log4perl qw(:easy);

use File::ShareDir qw(module_dir);
use Config::IniFiles;

my $log = Log::Log4perl->get_logger(__PACKAGE__);

use base qw(Sanger::CGP::Vagrent);

# constant values holding SO class identifiers
const my $SO_SUBSTITUTION_CLASS => 'SO:1000002:substitution';
const my $SO_INSERTION_CLASS => 'SO:0000667:insertion';
const my $SO_DELETION_CLASS => 'SO:0000159:deletion';
const my $SO_COMPLEXINDEL_CLASS => 'SO:1000032:indel';

# sequence variant classes
const my $SO_UNKNOWN_VARIANT_CLASS => 'SO:0001576:transcript_variant';
const my $SO_SYN_VARIANT_CLASS => 'SO:0001588:synonymous_codon';
const my $SO_NON_SYN_VARIANT_CLASS => 'SO:0001583:non_synonymous_codon';
const my $SO_STOP_GAINED_VARIANT_CLASS => 'SO:0001587:stop_gained';
const my $SO_STOP_LOST_VARIANT_CLASS => 'SO:0001578:stop_lost';
const my $SO_STOP_RETAINED_VARIANT_CLASS => 'SO:0001567:stop_retained_variant';
const my $SO_START_LOST_VARIANT_CLASS => 'SO:0001582:initiator_codon_change';
const my $SO_PREMATURE_START_GAINED_VARIANT_CLASS => 'SO:0001988:5_prime_UTR_premature_start_codon_gain_variant';
const my $SO_INTRON_VARIANT_CLASS => 'SO:0001627:intron_variant';
const my $SO_CODON_VARIANT_CLASS => 'SO:0001581:codon_variant';
const my $SO_5PRIME_UTR_VARIANT_CLASS => 'SO:0001623:5_prime_UTR_variant';
const my $SO_3PRIME_UTR_VARIANT_CLASS => 'SO:0001624:3_prime_UTR_variant';
const my $SO_NC_TRANS_VARIANT_CLASS => 'SO:0001619:nc_transcript_variant';
const my $SO_ESS_SPLICE_VARIANT_CLASS => 'SO:0001629:splice_site_variant';
const my $SO_SPLICE_REGION_VARIANT_CLASS => 'SO:0001995:extended_intronic_splice_region_variant';
const my $SO_FRAMESHIFT_VARIANT_CLASS => 'SO:0001589:frameshift_variant';
const my $SO_INFRAME_VARIANT_CLASS => 'SO:0001650:inframe_variant';
const my $SO_INFRAME_CODON_GAIN_VARIANT_CLASS => 'SO:0001651:inframe_codon_gain';
const my $SO_INFRAME_CODON_LOSS_VARIANT_CLASS => 'SO:0001652:inframe_codon_loss';
const my $SO_COMPLEX_CHANGE_VARIANT_CLASS => 'SO:0001577:complex_change_in_transcript';

const my $SO_5KB_UPSTREAM_VARIANT_CLASS => 'SO:0001635:5KB_upstream_variant';
const my $SO_2KB_UPSTREAM_VARIANT_CLASS => 'SO:0001636:2KB_upstream_variant';
const my $SO_5KB_DOWNSTREAM_VARIANT_CLASS => 'SO:0001633:5KB_downstream_variant';
const my $SO_500BP_DOWNSTREAM_VARIANT_CLASS => 'SO:0001634:500B_downstream_variant';

const my $SO_EXON_CLASS => 'SO:0000147:exon';
const my $SO_INTRON_CLASS => 'SO:0000188:intron';
const my $SO_SPLICE_REGION_CLASS => 'SO:0001996:extended_intronic_splice_region';
const my $SO_ESS_SPLICE_SITE_CLASS => 'SO:0001993:extended_cis_splice_site';

# mature transcript mRNA region classes
const my $SO_CDS_CLASS => 'SO:0000316:CDS';
const my $SO_5PRIME_UTR_CLASS => 'SO:0000204:five_prime_UTR';
const my $SO_3PRIME_UTR_CLASS => 'SO:0000205:three_prime_UTR';

# gene attribute classes
const my $SO_PROTEIN_CODING_CLASS => 'SO:0000010:protein_coding';
const my $SO_NON_PROTEIN_CODING_CLASS => 'SO:0000011:non_protein_coding';

const my $TERM_SUMMARY_INI => 'SequenceOntologySummary.ini';

sub DESTROY {
  my $self = shift;
  if(defined $self->{'_SOsum'}){
    foreach my $k( sort {$self->{'_notSummary'}->{$b} <=> $self->{'_notSummary'}->{$a}} keys %{$self->{'_notSummary'}}){
      print $self->{'_notSummary'}->{$k},' - ',$k,"\n" unless $self->{'_notSummary'}->{$k} == 1;
    }      
  }
}



sub _ontologyInit {
  my $self = shift;
  my %vars = @_;
  if(exists $vars{'ontologySymmary'} && defined $vars{'ontologySymmary'}){
    $self->{'_SOsum'} = $vars{'ontologySymmary'};
  } else {
    $self->_loadOntologySummaryIni();
  }
}

sub _loadOntologySummaryIni {
  my $self = shift;
  my $share_path = dirname(abs_path($0)).'/../share';
  $share_path = module_dir('Sanger::CGP::Vagrent') unless(-e $share_path);
  $self->{'_SOsum'} = new Config::IniFiles( -file => File::Spec->catfile($share_path,$TERM_SUMMARY_INI));
}

sub getOntologySummary {
  my ($self,$anno) = @_;
  my $mrna = $anno->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext);
	my $cds = $anno->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext);
	my $prot = $anno->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext);
  my @class = $anno->getClassifications;  
  my @terms;
  foreach my $a($anno,$mrna,$cds,$prot){
    next unless defined $a;
    foreach my $c($a->getClassifications){
      next if first {$c eq $_} @terms;
      push @terms, $c;
    }
  }
  my $thing = join ',', @terms;
 
  my $sum = $self->{'_SOsum'}->val('SUMMARY',$thing,$thing);  
  if($sum eq $thing){
     $self->{'_notSummary'}->{$thing}++;
  }

  return $sum;
}

sub getSummaryCache {
  return shift->{'_SOsum'};
}

sub classifyTranscript {
	my ($self,$tran) = @_;
	if($tran->isProteinCoding){
		return $self->getProteinCodingClass;
	} else {
		return $self->getNonProteinCodingClass;
	}
}

# variant type classes

sub getSubstitutionClass {
	my $self = shift;
	return $SO_SUBSTITUTION_CLASS;
}
sub getInsertionClass {
	my $self = shift;
	return $SO_INSERTION_CLASS;
}
sub getDeletionClass {
	my $self = shift;
	return $SO_DELETION_CLASS;
}
sub getComplexIndelClass {
	my $self = shift;
	return $SO_COMPLEXINDEL_CLASS;
}
# sequence variant classes

sub getUnknownVariantClass {
	my $self = shift;
	return $SO_UNKNOWN_VARIANT_CLASS;
}
sub getIntronVariantClass {
	my $self = shift;
	return $SO_INTRON_VARIANT_CLASS;
}
sub getSynonymousVariantClass {
	my $self = shift;
	return $SO_SYN_VARIANT_CLASS;
}
sub getNonSynonymousVariantClass {
	my $self = shift;
	return $SO_NON_SYN_VARIANT_CLASS;
}
sub getStopGainedVariantClass {
	my $self = shift;
	return $SO_STOP_GAINED_VARIANT_CLASS;
}
sub getStopLostVariantClass {
	my $self = shift;
	return $SO_STOP_LOST_VARIANT_CLASS;
}
sub getStopRetainedVariantClass {
	my $self = shift;
	return $SO_STOP_RETAINED_VARIANT_CLASS;
}
sub getStartLostVariantClass {
	my $self = shift;
	return $SO_START_LOST_VARIANT_CLASS;
}
sub getPrematureStartGainedVariantClass{
	my $self = shift;
	return $SO_PREMATURE_START_GAINED_VARIANT_CLASS;
}
sub getCodonVariantClass{
	my $self = shift;
	return $SO_CODON_VARIANT_CLASS;
}
sub get5PrimeUtrVariantClass{
	my $self = shift;
	return $SO_5PRIME_UTR_VARIANT_CLASS;
}
sub get3PrimeUtrVariantClass{
	my $self = shift;
	return $SO_3PRIME_UTR_VARIANT_CLASS;
}
sub getNonCodingTranscriptVariantClass{
	my $self = shift;
	return $SO_NC_TRANS_VARIANT_CLASS;
}
sub getEssentialSpliceSiteVariantClass{
	my $self = shift;
	return $SO_ESS_SPLICE_VARIANT_CLASS;
}
sub getSpliceRegionVariantClass{
	my $self = shift;
	return $SO_SPLICE_REGION_VARIANT_CLASS;
}
sub getFrameShiftVariantClass{
	my $self = shift;
	return $SO_FRAMESHIFT_VARIANT_CLASS;
}
sub getInFrameVariantClass{
	my $self = shift;
	return $SO_INFRAME_VARIANT_CLASS;
}
sub getInFrameCodonGainVariantClass {
	my $self = shift;
	return $SO_INFRAME_CODON_GAIN_VARIANT_CLASS;
}
sub getInFrameCodonLossVariantClass {
	my $self = shift;
	return $SO_INFRAME_CODON_LOSS_VARIANT_CLASS;
}
sub getComplexChangeVariantClass {
	my $self = shift;
	return $SO_COMPLEX_CHANGE_VARIANT_CLASS;
}
sub get5KBUpStreamVariantClass {
	my $self = shift;
	return $SO_5KB_UPSTREAM_VARIANT_CLASS;
}
sub get2KBUpStreamVariantClass {
	my $self = shift;
	return $SO_2KB_UPSTREAM_VARIANT_CLASS;
}
sub get5KBDownStreamVariantClass {
	my $self = shift;
	return $SO_5KB_DOWNSTREAM_VARIANT_CLASS;
}
sub get500BPDownStreamVariantClass {
	my $self = shift;
	return $SO_500BP_DOWNSTREAM_VARIANT_CLASS;
}
# variant transcript region classes

sub getExonClass {
	my $self = shift;
	return $SO_EXON_CLASS;
}

sub getIntronClass {
	my $self = shift;
	return $SO_INTRON_CLASS;
}

sub getSpliceRegionClass {
	my $self = shift;
	return $SO_SPLICE_REGION_CLASS;
}

sub getEssentialSpliceSiteClass {
	my $self = shift;
	return $SO_ESS_SPLICE_SITE_CLASS;
}

# mature transcript mRNA region classes

sub getCDSClass {
	my $self = shift;
	return $SO_CDS_CLASS;
}

sub get5PrimeUtrClass {
	my $self = shift;
	return $SO_5PRIME_UTR_CLASS;
}

sub get3PrimeUtrClass {
	my $self = shift;
	return $SO_3PRIME_UTR_CLASS;
}

# gene/transcript type classes

sub getProteinCodingClass {
	my $self = shift;
	return $SO_PROTEIN_CODING_CLASS;
}

sub getNonProteinCodingClass {
	my $self = shift;
	return $SO_NON_PROTEIN_CODING_CLASS;
}

1;
