package Sanger::CGP::Vagrent::Data::Transcript;

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
use Carp;
use Data::Dumper;
use Sanger::CGP::Vagrent qw($VERSION);
use Sort::Key qw(nkeysort);
use Const::Fast qw(const);

use base qw(Sanger::CGP::Vagrent);

const my $PROTEIN_CODING_TYPE => 'ProteinCoding';
const my $MIRCO_RNA_TYPE => 'miRNA';
const my $LINC_RNA_TYPE => 'lincRNA';
const my $SNO_RNA_TYPE => 'snoRNA';
const my $SN_RNA_TYPE => 'snRNA';
const my $R_RNA_TYPE => 'rRNA';

const my @ALL_GENE_TYPES => ($PROTEIN_CODING_TYPE,$MIRCO_RNA_TYPE,$LINC_RNA_TYPE,$SNO_RNA_TYPE,$SN_RNA_TYPE,$R_RNA_TYPE);

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
		if($k eq 'db'){
			$self->{_db} = $vars{db};
		} elsif($k eq 'dbversion'){
			$self->{_dbversion} = $vars{dbversion};
		} elsif($k eq 'exons' && ref($vars{exons}) eq 'ARRAY'){
			$self->{_exons} = $vars{exons};
		} elsif($k eq 'strand'){
			$self->{_strand} = $vars{strand};
		} elsif($k eq 'genomicminpos'){
			$self->{_genomicminpos} = $vars{genomicminpos};
		} elsif($k eq 'genomicmaxpos'){
			$self->{_genomicmaxpos} = $vars{genomicmaxpos};
		} elsif($k eq 'cdnaseq'){
			$self->{_cdnaseq} = $vars{cdnaseq};
		} elsif($k eq 'cdsphase'){
			$self->{_cdsphase} = $vars{cdsphase};
		} elsif($k eq 'cdsminpos'){
			$self->{_cdsminpos} = $vars{cdsminpos};
		} elsif($k eq 'cdsmaxpos'){
			$self->{_cdsmaxpos} = $vars{cdsmaxpos};
		} elsif($k eq 'acc'){
			$self->{_acc} = $vars{acc};
		} elsif($k eq 'accversion'){
			$self->{_accversion} = $vars{accversion};
		} elsif($k eq 'proteinacc'){
			$self->{_proteinacc} = $vars{proteinacc};
		} elsif($k eq 'proteinaccversion'){
			$self->{_proteinaccversion} = $vars{proteinaccversion};
		} elsif($k eq 'ccds'){
			$self->{_ccds} = $vars{ccds};
		} elsif($k eq 'genename'){
			$self->{_genename} = $vars{genename};
		} elsif($k eq 'genetype'){
			$self->{_genetype} = $vars{genetype};
			my $good = undef;
			foreach my $type(@ALL_GENE_TYPES){
				if($vars{genetype} eq $type){
					$good = $vars{genetype};
					last;
				}
			}
			if(defined($good)){
				$self->{_genetype} = $good;
			} else {
				croak('recieved unknown gene type: '.$vars{genetype});
			}
		}
	}
}

sub getAllGeneTypes {
  return @ALL_GENE_TYPES;
}

sub getGeneType {
	return shift->{_genetype};
}

sub getCCDS {
	return shift->{_ccds};
}

sub getGeneName {
	return shift->{_genename};
}

sub getAccession {
	return shift->{_acc};
}

sub getAccessionVersion {
	return shift->{_accversion};
}

sub getProteinAccession {
	return shift->{_proteinacc};
}

sub getProteinAccessionVersion {
	return shift->{_proteinaccversion};
}

sub getDatabase {
	return shift->{_db};
}

sub getDatabaseVersion {
	return shift->{_dbversion};
}

sub getExons {
  my $self = shift;
  return nkeysort {$_->getRnaMinPos} @{$self->{_exons}};
}

sub getExonsGenomicOrder {
  my $self = shift;
  return nkeysort {$_->getMinPos} @{$self->{_exons}};
}

sub getStrand {
	return shift->{_strand};
}

sub getcDNASeq {
	return shift->{_cdnaseq};
}

sub getCdsSeq {
	my $self = shift;
	return substr($self->{_cdnaseq},($self->{_cdsminpos} - 1),(($self->{_cdsmaxpos} - $self->{_cdsminpos}) + 1));
}

sub getCdsPhase {
	return shift->{_cdsphase};
}

sub getCdsMinPos {
	return shift->{_cdsminpos};
}

sub getCdsMaxPos {
	return shift->{_cdsmaxpos};
}

sub getCdsLength {
	my $self = shift;
	unless(defined($self->{_cdslength})){
		$self->{_cdslength} = ($self->getCdsMaxPos - $self->getCdsMinPos) + 1;
	}
	return $self->{_cdslength};
}

sub getmRNALength {
  my $self = shift;
  unless(defined $self->{_mrnalength}){
    my $size = 0;
    foreach my $e (@{$self->{_exons}}){
      $size += $e->getLength();
    }
    $self->{_mrnalength} = $size;
  }
  return $self->{_mrnalength};
}


sub getGenomicMinPos {
	return shift->{_genomicminpos};
}

sub getGenomicMaxPos {
	return shift->{_genomicmaxpos};
}

sub isProteinCoding {
	if(shift->getGeneType eq $PROTEIN_CODING_TYPE){
		return 1;
	} else {
		return 0;
	}
}

sub getProteinCodingType {
	return $PROTEIN_CODING_TYPE;
}

sub getMicroRnaType {
	return $MIRCO_RNA_TYPE;
}

sub getLincRnaType {
	return $LINC_RNA_TYPE;
}

sub getSnoRnaType {
	return $SNO_RNA_TYPE;
}

sub getSnRnaType {
	return $SN_RNA_TYPE;
}

sub getRRnaType {
	return $R_RNA_TYPE;
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::Data::Transcript - Data object representing a transcript

=head1 DESCRIPTION

This is a data class to hold details about a transcript. This allows the transcript information
to be abstracted away from its original source.

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $trans = Sanger::CGP::Vagrent::Data::Transcript->new(%params);

=item Function :

Builds a new Sanger::CGP::Vagrent::Data::Transcript object

=item Returns :

Sanger::CGP::Vagrent::Data::Transcript object initialized with parameter values

=item Params :

 db                => source DB the transcript came from
 dbversion         => version of the source DB
 exons             => array ref of Sanger::CGP::Vagrent::Data::Exon objects
 strand            => Genomic strand of the transcript (1 or -1)
 cdnaseq           => cDNA sequence string
 cdsminpos         => CDS minimum coordinate in the cDNA sequence
 cdsmaxpos         => CDS maximum coordinate in the cDNA sequence
 cdsphase          => The phase of translation on the CDS (ie the number of N's to add to the begining of the CDS to make it translate)
 acc               => Transcript accession number (typically the cDNA accession)
 accversion        => Version of the transcript accession number
 proteinacc        => Protein accession number (optional)
 proteinaccversion => Version of the protein accession number (optional)
 ccds              => CCDS id for the transcript (optional)
 genename          => Gene name
 genetype          => Gene type (defined by type constant)


=back

=head2 Attributes

=head3 getGeneName

=over

=item Usage :

 my $name = $trans->getGeneName;

=item Function

The gene name for the transcript

=item Returns

String

=back

=head3 getGeneType

=over

=item Usage :

 my $type = $trans->getGeneType;

=item Function

Returns the gene type, will be equal to one of the type constant values

=item Returns

String

=back

=head3 getCCDS

=over

=item Usage :

 my $ccds = $trans->getCCDS;

=item Function

The CCDS identifier for the transcript

=item Returns

String

=back

=head3 getAccession

=over

=item Usage :

 my $accession = $trans->getAccession;

=item Function

The accession of the sequence record for the transcript

=item Returns

String

=back

=head3 getAccessionVersion

=over

=item Usage :

 my $accessionVers = $trans->getAccessionVersion;

=item Function

The accession verson of the sequence record for the transcript

=item Returns

String

=back

=head3 getProteinAccession

=over

=item Usage :

 my $protAccession = $trans->getProteinAccession;

=item Function

The protein accession of the sequence record for the transcript

=item Returns

String

=back

=head3 getProteinAccessionVersion

=over

=item Usage :

 my $protAccessionVers = $trans->getProteinAccessionVersion;

=item Function

The protein accession verson of the sequence record for the transcript

=item Returns

String

=back

=head3 getDatabase

=over

=item Usage :

 my $sourceDB = $trans->getDatabase;

=item Function

The name of database the transcript originated from

=item Returns

String

=back

=head3 getDatabaseVersion

=over

=item Usage :

 my $sourceDBversion = $trans->getDatabaseVersion;

=item Function

The version of database the transcript originated from

=item Returns

String

=back

=head3 getGenomicMinPos

=over

=item Usage :

 my $genomicMin = $trans->getGenomicMinPos;

=item Function

The minimum position of the cDNA in the genome

=item Returns

Integer

=back

=head3 getGenomicMaxPos

=over

=item Usage :

 my $genomicMax = $trans->getGenomicMaxPos;

=item Function

The maximum position of the cDNA in the genome

=item Returns

Integer

=back

=head3 getStrand

=over

=item Usage :

 my $strand = $trans->getStrand;

=item Function

The genomic strand of the transcript

=item Returns

Integer (either 1 or -1)

=back

=head3 getCdsMinPos

=over

=item Usage :

 my $cdsMin = $trans->getCdsMinPos;

=item Function

The minimum position of the CDS within the cDNA sequence

=item Returns

Integer

=back

=head3 getCdsMaxPos

=over

=item Usage :

 my $cdsMax = $trans->getCdsMaxPos;

=item Function

The maximum position of the CDS within the cDNA sequence

=item Returns

Integer

=back

=head3 getcDNASeq

=over

=item Usage :

 my $cDNAseq = $trans->getcDNASeq;

=item Function

The cDNA sequence string

=item Returns

String - DNA sequence

=back

=head3 getCdsPhase

=over

=item Usage :

 my $phase = $trans->getCdsPhase;

=item Function

The CDS phase, ie the number of N's needed to append to the start of CDS sequence to create the right translation frame

=item Returns

Integer

=back

=head3 getExons

=over

=item Usage :

 my @exons = $trans->getExons;

=item Function

Returns the Sanger::CGP::Vagrent::Data::Exon objects that make up the transcript sorted into transcript order

=item Returns

Array of Sanger::CGP::Vagrent::Data::Exon objects

=back

=head2 Functions

=head3 isProteinCoding

=over

=item Usage :

 my $isCoding = $trans->isProteinCoding;

 if($trans->isProteinCoding){
   ......
 }

=item Function

Convenience function. Returns true if the gene type is protein coding, equivalent to

 if($trans->getType eq $trans->getProteinCodingType){
   ......
 }


=item Returns

Integer boolean (either 1 or 0)

=back

=head3 getCdsSeq

=over

=item Usage :

 my $CDSseq = $trans->getCdsSeq;

=item Function

CDS sequence string, ie the subsequence of the cDNA string defined by the CDS min and max coordinates

=item Returns

String - DNA sequence

=back

=head3 getCdsLength

=over

=item Usage :

 my $length = $trans->getCdsLength;

=item Function

The length of the CDS, based on the CDS min and max coordinates on the cDNA

=item Returns

Integer

=back

=head2 Constants

=head3 getProteinCodingType

=over

=item Usage :

 my $type = $trans->getProteinCodingType();
 my $type = Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType();

=item Function :

Constant lookup, returns the protein coding gene type

=item Returns :

String

=back

=head3 getMicroRnaType

=over

=item Usage :

 my $type = $trans->getMicroRnaType();
 my $type = Sanger::CGP::Vagrent::Data::Transcript::getMicroRnaType();

=item Function :

Constant lookup, returns the microRNA gene type

=item Returns :

String

=back

=head3 getLincRnaType

=over

=item Usage :

 my $type = $trans->getLincRnaType();
 my $type = Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType();

=item Function :

Constant lookup, returns the lincRNA gene type

=item Returns :

String

=back

=head3 getSnoRnaType

=over

=item Usage :

 my $type = $trans->getSnoRnaType();
 my $type = Sanger::CGP::Vagrent::Data::Transcript::getSnoRnaType();

=item Function :

Constant lookup, returns the snoRNA gene type

=item Returns :

String

=back

=head3 getSnRnaType

=over

=item Usage :

 my $type = $trans->getSnRnaType();
 my $type = Sanger::CGP::Vagrent::Data::Transcript::getSnRnaType();

=item Function :

Constant lookup, returns the snRNA gene type

=item Returns :

String

=back

=head3 getRRnaType

=over

=item Usage :

 my $type = $trans->getRRnaType();
 my $type = Sanger::CGP::Vagrent::Data::Transcript::getRRnaType();

=item Function :

Constant lookup, returns the rRNA gene type

=item Returns :

String

=back
