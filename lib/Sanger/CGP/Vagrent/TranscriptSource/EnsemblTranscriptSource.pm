package Sanger::CGP::Vagrent::TranscriptSource::EnsemblTranscriptSource;

use strict;

use Sanger::CGP::Vagrent::Data::Transcript;
use Sanger::CGP::Vagrent::Data::Exon;

use Carp qw(cluck);
use Log::Log4perl;
use Cwd qw(abs_path);
use Data::Dumper;
use Config::IniFiles;
use File::Spec;

use base qw(Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource);

my $log = Log::Log4perl->get_logger(__PACKAGE__);

1;

use constant SPECIES_INI_PARAMS => qw(transStatus transBiotype transNames protNames);

sub getTranscripts {
	my ($self,$gp) = @_;
	unless(defined($gp) && $gp->isa('Sanger::CGP::Vagrent::Data::AbstractGenomicPosition')){
		$log->error("Did not recieve a Sanger::CGP::Vagrent::Data::AbstractGenomicPosition object");
		return undef;
	}
	my $slice = $self->_getSlice($gp);
	my $rawTransList = $self->_getTranscriptAdaptor($gp)->fetch_all_by_Slice($slice);
	my $filteredTransList = $self->_filter($rawTransList,$gp);
	return $self->_convert($filteredTransList,$gp);
}

sub getTranscriptsForNextGeneInDumpRegion {
	my ($self) = @_;
	$self->throw('must define a dump region before trying to loop over the genes') unless defined $self->getDumpRegion;
	my $gr = $self->getDumpRegion;
	unless(defined $self->{_dumpInfo}->{_geneList}){
		my $sa = $self->_getSliceAdaptor($gr);
		my $ga = $self->_getGeneAdaptor($gr);
		my $slice;
		if($gr->getMinPos == $gr->getMaxPos && $gr->getMaxPos == 0){
			# full sequence region scan
			$slice = $sa->fetch_by_region('chromosome',$gr->getChr,undef,undef,1,$gr->getGenomeVersion);
			unless(defined $slice){
				$slice = $sa->fetch_by_region(undef,$gr->getChr,undef,undef,1,$gr->getGenomeVersion);
			}
			$self->{_dumpInfo}->{_fullSeq} = 1;
		} else {
			$slice = $sa->fetch_by_region('chromosome',$gr->getChr,$gr->getMinPos,$gr->getMaxPos,1,$gr->getGenomeVersion);
			unless(defined $slice){
				$slice = $sa->fetch_by_region(undef,$gr->getChr,$gr->getMinPos,$gr->getMaxPos,1,$gr->getGenomeVersion);
			}
			if($gr->getMinPos == 1){
				my $tmpSlice = $sa->fetch_by_region('chromosome',$gr->getChr,undef,undef,1,$gr->getGenomeVersion);
				unless(defined $slice){
					$tmpSlice = $sa->fetch_by_region(undef,$gr->getChr,undef,undef,1,$gr->getGenomeVersion);
				}
				unless(defined $tmpSlice){
					return;
				}
				if($tmpSlice->length == $gr->getMaxPos){
					$self->{_dumpInfo}->{_fullSeq} = 1;
				} else {
					$self->{_dumpInfo}->{_fullSeq} = 0;
				}
			} else {
				$self->{_dumpInfo}->{_fullSeq} = 0;
			}

		}
		$self->{_dumpInfo}->{_geneList} = $ga->fetch_all_by_Slice($slice);
		$self->{_dumpInfo}->{_counter} = 0;
	}

	for(my $i = $self->{_dumpInfo}->{_counter} ; $i <= scalar(@{$self->{_dumpInfo}->{_geneList}}) ; $i++){
		$self->{_dumpInfo}->{_counter}++;
		next unless(defined $self->{_dumpInfo}->{_geneList}->[$i]);
		my $filteredTransList = $self->_filter($self->{_dumpInfo}->{_geneList}->[$i]->get_all_Transcripts,$gr);
		next unless(defined $filteredTransList && scalar(@$filteredTransList) > 0);
		return $self->_convert($filteredTransList,$gr);
	}
	return;
}

sub _convert {
	my ($self,$tList,$gp) = @_;
	return undef unless(defined($tList) && scalar(@$tList) > 0);
	my @converted;
	foreach my $rawT(@$tList){
		my @exons;
# building Sanger::CGP::Vagrent::Data::Exon objects from ensembl exons
		foreach my $rawE (@{$rawT->get_all_Exons()}){
			my $convE = Sanger::CGP::Vagrent::Data::Exon->new(
								species => $gp->getSpecies,
								genomeVersion => $gp->getGenomeVersion,
								chr => $rawE->slice->seq_region_name,
								minpos => $rawE->seq_region_start,
								maxpos => $rawE->seq_region_end,
								rnaminpos => $rawE->cdna_start($rawT),
								rnamaxpos => $rawE->cdna_end($rawT),);
			push(@exons,$convE);
		}

		my $type = $self->_getGeneTypeForTranscript($rawT);
		my $protAcc = undef;
		my $protAccVers = undef;
		if($self->_hasProtein($type)){
			my $tslat = $rawT->translation;
			if(defined($tslat)){
				$protAcc = $tslat->stable_id;
				$protAccVers = 1;
			} else {
				# some Ensembl weirdness here, occasionally a protein coding gene doesn't have a translation!
				# skip it, its not right.
				next;
			}
		}

		my $cdsPhase = $self->_calculateCdsPhase($rawT);

# build an Sanger::CGP::Vagrent::Data::Transcript from the ensembl transcript

		my @sortedExons = sort {$a->getMinPos <=> $b->getMinPos} @exons;

		my $convT = Sanger::CGP::Vagrent::Data::Transcript->new(
																			db => 'Ensembl',
																			dbversion => $rawT->slice->adaptor->db->dbc->dbname(),
																			acc => $rawT->stable_id,
																			accversion => 1,
																			proteinacc => $protAcc,
																			proteinaccversion => $protAccVers,
																			ccds => $self->_getCCDSNameForTranscript($rawT),
																			genename => $self->_getGeneNameForTranscript($gp,$rawT,$type),
																			genetype => $type,
																			strand => $rawT->seq_region_strand,
																			cdnaseq => $rawT->seq->seq,
																			cdsminpos => $rawT->cdna_coding_start,
																			cdsmaxpos => $rawT->cdna_coding_end,
																			cdsphase => $cdsPhase,
																			genomicminpos => $sortedExons[0]->getMinPos,
																			genomicmaxpos => $sortedExons[-1]->getMaxPos,
																			exons => \@exons,);
		push(@converted,$convT);
	}
	return @converted;
}

sub _calculateCdsPhase {
	my ($self,$rawT) = @_;
	# check the transcript is translatable
	my $tmp = $rawT->translateable_seq;
	unless(defined($tmp) && length($tmp) > 0){
		return -1;
	}
	# count the number of N's at the start of the translateable_seq
	$tmp =~ m/^(n+)/i;
	my $match = $1;
	if(defined($match)){
		return length($match);
	} else {
		return 0;
	}
}

sub _hasProtein {
	my ($self,$type) = @_;
	if($type eq Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType()){
		return 1;
	}
	return 0;
}

sub _getGeneTypeForTranscript {
	my ($self,$t) = @_;
	if($t->biotype eq 'protein_coding'){
		return Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType();
	} elsif($t->biotype eq 'miRNA'){
		return Sanger::CGP::Vagrent::Data::Transcript::getMicroRnaType();
	} elsif($t->biotype eq 'lincRNA'){
		return Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType();
	} elsif($t->biotype eq 'snoRNA'){
		return Sanger::CGP::Vagrent::Data::Transcript::getSnoRnaType();
	} elsif($t->biotype eq 'snRNA'){
		return Sanger::CGP::Vagrent::Data::Transcript::getSnRnaType();
	} elsif($t->biotype eq 'rRNA'){
		return Sanger::CGP::Vagrent::Data::Transcript::getRRnaType();
	} else {
		$self->throw("Unhandled type ".$t->biotype);
	}
}

sub _getGeneNameForTranscript {
	my ($self,$gp,$t,$type) = @_;
	my $ga = $self->_getGeneAdaptor($gp);
	my $g = $ga->fetch_by_transcript_stable_id($t->stable_id);
	my @links = @{$g->get_all_DBEntries};

	my @protNames = $self->_getProteinNameTypesForSpecies($gp);
	my @rnaNames = $self->_getTranscriptNameTypesForSpecies($gp);

	if($type eq Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType()){
		foreach my $gn(@protNames){
			foreach my $l(@links){
				return $l->display_id if($l->dbname eq $gn);
			}
		}
	} elsif($type eq Sanger::CGP::Vagrent::Data::Transcript::getMicroRnaType() ||
					$type eq Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType() ||
					$type eq Sanger::CGP::Vagrent::Data::Transcript::getSnoRnaType() ||
					$type eq Sanger::CGP::Vagrent::Data::Transcript::getSnRnaType() ||
					$type eq Sanger::CGP::Vagrent::Data::Transcript::getRRnaType()) {
		foreach my $gn(@rnaNames){
			foreach my $l(@links){
				return $l->display_id if($l->dbname eq $gn);
			}
		}
	} else {
		foreach my $l(@links){
			print $l->dbname.' - '.$l->display_id."\n";
		}
		$self->throw("unhandled type ".$type);
	}
	return $g->stable_id;
}

sub _getCCDSNameForTranscript {
	my ($self,$t) = @_;
	my @links = @{$t->get_all_DBEntries};
	foreach my $l(@links){
		return $l->display_id if($l->dbname eq 'CCDS');
	}
	return '';
}

sub _filter {
	my ($self,$rawList,$gp) = @_;
	my $goodList;
	ENSTRANFILTER: foreach my $t(@$rawList){
		my $gdStatus = 0;
		foreach my $gst ($self->_getTranscriptStatusesForSpecies($gp)){
			if($t->status eq $gst){
				$gdStatus = 1;
				last;
			}
		}
		next unless($gdStatus == 1);
		foreach my $gbt ($self->_getTranscriptBiotypesForSpecies($gp)){
			if($t->biotype eq $gbt){
				# we want this one
				push(@$goodList,$t);
				last;
			}
		}
	}
	return $goodList;
}

sub _getSlice {
	my ($self,$gp) = @_;
	my $sa = $self->_getSliceAdaptor($gp);
	my $min = $gp->getMinPos - 10000;
	my $max = $gp->getMaxPos + 10000;
	my $s = $sa->fetch_by_region('chromosome',$gp->getChr,$min,$max,1,$gp->getGenomeVersion);
	unless(defined($s)){
		$s = $sa->fetch_by_region(undef,$gp->getChr,$min,$max,1,$gp->getGenomeVersion);
	}
	return $s;
}

sub _getTranscriptAdaptor {
	my ($self,$gp) = @_;
	return $self->_registry()->get_adaptor($self->_getSanitisedSpecies($gp),'core','transcript');
}

sub _getSliceAdaptor {
	my ($self,$gp) = @_;
	return $self->_registry()->get_adaptor($self->_getSanitisedSpecies($gp),'core','slice');
}

sub _getGeneAdaptor {
	my ($self,$gp) = @_;
	return $self->_registry()->get_adaptor($self->_getSanitisedSpecies($gp),'core','gene');
}

sub _registry {
	return shift->{_registry};
}

sub _getSanitisedSpecies {
	my ($self,$gp) = @_;
	$self->throw("Species ".$gp->getSpecies()." is not recognised, are you sure its expected") unless defined $self->{_speciesAlias}->{lc($gp->getSpecies())};
	return $self->{_speciesAlias}->{lc($gp->getSpecies())};
}

sub _getProteinNameTypesForSpecies {
	my ($self,$gp) = @_;
	return @{$self->{_speciesConfig}->{$self->_getSanitisedSpecies($gp)}->{'protNames'}};
}

sub _getTranscriptNameTypesForSpecies {
	my ($self,$gp) = @_;
	return @{$self->{_speciesConfig}->{$self->_getSanitisedSpecies($gp)}->{'transNames'}};
}

sub _getTranscriptStatusesForSpecies {
	my ($self,$gp) = @_;
	return @{$self->{_speciesConfig}->{$self->_getSanitisedSpecies($gp)}->{'transStatus'}};
}

sub _getTranscriptBiotypesForSpecies {
	my ($self,$gp) = @_;
	return @{$self->{_speciesConfig}->{$self->_getSanitisedSpecies($gp)}->{'transBiotype'}};
}

sub _init {
	my $self = shift;
  	my %vars = @_;
	foreach my $k(keys(%vars)){
		if($k eq 'registry' && $vars{'registry'}->isa('Bio::EnsEMBL::Registry')){
			$self->{_registry} = $vars{'registry'};
		}
	}
	$self->_parseIniFile();
}

sub _parseIniFile {
	my $self = shift;
	my $progPath = abs_path($0);
	my $mod = ref($self);
	$mod =~ s|::|/|g;
	$mod .= '.pm';
	my ($vol,$dirs,$file) = File::Spec->splitpath($INC{$mod});
	my @iniDirs;
	my $keep = 0;
	foreach my $d(reverse File::Spec->splitdir($dirs)){
		if($keep == 1){
			unshift(@iniDirs,$d);
		} elsif($d eq 'perl'){
			push(@iniDirs,$d,'config');
			$keep = 1;
		}
	}
	my $iniFile = File::Spec->catpath( $vol, File::Spec->catdir( @iniDirs ), 'EnsemblTranscriptSource.ini' );
	my $cfg = new Config::IniFiles( -file => $iniFile );
	my $tmpAlias;
	foreach my $key($cfg->Parameters('speciesAlias')){
		$self->{_speciesAlias}->{$key} = $cfg->val('speciesAlias',$key);
		$tmpAlias->{$cfg->val('speciesAlias',$key)} = 1;
	}

	foreach my $key (keys %$tmpAlias){
		$self->throw("$key section not found in ini file: $iniFile") unless($cfg->SectionExists($key));
		foreach my $param (SPECIES_INI_PARAMS){
			my $val = $cfg->val($key,$param);
			$self->throw("config undefined for species $key and parameter $param in ini file: $iniFile") unless defined $val;
			@{$self->{_speciesConfig}->{$key}->{$param}} = split(m/\s/,$val);
		}
	}
}

__END__

=head1 NAME

Sanger::CGP::Vagrent::TranscriptSource::EnsemblTranscriptSource - Ensembl backed transcript source

=head1 DESCRIPTION

This class will find Transcripts overlapping the specified genomic position from Ensembl using the specified Ensembl api and registry.
It will parse and filter the returned Ensembl transcripts according to the rules stored in a config file (found in perl/config/EnsemblTranscriptSource.ini).
Transcripts are filtered by status and biotype, and names are selected in priority order.
The ensembl data is that transfered into L<Transcript|Sanger::CGP::Vagrent::Data::Transcript> objects.

It inherits from L<Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource|Sanger::CGP::Vagrent::TranscriptSource::AbstractTranscriptSource>

=head1 METHODS

=head2 Constructor

=head3 new

=over

=item Usage :

 my $source = Sanger::CGP::Vagrent::TranscriptSource::EnsemblTranscriptSource->new(%params);

=item Function :

Builds a new Sanger::CGP::Vagrent::TranscriptSource::EnsemblTranscriptSource object with the supplied parameters

=item Returns :

Sanger::CGP::Vagrent::TranscriptSource::EnsemblTranscriptSource object initialized with parameter values

=item Params :

 registry => An Bio::EnsEMBL::Registry object representing the Ensembl database connection

=back

=head2 Functions

=head3 getTranscripts

=over

=item Usage :

 my @transList = $source->getTranscripts($genomicPosition);

=item Function :

Retrieves a list of L<Transcripts|Sanger::CGP::Vagrent::Data::Transcript> objects within 10Kb of the specified L<GenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition>

=item Returns :

An array of L<Sanger::CGP::Vagrent::Data::Transcript|Sanger::CGP::Vagrent::Data::Transcript> objects

=item Params :

Any L<Sanger::CGP::Vagrent::Data::AbstractGenomicPosition|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition> implementing object

=back

=head3 getTranscriptsForNextGeneInDumpRegion

=over

=item Usage :

 while (my @transList = $source->getTranscriptsForNextGeneInDumpRegion()){
 	....
 }

=item Function :

Retrieves a list of L<Transcripts|Sanger::CGP::Vagrent::Data::Transcript> objects belonging to the next gene inside the previously specified L<GenomicRegion|Sanger::CGP::Vagrent::Data::AbstractGenomicPosition>.
Requires $source->setDumpRegion($genomicPosition) to have been previously called.

=item Returns :

An array of L<Sanger::CGP::Vagrent::Data::Transcript|Sanger::CGP::Vagrent::Data::Transcript> objects

=item Params :

None

=back
