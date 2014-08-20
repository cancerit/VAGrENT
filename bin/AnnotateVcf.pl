#!/usr/bin/perl

use strict;
use English qw(-no_match_vars);
use warnings FATAL => 'all';
use Carp;
use Const::Fast qw(const);
use Try::Tiny;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

use List::Util qw(first);

use Vcf;

use Sanger::CGP::Vagrent;

# reference data storage
use Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource;

# bookmarkers
use Sanger::CGP::Vagrent::Bookmarkers::RepresentativeTranscriptBookmarker;
use Sanger::CGP::Vagrent::Bookmarkers::MostDeleteriousBookmarker;

# annotators
use Sanger::CGP::Vagrent::Annotators::AnnotatorCollection;
use Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator;
use Sanger::CGP::Vagrent::Annotators::InsertionAnnotator;
use Sanger::CGP::Vagrent::Annotators::DeletionAnnotator;
use Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator;

# input data structure
use Sanger::CGP::Vagrent::Data::Insertion;
use Sanger::CGP::Vagrent::Data::Deletion;
use Sanger::CGP::Vagrent::Data::Substitution;
use Sanger::CGP::Vagrent::Data::ComplexIndel;

const my @VCF_KEYS_EXACT => qw(SAMPLE contig vcfProcessLog reference INFO FORMAT);
const my @VCF_KEYS_STARTS_WITH => qw(source_ vcfProcessLog_);
const my @VCF_KEYS_EXCLUDE => qw(fileDate fileformat);

const my $TAB => "\t";
const my $NL => "\n";
const my $CHR_COL => 0;
const my $POS_COL => 1;
const my $REF_COL => 3;
const my $ALT_COL => 4;
const my $INFO_COL => 7;

const my $REPRE_BM => Sanger::CGP::Vagrent::Bookmarkers::RepresentativeTranscriptBookmarker->new();
const my $WORST_BM => Sanger::CGP::Vagrent::Bookmarkers::MostDeleteriousBookmarker->new();


my $header_already_parsed = 0;

eval {
  my $options = option_builder();
  Vcf::validate($options->{'input'});
  my $vcf_in = Vcf->new( file => $options->{'input'} );
  unless(defined $options->{'species'} && defined $options->{'assembly'}) {
    croak 'unable to determine species and assembly from VCF file, please specify on command line' unless find_species_in_vcf($vcf_in,$options);
  }
  open my $OUT_FH, '>', $options->{'output'} or croak 'Failed to create: '.$options->{'output'};
  my $annotator = get_annotator($options);
  
  process_data($vcf_in,$OUT_FH,$annotator,$options);
  close $OUT_FH or croak 'Failed to close: '.$options->{'output'};
  Vcf::validate($options->{'output'});
  1;
} or do {
  warn "EVAL_ERROR: $EVAL_ERROR\n" if($EVAL_ERROR);
  warn "CHILD_ERROR: $CHILD_ERROR\n" if($CHILD_ERROR);
  warn "OS_ERROR: $OS_ERROR\n" if($OS_ERROR);
  croak 'A problem occurred';
};

sub process_data {
  my ($in,$out,$anno,$opts) = @_;
  print $out generate_header($in,$opts);
  my $c = 0;
  while(my $record = $in->next_data_array) {
    $c++;
    generate_annotation($in,$anno,$opts,$record);
    print $out join $TAB, @{$record};
    print $out $NL;
  }
}

sub generate_annotation {
  my ($vcf,$annotator,$opts,$rec) = @_;
  my $var = parse_vcf_record($vcf,$opts,$rec); 
  $rec->[$INFO_COL] = $vcf->add_info_field($rec->[$INFO_COL], 'VT' => varType($var));

  my @groups = annotate($annotator,$var);
  
  if(scalar @groups > 0 && defined $groups[0]){
    my $default = getBookmarkedGroup($REPRE_BM,@groups);
    my $worst = getBookmarkedGroup($WORST_BM,@groups);
    if(defined $default){
      $rec->[$INFO_COL] = $vcf->add_info_field($rec->[$INFO_COL], 'VD' => stringifyAnnotation($default));
      $rec->[$INFO_COL] = $vcf->add_info_field($rec->[$INFO_COL], 'VC' => $annotator->getOntologySummary($default));
    } 
    if(defined $worst){
      $rec->[$INFO_COL] = $vcf->add_info_field($rec->[$INFO_COL], 'VW' => stringifyAnnotation($worst));
      if(!defined $default){
        $rec->[$INFO_COL] = $vcf->add_info_field($rec->[$INFO_COL], 'VC' => $annotator->getOntologySummary($worst));
      }
    }
  } 

#   my @annotation = annotate($annotator,$var);
#   if(defined $annotation[0] && $annotation[0] ne q{}){
#     $rec->[$INFO_COL] = $vcf->add_info_field($rec->[$INFO_COL], 'VD' => $annotation[0]);
#   	my $terms = (split('\|',$annotation[0]))[5];
# 		$rec->[$INFO_COL] = $vcf->add_info_field($rec->[$INFO_COL], 'VC' => $terms);
#   }  
#   if(defined $annotation[1] && $annotation[1] ne q{}){
#     $rec->[$INFO_COL] = $vcf->add_info_field($rec->[$INFO_COL], 'VW' => $annotation[1]);
#   }  


  return;
}

sub annotate {
  my ($annotator,$var) = @_;
  my @annotationGroups;
  try {
    @annotationGroups = $annotator->getAnnotation($var);
  } catch {
      warn "caught error: $_\n"; # not $@
  };  
#   my $repAnno = q{};
#   my $worstAnno = q{};
#   my $representative = undef;
# 	my $worst = undef;
#   if(scalar(@annotationGroups) > 0 && defined($annotationGroups[0])){	
#     my $sameAnno = 0;
#     foreach my $a (@annotationGroups){
#   		my $wasRep = 0;
#   		my $wasWor = 0;
#   		if(!defined($representative) && $a->hasBookmark($REPRE_BM)){
#   			$wasRep = 1;
#   			$representative = $a;
#   		}
#   		if(!defined($worst) && $a->hasBookmark($WORST_BM)){
#   			$wasWor = 1;
#   			$worst = $a;
#   		}
#       if($wasRep && $wasWor){
# 			  $sameAnno = 1;
# 		  }
# 		  last if(defined($representative) && defined($worst));
#     }  
#   }
#   return (stringifyAnnotation($representative),stringifyAnnotation($worst));
  return @annotationGroups;
}

sub getBookmarkedGroup {
  my ($bm, @groups) = @_;
  my $out = undef;
  foreach my $g(@groups){
    if($g->hasBookmark($bm)){
      $out = $g;
    }
  }
  return $out;
}

sub stringifyAnnotation {
  my $anno = shift;
  return undef unless defined $anno;
  my $mrna = $anno->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext);
	my $cds = $anno->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext);
	my $prot = $anno->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext);

	my $desc;
	$desc = $anno->getLabel.'|';
	if(defined($anno->getCCDS) && $anno->getCCDS ne ''){
		$desc .= $anno->getCCDS.'|';
	} else {
		$desc .= $anno->getAccession.'|';
	}

	$desc .= $mrna->getDescription.'|';
	if(defined($cds)){
		$desc .= $cds->getDescription.'|';
	} else {
		$desc .= '-|';
	}
	if(defined($prot)){
		$desc .= $prot->getDescription.'|';
	} else {
		$desc .= '-|';
	}

 	my @classIds;
 	my @classText;

 	foreach my $a($anno,$mrna,$cds,$prot){
 		next unless defined $a;

 		foreach my $term ($a->getClassifications){
 			my ($id,$text) = (split(':',$term))[1,2];
 			if((scalar @classIds) == 0){
 				push @classIds, $id;
 				push @classText, $text;
 			}
 			else {
 				next if(first {$_ eq $id} @classIds);
 				push @classIds, $id;
 				push @classText, $text;
 			}
 		}
 	}

	$desc .= join(':',@classText).'|';
	$desc .= 'SO:'.join(':SO:',@classIds);

	return $desc;
}

sub parse_vcf_record {
  my ($vcf,$opts,$record) = @_;
  my $var;
  my ($ref, $alt) = @{$record}[$REF_COL,$ALT_COL];
  my $ref_len = length $ref;
  my $alt_len = length $alt;

  my $pos = $record->[$POS_COL];

  # test in prevalence order ()based on genome example as no reason why all types cannot be in a single vcf file

  # a straight deletion
  if($ref_len > 1 && $alt_len == 1) {
    my $min = $pos + 1; # as vcf records base before for non-insert
    my $max = $min + ($vcf->get_info_field($record->[$INFO_COL], 'LEN') - 1);
    $var = Sanger::CGP::Vagrent::Data::Deletion->new(
							'species'				=> $opts->{'species'},
							'genomeVersion' => $opts->{'assembly'},
							'chr'	          => $record->[$CHR_COL],
							'minpos'        => $min,
							'maxpos'        => $max,
							'delseq'        => substr($ref, 1)); # as ref always includes the reference base prior to the change
  }
  # a straight insertion
  elsif($ref_len == 1 && $alt_len > 1) {
    $var = Sanger::CGP::Vagrent::Data::Insertion->new(
							'species'				=> $opts->{'species'},
							'genomeVersion' => $opts->{'assembly'},
							'chr'	          => $record->[$CHR_COL],
							'minpos'        => $pos,
							'maxpos'        => $pos+1,
							'insseq'        => substr($alt, 1)); # as alt always includes the reference base prior to the change
  }
  # this is a sub
  elsif($ref_len == 1 && $alt_len == 1) {
    $var = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> $opts->{'species'},
			'genomeVersion' => $opts->{'assembly'},
			'chr' 					=> $record->[$CHR_COL],
			'minpos'				=> $pos,
			'maxpos'				=> $pos,
			'wt' 						=> $ref,
			'mt'						=> $alt,);
  }
  # complex indel
  elsif($ref_len > 1 && $alt_len > 1) {
    my $min = $pos + 1; # as vcf records base before for non-insert
    my $max = $min + ($vcf->get_info_field($record->[$INFO_COL], 'LEN') - 1);
    $var = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
							'species'				=> $opts->{'species'},
							'genomeVersion' => $opts->{'assembly'},
							'chr'	          => $record->[$CHR_COL],
							'minpos'        => $min,
							'maxpos'        => $max,
							'delseq'        => substr($ref, 1),  # as ref always includes the reference base prior to the change
							'insseq'        => substr($alt, 1)); # as alt always includes the reference base prior to the change
  }
  else {
    croak "Unable to interpret this vcf result as Sub, Ins, Del or ComplexInDel";
  }
  return $var;
}

sub generate_header {
  my ($in,$opts) = @_;
  $in->parse_header unless $header_already_parsed;
  $in->add_header_line({'key'=>'source', 'value' => 'AnnotateVcf.pl'}, 'append' => 1);
  $in->add_header_line(make_process_log($opts), 'append' => 1);
  $in->add_header_line({key => 'cgpAnalysisProc',value => $opts->{'p'}}, 'append' => 1 ) if defined $opts->{'process'};
  foreach my $info (make_info_fields()){
    $in->add_header_line($info);
  }
  return sort_header($in);
}

sub sort_header {
  my @lines = split /\n/, shift->format_header;
  my $format_head = shift @lines;
  my $col_head = pop @lines;
  my %header_sets;
  for my $h_line(@lines) {
    my ($h_type) = $h_line =~ m/^\#\#([^=]+)/;
    unless(exists $header_sets{$h_type}) {
      $header_sets{$h_type} = [];
    }
    push @{$header_sets{$h_type}}, $h_line;
  }
  my @final_header;
  push @final_header, $format_head;
  for my $h_type(sort keys %header_sets) {
    push @final_header, @{$header_sets{$h_type}};
  }
  push @final_header, $col_head;
  return (join "\n", @final_header)."\n";
}

sub make_info_fields {
  my @out;

  push @out, {key => 'INFO',ID => 'VD',Number => 1,Type => 'String',Description => 'Vagrent Default Annotation'},
              {key => 'INFO',ID => 'VW',Number => 1,Type => 'String',Description => 'Vagrent Most Deleterious Annotation'},
              {key => 'INFO',ID => 'VT',Number => 1,Type => 'String',Description => 'Variant type based on the Vagrent Default Annotation'},
              {key => 'INFO',ID => 'VC',Number => 1,Type => 'String',Description => 'Variant consequence based on the Vagrent Default Annotation'};

  return @out;
}

sub make_process_log {
  my ($opts) = @_;
  my $params;
  foreach my $key (keys %$opts){
    $params->{$key} = $opts->{$key} if(defined $opts->{$key});
  }
  return {'key'=>'vcfProcessLog',
			InputVCF => $opts->{'i'},
			InputVCFSource => 'AnnotateVcf.pl',
			InputVCFVer => Sanger::CGP::Vagrent->VERSION,
			InputVCFParam => $params,
		};
}

sub get_annotator {
	my $options = shift;

  # creating an EnsemblTranscriptSource using the Ensembl registry
	my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => $options->{'cache'});

	# creating an AnnotatorCollection
	my $annotator = Sanger::CGP::Vagrent::Annotators::AnnotatorCollection->new(
							bookmarker => [$REPRE_BM,$WORST_BM],
							only_bookmarked => 1);

	# creating and adding a SimpleSubstitutionAnnotator to the AnnotatorCollection
	$annotator->addAnnotator(Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts, ontologySymmary => $annotator->getSummaryCache));
	$annotator->addAnnotator(Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts, ontologySymmary => $annotator->getSummaryCache));
	$annotator->addAnnotator(Sanger::CGP::Vagrent::Annotators::DeletionAnnotator->new(transcriptSource => $ts, ontologySymmary => $annotator->getSummaryCache));
	$annotator->addAnnotator(Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts, ontologySymmary => $annotator->getSummaryCache));

	return $annotator;
}

sub find_species_in_vcf {
  my ($vcf,$opts) = @_;
  my $out = 0;  
  $vcf->parse_header;
  $header_already_parsed = 1;
  
  my ($ctgs) = @{$vcf->get_header_line(key=>'contig')};
  if(defined $ctgs && ref $ctgs eq 'HASH'){
    my $ctg = $ctgs->{(keys %$ctgs)[0]};
    if(defined $ctg){
      if(exists $ctg->{'species'} && defined $ctg->{'species'} && exists $ctg->{'assembly'} && defined $ctg->{'assembly'} ){
        $opts->{'species'} = $ctg->{'species'};
        $opts->{'assembly'} = $ctg->{'assembly'};
        $out = 1;
      }
    }
  }
  return $out;
}

sub varType {
	my $var = shift;
	if($var->isa('Sanger::CGP::Vagrent::Data::Substitution')){
		return 'Sub';
	} elsif($var->isa('Sanger::CGP::Vagrent::Data::Deletion')){
	 	return 'Del';
	} elsif($var->isa('Sanger::CGP::Vagrent::Data::Insertion')){
	  return 'Ins';
	} elsif($var->isa('Sanger::CGP::Vagrent::Data::ComplexIndel')){
	  return 'Complex';
	} else {
		return 'Unknown'
	}
}

sub option_builder {
  my ($factory) = @_;

  my %opts = ();

  my $result = &GetOptions (
    'h|help' => \$opts{'help'},
    'v|version' => \$opts{'version'},
    'i|input=s' => \$opts{'input'},
    'o|output=s' => \$opts{'output'},
    'c|cache=s' => \$opts{'cache'},
    't|tabix=s' => \$opts{'tabix'},
    'p|process=n' => \$opts{'process'},
    'sp|species=s' => \$opts{'species'},
    'as|assembly=s' => \$opts{'assembly'},
    
  );

  pod2usage() if($opts{'help'});
  
  if($opts{'version'}){
    print 'Version: '.Sanger::CGP::Vagrent->VERSION."\n";
    exit;
  }  

  pod2usage(q{'-i' must be defined}) unless($opts{'input'});
  pod2usage(q{'-i' must exist}) unless(-e $opts{'input'});
  pod2usage(q{'-i' must be a file}) unless(-f $opts{'input'});
  pod2usage(q{'-i' is an empty file}) unless(-s $opts{'input'});

  pod2usage(q{'-c' must be defined}) unless($opts{'cache'});
  pod2usage(q{'-c' must exist}) unless(-e $opts{'cache'});
	pod2usage(q{'-c' must be a file}) unless(-f $opts{'cache'});
	pod2usage(q{'-c' is an empty file}) unless(-s $opts{'cache'});

	pod2usage(q{'-o' must be defined}) unless($opts{'output'});
	
  if($opts{'tabix'}) {
    pod2usage(q{'-t' must exist}) unless(-e $opts{'tabix'});
    pod2usage(q{'-t' must be a file}) unless(-f $opts{'tabix'});
    pod2usage(q{'-t' is an empty file}) unless(-s $opts{'tabix'});
  }
  return \%opts;
}

__END__

=head1 NAME

AnnotateVcf.pl - Annotate variants - Sub/Snp, Insertion, Deletion, ComplexInDel

=head1 SYNOPSIS

AnnotateVcf.pl [-h] -i <IN_FILE> -o <OUT_FILE>

  General Options:

    --help      (-h)      Brief documentation

    --input     (-i)      Input vcf file (expects *.bgz)

    --output    (-o)      Output vcf

    --cache     (-c)      Reference data cache file

  Conditional (can be specified if missing from the input VCF file)

    --species   (-sp)     Species

    --assembly  (-as)     Genome assembly version

  Optional

    --process   (-p)      ID_ANALYSIS_PROCESS that generated this file

    --tabix     (-t)      Gene footprint file, expects co-located *.tbi

=cut
