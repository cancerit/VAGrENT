#!/usr/bin/perl

use strict;

use lib qw(lib ../lib t ../t);

use warnings FATAL => 'all';
use Test::More;

use Data::Dumper;

use Sanger::CGP::Vagrent::Data::Substitution;
use Sanger::CGP::Vagrent::Data::Transcript;
use Sanger::CGP::Vagrent::Data::Exon;
use Sanger::CGP::Vagrent::Data::Annotation;
use Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator;

use Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource;

use AnnotationTestUtils;

testIntronic();
testSplice();
testExonic();
testUpStreamDownStream();

done_testing();

sub testUpStreamDownStream {
	testUpsteamMilesAway_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testUpsteam5001bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testUpsteam5000bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testUpsteam2001bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testUpsteam2000bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testUpsteam1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	testDownstream1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testDownstream500bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testDownstream501bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testDownstream5000bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testDownstream5001bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testDownstreamMilesAwaybp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	testUpsteamMilesAway_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testUpsteam5001bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testUpsteam5000bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testUpsteam2001bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testUpsteam2000bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testUpsteam1bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

	testDownsteam1bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testDownsteam500bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testDownsteam501bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testDownsteam5000bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testDownsteam5001bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testDownsteamMilesAway_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
}
sub testIntronic {


#		Testing 5 PRIME UTR INTRONIC
		test5PrimeUtrSpliceRegion7_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
 		test5PrimeUTRIntronic1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
 		test5PrimeUTRIntronic2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
 		test5PrimeUTRIntronic3_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
 		test5PrimeUTRIntronic4_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
 		test5PrimeUtrSpliceRegion8_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
 
 		test5PrimeUtrSpliceRegion7_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
 		test5PrimeUTRIntronic1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
 		test5PrimeUTRIntronic2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
 		test5PrimeUTRIntronic3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
 		test5PrimeUTRIntronic4_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
 		test5PrimeUtrSpliceRegion8_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);


#		Testing CDS INTRONIC
 		testIntronic1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
 		testIntronic2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
 		testSpliceRegion1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
 		testSpliceRegion15_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
 		testIntronic3_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
 		testIntronic4_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
 
 		testIntronic1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
 		testIntronic2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
 		testSpliceRegion1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
 		testSpliceRegion15_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
 		testIntronic3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
 		testIntronic4_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
 
 		testIntronic1_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
 		testIntronic2_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
 		testSpliceRegion1_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
 		testSpliceRegion15_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
 		testIntronic3_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
 		testIntronic4_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
 
 		testIntronic1_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
 		testIntronic2_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
 		testSpliceRegion1_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
 		testSpliceRegion15_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
 		testIntronic3_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
 		testIntronic4_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);


#		Testing 3 PRIME UTR INTRONIC
 		test3PrimeUtrSpliceRegion98_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
 		test3PrimeUTRIntronic96_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
 		test3PrimeUTRIntronic97_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
 		test3PrimeUTRIntronic98_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
 		test3PrimeUTRIntronic99_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
 		test3PrimeUtrSpliceRegion99_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
 
 		test3PrimeUtrSpliceRegion98_FAM54A(AnnotationTestUtils::FAM54A_TRANSCRIPT);
 		test3PrimeUTRIntronic96_FAM54A(AnnotationTestUtils::FAM54A_TRANSCRIPT);
 		test3PrimeUTRIntronic97_FAM54A(AnnotationTestUtils::FAM54A_TRANSCRIPT);
 		test3PrimeUTRIntronic98_FAM54A(AnnotationTestUtils::FAM54A_TRANSCRIPT);
 		test3PrimeUTRIntronic99_FAM54A(AnnotationTestUtils::FAM54A_TRANSCRIPT);
 		test3PrimeUtrSpliceRegion99_FAM54A(AnnotationTestUtils::FAM54A_TRANSCRIPT);

}
sub testSplice {

#		PROTIEN CODING + STRAND SPLICE
		# detailed + strand tests for all splice positions
		testIntronic2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion3_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion4_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion5_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion6_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion7_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion8_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testEssentialSplice1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testEssentialSplice2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testEssentialSplice3_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testEssentialSplice4_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion9_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion10_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testEssentialSplice5_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion11_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion12_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion13_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion14_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion15_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testIntronic3_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

#		PROTIEN CODING - STRAND SPLICE
		# detailed - strand tests for all splice positions
		testIntronic2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion4_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion5_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion6_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion7_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion8_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testEssentialSplice1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testEssentialSplice2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testEssentialSplice3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testEssentialSplice4_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion9_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion10_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testEssentialSplice5_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion11_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion12_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion13_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion14_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion15_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testIntronic3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

# handful of tests for non-coding genes
#		NON PROTIEN CODING + STRAND SPLICE
		testIntronic2_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
		testSpliceRegion1_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
		testEssentialSplice2_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
		testEssentialSplice5_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
		testSpliceRegion15_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
		testIntronic3_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);

#		NON PROTIEN CODING - STRAND SPLICE
		testIntronic2_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
		testSpliceRegion1_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
		testEssentialSplice1_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
		testEssentialSplice4_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
		testSpliceRegion15_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
		testIntronic3_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
}
sub testExonic {
	subtest 'Testing Exonic ' => sub {
		subtest 'PROTIEN CODING SILENT' => sub {
			testSilent1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
			testSilent1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
			done_testing();
		};
		subtest 'PROTIEN CODING MISSENSE' => sub {
			testMissense1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
			testMissense1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
			done_testing();
		};
		subtest 'PROTIEN CODING NONSENSE' => sub {
			testNonsense1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
			testNonsense1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
			done_testing();
		};
		subtest 'PROTIEN CODING STOP LOST' => sub {
			testStopLost1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
			testStopLost1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
			done_testing();
		};
		subtest 'PROTIEN CODING START LOST' => sub {
			testStartLost1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
			testStartLost1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
			done_testing();
		};
		subtest 'PROTIEN CODING 5 PRIME UTR' => sub {
			test5PrimeUTR1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
			test5PrimeUtr1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
			done_testing();
		};
		subtest 'PROTIEN CODING PREMATURE START GAINED' => sub {
			testPrematureStartGained_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
			testPrematureStartGained_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
			done_testing();
		};
		subtest 'PROTIEN CODING 3 PRIME UTR' => sub {
			test3PrimeUtr1_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
			test3PrimeUtr1_FAM54A(AnnotationTestUtils::FAM54A_TRANSCRIPT);
			done_testing();
		};
		subtest 'NON-CODING' => sub {
			testExon1_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
			testExon1_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
			done_testing();
		};
		done_testing();
	};
}


# CEP350 protein coding gene with 5 prime utr exons on + strand of genome

sub testUpsteamMilesAway_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Upstream Miles Away + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179913873,
			'maxpos'				=> 179913873,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testUpsteam5001bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Upstream 5001bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179918872,
			'maxpos'				=> 179918872,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testUpsteam5000bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Upstream 5000bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179918873,
			'maxpos'				=> 179918873,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get5KBUpStreamVariantClass);

		done_testing();
	};
}
sub testUpsteam2001bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Upstream 2001 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179921872,
			'maxpos'				=> 179921872,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get5KBUpStreamVariantClass);

		done_testing();
	};
}
sub testUpsteam2000bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Upstream 2000 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179921873,
			'maxpos'				=> 179921873,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get2KBUpStreamVariantClass);

		done_testing();
	};
}
sub testUpsteam1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Upstream 1 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179923872,
			'maxpos'				=> 179923872,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get2KBUpStreamVariantClass);

		done_testing();
	};
}
sub testDownstream1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Downstream 1 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084016,
			'maxpos'				=> 180084016,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get500BPDownStreamVariantClass);

		done_testing();
	};
}
sub testDownstream500bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Downstream 500 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084515,
			'maxpos'				=> 180084515,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get500BPDownStreamVariantClass);

		done_testing();
	};
}
sub testDownstream501bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Downstream 501 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084516,
			'maxpos'				=> 180084516,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get5KBDownStreamVariantClass);
		done_testing();
	};
}
sub testDownstream5000bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Downstream 5000 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180089015,
			'maxpos'				=> 180089015,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get5KBDownStreamVariantClass);
		done_testing();
	};
}
sub testDownstream5001bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Downstream 5001 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180089016,
			'maxpos'				=> 180089016,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');
		done_testing();
	};
}
sub testDownstreamMilesAwaybp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Downstream Miles Away + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180099016,
			'maxpos'				=> 180099016,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');
		done_testing();
	};
}

sub test5PrimeUTR1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 prime UTR + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179923875,
			'maxpos'				=> 179923875,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getExonClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									3,0,3,0,'C','U','r.3c>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get5PrimeUtrVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','c.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testPrematureStartGained_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Premature Start Gained + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955310,
			'maxpos'				=> 179955310,
			'wt' 						=> 'T',
			'mt'						=> 'G',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getExonClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									412,0,412,0,'U','G','r.412u>g',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getPrematureStartGainedVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','c.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub test5PrimeUtrSpliceRegion7_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 prime UTR Splice Region + strand 1' => sub {
    my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   

    print Dumper($ts);
  
		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924287,
			'maxpos'				=> 179924287,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									405,10,405,10,'A','U','r.405+10a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass,$a->get5PrimeUtrVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','c.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Intronic + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924288,
			'maxpos'				=> 179924288,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Intronic + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924388,
			'maxpos'				=> 179924388,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic3_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Intronic + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955204,
			'maxpos'				=> 179955204,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic4_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Intronic + strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955293,
			'maxpos'				=> 179955293,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUtrSpliceRegion8_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 prime UTR Splice Region + strand 8' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955294,
			'maxpos'				=> 179955294,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									406,-10,406,-10,'A','U','r.406-10a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass,$a->get5PrimeUtrVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','c.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testStartLost1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Start Lost + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955317,
			'maxpos'				=> 179955317,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									419,0,419,0,'A','U','r.419a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									1,0,1,0,'A','T','c.1A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									1,0,1,0,'M','L','p.M1L',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getStartLostVariantClass);
		done_testing();
	};
}
sub testNonsense1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Nonsense + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179959753,
			'maxpos'				=> 179959753,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									650,0,650,0,'A','U','r.650a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									232,0,232,0,'A','T','c.232A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									78,0,78,0,'K','*','p.K78*',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getStopGainedVariantClass);
		done_testing();
	};
}
sub testIntronic1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Intronic + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961097,
			'maxpos'				=> 179961097,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Inronic + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961186,
			'maxpos'				=> 179961186,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testSpliceRegion1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961187,
			'maxpos'				=> 179961187,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-10,654,-10,'A','U','r.654-10a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-10,236,-10,'A','T','c.236-10A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961188,
			'maxpos'				=> 179961188,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-9,654,-9,'A','U','r.654-9a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-9,236,-9,'A','T','c.236-9A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion3_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961189,
			'maxpos'				=> 179961189,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-8,654,-8,'A','U','r.654-8a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-8,236,-8,'A','T','c.236-8A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion4_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region + strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961190,
			'maxpos'				=> 179961190,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-7,654,-7,'A','U','r.654-7a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-7,236,-7,'A','T','c.236-7A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion5_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region + strand 5' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961191,
			'maxpos'				=> 179961191,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-6,654,-6,'A','U','r.654-6a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-6,236,-6,'A','T','c.236-6A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion6_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region + strand 6' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961192,
			'maxpos'				=> 179961192,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-5,654,-5,'A','U','r.654-5a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-5,236,-5,'A','T','c.236-5A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion7_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region + strand 7' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961193,
			'maxpos'				=> 179961193,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-4,654,-4,'A','U','r.654-4a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-4,236,-4,'A','T','c.236-4A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion8_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region + strand 8' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961194,
			'maxpos'				=> 179961194,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-3,654,-3,'A','U','r.654-3a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-3,236,-3,'A','T','c.236-3A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testEssentialSplice1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961195,
			'maxpos'				=> 179961195,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-2,654,-2,'A','U','r.654-2a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-2,236,-2,'A','T','c.236-2A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testEssentialSplice2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961196,
			'maxpos'				=> 179961196,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-1,654,-1,'A','U','r.654-1a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-1,236,-1,'A','T','c.236-1A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSilent1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Silent + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961198,
			'maxpos'				=> 179961198,
			'wt' 						=> 'T',
			'mt'						=> 'C',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									655,0,655,0,'U','C','r.655u>c',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									237,0,237,0,'T','C','c.237T>C',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									79,0,79,0,'D','D','p.D79D',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSynonymousVariantClass);
		done_testing();
	};
}
sub testMissense1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Missense + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961199,
			'maxpos'				=> 179961199,
			'wt' 						=> 'G',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									656,0,656,0,'G','U','r.656g>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									238,0,238,0,'G','T','c.238G>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									80,0,80,0,'G','C','p.G80C',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getNonSynonymousVariantClass);
		done_testing();
	};
}
sub testEssentialSplice3_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961357,
			'maxpos'				=> 179961357,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getEssentialSpliceSiteClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,1,813,1,'A','U','r.813+1a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,1,395,1,'A','T','c.395+1A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testEssentialSplice4_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice + strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961358,
			'maxpos'				=> 179961358,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getEssentialSpliceSiteClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,2,813,2,'A','U','r.813+2a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,2,395,2,'A','T','c.395+2A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion9_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region + strand 9' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961359,
			'maxpos'				=> 179961359,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,3,813,3,'A','U','r.813+3a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,3,395,3,'A','T','c.395+3A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion10_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region + strand 10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961360,
			'maxpos'				=> 179961360,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,4,813,4,'A','U','r.813+4a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,4,395,4,'A','T','c.395+4A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testEssentialSplice5_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice + strand 5' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961361,
			'maxpos'				=> 179961361,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getEssentialSpliceSiteClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,5,813,5,'A','U','r.813+5a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,5,395,5,'A','T','c.395+5A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion11_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region + strand 11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961362,
			'maxpos'				=> 179961362,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,6,813,6,'A','U','r.813+6a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,6,395,6,'A','T','c.395+6A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion12_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region + strand 12' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961363,
			'maxpos'				=> 179961363,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,7,813,7,'A','U','r.813+7a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,7,395,7,'A','T','c.395+7A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion13_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region + strand 13' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961364,
			'maxpos'				=> 179961364,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,8,813,8,'A','U','r.813+8a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,8,395,8,'A','T','c.395+8A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion14_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region + strand 14' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961365,
			'maxpos'				=> 179961365,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,9,813,9,'A','U','r.813+9a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,9,395,9,'A','T','c.395+9A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion15_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region + strand 15' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961366,
			'maxpos'				=> 179961366,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,10,813,10,'A','U','r.813+10a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,10,395,10,'A','T','c.395+10A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testIntronic3_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Intronic + strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961367,
			'maxpos'				=> 179961367,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic4_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Intronic + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961456,
			'maxpos'				=> 179961456,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testStopLost1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 StopLost + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180080294,
			'maxpos'				=> 180080294,
			'wt' 						=> 'T',
			'mt'						=> 'G',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getExonClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									9770,0,9770,0,'U','G','r.9770u>g',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									9352,0,9352,0,'T','G','c.9352T>G',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									3118,0,3118,0,'*','G','p.*3118G',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getStopLostVariantClass);
		done_testing();
	};
}

# TOR1AIP2 protein coding gene with 5 prime utr exons on - strand of genome

sub testUpsteamMilesAway_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Upstream Miles Away - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179856934,
			'maxpos'				=> 179856934,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');
		done_testing();
	};
}
sub testUpsteam5001bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Upstream 5001bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179851935,
			'maxpos'				=> 179851935,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');
		done_testing();
	};
}
sub testUpsteam5000bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Upstream 5000bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179851934,
			'maxpos'				=> 179851934,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get5KBUpStreamVariantClass);
		done_testing();
	};
}
sub testUpsteam2001bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Upstream 2001bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179848935,
			'maxpos'				=> 179848935,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get5KBUpStreamVariantClass);
		done_testing();
	};
}
sub testUpsteam2000bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Upstream 2000bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179848934,
			'maxpos'				=> 179848934,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get2KBUpStreamVariantClass);
		done_testing();
	};
}
sub testUpsteam1bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Upstream 1bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846935,
			'maxpos'				=> 179846935,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get2KBUpStreamVariantClass);
		done_testing();
	};
}
sub testDownsteam1bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Downstream 1bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179809101,
			'maxpos'				=> 179809101,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get500BPDownStreamVariantClass);
		done_testing();
	};
}
sub testDownsteam500bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Downstream 500bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179808602,
			'maxpos'				=> 179808602,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get500BPDownStreamVariantClass);
		done_testing();
	};
}
sub testDownsteam501bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Downstream 501bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179808601,
			'maxpos'				=> 179808601,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get5KBDownStreamVariantClass);
		done_testing();
	};
}
sub testDownsteam5000bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Downstream 5000bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179804102,
			'maxpos'				=> 179804102,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get5KBDownStreamVariantClass);
		done_testing();
	};
}
sub testDownsteam5001bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Downstream 5001bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179804101,
			'maxpos'				=> 179804101,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');
		done_testing();
	};
}
sub testDownsteamMilesAway_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Downstream Miles Away - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179704101,
			'maxpos'				=> 179704101,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');
		done_testing();
	};
}

sub test5PrimeUtrSpliceRegion7_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 prime UTR Splice Region - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846364,
			'maxpos'				=> 179846364,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->get5PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									242,10,242,10,'U','A','r.242+10u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass,$a->get5PrimeUtrVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','c.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Intronic - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846363,
			'maxpos'				=> 179846363,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->get5PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Intronic - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846274,
			'maxpos'				=> 179846274,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->get5PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic3_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Intronic - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179822046,
			'maxpos'				=> 179822046,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->get5PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic4_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Intronic - strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821957,
			'maxpos'				=> 179821957,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->get5PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUtrSpliceRegion8_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 prime UTR Splice Region - strand 8' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821956,
			'maxpos'				=> 179821956,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->get5PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									243,-10,243,-10,'U','A','r.243-10u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass,$a->get5PrimeUtrVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','c.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testIntronic1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Intronic - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820598,
			'maxpos'				=> 179820598,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Intronic - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820509,
			'maxpos'				=> 179820509,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testSpliceRegion1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820508,
			'maxpos'				=> 179820508,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-10,423,-10,'U','A','r.423-10u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-10,35,-10,'T','A','c.35-10T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820507,
			'maxpos'				=> 179820507,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-9,423,-9,'U','A','r.423-9u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-9,35,-9,'T','A','c.35-9T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion3_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820506,
			'maxpos'				=> 179820506,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-8,423,-8,'U','A','r.423-8u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-8,35,-8,'T','A','c.35-8T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion4_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region - strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820505,
			'maxpos'				=> 179820505,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-7,423,-7,'U','A','r.423-7u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-7,35,-7,'T','A','c.35-7T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion5_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region - strand 5' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820504,
			'maxpos'				=> 179820504,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-6,423,-6,'U','A','r.423-6u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-6,35,-6,'T','A','c.35-6T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion6_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region - strand 6' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820503,
			'maxpos'				=> 179820503,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-5,423,-5,'U','A','r.423-5u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-5,35,-5,'T','A','c.35-5T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion7_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region - strand 7' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820502,
			'maxpos'				=> 179820502,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-4,423,-4,'U','A','r.423-4u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-4,35,-4,'T','A','c.35-4T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion8_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region - strand 8' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820501,
			'maxpos'				=> 179820501,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-3,423,-3,'U','A','r.423-3u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-3,35,-3,'T','A','c.35-3T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testEssentialSplice1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820500,
			'maxpos'				=> 179820500,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getEssentialSpliceSiteClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-2,423,-2,'U','A','r.423-2u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-2,35,-2,'T','A','c.35-2T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testEssentialSplice2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820499,
			'maxpos'				=> 179820499,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getEssentialSpliceSiteClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-1,423,-1,'U','A','r.423-1u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-1,35,-1,'T','A','c.35-1T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testEssentialSplice3_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819979,
			'maxpos'				=> 179819979,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getEssentialSpliceSiteClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,1,941,1,'U','A','r.941+1u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,1,553,1,'T','A','c.553+1T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testEssentialSplice4_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice - strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819978,
			'maxpos'				=> 179819978,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getEssentialSpliceSiteClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,2,941,2,'U','A','r.941+2u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,2,553,2,'T','A','c.553+2T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion9_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region - strand 9' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819977,
			'maxpos'				=> 179819977,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,3,941,3,'U','A','r.941+3u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,3,553,3,'T','A','c.553+3T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion10_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region - strand 10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819976,
			'maxpos'				=> 179819976,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,4,941,4,'U','A','r.941+4u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,4,553,4,'T','A','c.553+4T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testEssentialSplice5_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice - strand 5' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819975,
			'maxpos'				=> 179819975,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getEssentialSpliceSiteClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,5,941,5,'U','A','r.941+5u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,5,553,5,'T','A','c.553+5T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion11_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region - strand 11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819974,
			'maxpos'				=> 179819974,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,6,941,6,'U','A','r.941+6u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,6,553,6,'T','A','c.553+6T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion12_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region - strand 12' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819973,
			'maxpos'				=> 179819973,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,7,941,7,'U','A','r.941+7u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,7,553,7,'T','A','c.553+7T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion13_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region - strand 13' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819972,
			'maxpos'				=> 179819972,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,8,941,8,'U','A','r.941+8u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,8,553,8,'T','A','c.553+8T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion14_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region - strand 14' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819971,
			'maxpos'				=> 179819971,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,9,941,9,'U','A','r.941+9u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,9,553,9,'T','A','c.553+9T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testSpliceRegion15_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region - strand 15' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819970,
			'maxpos'				=> 179819970,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,10,941,10,'U','A','r.941+10u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,10,553,10,'T','A','c.553+10T>A',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testIntronic3_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Intronic - strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819969,
			'maxpos'				=> 179819969,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic4_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Intronic - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819880,
			'maxpos'				=> 179819880,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testPrematureStartGained_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Premature Start Gained - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846919,
			'maxpos'				=> 179846919,
			'wt' 						=> 'G',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getExonClass,$a->get5PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									16,0,16,0,'C','A','r.16c>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getPrematureStartGainedVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','c.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub test5PrimeUtr1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 prime UTR - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821801,
			'maxpos'				=> 179821801,
			'wt' 						=> 'G',
			'mt'						=> 'C',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getExonClass,$a->get5PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									388,0,388,0,'C','G','r.388c>g',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get5PrimeUtrVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','c.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub testStartLost1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Start Lost - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821800,
			'maxpos'				=> 179821800,
			'wt' 						=> 'T',
			'mt'						=> 'G',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getExonClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									389,0,389,0,'A','C','r.389a>c',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									1,0,1,0,'A','C','c.1A>C',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									1,0,1,0,'M','L','p.M1L',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getStartLostVariantClass);
		done_testing();
	};
}
sub testSilent1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Silent - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820497,
			'maxpos'				=> 179820497,
			'wt' 						=> 'G',
			'mt'						=> 'A',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getExonClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									424,0,424,0,'C','U','r.424c>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									36,0,36,0,'C','T','c.36C>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									12,0,12,0,'D','D','p.D12D',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSynonymousVariantClass);
		done_testing();
	};
}
sub testMissense1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Missense - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820496,
			'maxpos'				=> 179820496,
			'wt' 						=> 'A',
			'mt'						=> 'G',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getExonClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									425,0,425,0,'U','C','r.425u>c',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									37,0,37,0,'T','C','c.37T>C',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									13,0,13,0,'S','P','p.S13P',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getNonSynonymousVariantClass);
		done_testing();
	};
}
sub testNonsense1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Nonsense - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820382,
			'maxpos'				=> 179820382,
			'wt' 						=> 'T',
			'mt'						=> 'A',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getExonClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									539,0,539,0,'A','U','r.539a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									151,0,151,0,'A','T','c.151A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									51,0,51,0,'K','*','p.K51*',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getStopGainedVariantClass);
		done_testing();
	};
}
sub testStopLost1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Stop Lost - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179815206,
			'maxpos'				=> 179815206,
			'wt' 						=> 'T',
			'mt'						=> 'A',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getExonClass,$a->getCDSClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									1801,0,1801,0,'A','U','r.1801a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									1413,0,1413,0,'A','T','c.1413A>T',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getCodonVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									471,0,471,0,'*','Y','p.*471Y',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getStopLostVariantClass);
		done_testing();
	};
}

# MYB protein coding gene with 3 prime utr exons on + strand of genome

sub test3PrimeUtrSpliceRegion98_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Splice Region + strand 98' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135524472,
			'maxpos'				=> 135524472,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->get3PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1933,10,1933,10,'A','U','r.1933+10a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass,$a->get3PrimeUtrVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','c.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic96_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Intronic + strand 96' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135524473,
			'maxpos'				=> 135524473,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->get3PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic97_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Intronic + strand 97' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135524562,
			'maxpos'				=> 135524562,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->get3PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic98_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Intronic + strand 98' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135538902,
			'maxpos'				=> 135538902,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->get3PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic99_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Intronic + strand 99' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135538991,
			'maxpos'				=> 135538991,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->get3PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUtrSpliceRegion99_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Splice Region + strand 99' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135538992,
			'maxpos'				=> 135538992,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->get3PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1934,-10,1934,-10,'A','U','r.1934-10a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass,$a->get3PrimeUtrVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','c.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub test3PrimeUtr1_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135522777,
			'maxpos'				=> 135522777,
			'wt' 						=> 'G',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getExonClass,$a->get3PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									1715,0,1715,0,'G','U','r.1715g>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get3PrimeUtrVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','c.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}

# FAM54A protein coding gene with 3 prime utr exons on - strand of genome

sub test3PrimeUtrSpliceRegion98_FAM54A {
	my $file = shift;

	subtest 'Testing FAM54A 3 prime UTR Splice Region - strand 98' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136554453,
			'maxpos'				=> 136554453,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->get3PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1128,10,1128,10,'U','A','r.1128+10u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass,$a->get3PrimeUtrVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','c.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic96_FAM54A {
	my $file = shift;

	subtest 'Testing FAM54A 3 prime UTR Intronic - strand 96' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136554452,
			'maxpos'				=> 136554452,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->get3PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic97_FAM54A {
	my $file = shift;

	subtest 'Testing FAM54A 3 prime UTR Intronic - strand 97' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136554363,
			'maxpos'				=> 136554363,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->get3PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic98_FAM54A {
	my $file = shift;

	subtest 'Testing FAM54A 3 prime UTR Intronic - strand 98' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136552626,
			'maxpos'				=> 136552626,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->get3PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic99_FAM54A {
	my $file = shift;

	subtest 'Testing FAM54A 3 prime UTR Intronic - strand 99' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136552537,
			'maxpos'				=> 136552537,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getIntronClass,$a->get3PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUtrSpliceRegion99_FAM54A {
	my $file = shift;

	subtest 'Testing FAM54A 3 prime UTR Splice Region - strand 99' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136552536,
			'maxpos'				=> 136552536,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getSpliceRegionClass,$a->get3PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1129,-10,1129,-10,'U','A','r.1129-10u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass,$a->get3PrimeUtrVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','c.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}
sub test3PrimeUtr1_FAM54A {
	my $file = shift;

	subtest 'Testing FAM54A 3 prime UTR - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136554637,
			'maxpos'				=> 136554637,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getExonClass,$a->get3PrimeUtrClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									954,0,954,0,'G','A','r.954g>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->get3PrimeUtrVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','c.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass);
		done_testing();
	};
}

# AC068831.2 lincRNA gene on + strand of genome

sub testExon1_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Exon + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91413150,
			'maxpos'				=> 91413150,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getExonClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									3,0,3,0,'C','U','r.3c>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getNonCodingTranscriptVariantClass);
		done_testing();
	};
}
sub testIntronic1_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Intronic + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573253,
			'maxpos'				=> 91573253,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic2_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Intronic + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573342,
			'maxpos'				=> 91573342,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testSpliceRegion1_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Splice Region + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573343,
			'maxpos'				=> 91573343,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									25,-10,25,-10,'A','U','r.25-10a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);
		done_testing();
	};
}
sub testEssentialSplice2_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Essential Splice + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573352,
			'maxpos'				=> 91573352,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									25,-1,25,-1,'A','U','r.25-1a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);
		done_testing();
	};
}
sub testEssentialSplice5_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Essential Splice + strand 5' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573476,
			'maxpos'				=> 91573476,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									143,5,143,5,'A','U','r.143+5a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);
		done_testing();
	};
}
sub testSpliceRegion15_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Splice Region + strand 15' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573481,
			'maxpos'				=> 91573481,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									143,10,143,10,'A','U','r.143+10a>u',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);
		done_testing();
	};
}
sub testIntronic3_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Intronic + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573482,
			'maxpos'				=> 91573482,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic4_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Intronic + strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573571,
			'maxpos'				=> 91573571,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}

# AC011503.1 lincRNA gene on - strand of genome

sub testExon1_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Exon - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24140728,
			'maxpos'				=> 24140728,
			'wt' 						=> 'C',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getExonClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									5,0,5,0,'G','A','r.5g>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getNonCodingTranscriptVariantClass);
		done_testing();
	};
}
sub testIntronic1_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Intronic - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24135338,
			'maxpos'				=> 24135338,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic2_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Intronic - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24135249,
			'maxpos'				=> 24135249,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testSpliceRegion1_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Splice Region - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24135248,
			'maxpos'				=> 24135248,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									201,-10,201,-10,'U','A','r.201-10u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);
		done_testing();
	};
}
sub testEssentialSplice1_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Essential Splice - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24135240,
			'maxpos'				=> 24135240,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									201,-2,201,-2,'U','A','r.201-2u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);
		done_testing();
	};
}
sub testEssentialSplice4_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Essential Splice - strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134970,
			'maxpos'				=> 24134970,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									467,2,467,2,'U','A','r.467+2u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getEssentialSpliceSiteVariantClass);
		done_testing();
	};
}
sub testSpliceRegion15_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Splice Region - strand 15' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134962,
			'maxpos'				=> 24134962,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									467,10,467,10,'U','A','r.467+10u>a',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getSpliceRegionVariantClass);
		done_testing();
	};
}
sub testIntronic3_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Intronic - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134961,
			'maxpos'				=> 24134961,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic4_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Intronic - strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Substitution->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134872,
			'maxpos'				=> 24134872,
			'wt' 						=> 'A',
			'mt'						=> 'T',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::SimpleSubstitutionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - lincRNA');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'doesnt have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'doesnt have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getNonProteinCodingClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
