#!/usr/bin/perl

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

use lib qw(lib ../lib t ../t);

use Test::More;

use Data::Dumper;

use Sanger::CGP::Vagrent::Data::Insertion;
use Sanger::CGP::Vagrent::Data::Transcript;
use Sanger::CGP::Vagrent::Data::Exon;
use Sanger::CGP::Vagrent::Data::Annotation;
use Sanger::CGP::Vagrent::Annotators::InsertionAnnotator;

use Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource;

use AnnotationTestUtils;

testIntronic();
testSplice();
testExonic();
testUpStreamDownStream();
testStrangeCases();
done_testing();

sub testUpStreamDownStream {
	testUpsteamMilesAway_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam5001bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam5000bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam4999bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam2000bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam1999bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam0bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	testStartsDownstream0bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownstream1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownstream499bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownstream500bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownstream4999bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownstream5000bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownstream5001bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownstreamMilesAway_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	testUpsteamMilesAway_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEndsUpsteam5001bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEndsUpsteam5000bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEndsUpsteam4999bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEndsUpsteam2000bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEndsUpsteam1999bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEndsUpsteam1bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEndsUpsteam0bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

	testStartsDownstream0bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownstream1bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownstream499bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownstream500bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownstream4999bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownstream5000bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownstream5001bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownstreamMilesAway_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

}

sub testStrangeCases {
#CENTRE OF INTRONS ' => sub {
		testIntronic_DeadCenterOfEvenSizedIntron_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testIntronic_StartingDeadCenterOfOddSizedIntron_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testIntronic_EndingDeadCenterOfOddSizedIntron_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

		testCDSStartAdjacent_1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testCDSStartAdjacent_3bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testCDSEndAdjacent_1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testCDSEndAdjacent_3bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

		testCDSStartAdjacent_1bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testCDSStartAdjacent_3bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testCDSEndAdjacent_1bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

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
		testSpliceRegion13_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testIntronic3_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testIntronic4_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

		testIntronic1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testIntronic2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion13_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testIntronic3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testIntronic4_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

		testIntronic1_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
		testIntronic2_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
		testSpliceRegion1_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
		testSpliceRegion13_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
		testIntronic3_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
		testIntronic4_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);

		testIntronic1_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
		testIntronic2_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
		testSpliceRegion1_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
		testSpliceRegion13_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
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
		testEssentialSplice5_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion9_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion10_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion11_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion12_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testSpliceRegion13_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
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
		testEssentialSplice5_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion9_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion10_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion11_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion12_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testSpliceRegion13_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testIntronic3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

		# handful of tests for non-coding genes
#		NON PROTIEN CODING + STRAND SPLICE' => sub {
		testIntronic2_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
		testSpliceRegion1_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
		testEssentialSplice1_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
		testEssentialSplice5_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
		testSpliceRegion13_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
		testIntronic3_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);

#		NON PROTIEN CODING - STRAND SPLICE
		testIntronic2_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
		testSpliceRegion1_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
		testEssentialSplice1_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
		testEssentialSplice3_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
		testSpliceRegion13_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
		testIntronic3_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);

}
sub testExonic {

	#5 PRIME UTR 1bp INS
		test5PrimeUTR_1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		test5PrimeUTR_1bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);


	#5 PRIME UTR 3bp INS
		test5PrimeUTR_3bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		test5PrimeUTR_3bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);


	#5 PRIME UTR 5bp INS
		test5PrimeUTR_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		test5PrimeUTR_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);


	#5 PRIME UTR 9bp INS
		test5PrimeUTR_9bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		test5PrimeUTR_9bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);


	#3 PRIME UTR 1bp INS
		test3PrimeUtr_1bp_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
		test3PrimeUtr_1bp_FAM54A(AnnotationTestUtils::FAM54A_TRANSCRIPT);


	#3 PRIME UTR 3bp INS
		test3PrimeUtr_3bp_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
		test3PrimeUtr_3bp_FAM54A(AnnotationTestUtils::FAM54A_TRANSCRIPT);


	#3 PRIME UTR 5bp INS
		test3PrimeUtr_5bp_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
		test3PrimeUtr_5bp_FAM54A(AnnotationTestUtils::FAM54A_TRANSCRIPT);


	#3 PRIME UTR 9bp INS
		test3PrimeUtr_9bp_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
		test3PrimeUtr_9bp_FAM54A(AnnotationTestUtils::FAM54A_TRANSCRIPT);


	#EXON BOUNDRY
		testEssentialSplice1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testCodingExonBoundry_1bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testCodingExonBoundry_5bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testCodingExonBoundry_1bp_2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testCodingExonBoundry_5bp_2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testEssentialSplice2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

		testEssentialSplice1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testCodingExonBoundry_1bp_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testCodingExonBoundry_5bp_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testCodingExonBoundry_1bp_2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testCodingExonBoundry_5bp_2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
		testEssentialSplice2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

	#CODING EXON 5bp INS
		testCodingExon_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	#CODING EXON IN PHASE 9bp INS
		testCodingExonInPhase_9bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testCodingExonInPhase_9bp_2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);


	#CODING EXON OUT OF PHASE 3bp INS
		testCodingExonOutOfPhase_3bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testCodingExonOutOfPhase_3bp_2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	#CODING EXON 1bp INS
		testCodingExon_1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testCodingExon_1bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

	#CODING EXON IN PHASE 3bp INS
		testCodingExonInPhase_3bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testCodingExonInPhase_3bp_2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

		testCodingExonInPhase_3bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

	#CODING EXON IN FRAME 15bp INS STOP
		testCodingExonInPhase_InsStop_15bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testCodingExonOutOfPhase_InsStop_15bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
		testCodingExonOutOfPhase_InsStop_15bp_2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

}


# CEP350 protein coding gene with 5 prime utr exons on + strand of genome

sub testUpsteamMilesAway_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Upstream Miles Away + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179913873,
			'maxpos'				=> 179913874,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testEndsUpsteam5001bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 5001 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179918871,
			'maxpos'				=> 179918872,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testEndsUpsteam5000bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 5000 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179918872,
			'maxpos'				=> 179918873,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testEndsUpsteam4999bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 4999 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179918873,
			'maxpos'				=> 179918874,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get5KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam2000bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 2000 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179921872,
			'maxpos'				=> 179921873,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get5KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam1999bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 1999 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179921873,
			'maxpos'				=> 179921874,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get2KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 1 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179923871,
			'maxpos'				=> 179923872,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get2KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam0bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 0 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179923872,
			'maxpos'				=> 179923873,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get2KBUpStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownstream0bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 0 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084015,
			'maxpos'				=> 180084016,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get500BPDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownstream1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 1 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084016,
			'maxpos'				=> 180084017,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get500BPDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownstream499bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 499 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084514,
			'maxpos'				=> 180084515,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get500BPDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownstream500bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 500 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084515,
			'maxpos'				=> 180084516,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownstream4999bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 4999 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180089014,
			'maxpos'				=> 180089015,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownstream5000bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 5000 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180089015,
			'maxpos'				=> 180089016,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testStartsDownstream5001bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 5001 + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180089016,
			'maxpos'				=> 180089017,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testStartsDownstreamMilesAway_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream Miles Away + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180094015,
			'maxpos'				=> 180094016,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}

sub test5PrimeUtrSpliceRegion7_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 prime UTR Splice Region + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924286,
			'maxpos'				=> 179924287,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									405,9,405,10,'-','A','r.405+9_405+10insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass,$a->get5PrimeUtrVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924287,
			'maxpos'				=> 179924288,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Intronic + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924388,
			'maxpos'				=> 179924389,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic3_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Intronic + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955204,
			'maxpos'				=> 179955205,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic4_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Intronic + strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955293,
			'maxpos'				=> 179955294,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUtrSpliceRegion8_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 prime UTR Splice Region + strand 8' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955294,
			'maxpos'				=> 179955295,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									406,-10,406,-9,'-','A','r.406-10_406-9insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass,$a->get5PrimeUtrVariantClass);

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

sub test5PrimeUTR_1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 prime UTR 1bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179923875,
			'maxpos'				=> 179923876,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									3,0,4,0,'-','A','r.3_4insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTR_3bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 prime UTR 3bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179923875,
			'maxpos'				=> 179923876,
			'insseq' 				=> 'ATA');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									3,0,4,0,'-','AUA','r.3_4insaua',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTR_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 prime UTR 5bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179923875,
			'maxpos'				=> 179923876,
			'insseq' 				=> 'ATACC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									3,0,4,0,'-','AUACC','r.3_4insauacc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTR_9bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 prime UTR 9bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179923875,
			'maxpos'				=> 179923876,
			'insseq' 				=> 'ATACCGGGG');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									3,0,4,0,'-','AUACCGGGG','r.3_4insauaccgggg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get5PrimeUtrVariantClass);

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

sub testCDSStartAdjacent_1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 prime UTR 1bp + strand CDS start adjacent' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955316,
			'maxpos'				=> 179955317,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									418,0,419,0,'-','A','r.418_419insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get5PrimeUtrVariantClass);

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
sub testCDSStartAdjacent_3bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 prime UTR 3bp + strand CDS start adjacent' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955316,
			'maxpos'				=> 179955317,
			'insseq' 				=> 'ATA');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									418,0,419,0,'-','AUA','r.418_419insaua',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get5PrimeUtrVariantClass);

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

sub testCDSEndAdjacent_1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 prime UTR 1bp + strand CDS end adjacent' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180080296,
			'maxpos'				=> 180080297,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get3PrimeUtrClass,$a->getExonClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									9772,0,9773,0,'-','A','r.9772_9773insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get3PrimeUtrVariantClass);

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
sub testCDSEndAdjacent_3bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 prime UTR 3bp + strand CDS end adjacent' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180080296,
			'maxpos'				=> 180080297,
			'insseq' 				=> 'ACC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get3PrimeUtrClass,$a->getExonClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									9772,0,9773,0,'-','ACC','r.9772_9773insacc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get3PrimeUtrVariantClass);

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

sub testIntronic1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Intronic + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961097,
			'maxpos'				=> 179961098,
			'insseq' 						=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Inronic + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961186,
			'maxpos'				=> 179961187,
			'insseq' 						=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testSpliceRegion1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961187,
			'maxpos'				=> 179961188,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-10,654,-9,'-','A','r.654-10_654-9insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-10,236,-9,'-','A','c.236-10_236-9insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961188,
			'maxpos'				=> 179961189,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-9,654,-8,'-','A','r.654-9_654-8insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-9,236,-8,'-','A','c.236-9_236-8insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961189,
			'maxpos'				=> 179961190,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-8,654,-7,'-','A','r.654-8_654-7insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-8,236,-7,'-','A','c.236-8_236-7insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961190,
			'maxpos'				=> 179961191,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-7,654,-6,'-','A','r.654-7_654-6insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-7,236,-6,'-','A','c.236-7_236-6insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961191,
			'maxpos'				=> 179961192,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-6,654,-5,'-','A','r.654-6_654-5insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-6,236,-5,'-','A','c.236-6_236-5insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961192,
			'maxpos'				=> 179961193,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-5,654,-4,'-','A','r.654-5_654-4insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-5,236,-4,'-','A','c.236-5_236-4insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961193,
			'maxpos'				=> 179961194,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-4,654,-3,'-','A','r.654-4_654-3insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-4,236,-3,'-','A','c.236-4_236-3insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961194,
			'maxpos'				=> 179961195,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-3,654,-2,'-','A','r.654-3_654-2insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-3,236,-2,'-','A','c.236-3_236-2insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961195,
			'maxpos'				=> 179961196,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-2,654,-1,'-','A','r.654-2_654-1insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-2,236,-1,'-','A','c.236-2_236-1insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

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

sub testCodingExonBoundry_1bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Coding Exon Boundry 1bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961196,
			'maxpos'				=> 179961197,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-1,654,0,'-','A','r.654-1_654insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-1,236,0,'-','A','c.236-1_236insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									79,0,79,0,'D','EW*','p.D79fs*3',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getFrameShiftVariantClass);
		done_testing();
	};
}
sub testCodingExonBoundry_5bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Coding Exon Boundry 5bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961196,
			'maxpos'				=> 179961197,
			'insseq' 				=> 'CCGGT');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-1,654,0,'-','CCGGU','r.654-1_654insccggu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-1,236,0,'-','CCGGT','c.236-1_236insCCGGT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									79,0,79,0,'D','AGMVDTWMILGLMLQSPNPLNHEKRNLVVLSGPPPWRVM*','p.D79fs*40',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getFrameShiftVariantClass);
		done_testing();
	};
}

sub testCodingExon_1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Coding Exon 1bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961198,
			'maxpos'				=> 179961199,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									655,0,656,0,'-','A','r.655_656insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									237,0,238,0,'-','A','c.237_238insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									80,0,80,0,'G','R*','p.G80fs*2',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getFrameShiftVariantClass);
		done_testing();
	};
}

sub testCodingExonInPhase_3bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Coding Exon In Phase New Codon 3bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961198,
			'maxpos'				=> 179961199,
			'insseq' 				=> 'AAA');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									655,0,656,0,'-','AAA','r.655_656insaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									237,0,238,0,'-','AAA','c.237_238insAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									79,0,80,0,'-','K','p.D79_G80insK',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInFrameCodonGainVariantClass);
		done_testing();
	};
}
sub testCodingExonInPhase_3bp_2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Coding Exon In Phase Duplicate Codon 3bp + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961198,
			'maxpos'				=> 179961199,
			'insseq' 				=> 'GGA');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									655,0,656,0,'-','GGA','r.655_656insgga',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									237,0,238,0,'-','GGA','c.237_238insGGA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									80,0,81,0,'-','G','p.G80_R81insG',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInFrameCodonGainVariantClass);
		done_testing();
	};
}

sub testCodingExonOutOfPhase_3bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Coding Exon Out Of Phase 3bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961199,
			'maxpos'				=> 179961200,
			'insseq' 				=> 'AAA');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									656,0,657,0,'-','AAA','r.656_657insaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									238,0,239,0,'-','AAA','c.238_239insAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									80,0,80,0,'G','ES','p.G80delinsES',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass);
		done_testing();
	};
}
sub testCodingExonOutOfPhase_3bp_2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Coding Exon Out Of Phase 3bp + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961200,
			'maxpos'				=> 179961201,
			'insseq' 				=> 'AAA');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									657,0,658,0,'-','AAA','r.657_658insaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									239,0,240,0,'-','AAA','c.239_240insAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									80,0,81,0,'-','N','p.G80_R81insN',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInFrameCodonGainVariantClass);
		done_testing();
	};
}

sub testCodingExon_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Coding Exon 5bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961198,
			'maxpos'				=> 179961199,
			'insseq' 				=> 'ATTGG');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									655,0,656,0,'-','AUUGG','r.655_656insauugg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									237,0,238,0,'-','ATTGG','c.237_238insATTGG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									80,0,80,0,'G','IGVDTWMILGLMLQSPNPLNHEKRNLVVLSGPPPWRVM*','p.G80fs*39',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getFrameShiftVariantClass);
		done_testing();
	};
}

sub testCodingExonInPhase_9bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Coding Exon In Phase New Codon 9bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961201,
			'maxpos'				=> 179961202,
			'insseq' 				=> 'AAATTTAAA');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									658,0,659,0,'-','AAAUUUAAA','r.658_659insaaauuuaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									240,0,241,0,'-','AAATTTAAA','c.240_241insAAATTTAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									80,0,81,0,'-','KFK','p.G80_R81insKFK',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInFrameCodonGainVariantClass);
		done_testing();
	};
}
sub testCodingExonInPhase_9bp_2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Coding Exon In Phase Duplicate Codon 9bp + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961201,
			'maxpos'				=> 179961202,
			'insseq' 				=> 'GGATTTAAA');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									658,0,659,0,'-','GGAUUUAAA','r.658_659insggauuuaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									240,0,241,0,'-','GGATTTAAA','c.240_241insGGATTTAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									80,0,81,0,'-','GFK','p.G80_R81insGFK',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInFrameCodonGainVariantClass);
		done_testing();
	};
}

sub testCodingExonInPhase_InsStop_15bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Coding Exon In Phase Ins Stop 15bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961198,
			'maxpos'				=> 179961199,
			'insseq' 				=> 'GGAATTTAAATTGGA');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									655,0,656,0,'-','GGAAUUUAAAUUGGA','r.655_656insggaauuuaaauugga',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									237,0,238,0,'-','GGAATTTAAATTGGA','c.237_238insGGAATTTAAATTGGA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									80,0,81,0,'-','I*','p.G80_R81insI*',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInFrameCodonGainVariantClass,$a->getStopGainedVariantClass);
		done_testing();
	};
}
sub testCodingExonOutOfPhase_InsStop_15bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Coding Exon Out of Phase Ins Stop 15bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961197,
			'maxpos'				=> 179961198,
			'insseq' 				=> 'CCCCATTTAAATTCC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									654,0,655,0,'-','CCCCAUUUAAAUUCC','r.654_655insccccauuuaaauucc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									236,0,237,0,'-','CCCCATTTAAATTCC','c.236_237insCCCCATTTAAATTCC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									79,0,80,0,'-','PI*','p.D79_G80insPI*',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInFrameCodonGainVariantClass,$a->getStopGainedVariantClass);
		done_testing();
	};
}
sub testCodingExonOutOfPhase_InsStop_15bp_2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Coding Exon Out of Phase Ins Stop 15bp + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961196,
			'maxpos'				=> 179961197,
			'insseq' 				=> 'CCCCCATTTAAATTC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,-1,654,0,'-','CCCCCAUUUAAAUUC','r.654-1_654inscccccauuuaaauuc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									236,-1,236,0,'-','CCCCCATTTAAATTC','c.236-1_236insCCCCCATTTAAATTC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									79,0,79,0,'D','API*','p.D79delinsAPI*',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getStopGainedVariantClass);
		done_testing();
	};
}

sub testCodingExonBoundry_1bp_2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Coding Exon Boundry 1bp + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961356,
			'maxpos'				=> 179961357,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,0,813,1,'-','A','r.813_813+1insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,0,395,1,'-','A','c.395_395+1insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									133,0,133,0,'E','GNPWCTFQFQFQPSGIKARILCRC*','p.E133fs*25',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getFrameShiftVariantClass);
		done_testing();
	};
}
sub testCodingExonBoundry_5bp_2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Coding Exon Boundry 5bp + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961356,
			'maxpos'				=> 179961357,
			'insseq' 				=> 'ATGCC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,0,813,1,'-','AUGCC','r.813_813+1insaugcc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,0,395,1,'-','ATGCC','c.395_395+1insATGCC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									133,0,133,0,'E','CRKSMVHLPISVPAIWNQSTYTV*','p.E133fs*24',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getFrameShiftVariantClass);
		done_testing();
	};
}

sub testEssentialSplice2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961357,
			'maxpos'				=> 179961358,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,1,813,2,'-','A','r.813+1_813+2insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,1,395,2,'-','A','c.395+1_395+2insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

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
sub testEssentialSplice3_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961358,
			'maxpos'				=> 179961359,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,2,813,3,'-','A','r.813+2_813+3insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,2,395,3,'-','A','c.395+2_395+3insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961359,
			'maxpos'				=> 179961360,
			'insseq' 				=> 'AA');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,3,813,4,'-','AA','r.813+3_813+4insaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,3,395,4,'-','AA','c.395+3_395+4insAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961360,
			'maxpos'				=> 179961361,
			'insseq' 				=> 'AA');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,4,813,5,'-','AA','r.813+4_813+5insaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,4,395,5,'-','AA','c.395+4_395+5insAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961361,
			'maxpos'				=> 179961362,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,5,813,6,'-','A','r.813+5_813+6insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,5,395,6,'-','A','c.395+5_395+6insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961362,
			'maxpos'				=> 179961363,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,6,813,7,'-','A','r.813+6_813+7insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,6,395,7,'-','A','c.395+6_395+7insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961363,
			'maxpos'				=> 179961364,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,7,813,8,'-','A','r.813+7_813+8insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,7,395,8,'-','A','c.395+7_395+8insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961364,
			'maxpos'				=> 179961365,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,8,813,9,'-','A','r.813+8_813+9insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,8,395,9,'-','A','c.395+8_395+9insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961365,
			'maxpos'				=> 179961366,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									813,9,813,10,'-','U','r.813+9_813+10insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									395,9,395,10,'-','T','c.395+9_395+10insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961366,
			'maxpos'				=> 179961367,
			'insseq' 						=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic4_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Intronic + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961456,
			'maxpos'				=> 179961457,
			'insseq' 						=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}

sub testIntronic_DeadCenterOfEvenSizedIntron_CEP350{
my $file = shift;

	subtest 'Testing CEP350 Dead Center of Even Sized Intron + strand' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179960476,
			'maxpos'				=> 179960477,
			'insseq' 						=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic_StartingDeadCenterOfOddSizedIntron_CEP350{
my $file = shift;

	subtest 'Testing CEP350 Starting Dead Center of Odd Sized Intron + strand' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179963522,
			'maxpos'				=> 179963523,
			'insseq' 						=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic_EndingDeadCenterOfOddSizedIntron_CEP350{
my $file = shift;

	subtest 'Testing CEP350 Ending Dead Center of Odd Sized Intron + strand' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179963521,
			'maxpos'				=> 179963522,
			'insseq' 						=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}

# TOR1AIP2 protein coding gene with 5 prime utr exons on - strand of genome

sub testUpsteamMilesAway_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Upstream Miles Away - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179856933,
			'maxpos'				=> 179856934,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testEndsUpsteam5001bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Ends Upstream 5001 - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179851935,
			'maxpos'				=> 179851936,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testEndsUpsteam5000bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Ends Upstream 5000 - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179851934,
			'maxpos'				=> 179851935,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testEndsUpsteam4999bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Ends Upstream 4999 - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179851933,
			'maxpos'				=> 179851934,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get5KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam2000bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Ends Upstream 2000 - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179848934,
			'maxpos'				=> 179848935,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get5KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam1999bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Ends Upstream 1999 - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179848933,
			'maxpos'				=> 179848934,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get2KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam1bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Ends Upstream 1 - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846935,
			'maxpos'				=> 179846936,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get2KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam0bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Ends Upstream 0 - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846934,
			'maxpos'				=> 179846935,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get2KBUpStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownstream0bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream 0 - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179809101,
			'maxpos'				=> 179809102,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get500BPDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownstream1bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream 1 - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179809100,
			'maxpos'				=> 179809101,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get500BPDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownstream499bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream 499 - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179808602,
			'maxpos'				=> 179808603,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get500BPDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownstream500bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream 500 - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179808601,
			'maxpos'				=> 179808602,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownstream4999bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream 4999 - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179804102,
			'maxpos'				=> 179804103,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownstream5000bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream 5000 - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179804101,
			'maxpos'				=> 179804102,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testStartsDownstream5001bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream 5001 - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179804100,
			'maxpos'				=> 179804101,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testStartsDownstreamMilesAway_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream Miles Away - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179704100,
			'maxpos'				=> 179704101,
			'insseq' 						=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846364,
			'maxpos'				=> 179846365,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									242,9,242,10,'-','U','r.242+9_242+10insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass,$a->get5PrimeUtrVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846363,
			'maxpos'				=> 179846364,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Intronic - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846273,
			'maxpos'				=> 179846274,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic3_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Intronic - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179822046,
			'maxpos'				=> 179822047,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic4_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Intronic - strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821956,
			'maxpos'				=> 179821957,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUtrSpliceRegion8_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 prime UTR Splice Region - strand 8' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821955,
			'maxpos'				=> 179821956,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		print join("\n",$a->getMessages);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									243,-10,243,-9,'-','U','r.243-10_243-9insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass,$a->get5PrimeUtrVariantClass);

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

sub test5PrimeUTR_1bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 prime UTR 1bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821803,
			'maxpos'				=> 179821804,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									385,0,386,0,'-','U','r.385_386insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTR_3bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 prime UTR 3bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821803,
			'maxpos'				=> 179821804,
			'insseq' 				=> 'ATA');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									385,0,386,0,'-','UAU','r.385_386insuau',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTR_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 prime UTR 5bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821803,
			'maxpos'				=> 179821804,
			'insseq' 				=> 'ATACC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									385,0,386,0,'-','GGUAU','r.385_386insgguau',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTR_9bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 prime UTR 9bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821803,
			'maxpos'				=> 179821804,
			'insseq' 				=> 'ATACCGGGG');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									385,0,386,0,'-','CCCCGGUAU','r.385_386insccccgguau',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get5PrimeUtrVariantClass);

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

sub testCDSStartAdjacent_1bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 prime UTR 1bp - strand CDS start adjacent' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821800,
			'maxpos'				=> 179821801,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									388,0,389,0,'-','U','r.388_389insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get5PrimeUtrVariantClass);

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
sub testCDSStartAdjacent_3bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 prime UTR 1bp - strand CDS start adjacent' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821800,
			'maxpos'				=> 179821801,
			'insseq' 				=> 'AGG');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									388,0,389,0,'-','CCU','r.388_389insccu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get5PrimeUtrVariantClass);

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

sub testCDSEndAdjacent_1bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 prime UTR 1bp - strand CDS end adjacent' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179815205,
			'maxpos'				=> 179815206,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get3PrimeUtrClass,$a->getExonClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									1801,0,1802,0,'-','U','r.1801_1802insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get3PrimeUtrVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820597,
			'maxpos'				=> 179820598,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Intronic - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820508,
			'maxpos'				=> 179820509,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testSpliceRegion1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820507,
			'maxpos'				=> 179820508,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-10,423,-9,'-','U','r.423-10_423-9insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-10,35,-9,'-','T','c.35-10_35-9insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820506,
			'maxpos'				=> 179820507,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-9,423,-8,'-','U','r.423-9_423-8insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-9,35,-8,'-','T','c.35-9_35-8insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820505,
			'maxpos'				=> 179820506,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-8,423,-7,'-','U','r.423-8_423-7insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-8,35,-7,'-','T','c.35-8_35-7insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820504,
			'maxpos'				=> 179820505,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-7,423,-6,'-','U','r.423-7_423-6insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-7,35,-6,'-','T','c.35-7_35-6insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820503,
			'maxpos'				=> 179820504,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-6,423,-5,'-','U','r.423-6_423-5insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-6,35,-5,'-','T','c.35-6_35-5insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820502,
			'maxpos'				=> 179820503,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-5,423,-4,'-','U','r.423-5_423-4insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-5,35,-4,'-','T','c.35-5_35-4insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820501,
			'maxpos'				=> 179820502,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-4,423,-3,'-','U','r.423-4_423-3insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-4,35,-3,'-','T','c.35-4_35-3insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820500,
			'maxpos'				=> 179820501,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-3,423,-2,'-','U','r.423-3_423-2insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-3,35,-2,'-','T','c.35-3_35-2insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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

	subtest 'Testing TOR1AIP2 Essential Splice Region - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820499,
			'maxpos'				=> 179820500,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-2,423,-1,'-','U','r.423-2_423-1insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-2,35,-1,'-','T','c.35-2_35-1insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

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

sub testCodingExonBoundry_1bp_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Coding Exon Boundry 1bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820498,
			'maxpos'				=> 179820499,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-1,423,0,'-','U','r.423-1_423insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-1,35,0,'-','T','c.35-1_35insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									12,0,12,0,'D','VLSKGFGK*','p.D12fs*9',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getFrameShiftVariantClass);
		done_testing();
	};
}
sub testCodingExonBoundry_5bp_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Coding Exon Boundry 5bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820498,
			'maxpos'				=> 179820499,
			'insseq' 				=> 'AGATT');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-1,423,0,'-','AAUCU','r.423-1_423insaaucu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-1,35,0,'-','AATCT','c.35-1_35insAATCT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									12,0,12,0,'D','ESTLKRIWKMIHQ*','p.D12fs*14',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getFrameShiftVariantClass);
		done_testing();
	};
}

sub testCodingExon_1bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Coding Exon 1bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820492,
			'maxpos'				=> 179820493,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									428,0,429,0,'-','U','r.428_429insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									40,0,41,0,'-','T','c.40_41insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									14,0,14,0,'Q','LKGFGK*','p.Q14fs*7',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getFrameShiftVariantClass);
		done_testing();
	};
}

sub testCodingExonInPhase_3bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Coding Exon In Phase New Codon 3bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820493,
			'maxpos'				=> 179820494,
			'insseq' 				=> 'ATC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									427,0,428,0,'-','GAU','r.427_428insgau',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									39,0,40,0,'-','GAT','c.39_40insGAT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									13,0,14,0,'-','D','p.S13_Q14insD',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInFrameCodonGainVariantClass);
		done_testing();
	};
}

sub testCodingExonBoundry_1bp_2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Coding Exon Boundry 1bp - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819979,
			'maxpos'				=> 179819980,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,0,941,1,'-','U','r.941_941+1insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,0,553,1,'-','T','c.553_553+1insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									185,0,185,0,'E','VGWKSSTTNTKIGRN*','p.E185fs*16',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getFrameShiftVariantClass);
		done_testing();
	};
}
sub testCodingExonBoundry_5bp_2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Coding Exon Boundry 5bp - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819979,
			'maxpos'				=> 179819980,
			'insseq' 				=> 'ACCGT');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,0,941,1,'-','ACGGU','r.941_941+1insacggu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,0,553,1,'-','ACGGT','c.553_553+1insACGGT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									185,0,185,0,'E','DGRLEVIHNKHKNWKKLKRMHRTP*','p.E185fs*25',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getFrameShiftVariantClass);
		done_testing();
	};
}

sub testEssentialSplice2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice Region - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819978,
			'maxpos'				=> 179819979,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,1,941,2,'-','U','r.941+1_941+2insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,1,553,2,'-','T','c.553+1_553+2insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

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

	subtest 'Testing TOR1AIP2 Essential Splice Region - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819977,
			'maxpos'				=> 179819978,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,2,941,3,'-','U','r.941+2_941+3insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,2,553,3,'-','T','c.553+2_553+3insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

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

	subtest 'Testing TOR1AIP2 Essential Splice Region - strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819976,
			'maxpos'				=> 179819977,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,3,941,4,'-','U','r.941+3_941+4insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,3,553,4,'-','T','c.553+3_553+4insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

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

	subtest 'Testing TOR1AIP2 Essential Splice Region - strand 5' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819975,
			'maxpos'				=> 179819976,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,4,941,5,'-','U','r.941+4_941+5insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,4,553,5,'-','T','c.553+4_553+5insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819974,
			'maxpos'				=> 179819975,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,5,941,6,'-','U','r.941+5_941+6insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,5,553,6,'-','T','c.553+5_553+6insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819973,
			'maxpos'				=> 179819974,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,6,941,7,'-','U','r.941+6_941+7insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,6,553,7,'-','T','c.553+6_553+7insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819972,
			'maxpos'				=> 179819973,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,7,941,8,'-','U','r.941+7_941+8insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,7,553,8,'-','T','c.553+7_553+8insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819971,
			'maxpos'				=> 179819972,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,8,941,9,'-','U','r.941+8_941+9insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,8,553,9,'-','T','c.553+8_553+9insT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819970,
			'maxpos'				=> 179819971,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									941,9,941,10,'-','A','r.941+9_941+10insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									553,9,553,10,'-','A','c.553+9_553+10insA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819969,
			'maxpos'				=> 179819970,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic4_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Intronic - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179819880,
			'maxpos'				=> 179819881,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}

# MYB protein coding gene with 3 prime utr exons on + strand of genome

sub test3PrimeUtrSpliceRegion98_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Splice Region + strand 98' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135524471,
			'maxpos'				=> 135524472,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1933,9,1933,10,'-','U','r.1933+9_1933+10insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass,$a->get3PrimeUtrVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135524472,
			'maxpos'				=> 135524473,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic97_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Intronic + strand 97' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135524561,
			'maxpos'				=> 135524562,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic98_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Intronic + strand 98' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135538902,
			'maxpos'				=> 135538903,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic99_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Intronic + strand 99' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135538991,
			'maxpos'				=> 135538992,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUtrSpliceRegion99_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Splice Region + strand 99' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135538992,
			'maxpos'				=> 135538993,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1934,-10,1934,-9,'-','U','r.1934-10_1934-9insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass,$a->get3PrimeUtrVariantClass);

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

sub test3PrimeUtr_1bp_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR 1bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135522777,
			'maxpos'				=> 135522778,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									1715,0,1716,0,'-','A','r.1715_1716insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUtr_3bp_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR 3bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135522777,
			'maxpos'				=> 135522778,
			'insseq' 				=> 'AGG');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									1715,0,1716,0,'-','AGG','r.1715_1716insagg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUtr_5bp_MYB {
	my $file = shift;

	subtest 'Testing MYB 5 prime UTR 3bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135522777,
			'maxpos'				=> 135522778,
			'insseq' 				=> 'AGGTT');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									1715,0,1716,0,'-','AGGUU','r.1715_1716insagguu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUtr_9bp_MYB {
	my $file = shift;

	subtest 'Testing MYB 5 prime UTR 9bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135522777,
			'maxpos'				=> 135522778,
			'insseq' 				=> 'AGGTTATCG');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									1715,0,1716,0,'-','AGGUUAUCG','r.1715_1716insagguuaucg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get3PrimeUtrVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136554453,
			'maxpos'				=> 136554454,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1128,9,1128,10,'-','A','r.1128+9_1128+10insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass,$a->get3PrimeUtrVariantClass);

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


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136554452,
			'maxpos'				=> 136554453,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic97_FAM54A {
	my $file = shift;

	subtest 'Testing FAM54A 3 prime UTR Intronic - strand 97' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136554363,
			'maxpos'				=> 136554364,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic98_FAM54A {
	my $file = shift;

	subtest 'Testing FAM54A 3 prime UTR Intronic - strand 98' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136552626,
			'maxpos'				=> 136552627,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic99_FAM54A {
	my $file = shift;

	subtest 'Testing FAM54A 3 prime UTR Intronic - strand 99' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136552536,
			'maxpos'				=> 136552537,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUtrSpliceRegion99_FAM54A {
	my $file = shift;

	subtest 'Testing FAM54A 3 prime UTR Splice Region - strand 99' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136552535,
			'maxpos'				=> 136552536,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1129,-10,1129,-9,'-','A','r.1129-10_1129-9insa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass,$a->get3PrimeUtrVariantClass);

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

sub test3PrimeUtr_1bp_FAM54A {
	my $file = shift;

	subtest 'Testing FAM54A 3 prime UTR 1bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136554636,
			'maxpos'				=> 136554637,
			'insseq' 				=> 'A');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									954,0,955,0,'-','U','r.954_955insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUtr_3bp_FAM54A {
	my $file = shift;

	subtest 'Testing FAM54A 3 prime UTR 3bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136554636,
			'maxpos'				=> 136554637,
			'insseq' 				=> 'ACG');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									954,0,955,0,'-','CGU','r.954_955inscgu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUtr_5bp_FAM54A {
	my $file = shift;

	subtest 'Testing FAM54A 3 prime UTR 5bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136554636,
			'maxpos'				=> 136554637,
			'insseq' 				=> 'ACGCC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									954,0,955,0,'-','GGCGU','r.954_955insggcgu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUtr_9bp_FAM54A {
	my $file = shift;

	subtest 'Testing FAM54A 3 prime UTR 9bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 136554636,
			'maxpos'				=> 136554637,
			'insseq' 				=> 'ACGCCTTTT');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									954,0,955,0,'-','AAAAGGCGU','r.954_955insaaaaggcgu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->get3PrimeUtrVariantClass);

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

sub testIntronic1_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Intronic + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573252,
			'maxpos'				=> 91573253,
			'insseq' 				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic2_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Intronic + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573342,
			'maxpos'				=> 91573343,
			'insseq' 				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testSpliceRegion1_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Splice Region + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573343,
			'maxpos'				=> 91573344,
			'insseq' 				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									25,-10,25,-9,'-','C','r.25-10_25-9insc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);
		done_testing();
	};
}
sub testEssentialSplice1_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Essential Splice + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573351,
			'maxpos'				=> 91573352,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									25,-2,25,-1,'-','U','r.25-2_25-1insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);
		done_testing();
	};
}
sub testEssentialSplice5_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Essential Splice + strand 5' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573475,
			'maxpos'				=> 91573476,
			'insseq' 				=> 'T');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									143,4,143,5,'-','U','r.143+4_143+5insu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);
		done_testing();
	};
}
sub testSpliceRegion13_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Splice Region + strand 13' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573480,
			'maxpos'				=> 91573481,
			'insseq' 				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									143,9,143,10,'-','C','r.143+9_143+10insc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);
		done_testing();
	};
}
sub testIntronic3_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Intronic + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573481,
			'maxpos'				=> 91573482,
			'insseq' 				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic4_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Intronic + strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573570,
			'maxpos'				=> 91573571,
			'insseq' 				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}

# AC011503.1 lincRNA gene on - strand of genome

sub testIntronic1_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Intronic - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24135338,
			'maxpos'				=> 24135339,
			'insseq' 				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic2_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Intronic - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24135248,
			'maxpos'				=> 24135249,
			'insseq' 				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testSpliceRegion1_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Splice Region - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24135247,
			'maxpos'				=> 24135248,
			'insseq' 				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									201,-10,201,-9,'-','G','r.201-10_201-9insg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);
		done_testing();
	};
}
sub testEssentialSplice1_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Essential Splice - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24135239,
			'maxpos'				=> 24135240,
			'insseq' 				=> 'CCC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									201,-2,201,-1,'-','GGG','r.201-2_201-1insggg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);
		done_testing();
	};
}
sub testEssentialSplice3_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Essential Splice - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134969,
			'maxpos'				=> 24134970,
			'insseq' 				=> 'CCC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									467,2,467,3,'-','GGG','r.467+2_467+3insggg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getEssentialSpliceSiteVariantClass);
		done_testing();
	};
}
sub testSpliceRegion13_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Splice Region - strand 13' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134962,
			'maxpos'				=> 24134963,
			'insseq' 				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getInsertionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									467,9,467,10,'-','G','r.467+9_467+10insg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getInsertionClass,$a->getSpliceRegionVariantClass);
		done_testing();
	};
}
sub testIntronic3_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Intronic - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134961,
			'maxpos'				=> 24134962,
			'insseq' 				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic4_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Intronic - strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::Insertion->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134871,
			'maxpos'				=> 24134872,
			'insseq' 				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::InsertionAnnotator->new(transcriptSource => $ts);

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
									$a->getInsertionClass,$a->getIntronVariantClass);
		done_testing();
	};
}
