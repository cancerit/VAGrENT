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

use warnings FATAL => 'all';
use Test::More;

use Data::Dumper;

use Sanger::CGP::Vagrent::Data::ComplexIndel;
use Sanger::CGP::Vagrent::Data::Transcript;
use Sanger::CGP::Vagrent::Data::Exon;
use Sanger::CGP::Vagrent::Data::Annotation;
use Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator;

use Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource;

use AnnotationTestUtils;

testUpStreamDownStream();
testIntronic();
testSplice();
testNonCodingExon();
testCodingExon();
testBoundyCrossing();
testLargerChanges();
done_testing();

sub testIntronic {
	#Testing 5 PRIME UTR INTRONIC
	test5PrimeUTRIntronic_1bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTRIntronic_1bp_2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTRIntronic_1bp_3_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTRIntronic_10bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTRIntronic_10bp_2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	test5PrimeUTRIntronic_1bp_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUTRIntronic_1bp_2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUTRIntronic_1bp_3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUTRIntronic_10bp_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUTRIntronic_10bp_2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

	#Testing CDS INTRONIC

	testCDSIntronic_1bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testCDSIntronic_10bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testCDSIntronic_1bp_2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testCDSIntronic_10bp_2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testCDSIntronic_1bp_3_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testCDSIntronic_10bp_3_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	testCDSIntronic_1bp_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testCDSIntronic_10bp_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testCDSIntronic_1bp_2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testCDSIntronic_10bp_2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testCDSIntronic_1bp_3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testCDSIntronic_10bp_3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

	#Testing 3 PRIME UTR INTRONIC

	test3PrimeUTRIntronic_1bp_96_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
	test3PrimeUTRIntronic_10bp_96_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);

	test3PrimeUTRIntronic_1bp_99_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
	test3PrimeUTRIntronic_10bp_99_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);

	#Testing NON-CODING TRANSCRIPT INTRONIC

	testIntronic_1bp_1_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
	testIntronic_1bp_2_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
	testIntronic_1bp_3_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);
	testIntronic_1bp_4_AC068831(AnnotationTestUtils::AC068831_TRANSCRIPT);

}
sub testUpStreamDownStream {
	testUpsteamMilesAway_1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam5001bp_1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam5000bp_1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam2001bp_1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam2000bp_1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam1bp_1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownstream1bp_1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownstream500bp_1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownstream501bp_1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownstream5000bp_1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownstream5001bp_1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownstreamMilesAwaybp_1bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	testUpsteamMilesAway_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam5001bp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam5000bp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam2001bp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam2000bp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam1997bp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam1996bp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEndsUpsteam1bp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownsteam1bp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownsteam496bp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownsteam497bp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownsteam498bp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownsteam499bp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownsteam500bp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownsteam501bp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownsteam5000bp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownsteam5001bp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testStartsDownsteamMilesAwaybp_5bp_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	testUpsteamMilesAway_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEndsUpsteam5001bp_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEndsUpsteam5000bp_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEndsUpsteam2001bp_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEndsUpsteam2000bp_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEndsUpsteam1997bp_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEndsUpsteam1996bp_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEndsUpsteam1bp_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownsteam1bp_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownsteam496bp_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownsteam497bp_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownsteam498bp_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownsteam499bp_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownsteam500bp_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownsteam501bp_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownsteam5000bp_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownsteam5001bp_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testStartsDownsteamMilesAway_5bp_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
}
sub testSplice {
	#5 PRIME UTR SPLICE
	test5PrimeUTREssentialSpliceSite_1bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTREssentialSpliceSite_1bp_3_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTREssentialSpliceSite_1bp_5_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTRSpliceRegion_1bp_6_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTRSpliceRegion_1bp_10_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTRIntronic_1bp_11_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTRIntronic_1bp_n11_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTRSpliceRegion_1bp_n10_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTRSpliceRegion_1bp_n3_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTREssentialSpliceSite_1bp_n2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTREssentialSpliceSite_1bp_n1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	test5PrimeUtrEssentialSpliceSite_1bp_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUtrEssentialSpliceSite_1bp_3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUtrEssentialSpliceSite_1bp_5_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUtrSpliceRegion_1bp_6_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUtrSpliceRegion_1bp_10_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUTRIntronic_1bp_11_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUTRIntronic_1bp_n11_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUtrSpliceRegion_1bp_n10_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUtrSpliceRegion_1bp_n3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUtrEssentialSpliceSite_1bp_n2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUtrEssentialSpliceSite_1bp_n1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

	#3 PRIME UTR SPLICE
	test3PrimeUTREssentialSpliceSite_1bp_1_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
	test3PrimeUTREssentialSpliceSite_1bp_5_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
	test3PrimeUTRSpliceRegion_1bp_6_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
	test3PrimeUTRSpliceRegion_1bp_10_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
	test3PrimeUTRIntronic_1bp_11_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
	test3PrimeUTRIntronic_1bp_n11_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
	test3PrimeUTRSpliceRegion_1bp_n10_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
	test3PrimeUTRSpliceRegion_1bp_n3_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
	test3PrimeUTREssentialSpliceSite_1bp_n2_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
	test3PrimeUTREssentialSpliceSite_1bp_n1_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);

	# NON-CODING GENE SPLICE
	testIntronic_1bp_n11_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
	testSpliceRegion_1bp_n10_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
	testSpliceRegion_1bp_n3_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
	testEssentialSpliceSite_1bp_n2_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
	testEssentialSpliceSite_1bp_n1_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
	testEssentialSpliceSite_1bp_1_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
	testEssentialSpliceSite_1bp_2_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
	testEssentialSpliceSite_1bp_3_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
	testEssentialSpliceSite_1bp_4_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
	testEssentialSpliceSite_1bp_5_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
	testSpliceRegion_1bp_6_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
	testSpliceRegion_1bp_10_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
	testIntronic_1bp_11_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);

	# CDS SPLICE

	testEssentialSpliceSite_1bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEssentialSpliceSite_1bp_3_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEssentialSpliceSite_1bp_5_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testSpliceRegion_1bp_6_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testSpliceRegion_1bp_10_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testIntronic_1bp_11_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testIntronic_1bp_n11_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testSpliceRegion_1bp_n10_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testSpliceRegion_1bp_n3_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEssentialSpliceSite_1bp_n2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEssentialSpliceSite_1bp_n1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	testEssentialSpliceSite_3bp_1t3_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEssentialSpliceSite_3bp_3t5_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEssentialSpliceSite_3bp_4t6_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEssentialSpliceSite_3bp_5t7_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testSpliceRegion_3bp_6t8_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testSpliceRegion_3bp_8t10_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testSpliceRegion_3bp_9t11_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testSpliceRegion_3bp_10t12_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testIntronic_3bp_11t13_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testIntronic_3bp_n13tn11_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testSpliceRegion_3bp_n12tn10_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testSpliceRegion_3bp_n11tn9_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testSpliceRegion_3bp_n10tn8_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testSpliceRegion_3bp_n5tn3_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEssentialSpliceSite_3bp_n4tn2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEssentialSpliceSite_3bp_n3tn1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	testEssentialSpliceSite_1bp_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEssentialSpliceSite_1bp_3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEssentialSpliceSite_1bp_5_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testSpliceRegion_1bp_6_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testSpliceRegion_1bp_10_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testIntronic_1bp_11_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testIntronic_1bp_n11_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testSpliceRegion_1bp_n10_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testSpliceRegion_1bp_n3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEssentialSpliceSite_1bp_n2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEssentialSpliceSite_1bp_n1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

	testEssentialSpliceSite_3bp_1t3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEssentialSpliceSite_3bp_3t5_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEssentialSpliceSite_3bp_4t6_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEssentialSpliceSite_3bp_5t7_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testSpliceRegion_3bp_6t8_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testSpliceRegion_3bp_8t10_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testSpliceRegion_3bp_9t11_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testSpliceRegion_3bp_10t12_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testIntronic_3bp_11t13_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testIntronic_3bp_n13tn11_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testSpliceRegion_3bp_n12tn10_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testSpliceRegion_3bp_n11tn9_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testSpliceRegion_3bp_n10tn8_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testSpliceRegion_3bp_n5tn3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEssentialSpliceSite_3bp_n4tn2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEssentialSpliceSite_3bp_n3tn1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

}
sub testNonCodingExon {
	#5 PRIME UTR EXON
	test5PrimeUTRExon_1bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTRExon_1bp_2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTRExon_1bp_3_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTRExon_1bp_4_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	test5PrimeUTRExon_2bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	test5PrimeUTRExon_2bp_2_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	test5PrimeUtrExon_1bp_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUtrExon_1bp_2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUtrExon_1bp_3_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUtrExon_1bp_4_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

	test5PrimeUtrExon_2bp_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	test5PrimeUtrExon_2bp_2_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

	#3 PRIME UTR EXON

	test3PrimeUTRExon_1bp_1_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
	test3PrimeUTRExon_1bp_2_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
	test3PrimeUTRExon_1bp_3_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);
	test3PrimeUTRExon_1bp_4_MYB(AnnotationTestUtils::MYB_TRANSCRIPT);

	#NON-CODING EXON

	testExon_1bp_1_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
	testExon_1bp_2_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
	testExon_1bp_3_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);
	testExon_1bp_4_AC011503(AnnotationTestUtils::AC011503_TRANSCRIPT);



}
sub testCodingExon {

	# + STRAND

	testInFrameInPhasePreserveLength_3bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testInFrameOutOfPhasePreserveLength_2bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testInFrameOutOfPhasePreserveLength_4bp_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	testInFrameInPhaseShorten_9to3_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testInFrameInPhaseLengthen_3to9_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	testInFrameOutOfPhaseShorten_8to2_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testInFrameOutOfPhaseLengthen_2to8_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	testFrameShiftShorten_9to1_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testFrameShiftLengthen_1to9_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	# - STRAND

	testInFrameInPhasePreserveLength_3bp_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testInFrameOutOfPhasePreserveLength_2bp_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testInFrameOutOfPhasePreserveLength_4bp_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

	testInFrameInPhaseShorten_9to3_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testInFrameInPhaseLengthen_3to9_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

	testInFrameOutOfPhaseShorten_8to2_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testInFrameOutOfPhaseLengthen_2to8_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

	testFrameShiftShorten_9to1_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testFrameShiftLengthen_1to9_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
}
sub testBoundyCrossing {
	# + STRAND
	testIntronToExonBoundry_15to1_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testSpliceRegionToExonBoundry_6to1_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testEssSpliceToExonBoundry_4to1_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	testExonToEssSpliceBoundry_4to1_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testExonToSpliceRegionBoundry_9to1_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testExonToIntronBoundry_20to1_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);

	# - strand
	testIntronToExonBoundry_15to1_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testSpliceRegionToExonBoundry_10to1_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testEssSpliceToExonBoundry_4to1_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

	testExonToEssSpliceBoundry_4to1_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testExonToSpliceRegionBoundry_10to1_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);
	testExonToIntronBoundry_16to1_1_TOR1AIP2(AnnotationTestUtils::TOR1AIP2_TRANSCRIPT);

}
sub testLargerChanges {
	# + strand
	testExonRemovalWithIntron_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testFrontHalfTranscriptRemoval_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);
	testRearHalfTranscriptRemoval_1_CEP350(AnnotationTestUtils::CEP350_TRANSCRIPT);


}


# CEP350 protein coding gene with 5 prime utr exons on + strand of genome

sub testUpsteamMilesAway_1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Upstream Miles Away 1bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' 		=> 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179913873,
			'maxpos'				=> 179913873,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testEndsUpsteam5001bp_1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 5001 1bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179918872,
			'maxpos'				=> 179918872,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testEndsUpsteam5000bp_1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 5000 1bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179918873,
			'maxpos'				=> 179918873,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get5KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam2001bp_1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 2001 1bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179921872,
			'maxpos'				=> 179921872,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get5KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam2000bp_1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 2000 1bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179921873,
			'maxpos'				=> 179921873,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get2KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam1bp_1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 1 1bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179923872,
			'maxpos'				=> 179923872,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get2KBUpStreamVariantClass);

		done_testing();
	};
}

sub test5PrimeUTRIntronic_1bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Intronic 1bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924289,
			'maxpos'				=> 179924289,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic_1bp_2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Intronic 1bp + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179934289,
			'maxpos'				=> 179934289,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic_1bp_3_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Intronic 1bp + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955293,
			'maxpos'				=> 179955293,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic_10bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Intronic 10bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924289,
			'maxpos'				=> 179924298,
			'delseq' 				=> 'ATATATATAT',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic_10bp_2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Intronic 10bp + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955284,
			'maxpos'				=> 179955293,
			'delseq' 				=> 'ATATATATAT',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}

sub testCDSIntronic_1bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 CDS Intronic 1bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961186,
			'maxpos'				=> 179961186,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testCDSIntronic_10bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 CDS Intronic 10bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961177,
			'maxpos'				=> 179961186,
			'delseq' 				=> 'ATATATATAT',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testCDSIntronic_1bp_2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 CDS Intronic 1bp + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961367,
			'maxpos'				=> 179961367,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testCDSIntronic_10bp_2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 CDS Intronic 10bp + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179961367,
			'maxpos'				=> 179961376,
			'delseq' 				=> 'ATATATATAT',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testCDSIntronic_1bp_3_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 CDS Intronic 1bp + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179963367,
			'maxpos'				=> 179963367,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testCDSIntronic_10bp_3_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 CDS Intronic 10bp + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179963367,
			'maxpos'				=> 179963376,
			'delseq' 				=> 'ATATATATAT',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}

sub testStartsDownstream1bp_1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 1 1bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084016,
			'maxpos'				=> 180084016,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get500BPDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownstream500bp_1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 500 1bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084515,
			'maxpos'				=> 180084515,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get500BPDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownstream501bp_1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 501 1bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084516,
			'maxpos'				=> 180084516,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownstream5000bp_1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 5000 1bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180089015,
			'maxpos'				=> 180089015,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownstream5001bp_1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 5001 1bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180089016,
			'maxpos'				=> 180089016,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testStartsDownstreamMilesAwaybp_1bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream MilesAway 1bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 190089015,
			'maxpos'				=> 190089015,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}

sub testUpsteamMilesAway_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Upstream Miles Away 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' 		=> 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179913873,
			'maxpos'				=> 179913877,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testEndsUpsteam5001bp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 5001 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179918868,
			'maxpos'				=> 179918872,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testEndsUpsteam5000bp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 5000 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179918869,
			'maxpos'				=> 179918873,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get5KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam2001bp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 2001 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179921868,
			'maxpos'				=> 179921872,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get5KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam2000bp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 2000 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179921869,
			'maxpos'				=> 179921873,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get2KBUpStreamVariantClass,$a->get5KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam1997bp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 1997 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179921872,
			'maxpos'				=> 179921876,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get2KBUpStreamVariantClass,$a->get5KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam1996bp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 1996 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179921873,
			'maxpos'				=> 179921877,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get2KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam1bp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ends Upstream 1 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179923868,
			'maxpos'				=> 179923872,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get2KBUpStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam1bp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 1 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084016,
			'maxpos'				=> 180084020,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get500BPDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam496bp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 496 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084511,
			'maxpos'				=> 180084515,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get500BPDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam497bp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 497 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084512,
			'maxpos'				=> 180084516,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get500BPDownStreamVariantClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam498bp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 498 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084513,
			'maxpos'				=> 180084517,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get500BPDownStreamVariantClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam499bp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 499 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084514,
			'maxpos'				=> 180084518,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get500BPDownStreamVariantClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam500bp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 500 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084515,
			'maxpos'				=> 180084519,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get500BPDownStreamVariantClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam501bp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 501 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180084516,
			'maxpos'				=> 180084520,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam5000bp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 5000 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180089015,
			'maxpos'				=> 180089019,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam5001bp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 5001 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180089016,
			'maxpos'				=> 180089020,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testStartsDownsteamMilesAwaybp_5bp_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Starts Downstream 5001 5bp long + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 190089016,
			'maxpos'				=> 190089020,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}

sub test5PrimeUTRExon_1bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Exon 1bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924276,
			'maxpos'				=> 179924276,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									404,0,404,0,'A','UUUU','r.404delAinsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTRExon_1bp_2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Exon 1bp + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924277,
			'maxpos'				=> 179924277,
			'delseq' 				=> 'G',
			'insseq'				=> 'AAAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									405,0,405,0,'G','AAAA','r.405delGinsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTREssentialSpliceSite_1bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Essential Splice Site 1bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924278,
			'maxpos'				=> 179924278,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									405,1,405,1,'G','UUUU','r.405+1delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTREssentialSpliceSite_1bp_3_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Essential Splice Site 1bp + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924280,
			'maxpos'				=> 179924280,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									405,3,405,3,'G','UUUU','r.405+3delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTREssentialSpliceSite_1bp_5_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Essential Splice Site 1bp + strand 5' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924282,
			'maxpos'				=> 179924282,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									405,5,405,5,'G','UUUU','r.405+5delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTRSpliceRegion_1bp_6_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Splice Region 1bp + strand 6' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924283,
			'maxpos'				=> 179924283,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									405,6,405,6,'G','UUUU','r.405+6delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTRSpliceRegion_1bp_10_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Splice Region 1bp + strand 10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924287,
			'maxpos'				=> 179924287,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									405,10,405,10,'G','UUUU','r.405+10delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTRIntronic_1bp_11_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Intronic 1bp + strand 11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924288,
			'maxpos'				=> 179924288,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic_1bp_n11_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Intronic 1bp + strand -11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955293,
			'maxpos'				=> 179955293,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRSpliceRegion_1bp_n10_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Splice Region 1bp + strand -10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955294,
			'maxpos'				=> 179955294,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									406,-10,406,-10,'G','UUUU','r.406-10delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTRSpliceRegion_1bp_n3_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Splice Region 1bp + strand -3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955301,
			'maxpos'				=> 179955301,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									406,-3,406,-3,'G','UUUU','r.406-3delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTREssentialSpliceSite_1bp_n2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Essential Splice Site 1bp + strand -2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955302,
			'maxpos'				=> 179955302,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									406,-2,406,-2,'G','UUUU','r.406-2delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTREssentialSpliceSite_1bp_n1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Essential Splice Site 1bp + strand -1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955303,
			'maxpos'				=> 179955303,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									406,-1,406,-1,'G','UUUU','r.406-1delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTRExon_1bp_3_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Exon 1bp + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955304,
			'maxpos'				=> 179955304,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									406,0,406,0,'G','UUUU','r.406delGinsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTRExon_1bp_4_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Exon 1bp + strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955305,
			'maxpos'				=> 179955305,
			'delseq' 				=> 'T',
			'insseq'				=> 'AAAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									407,0,407,0,'U','AAAA','r.407delUinsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get5PrimeUtrVariantClass);

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

sub test5PrimeUTRExon_2bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Exon 2bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179924276,
			'maxpos'				=> 179924277,
			'delseq' 				=> 'AG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									404,0,405,0,'AG','UUUU','r.404_405delAGinsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTRExon_2bp_2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 5 Prime UTR Exon 2bp + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955304,
			'maxpos'				=> 179955305,
			'delseq' 				=> 'GT',
			'insseq'				=> 'AAAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									406,0,407,0,'GU','AAAA','r.406_407delGUinsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get5PrimeUtrVariantClass);

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

sub testEssentialSpliceSite_1bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice Site 1bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179990145,
			'maxpos'				=> 179990145,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3653,1,3653,1,'G','UUUU','r.3653+1delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3235,1,3235,1,'G','TTTT','c.3235+1delGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

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
sub testEssentialSpliceSite_1bp_3_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice Site 1bp + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179990147,
			'maxpos'				=> 179990147,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3653,3,3653,3,'G','UUUU','r.3653+3delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3235,3,3235,3,'G','TTTT','c.3235+3delGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

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
sub testEssentialSpliceSite_1bp_5_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice Site 1bp + strand 5' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179990149,
			'maxpos'				=> 179990149,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3653,5,3653,5,'G','UUUU','r.3653+5delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3235,5,3235,5,'G','TTTT','c.3235+5delGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

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
sub testSpliceRegion_1bp_6_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region 1bp + strand 6' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179990150,
			'maxpos'				=> 179990150,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3653,6,3653,6,'G','UUUU','r.3653+6delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3235,6,3235,6,'G','TTTT','c.3235+6delGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

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
sub testSpliceRegion_1bp_10_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region 1bp + strand 10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179990154,
			'maxpos'				=> 179990154,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3653,10,3653,10,'G','UUUU','r.3653+10delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3235,10,3235,10,'G','TTTT','c.3235+10delGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

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
sub testIntronic_1bp_11_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Intronic 1bp + strand 11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179990155,
			'maxpos'				=> 179990155,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic_1bp_n11_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Intronic 1bp + strand -11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179991822,
			'maxpos'				=> 179991822,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testSpliceRegion_1bp_n10_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region 1bp + strand -10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179991823,
			'maxpos'				=> 179991823,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3654,-10,3654,-10,'G','UUUU','r.3654-10delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3236,-10,3236,-10,'G','TTTT','c.3236-10delGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

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
sub testSpliceRegion_1bp_n3_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region 1bp + strand -3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179991830,
			'maxpos'				=> 179991830,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3654,-3,3654,-3,'G','UUUU','r.3654-3delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3236,-3,3236,-3,'G','TTTT','c.3236-3delGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

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
sub testEssentialSpliceSite_1bp_n2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice Site 1bp + strand -2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179991831,
			'maxpos'				=> 179991831,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3654,-2,3654,-2,'G','UUUU','r.3654-2delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3236,-2,3236,-2,'G','TTTT','c.3236-2delGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

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
sub testEssentialSpliceSite_1bp_n1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice Site 1bp + strand -1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179991832,
			'maxpos'				=> 179991832,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3654,-1,3654,-1,'G','UUUU','r.3654-1delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3236,-1,3236,-1,'G','TTTT','c.3236-1delGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

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

sub testEssentialSpliceSite_3bp_1t3_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice Site 3bp + strand 1 to 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179990145,
			'maxpos'				=> 179990147,
			'delseq' 				=> 'GGG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3653,1,3653,3,'GGG','UUUU','r.3653+1_3653+3delggginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3235,1,3235,3,'GGG','TTTT','c.3235+1_3235+3delGGGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

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
sub testEssentialSpliceSite_3bp_3t5_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice Site 3bp + strand 3 to 5' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179990147,
			'maxpos'				=> 179990149,
			'delseq' 				=> 'GGG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3653,3,3653,5,'GGG','UUUU','r.3653+3_3653+5delggginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3235,3,3235,5,'GGG','TTTT','c.3235+3_3235+5delGGGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

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
sub testEssentialSpliceSite_3bp_4t6_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice Site 3bp + strand 4 to 6' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179990148,
			'maxpos'				=> 179990150,
			'delseq' 				=> 'GGG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getEssentialSpliceSiteClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3653,4,3653,6,'GGG','UUUU','r.3653+4_3653+6delggginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3235,4,3235,6,'GGG','TTTT','c.3235+4_3235+6delGGGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->getSpliceRegionVariantClass);

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
sub testEssentialSpliceSite_3bp_5t7_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice Site 3bp + strand 5 to 7' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179990149,
			'maxpos'				=> 179990151,
			'delseq' 				=> 'GGG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getEssentialSpliceSiteClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3653,5,3653,7,'GGG','UUUU','r.3653+5_3653+7delggginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3235,5,3235,7,'GGG','TTTT','c.3235+5_3235+7delGGGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->getSpliceRegionVariantClass);

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
sub testSpliceRegion_3bp_6t8_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region 3bp + strand 6 to 8' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179990150,
			'maxpos'				=> 179990152,
			'delseq' 				=> 'GGG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3653,6,3653,8,'GGG','UUUU','r.3653+6_3653+8delggginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3235,6,3235,8,'GGG','TTTT','c.3235+6_3235+8delGGGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

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
sub testSpliceRegion_3bp_8t10_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region 3bp + strand 8 to 10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179990152,
			'maxpos'				=> 179990154,
			'delseq' 				=> 'GGG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3653,8,3653,10,'GGG','UUUU','r.3653+8_3653+10delggginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3235,8,3235,10,'GGG','TTTT','c.3235+8_3235+10delGGGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

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
sub testSpliceRegion_3bp_9t11_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region 3bp + strand 9 to 11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179990153,
			'maxpos'				=> 179990155,
			'delseq' 				=> 'GGG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3653,9,3653,11,'GGG','UUUU','r.3653+9_3653+11delggginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getIntronVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3235,9,3235,11,'GGG','TTTT','c.3235+9_3235+11delGGGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getIntronVariantClass);

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
sub testSpliceRegion_3bp_10t12_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region 3bp + strand 10 to 12' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179990154,
			'maxpos'				=> 179990156,
			'delseq' 				=> 'GGG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3653,10,3653,12,'GGG','UUUU','r.3653+10_3653+12delggginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getIntronVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3235,10,3235,12,'GGG','TTTT','c.3235+10_3235+12delGGGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getIntronVariantClass);

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
sub testIntronic_3bp_11t13_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Intronic 3bp + strand 11 to 13' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179990155,
			'maxpos'				=> 179990157,
			'delseq' 				=> 'GGG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic_3bp_n13tn11_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Intronic 3bp + strand -13 to -11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179991820,
			'maxpos'				=> 179991822,
			'delseq' 				=> 'GGG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testSpliceRegion_3bp_n12tn10_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region 3bp + strand -12 to -10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179991821,
			'maxpos'				=> 179991823,
			'delseq' 				=> 'GGG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3654,-12,3654,-10,'GGG','UUUU','r.3654-12_3654-10delggginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getIntronVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3236,-12,3236,-10,'GGG','TTTT','c.3236-12_3236-10delGGGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getIntronVariantClass);

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
sub testSpliceRegion_3bp_n11tn9_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region 3bp + strand -11 to -9' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179991822,
			'maxpos'				=> 179991824,
			'delseq' 				=> 'GGG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3654,-11,3654,-9,'GGG','UUUU','r.3654-11_3654-9delggginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getIntronVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3236,-11,3236,-9,'GGG','TTTT','c.3236-11_3236-9delGGGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getIntronVariantClass);

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
sub testSpliceRegion_3bp_n10tn8_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region 3bp + strand -10 to -8' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179991823,
			'maxpos'				=> 179991825,
			'delseq' 				=> 'GGG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3654,-10,3654,-8,'GGG','UUUU','r.3654-10_3654-8delggginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3236,-10,3236,-8,'GGG','TTTT','c.3236-10_3236-8delGGGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

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
sub testSpliceRegion_3bp_n5tn3_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region 3bp + strand -5 to -3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179991828,
			'maxpos'				=> 179991830,
			'delseq' 				=> 'GGG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3654,-5,3654,-3,'GGG','UUUU','r.3654-5_3654-3delggginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3236,-5,3236,-3,'GGG','TTTT','c.3236-5_3236-3delGGGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

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
sub testEssentialSpliceSite_3bp_n4tn2_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice Site 3bp + strand -4 to -2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179991829,
			'maxpos'				=> 179991831,
			'delseq' 				=> 'GGG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3654,-4,3654,-2,'GGG','UUUU','r.3654-4_3654-2delggginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3236,-4,3236,-2,'GGG','TTTT','c.3236-4_3236-2delGGGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getEssentialSpliceSiteVariantClass);

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
sub testEssentialSpliceSite_3bp_n3tn1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Essential Splice Site 3bp + strand -3 to -1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179991830,
			'maxpos'				=> 179991832,
			'delseq' 				=> 'GGG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3654,-3,3654,-1,'GGG','UUUU','r.3654-3_3654-1delggginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									3236,-3,3236,-1,'GGG','TTTT','c.3236-3_3236-1delGGGinsTTTT',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getEssentialSpliceSiteVariantClass);

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

sub testInFrameInPhasePreserveLength_3bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 In Frame In Phase Preserve Length + strand 3bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955320,
			'maxpos'				=> 179955322,
			'delseq' 				=> 'AGG',
			'insseq'				=> 'CCC',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									422,0,424,0,'AGG','CCC','r.422_424delAGGinsccc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									4,0,6,0,'AGG','CCC','c.4_6delAGGinsCCC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									2,0,2,0,'R','P','p.R2P',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getNonSynonymousVariantClass);
		done_testing();
	};


}
sub testInFrameOutOfPhasePreserveLength_2bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 In Frame Out Of Phase Preserve Length + strand 2bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955320,
			'maxpos'				=> 179955321,
			'delseq' 				=> 'AG',
			'insseq'				=> 'CC',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									422,0,423,0,'AG','CC','r.422_423delAGinscc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									4,0,5,0,'AG','CC','c.4_5delAGinsCC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									2,0,2,0,'R','P','p.R2P',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getNonSynonymousVariantClass);
		done_testing();
	};


}
sub testInFrameOutOfPhasePreserveLength_4bp_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 In Frame Out Of Phase Preserve Length + strand 4bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955320,
			'maxpos'				=> 179955323,
			'delseq' 				=> 'AGGA',
			'insseq'				=> 'CCCC',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									422,0,425,0,'AGGA','CCCC','r.422_425delAGGAinscccc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									4,0,7,0,'AGGA','CCCC','c.4_7delAGGAinsCCCC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									2,0,3,0,'RS','PR','p.R2_S3delinsPR',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass);
		done_testing();
	};


}

sub testInFrameInPhaseShorten_9to3_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 In Frame In Phase Shorten + strand 9bp to 3bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955320,
			'maxpos'				=> 179955328,
			'delseq' 				=> 'AGGAGCAGC',
			'insseq'				=> 'CCC',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									422,0,430,0,'AGGAGCAGC','CCC','r.422_430delAGGAGCAGCinsccc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									4,0,12,0,'AGGAGCAGC','CCC','c.4_12delAGGAGCAGCinsCCC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									2,0,4,0,'RSS','P','p.R2_S4delinsP',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass);
		done_testing();
	};


}
sub testInFrameInPhaseLengthen_3to9_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 In Frame In Phase Lengthen + strand 3bp to 9bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955320,
			'maxpos'				=> 179955322,
			'delseq' 				=> 'AGG',
			'insseq'				=> 'CCCCCCCCC',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									422,0,424,0,'AGG','CCCCCCCCC','r.422_424delAGGinsccccccccc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									4,0,6,0,'AGG','CCCCCCCCC','c.4_6delAGGinsCCCCCCCCC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									2,0,2,0,'R','PPP','p.R2delinsPPP',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass);
		done_testing();
	};


}

sub testInFrameOutOfPhaseShorten_8to2_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 In Frame Out Of Phase Shorten + strand 8bp to 2bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955320,
			'maxpos'				=> 179955327,
			'delseq' 				=> 'AGGAGCAG',
			'insseq'				=> 'CC',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									422,0,429,0,'AGGAGCAG','CC','r.422_429delAGGAGCAGinscc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									4,0,11,0,'AGGAGCAG','CC','c.4_11delAGGAGCAGinsCC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									2,0,4,0,'RSS','P','p.R2_S4delinsP',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass);
		done_testing();
	};


}
sub testInFrameOutOfPhaseLengthen_2to8_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 In Frame Out Of Phase Lengthen + strand 2bp to 8bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955320,
			'maxpos'				=> 179955321,
			'delseq' 				=> 'AG',
			'insseq'				=> 'CCCCCCCC',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									422,0,423,0,'AG','CCCCCCCC','r.422_423delAGinscccccccc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									4,0,5,0,'AG','CCCCCCCC','c.4_5delAGinsCCCCCCCC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									2,0,2,0,'R','PPP','p.R2delinsPPP',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass);
		done_testing();
	};


}

sub testFrameShiftShorten_9to1_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Frame Shift Shorten + strand 9bp to 1bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955320,
			'maxpos'				=> 179955328,
			'delseq' 				=> 'AGGAGCAGC',
			'insseq'				=> 'C',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									422,0,430,0,'AGGAGCAGC','C','r.422_430delAGGAGCAGCinsc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									4,0,12,0,'AGGAGCAGC','C','c.4_12delAGGAGCAGCinsC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									2,0,2,0,'R','QIKRGAFTKSKELSKQGYCSSRYNHIVGCTFSNQGCSETH*','p.R2fs*41',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getFrameShiftVariantClass);
		done_testing();
	};


}
sub testFrameShiftLengthen_1to9_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Frame Shift Lengthen + strand 1bp to 9bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179955320,
			'maxpos'				=> 179955320,
			'delseq' 				=> 'A',
			'insseq'				=> 'CCCCCCCCC',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									422,0,422,0,'A','CCCCCCCCC','r.422delAinsccccccccc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									4,0,4,0,'A','CCCCCCCCC','c.4delAinsCCCCCCCCC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									2,0,2,0,'R','PPPGAANQKRCLYQIQGTLKARILFKQI*','p.R2fs*29',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getFrameShiftVariantClass);
		done_testing();
	};


}

sub testIntronToExonBoundry_15to1_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Intron To Exon Boundry + strand 15bp to 1bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179965675,
			'maxpos'				=> 179965689,
			'delseq' 				=> 'TTTTTTATTGTAGGG',
			'insseq'				=> 'C',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass,$a->getEssentialSpliceSiteClass,$a->getSpliceRegionClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									814,-13,815,0,'UUUUUUAUUGUAGGG','C','r.814-13_815deluuuuuuauuguagGGinsc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									396,-13,397,0,'TTTTTTATTGTAGGG','C','c.396-13_397delTTTTTTATTGTAGGGinsC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

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
sub testSpliceRegionToExonBoundry_6to1_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Splice Region To Exon Boundry + strand 6bp to 1bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179965684,
			'maxpos'				=> 179965689,
			'delseq' 				=> 'GTAGGG',
			'insseq'				=> 'C',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass,$a->getEssentialSpliceSiteClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									814,-4,815,0,'GUAGGG','C','r.814-4_815delguagGGinsc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									396,-4,397,0,'GTAGGG','C','c.396-4_397delGTAGGGinsC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

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
sub testEssSpliceToExonBoundry_4to1_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Ess Splice To Exon Boundry + strand 4bp to 1bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179965686,
			'maxpos'				=> 179965689,
			'delseq' 				=> 'AGGG',
			'insseq'				=> 'C',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									814,-2,815,0,'AGGG','C','r.814-2_815delagGGinsc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									396,-2,397,0,'AGGG','C','c.396-2_397delAGGGinsC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

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

sub testExonToEssSpliceBoundry_4to1_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Exon To Ess Splice Boundry + strand 4bp to 1bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179966309,
			'maxpos'				=> 179966312,
			'delseq' 				=> 'AGGT',
			'insseq'				=> 'C',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1435,0,1436,2,'AGGU','C','r.1435_1436+2delAGguinsc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1017,0,1018,2,'AGGT','C','c.1017_1018+2delAGGTinsC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

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
sub testExonToSpliceRegionBoundry_9to1_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Exon to Splice Region Boundry + strand 9bp to 1bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179966309,
			'maxpos'				=> 179966317,
			'delseq' 				=> 'AGGTTTGTA',
			'insseq'				=> 'C',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass,$a->getEssentialSpliceSiteClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1435,0,1436,7,'AGGUUUGUA','C','r.1435_1436+7delAGguuuguainsc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1017,0,1018,7,'AGGTTTGTA','C','c.1017_1018+7delAGGTTTGTAinsC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

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
sub testExonToIntronBoundry_20to1_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Exon to Intron Boundry + strand 20bp to 1bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179966309,
			'maxpos'				=> 179966328,
			'delseq' 				=> 'AGGTTTGTAATCAGTCTTTA',
			'insseq'				=> 'C',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass,$a->getEssentialSpliceSiteClass,$a->getSpliceRegionClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1435,0,1436,18,'AGGUUUGUAAUCAGUCUUUA','C','r.1435_1436+18del20insc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1017,0,1018,18,'AGGTTTGTAATCAGTCTTTA','C','c.1017_1018+18del20insC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

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

sub testExonRemovalWithIntron_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Exon Removed with Intron + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179972259,
			'maxpos'				=> 179972472,
			'delseq' 				=> 'CTTATAAAGCAACTTTGATAAGAATATTCTGTGTTACTGGTATGTATTAGGTTTCAACCCTTCAGAGACCAAGATTCGAACACCTGATGGGAAAGTGTGGCAGGAGGCTGAGTTTCAAAACATGAGTAGAGAACTGTATCGAGATTTAGCACTTCACTTTGCAGGTGAGAATAGAGCTTTCTTATGTTTTGAGTATTTTGTCATACATGCTATC',
			'insseq'				=> 'G',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass,$a->getEssentialSpliceSiteClass,$a->getSpliceRegionClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1437,-50,1550,50,'CUUAUAAAGCAACUUUGAUAAGAAUAUUCUGUGUUACUGGUAUGUAUUAGGUUUCAACCCUUCAGAGACCAAGAUUCGAACACCUGAUGGGAAAGUGUGGCAGGAGGCUGAGUUUCAAAACAUGAGUAGAGAACUGUAUCGAGAUUUAGCACUUCACUUUGCAGGUGAGAAUAGAGCUUUCUUAUGUUUUGAGUAUUUUGUCAUACAUGCUAUC','G',
									'r.1437-50_1550+50del214insg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1019,-50,1132,50,'CTTATAAAGCAACTTTGATAAGAATATTCTGTGTTACTGGTATGTATTAGGTTTCAACCCTTCAGAGACCAAGATTCGAACACCTGATGGGAAAGTGTGGCAGGAGGCTGAGTTTCAAAACATGAGTAGAGAACTGTATCGAGATTTAGCACTTCACTTTGCAGGTGAGAATAGAGCTTTCTTATGTTTTGAGTATTTTGTCATACATGCTATC','G',
									'c.1019-50_1132+50del214insG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

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
sub testFrontHalfTranscriptRemoval_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Front Half Transcript Removal + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179922873,
			'maxpos'				=> 179961456,
			'delseq' 				=> 'CCATACGATTACCCACATTTCCTGTCCCCTCCAGTGCAGCTTCCTGGAGACACTTCACACGTATGCAGGACTTCTGTACATCCTACTCAGGCACAGCTGTAGGCCCTACAAAGAGGGACCATTTACTGCCTCCTGCTGTCATTGATATTCGACAAACATTTATTGAGCACCTACAATGTCCTAGTTACTGTGCTAGGCCATTTTGGGGTATAGAGAGATAAAACAGATCTTAAAGTCTTCACGACACTTAAAGATGAGAAACAGTCTCTGGGTAAATAACTATCGGAGTTCAGAAGGTGGACAGAAGTGGTGGGCGGCGAGTCTGTGTGCAGTGAGCACCTAGCGGAAGATGGCATTATTTGCAGAAGAGAGGTCCTTTAAGCCAAGTCTTGATAGATGATTTCACCACGCGGAGATTGGGAAACGAGGACAAAGATATGAGAAGATATGTTGGCCACAGAGAGCCTCTCAAAGCACTGTTAAAGCAAAGCCAATACTTAATGGGAGGAAACTCCGTTTTACTTGCGTGGAGTCTTAAGTCATTAGATGCATGCTAAACAGAACTGTAACCCTCGCGAATAAAACGGCCACAGTGACAGCAGTGGGTAACGGCAGAGCTGCAACTACTCTGACAGGGGTGACGACTTATTTGGACGACAGTTGGTGAAACATTGTAAACTCTCAGAAGCTGCCGTTTCCAGACTTTGAGTTCATTCTTGTCGAATCCCCAGCTTGTTTTCTCAGCCCGGGGTGAAATGCGCAGCCAAAAGGCCCAATGAGCGTCAAAGCGAAAGGGCAGGGGTGGTCTCTTCTGCGGATATACCTTGGCTCCTCCCAAGTCCTACGAGTTTAAACTTCGACGCACCAAACATGGCTTCTTCGCGGTTCTGCTCCTTGTTTTATTTCGCTCCCACCTTTCGGCCACACTCAACATGGCCCTAACCGAAGGTTGGCATGCCTGTGGCACTGTCTCGGGACAGTTCTCGCGAGGTTTGCCTGGGCCCGACTGCTCCCGCCCTGCCTCGCGCCGCCGCCACAGCCGCTGCCGCTGGCGCCGCTCCTCCTCCGTGTCAGTTGTTGGGCTGTAATGGCGACTGGGCCGCCCCTGACGAAGTGACTCCCGGGCCGGGGAGCGGGGCCAGACCTGCGCCAGAGAGAACTGCAGGGAGCCGCAGCTCGGGGGGTGGCTTGCCCTGAGGGAGGGGAGGCAGCCTTTCCGCCTTGTCTTCCTTCCCAGCGGACCGGCGGATCCCCGGAGCCGGTGCGAGGAGGGCACCCGGTGCGTCCCCGGAGCGGGGAGGCCAGGCCGGGCAGCCCTGGGGCCGGTCGGGGCGGCGTCACTGCACCCTCCGCCAGGCTCCGCGGGATGCACCGTGGTAGCCGAGGGCGGAGGCGACACTCTCAGGTGAGCTCCTGTAGGACCTGGCTGGGACCGCGGAGTCTCGCTCAGCCTCAACTGTCTCTGTCCGCGGTCCCGGGTCCCCGGTGGCGGCAAGAGAGGCGGGGCCGGGCGCGGAGAGTCTTGGCCCTGCTTGGCGTCCCCCTCGGGGAACCTCCCTGGGTTCTGAGGCTTTTCGGGGACCTAGGTTTCGGGTCCCACCAGTGCCTGCTGAAGACGGAATTTCCTTCCTTACGTGCCAGTGTCCTCCCGCCTCCTCCACCCCTCCTCTAGCCTGCAGAGGTAAGATACCTCTAGCCTTTGCTGCCCGGAGAAGTGTCCCTTCTTCATCTGGCGATCAGGCGTATTTGGGGAGATTCCGGGGACTGGCTCCAGTGGCAGGGATGATCAGTTTCTGTTCCTTTTCAGTGGGGATCACATTCTCAAGTATAGGGAAAACTTTGATCTACTCTAGTCTAGCGGCTGCTGAAGGTAGTTTATTCACAGATAATTTGGGGAAAGATTTAAAAGGCTTTATGTTTTGGAAGTTGTCAGTATTTTCCTCTGCCTTTTAAATGCAGCTCAACATGCTAAGTAATATTTAGCTAATGTCACAAGATATGAATTGGAATTTGCCTGGATTTTAAATATTCAGATAGAAAGTTGCACATAGGCATAGATTTTAATGTTTTAATGCAGTAGCAAATCCTTCATATGGGAAAAAACAGGTTGTCAGAAATAGTGATCAATATCCATGTACTCCACCTGATTGTTTTGTTTTGATGGTGATCCTTGTTTTGAATTTGGGGAAGTCGTTATTGTTGCCATGTTAACGATACCTCCTGCCCTGTATTTGAAGAGTTTTAAAAGAAAGACTATCTTTGGCAGAGCTAGTAACAATCTGAAAAGTTCTGTAGAAACAATGATAAAATACTAAAAGTTGCCCTCCATTCTGAAAATTTCTTCGTTTAATATATTGTTTTTAAACAAGGAAATGCACTGCTAATGACATGGAATTGTGACTTTGAGGGATTGTTAGACCCATGCTGGAAGTTGCTTCAAAGCCATGGGCGTATAATAATTAGTTTCGGTGTGCTAAATTGACAAATTTTATTTAATGTAATGAAATAAATCTGAATCATTTGAACTGTCAGTGCATTAAATTTAAGATACTGTAATCTGAAAATAAACTTCAAGCTAGAATACATTGTATATGTGGATAAATACATGTGTATCTCTGTATCTCTTCATCCTCAGTAGCAAAATTTTGTGTATGTAGAGATGGCAAGTTGTATATGTAAAAAATTTCATACGTATGTGAACATATTTCTTTGATTTTCTGAGATAAATTGGTGTGTGTGTCATGTAAGAGCGTTAATACTTAAACAGCCAACCTTCTGTGGGGCAATCAATTTATGAATTATAAAGTTTTTGGATTTTTTCCTTCTGCTTCTGATTGTCAGTCTTTTGCCCATTTTATAACACTTCTACCATATATGAGCAGAAGAAAGGAGATGAAGTACTGCTATAATCAATTTTATTTGTTAAATTGCTATAATTTATAGTGGGAGCAAAAAAGAAGTCTAAAATAACTGACATCAGTTTTTGGTGTGCTTCTATATACATTATTTTATTTAGTTTTCACAGTGATCTTGTGAAGTAGGTATGGCAGATATCCTTATTTTGCAGTCGTGGAAACTGACACACAAGGAAGATCATTTATGTGAATGCCTCTCATTTCTGCAGATAATGTGGAATTATTTTTATTTTAAACTAAATATTAATATAGGTGCTCTGGCATTTGGTGTTATTTTAAAATATTAAAAAATTATTCTATTTTTGTTCCTTTTTACATGCCACTTGCCCCTTATGCTGATGATCGTCAGTTCTTTAAAAACATTAGTATTTTTATAACGTTTTGTTAATTATTGACAAGTCATGAATTATGATGGATCACTTACTTCTTTATGTATATACTTTTTGCTTGAACAAGATTATTAAAGTGTAGTTCATTTATTATCTTTCACTTTGGATACTCAGTTCTGTCATCTTTATTCTTTCAGTCATTTATACTCTTAGTTATGTAAAACAAACATTCAAGCCACATGATCAGTAACATAACTGAAAAGTAGTAATATCTTGATGATTTTATTGTTGAGACAACTGATGATTTAAAAAGTTAATTCTCTGCATATAATTGTGATTGTTATTTTTGCAACATGAGGCAAATTTTAGTGTTATAGTGTAAAAAAATCATCAGATTGCATAGTATTTTTTTCTGTATTGATTTGTTGTTGTTGTTGTTCATTGTTATCATAAGTGTGATTTTACATTTACTTGTGGTTATTTTGAGTTATTCTTTGCCTTTCTCATCTGGTCTGTAAGTTCCATGATGGCAGAGATCCTTTAATCACTATTATGCACCCAGTGCCTAGCATATGGTAACATAGTCAATAAATACTAATTGAAGAAATTGACTAAAGAATAGTCCCAGAATAATATTTTAAATGCAAAACAAAACAAAAGACTAGTTCTTAGTTCCTAAGCACACCTCAACTCTTGGGAATATATTCCAGAGATATGCTCAGATTTGCTGGGTAGCAAGGGAATGACAGCCTGGGATTAAAAAGCCACTGGGAATAAGCCACAGTGGAAACTTCATTTACTATTTTAATGATTTGCCTAGGTGCTCTGGCAACAGGTTACGGTAATTTAAAGAAAATATAAACTTTAAGACCCTAACCTGTTATACTAGTACTATAATTATAGCCTTGTAATAACTTTAAGTCACTTCTTAATCTCAGTATTGGGAGTTTAAATTCTTCAAGGCCAGCTAGTCCAATTCCCATATTTTATACATAAGAAAACAGAAGGTTGTATAGGGGCGATATCTCATGCATAGTTGATTTCATAACTGGGAGAAGAATCCATTGCACAGCTGCTAAGTCTACATTTGTTTAATATGTCAGATAATTATAAAATCAGAACATAATTGTGATCAACTTGGCAGTCTATTCAGTATTAATTAAAAACTCAAAATGGGCCTTTTTTTTCATGAATCCAAAGTGTTAACTTAAAAATGAATCTGAACCCTTTCCTCTGTAGTTGGGACTATAAATTGATACCATTGGAGCAACCTTCTGGAAGGAAATTTGGCATTCCCTATCAAGATTAAAAGTATATATATTGTTTGATAGAATCATTTCCTGCTTAAAATTTATGTTACATAAATTTTAACAGCATAGGAATTGTATAAGAATGTTTATATTAGCATTATTTCAAGATGTCAATAACATTTATAAATAATTTAAAAATGTCCTAAATGTCTGTTGGTAGGCTAATAGTTAAAAGTCTAGTATTTTTATGCTGTGAATTACTGTGCAGTTGTTAACAGGTTTCAGGGTATCTATATTTACTTCATGAGCAGCAGCTGAGCATATCCTTCAGCTCCTGGAACCATGCGGTGGACCTTCTGGATTCAAAACCTAACTGCACTACTTTCTAGCTGTGACTCTTGGCGAATTACTTAAAGTCAACTATAAATTGGATATACTAGTAGTTCTCATGTTACAGGGGTACCATGAGAATTAAATGAGTTAACAGTTAGGCCCTAAATAAGCATTAATAAGCGCCACATAAGTGTTTGCTGTTACAATTATGTATTGACTTGGAAGAATGTTCATGATATAAAATTAAGTGAAAAAATAAGTTGCACTATACTATAGAATATGATTTATTTAAAAAACACATAGACAAATGTAGTTTTGTTTTCATTTTGCATGTGAATAGAAAAAGGTCTATAAGGATAGACATTAAACATGTTATTCCTGAGGGAGTGGATGACTTTTCTTTTTTACTTTAGACATTTTTGTATGGTTTTATTGTAAGAATCACTTTTATAATTTTTGAAACTTTTTTCTAAAAAAAGGAAGTAAATTTTTTCTTCTATTAAATAGTGATATAAATGTATTTTTGAATGTGATTTTTAGAGATCTACACATTTCATCAACAGAGTGTCTCATAATTAGATCCTGTGTCAAAACATTTGTTGATGGAAGGGCTGGGACTTGGTGGCTCACGCCTGTAATCCCAGCTCTTAGAGAGGCCGAGGCAAGTGGGTCACTTGAGGTCAGGAGTTCGAGACCGGCCTGGCCAACAAGGCGAAACCCATCTCTACTAAAAATACGAAAACTAGCAAGGCGTGGTGGCAGGCCTGTAATCCCAGCTACTGGGGAGGCCAAGGCAGGAGAATTGCTTGAACCAGGGAGGTGGAGGTTGCAGTGAGCCAAGACTGCCACTGCACTCCAGCCTGGGCAACAGAGCAAGACTCCATCTCAAAACAAACAAAACTTTTGTTGGCTGGGTGCGGTGTCTCATGCCTATAATCCCAGCACTTTGGGAGGCTGAGGCGGGCAGATCACCTGAGATTAGGAGCTTCAGACCAGCCTGGCCAACATAAGGAAACCCCATCTCTACTAAAAGTAAAAAACTAGCTGGGCATGGTGGCAGGCTCCTATAATCCCAGCTGCTGGGGAAGCTGAGGCAGGAGAATTGCTTGAACCTGAGTGGCAGAGGTTGCAGTGAGCCAAGATTGCACCACTGCACTCCAGCCTAGGTGACAGAGGAAGACTCTGGCTCAAAACAGACAAACAAACAAAAAAACTTTTGTTGGCAGAAGGTTGAAATTTAGGAAAATCTATAATAAATATATTTGTTTTTCTCCTAGCTTTGTAGCCTTAAGTCTTATTTCTTACTTTACAAAGTAAACCCGTAGGCTAGATCAGTGTTTCTTAAATCCAGTATTCCTAAACTGTGCACTGAAGTAAATCCTCACTTAATGTTGTCCATAAGTACTTGGAAACTGCAACTTCTAAGACCCCAAACACCTCTAATATTAAACATTGAAATAAATGGGTCCTATGCATACATTTAAGAAATACTAAAAACAAGATAATTACCCACTTATTTCAGTTCAGGATGGCAGGTGGCTGGAGCCTGTCCCGGCAAGTAGGGCACAAGATAGGCTGCCATTCCATTGCAGGTAGCATTTACAGACACACACACACCCACCCACCCACCCACATTCACTCAGACTGGGACCATATAGACATGCCAATTCATGTCAGTTCACCTAATGCGTGCACCTCTCTGATGTGGGAGGAAACTGGAGCATACGGAGAAAAGCCACACAGACATGGGGAGAACATGCAGGCTCCACATAGATGGTAGCCCTGGTGGGGAAATCAAAATTTTTTTCTCATCAAAATTGTGAAGAAATGATATTGAATGAAACTATGTTATGGAGGACCTGCCGTACTGTGGAGCACCATAGTAACTCACCTGGACACTATGCATTATTTGAGAAATTTTAATGGAAACAGCAATATTGGATAGCTGTGGGATATTTGACATCTATAGGACACTACTCAAACTACCAGTTTGAGATACCACACAAGAGATTGTGTTCCATTTTTTTTGATGACACTGTGTCTTCCAAACGCTAGGTTTTGGGAGTTTGCTGTGATAAAAAGCAAGTACTATGCAAAAATCAATTTGGAACAGTAAATGAAGGTGGTGGTGTCCAAAATGATTTCAAGATTTGAGAACTTGTGTAATATCTAGTAGGTACACACATCCTATTAAGTAGTTGTGGTAATTTAGTAATAAAAATATTTATTTCAATTTATACATATTTATTTTTCGATCTGCTACTAAGTTGGTAGGCTATAAATACTTTTTAGGTTGTGGACTGAACTACTTAGTAAATGGAACTGTTAGGTATTTTTTTGGCCTAGGGGTTCCATGAAAAAAATAACTGGCATACTGAGGGGTCTGTGAACTGAGAAAGTTTAGGAACCTCTTGAACTTTTGACCAGGACACCTATATATATTTTCCATGTTCTGTGCTCCCAAAGGTTGATTTTTTACATTTATCTTTATTTTTTATTTTTTTAAATAAAAAGAATAGATGTGTGTTATACGTGTATGTGTATACATATAATGAATCAAAAGTTTCACTTTTCAAAGTTTTACTTTAAACAGCCTCATGGCCAGGCGCGGTGGCTCATGCCTATAATGCCAGCACTTTGGAAGGCCGAGGTAGGCAGATCATTTGAGGTCAGGAGTTCAAGACCAACCTGGCCAACGTGGTGAAACACCATCTCTACTAAAAATACAAAAAAATTAGCTCGGCGTGGTGGCGCGCACCTGTAATCCCAGTTACTCAGTAGGCTAAGGCAGGAGAATCACTTGAACCTGGGAGGTGTAGGTTGCAGTGAGCCGAGATCGCGTCACTGCACTCCAGCCAGGGCGACAGAGAGACTGTGTCTCAAAAATAAAATAAAATAAAATAAAGCGCTCACTACCTTCAATACACTTGATATTTGCAAATCAATCCTATTTAGTTCCATTAAAAAAATGCTAGCTCCAATCCACGATATTGTTTTGTTTTGTTTTAATCTACTAACGGGTCATGACTACACACTTGATAATATGATGCTGGTTGACTCCTAACGTCCAGCTTTAATACGGTTTGTTTCTCTGGTTAGTTATAAATGAAGTTAGACTACGCTTTAGCTTTTATTTTAGATTTGGGTTGTGGGGGATACACGTGCAGGTTCGTTACATGGGTATATTGCATCATGCTGAGGTTTGGGCTTCTGTTGAACCTATCACTGAAATAGCGAACATAGTACCCAGTAGGTTGTTTTTCAACCTTTACCCCTTACCTGTCACTTCCCACCTTTTTTTTTGAGATGGAGTCTCTGTCACCCAGGCTGGAGTGCAGTGGCTTGATCTAGGCTCATTGCAACTTCTGCCTCCCAGGTTTAAGCGATTCTCCTGCCTCAGCCTCCCGAGTAGCTTGGACTACAGGTACGAGCCACCATGCCCGGCTAATTTTTTGTATTTTTAGGAGAGATGGGGTTTCACCATGTTGGCCAGGCTGGTCTCAAACTTGTGAGCTCAAGTGGTCTGCCCACCTTGGCCTCCCAAGGTACTGGTATTACAGGCGTGAGCCACCACGCCCGGCCCTTTCCCACTTTTGGAATCCACAGTATCTTATTTTCATCTTTATGTCTATATTTACCTATTGTTTAGACCCCACTTACAAGTGAGAGAAGGTGTGATATTTGATTTTCTGTTTCTGCATTAATTTGCTTAGGATAATGACCTCCAGCTGCATCCATGTTGCTCCAAAGGACATGATTTCATTCATTTTTATGGCTACGTAGTATTCCACGGTGTATATGTACCACATTTTCTTTATCTAATCTACTGTTTTTTGTTTGTTTGTTTGTTTTGTTTTGTTTTTGAGGCAGAGTCTCGCTCTGTTGCTCACGCTGGAGCGCAGTGACGTGATCTTGGCTCACTGCAACTTCCGCCTCCCGGGTCCAAGCAATTCTCCTGCCTCAGCCTCCTGAGTAGCTGGGATTACAGGCACATGCCACCATGCCTGGCTGATTTTTGTATTTTTAGTAGAGATGGAGTTTCACTATGTTGGTCAGGCTGGTCTCAAACTCCTGACCTCATGATCTGCCCACCTCGGCCTCCCAAAGTGCTGGGATTACAGGCATGAGCCACTGCTCCCGGCCTATCTGATCTACTGTTGATAGGTACTTAGTTTGATTCAATGACTTTGCTGTTGTGAATAGTGCTGTGATAAACATACTAGTGCAGTGTCTTTTTGATAGAATGATTTCTTTTCCTTTGGATAGATACCCAGTAGTGGGAATGCTGGGTCTAATGATAGTTCTAAGTTCTTTGAGAAATCCCTCTACTCTTTTCCACAGGGGTTGAACTAATTTACATTCCCACCAGCAGTGTATGATACTCCCTTTTCTCTGCATCCTTGTCATCTGTTATTTTTTAACTCTTTAATAATAGCAATTCTGACTGGTGTGAGATGGTATCAGTGTGGTTTTAATTTGTATTTTATTGATGATTAGTGATTCTGAGCATTTTTTTCATATGTTTGTTGGCCACTTGTATGTCTTCTTTTGAGAAGTATCCGTTTATGTCCATTGCTCACTTTTTAATGGAGTTATTTTTTTTTTCTTGTTCATTTGTTTGAGATTCCTTGTAGATTCTGGATCCTAGTCCTTTGTCAGATGTACAGTTTGCAAATATTTTTTTCCCATTCTTTAGGTTACCTCTTTACTTTGTTGATTCTTTTGCTATGTAGAAGCTTTTTAGTTTAATTAAATCTTGTTTATGAATTTTTGTTTTGTTGCATTTGCCTTTGAGGTTTTAGTCATAAATTCTTTGCCTAGGCCCATGTGCAGAAGAGTTTTTCCTAGGTTTTTTAAAATAGGATCTTTATAGTTTGAGGTTTTATGTTTAAGTCTTTAATACGTCTTGAGTTAATTTTTGTATGTGGTGAAAGGTAGGGGGTCCAGTTTCATTCTTCTGCACATGCTAGCCAGTTTTTTTGTGTGTTTTTTGGGTTTTTTTTTCTTTTTAGCACTGTTTATTGAATAGGTGTCCTTTCCCCATTTTTGCTTTTGTCAGCTTTGTTCAAGATCAACTGGTTATAGGTGTATGGCTTTATTTCTGGGTTCTCTATTCTGTTACATTGATACAGATGTCTGTTTTTGTATTGGTAACATGCTGTTTTGGTTAGGATATCCTTGTAGTATAGTTTGAAGTTGGGTAATGTGATGCCTCCTGCTGTGTTCTTTTTGCTTGGCTATTCAGGCTCTTTTTTGGTTCCATATGAATTTTAGAATTGTTTTTTCTAATTTGTGAAAAATGATGTTGCTAATTTGATGGGAATAGAATTTTTGGATTGTTTTGGGCAGTATGGTCATTTAAATTTTGATTCTTCCAATCCATGAGCATGGAATGTTTTTCCATTTGTTTGTGTCATCTGTTATTTCCTTCATCAGTGTTTTGTAGTTCTCCTTGTAGAAATCTTTCACCTCCCTGGTTAAATTTATTCCTGGGTATTTTATTTTTCATGTGGCTATTATAGATGGGATTGCATTTTTGATTTGGTTCTCAGCTTGAACATTATTGGCATATAGAAATGCCACTGACTGTTGTAAGTTGATTTTGCATCCTGAAACTTTACTGAAGTCTTTTTCTTAGGTCTAGGAGTTTTTTGAAGGACTCCTTAGGGTTGAACAGAGATAATTTGAGTTTCTCTTTTCCTATATGGATGCCTTTTTTTAAATTTTAATTTTCTGCCTGATTCCTATGGGTAGGATGGAAAGTATTATGTTGAATAGGAGTCGTAAGAGTGCACATCCTTGTCTTTTTTCAGTTCATAGAGGGGATGCTTCTAACTTTTGTCCATTCAGTGTGATCTTGATTTTTTGGAGTAGTTTCAGTAAGATTGATACCAGCTCTTCTTTGTTTGTCTGGTATAATTTGGCTGTGAATGCATCTGATCCTGGTCTTTTTTTGGTTTGTAGATTTGTTTTTATAACTGATTCAATTTTGTAAATTGTTATGATCTGTTCAGGATTTCAGTTTCTTCCTTGTTCAATCTTGGGAGGTTGTGTGTTTCCAAGAATTTATCCATTTCCTCTGGTTTTTCTAGTTTGTGCACATAAAGATGTTCATAGTAGTCTCTGAGGCTCTTTTGTATTTTTGTAGTATGAGTTCTGTCACCTTTGTCATTTCTGATTATGCTTCTTAGGATCTTCTCACTTTTTTTTTTTTTAATTTACTTAGAGGTTTATCTATCTTGCTTATCCCTTCAAAGAACTAACTTTGTTCTGTTGATCCTTTGTATGTTTTTTTTGTTGTTGTTTGTTTGTTTGTTTGTTTTGGTCTCAATTAATTTAGTTCTGCTCTGATCTTAGTTACTTCTTTTTTTTCTGCTAGCTTTGGGTTTAGTTTGTTCTTTTTCTGGTTCCTTTAAGTTCAACTTTAAGTTGTCAATTTGATATCTTACTGTCTTTTTGATGTAGGCATTTAGTGCTGTAAACTTTCTCCTTAACATTGCTTTTGCTGTATCCCAGAGGTTTTGGTATGTTTGTTTCTCTATTTTCATTTGTTTCAAAGAATTTTTTTATTTCTGCCTTAATTTGTTGTTTAACCAAGTTAGTTAGGAGCAAGTTGTTTAGTTTCCATGTATTTGTGGTTTTGAGAGTTCCTTTTGGTATTGATTTCTAATTTTATTCTACTGTAGTCTGATAATATGCTTGACATGATTTCATTTTTTAAGAATGTATTGAGACTTGCTTTATGACTGAGCTTGTGGTCCATCTTAGAGAATGTTCCATGTGCAGATGAGAGTGTATATTCTGTGATTTTTGGGTACAGTATTCTGTAGAAGTCTGTTAGGTCCGTTTGGTCAAGTGTCTAAGTCCAGAGTTTCTTTGTTAGTTTTATGCTTTGACGTTCTGTTTCTAATGTATATAGGCACATACATTTACAAATGGAAATCTCCTTTATAGATGTAAATTTCTTTTACAGAAGGATTTCAGAATAGCCAGCTAAATGCCAGAAAGGTTTATTTTAATTGGACAAGTGATCTTTTTTCTTTTTCTTTTCTTTTTTTTTTTTTTTTTTTTTTGAGATAGAGTCTTACTCTTGTCACCTAGGCTGGAGTGCAGTTGTGCAATCTCGGCTTACTGCAACCTCCGCCTCCTGGGTTCAGGCGATACTCCTGTCTCAGCCTATTGAATAGCTGGGATTACAGGTGCCCGCCACCACACCTGGCTAATTTTTTTTTTGTACTTTTAGTAGAGACGGGGTTTCACCATGTTGGCCAGGCTGGTCTTGAACTCCTGACCTCAGGTGATCCACTTTCCTTGGCCTTCCAAAGTGCTGGGATTACAGGCATGAGCCACTGTGCTCGGCCTACAGGTGATCTTTTTAACCTAGCTACTATTTCTTAGCTGAAATTACTAAGTTCAGGATGGAGCCCATTAAGGAATAGGGCAAAGAAAGCCTTCTCTGTACCTGGACTCAGCAAAAATAGATCTGAAAAAGGAGGAAGCTTACTTTACCTGGCAACCTACCTTTTATAAACGCTATCTAGGATAACTTTCTTTTCACCTTCAGGGAAGAGTAGTAACGGAGCTGAAAGATTAGCAGATTTAATTTTTCTTAACAATTAGTTGCTTAAGCTTTTTATTTGCCTTCTATAAAGAGTCTTTTTATAAAAAGCAATAAAAATATTGAAATCTTTTTAGAAGCTTCTGCACACCAACAGGCATCCCTAGATGAGACTGATTTGGGAGCCCTCATTTTCAAATGTACTTCTTAAAGTGCAGTGTTCATTTGGAATGTACCGTTGTGATTTTAAATTACCTTTAGTAAGATTTCATCATTTTTGTAAGCATTTGCTGCTTCTGAGGCCTAATAATAATGCATGTATAAGCTGGAAGGTGTAGTACTCAGGTCTTCAGAAATTAAGGATCCCATTTTTACCTTGAATAGCAGCTTTGGTCCTTGGATCCCTTTGATCAACTTAGCCAATGATTTTTCCCTCCCTAAGCATGCAAGAAAAATAAACAAAGGGGGCTAGAACACAAAAATCTATGCAAATTTCTGAAAGCTGAAATTGGTACCCCCTTCAGTATTATCATTTATACTAATTTCTGTCTAAACCAGCCAGACATATGAGGCCTCTAACTAGATCCAAGACAGTTAGATACAGTTGGATCCTGGACCCAGTCCAATTTCTGTTGTAACTTCCAAACCTAGTTTGGATCAGAAATTTGCTCAAACTCTGATAGCTCAAATGCACACTGATGGAGCTTCAGAATCTGAGAGAGAACTTACCCACAATCTCCAGTTGCTTTGAGAGAGCAGTGGACAAAATGAGTCCAGTGGGTACCTTGCTTTGTCACTTGGTGTTCCTGGGAGTTGCTAGAAGCTTTAGTTTGGATCCCACTTCTGACACCATCTATTAAAATTAAAAACCTTAGACAAATTAAATTTAACAGAATTTGATTGAGCAAAGAAGAGTTCATAAATCAGCCAGCTTCCAAACCAGAATAGGTTCAGAGAGGCTTCTTAGACTATACTTTATACAATAGATCATCTCTGAAAATGTGTACAGTTGCTCTATTTTTGTATGATAAGTTATTAAATGTTATAGAGCAATAATGATGTAAAGAAGACCTTTGTACAAATAGAAGGATATTTTTCCAAAGAAAATCTGATTTATCTTTGATTAGTTTGGGTTTTGGTGATTTTTTTAAAAAGCCATTTATAGTTTTTGTTGTTGTTAGTGGGCTTGTAGAAGTTGTTTCTGGAATTGTGGCACCCTCCTTTTTAAAAACTATATCTAAAAATCAGTGATAAGTAATAGTAATAGCACTGGCTATAACTTAATGTTTCATTTCTTTTCTTTTTTTGAGACTGAGTCTCACTCTGTCGCCCAGGAGGGAGTACAGTAGCGTGATCTTGGCTCACTGCAACCTCCGCCTTCCGGGTTCAAGCAGTTCTCCTGTCTTAGCCTCCCAAGTAGCTGGGATTACAGACATGCACCACCACGTCCAGCCAATTTTTGTATTTTTAATAGAGATGGGGTTTCACCATGTTGGCCAGGCTGGTCTCCAACTCCTGACCCCAGGTGATCGCCCACCTCGGCCTCCCAAAGTGCTGGTATTACAGGTGTAAGCCACCGGGCCCAGCTGGCTTCATTTCATTCATTATTTCATGTCTAAATTTTAAAAGTGAAAATATTTTATGAAAGAATTTATTTTTAATTTTGAGATTGTTTTAAAATATAAAATATGTATGCACACGTATACTTATTTTTTTACGTGTTACCGTGTAAGATATTTTCTGCTTAGGAATGTCTAACCTTTTCATGTTATATGAAATACCTGGCTTCATATGTGTATTTGTGCCTTTTTGCATTTAAGAATTGGATTCTGAGACTTATCAATGGTACTGTCATTTGGCAAACAGGTCTGGAGCACTCAATGTGGCTATATACTGTGCGTGAACTAGGGATATAAAGAGGAATAAGAGCAGGTGCTGTGGCTCAGGCCTATAGTCCCAGCACTCTGGGAGGCCAAGGCGGGTGCATCACCTGAGGTCAGGAGATCGAGACCAGCCTGGCCAACATGGCGAAACCCCATCTCTACTAAAAATATAAAAATTAGCTGGGCATGATGGTATGTGCCTGTAATCCCAGCTACTAGGGAGGCTGAGGCAAGAGAATCGCTTGATCCCGGGAGGTAGAGGTTGCAGTGAGCCGAGATTGCGCCACTGCACTCCCGCCTGGATGACAAGAGCAAGACTCTGTCTCCAAAAAAAAAAAAAAAAAGAGGAATAAAATATCTCTGCATCTCTGCAGTTTTATGTAAGAGACAGATATGTAAATAATAGTACAGAGGGGTAAAGAGAATTGATAACATTCTGTGGATGGATTCATAAAGAGGAAGCAAGTATACTGCATTAAGGAGTCATTAAAGGCTTCCAATGAAAGAAAAATTGTTTCCATTTCCTTAAACAGCCTCTCAGTTGACTAAGTTATTGAACGTGGAAGGAGACTTTTGAAAATTATCTTGCTTTCCACGCTACCCAGAGAGGGGTGCATACAGTGTTGTTCTGGATTCCCATTGTAACTTAAAGGGAAACTTTTACAATATCCAGAGGCCTTGATGTCCTTAAGTTCCTTGCAGGAGGAACCCACTTAGATGGCACCAGCCTTGACTTCCAAATGGAACAGTACATCTATAAAAGGAAAAGTAATGGCATCTACATCATAAATCTGAAGAGGACCTGAAAGAAGCTTCTGCTGGCAGCTCGTGCCATTGTTACCTTTGAAAACCCCGCTGATGTCAGTGTCATATTCTGCAGAAATACTGGCCAGAAGGCTGTGTTGAAGTTTGCTGCTGCCACTGGAGCCACTCCAATTGCTGACTGCTTTACTCCTGGAACCTTCACTAACCAGATCCAGGCAGCCTTCCGGGAGCCATGGCTTGTGGTGGTTTGACCACCAGCCTCAAACTGAAGCATCTTAGGTTAACCTACCTACCATTGCTCTGTGTAACACAGAGTCTCCTTTGCGCTGTGTGGACATTGCCATCCCGTGCAACAATAAGGAAGCTCACTCAGTGGGTTTGATATGGTGGGTGCTGGCCGGGGAAGTTCTGTGTGTGTGTGCCACCATTTCCTGTGAACACCCGTGGGAAGTCATGCCTGATCTCTACTTCTACAGAGATCCTGAAGAGACTGAAAAAGCAGAGCAGGCTGCTGCTGAAAACGCTGTGACCAAGGAAGAATTTTAGGGTGAATGGACTGTGTCAGCTCTGGAGTTTAACACTACTCAGCCTGAGGTTGCAGACTGGTCTGAGGGCGTGCGGGTGCCCTCTGTGCCTATTCAGCAGTTTCCTACTAAAGACTGGAGTGCTCAGAAGACTGGTCTACAGCTCCCACTGCTCAGGCCACTGAATGGGTAGGAATAAGCACTGAATTGGCCTTAAGCTGTTCTTGCATGGACTCTTAAGCAACATGGAAACAGGGTTGATGGAAAATAAGCAGCAGTTTCTAAAAAAGACAAAAACAAATTATTTCTTTTTTCCTTGAAGTTTTCCTCAAAAAGTGATCAGTGTGTCTGCTTTAAAAAAAGAATTCAGAAACAGAAGCTGGGTGTGGTGGTGTGTGCTGTAATCCCAGCTACTTGGGAGGCTGAAGTGGGAGGGTTGCTTGAGCCCAGGAATTCAGCAAGACTAGCCTGGGCAACAATGAGACCCTGTCTTCTAAAAAATAAAAAAAATATACACTGAAAAATCTACTTCCTATCCTGGACCCTCACTTCCTCTCCTTGGAGTTAACTAATAGTAACACTTTCTTATGTATAAGGGAAAATTTTATAGTGAAACTTTTATACTTTTAGGCAAAAAATAATCCTGACTTGCAAACTGTAAGACACCTTTTTGTAGATATATTCCTTGAAATTAATTTAGAAAAGTACTCAGGATAACTAGATAGTGAGTAAGCAAAAGTGTAGGCCTATTTTTAGGGCTCCCCCCTTCTTTTGGTTATTTCTGGTCCTGTTTGCTGCCCTCAGCTCTCCTAAAATAAAAAACTATGGGTGTCTTGTTTACGAAGTTTAGCTCTTTATTTCATCAGTGTGGGACATTGTGTTCTCTTTTATCCTTACTCTACTTCCTATGCTGCAAAGGCCAAAATAGTACTAATGTTACATAGTATCTTCTGGATTTCCTGGATTACAGTAAGTCACTCAGTTTTTGTAAAAGTAGTTTTTGAGGTTTATATTTTATAAGTAGAAACAGGTGAGACTAAATTCTCAGATTAGTGGTATATTAGAATATGTTATAGGACTGTTGGTTATTCTGGAATTTACTTGGAACAGTGTAATAGATCTGAGTTATAACTGAGTATTTTCCACTTTTTTGTCTGCTGGTCCTTACTAAATTTACCATTTGTCCATTTTATGTTGAGGATAGATAGTAGGAATATTGCAGAGTTGTCTTTTCTTAACCTTTCAGTTAACTTTTTCTTAACCTATCATAGATCTGATATAGATCAGAGAACTGAAATGTAGCAAAGGTGTGGAAAATATTCTAAGAGCCAGAATGAAGTCATTTTTTTCATCTTTTGTCCCTTAATGCCTGGCACTTAGTAGGAGGCCAATAAATATTGAGCAAAATGAAACTATGTGATTGGATGAAAACACTGTCAATACATGTGAAAGAAAAATTTGCTGTTAGAAGATAATTGAGAATTGGCTATTTCACAATGTCATCTTATATCTTACTTGTGTATACACACCACAGGTATAGATCATTGTTTCACAGATGTTTTGGTATAAGACTGCCTCATGCAGTTGACTTATGGAAGACTCAAAAAGATTTTGTTCATGTGGGATTTATCTGTTGATGTTTACTGTATTTGAAATTTTCAGTGTACTTGTTAGTTTATTTAAAATAATGATAAACCTATTACATGTTAAAATAATAACTTTTTGTTGTTGTTGTTTGTTTTTTTTTTTTTTGAGACAGAGTCTTGCTCTGTCGCTCAGGCTGGAGTGCAGCAGTGCAATCTCGGCTCATTGCAACCTCTGCCTCCCGGGTTTAAGCAATTCTCCTGCTTTAGCCTCCCAAGTAGCTGGGACTACAGGCACATGCCACCACACCCAGCTAATTTTTGTATTTTTAGTAGAGACTGGGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACTGCAGGTGATCTGCCTGCCTCGACCTCCCAAAGTGCTGAGATTATAGGCGTGAGCCACTGCGCCCTGCCTTGTTTTGTTTTTAAGACAGGGTCTCACCCCATTGCCCAGGCTGGAGTGTAGTAGTGCAATTATAGCTCACTGCAGTCTTGACGTTTAAGGCTCAAGTGATCCTCCAACCTCAATCTCCTGAGTAGTGGGGACTACAGGTGCATGCCACCACACCTGGCTAATTAAAAATTTTTTTTTGTGGACACTGGGTCTCACTTTGTTGTCCAGGCTGGTCTTGAACTCTTGGACTCAAGCAATCTTCCTATCTTGGCCTCCCAAAGCGCTGGGATTACAGGTGTGACCCACTGCACCTGGGCTTAACATTTTAAAATGGAAAATAACTAACTGTATCTCAACAATAAAAAGCAGTGAGTGGCATTGTTTGACATTTTTCAACTCTTTTTAGTGTCTGGTTTAATAGAAGACATCTGGATTCTCATATCTGCTTCTGTTTGCATTCAGTCTGCTGTGATATCATATGTCATGTGGCCTTTGGAAAACTCCATTGTACATTCATAAGAGAATGAGAATAAAAAAAGACAAAGTCTTAATATTATTATGGAAATAGTTTTGATTTCAAAGATGCTTCCTGAAAGAGTGTGAGAGACCCTCAAAGATCCCTGAACCATACTTTGAGAACTGGTAGTATAGATTATAACTTCAGAGTAGAAAAGAACCTTACAGGTCATATGGTCATCTCTTCTTGAGGTCACCTTTTCTTGTAGTTGTTTCTGGTTTATGCAGCTGGATGACTGAGATAAGGAATACCTCGTTTGTTTTTTGTGCGATAATGAGCTGGCTCTTGGAGATGGTGAGTTTAATGTGCCTTTCAGATTTCTAAATAAGCAGTTTGATATGTGGGTTAGGATATTAGAGGTAAAATCTGATGTTTTGGTAGTCATTGGCATGTAGGTGGTTTTGAAAGCCTTTTGCAAGAGTGAGGTCGCCCAGCGAGAGAGTATAGTGTTAGAAGCAAAAAAGCTGAAGGTCAAATCCTCCTGAACTTCAAGGTTTAAAAGGCATATGGAGAAAGACTAGCTAGTAAAGGAGACTGAGAAAGAGAACTAGAGAAGTGGGAAGATGTATGTATTGTATGTATGGAGTGTTGTGTATAGAAGTCAAGGGAAAAGAGTATGTATGGAAAATGTGTGTATGGAGTGTAGTGTATAGAAGTCAAGGGAAAAGAGTGCTGTTTATTTTATTTTTTATTTTTTGAGAGGGAGTCTTGCTCTGTCACCCAGGCTGGAGTGCAGTGGTGCAATATCAGCTCACTGCAACCTCCGACTCCCTGGTTCAAGCAATTCTCCTGCCTCAGCCTCCTGAGTAGCTGGGATTACAGGCACACGCCACCACACCCTGCCGAAAAGAGTGTTTTTAGAGAGTGGTCAATTGTAATGAAGCTGCCAGGAGGTCAAGTAGGTTGAGGATTGAAAATACCCATTAGATTTGGAGATGAAGAGGTTTTTGGTGATTTTAATGAGAGTTATTCAGAAATGGGGATGGCAACCACATTGGAGTTGCTTGAAGAGTCAGAGTCTAAATGGACACTGAGTGTATTTAACTTTTTTTTTTTTTTTTGAAACAGAGTCTCATTCTGTCGCCAGGCTGGAGTGCAGTGGCATGATCTCGGCTCACTGCAAGCTCCGCCTCCGGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCACCACCACACCCCTTTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCATGTTAGCCAGGATGGTCTTGATCTCCTGACCTCGTGATCCATCTGCCTCAGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCTCTCCCGGCCGCAATAGGTACAACATTTTTAACAAGGCAAATATGTTTTCTAAGTTGTAGCCTAATAATCAGCAGTGGTACCCCTCTAGTAATTCCTGCAAATAATCTATGTACATCATTCCCTTTTAAGAAACACTAACATTTTTTGTGAGGTTGAAGCCATAGACCTAGTAGTCAATCTTTGCTTTGCATGGGTAGCTGAGAGAAATGTCTTCCCTTCTGAGTAGTTTCTGCTGTTCAGTGAATCTAGATGGCAGTGTTTTCATTGTTCCAGTTTAACTGTTAATCACCCTGTTCCTGTCAGTGAGTGTATAAATGAGTATCCAGTAACGTTTATATGTATTTCTTATACTTTTTTTTAGTGTAATCATAATTAGGCTTTTACTTGGTGATAGAGCAGTAGTATTATCTTTTACTCCAGTTTTCTGATTTATTGCACATTAATGATGACATGCCTTTATCAAACGGCTCCCATTATTTCCATGACATCATATATTTAGCTTATTATGTGACTACTTGTAGCGCTTTTCTTTGCCCTGTTATCAGCGGGTCAAATTTTGTGTGTGTGTTCCTTGTAGTTTATTTTACTATGTATTCTTACAGGCTACTTGCAAAATAGTATGTGAATTTCATTTCCTTATTTTTCATTTCTTTTCATTTTCTATTTTTTCAATTTCTTTTTCTCCCCCTCATTTACTGCCGTAGTGATAAATATGTTAATTGTTCCTCTTTCTTTCTTATGAGTGGTTTTGCAGCCTGGCCCTTTAACTACCTCTTCTGCTTCTGAGATCATGTTTCTTAACGCTTCTGTTGACCCATCATACCTCCTTCTTCACCATTCATCCTTCCACTGGCAGCTGTGCCTTTGATTCTCCTTTATTTTTTTCTGTTTTTTTTTTGTTGTTGTTGCTGGTGTTTGTTTGTTTGTTTTGTTTTGTTTTGAGATGGAGTCTCGCTCTGTTGCCCAGACTGGAGTGCAGTGGCGTGATCTCGGCTCACTGCAAGCTCTGCCTCTCGGGTTCATGCCATTCTCCTGCCTCAGCCTCCCACGTAGCTGGGACCACAGGTGCCTGCCACCACGTCTGGCTAATTTTTTGGTATTTTTAGTAGAGACAGAGTTTCACCATGTTAGCCAGGATGGTCTCGATCTCCTGACCTTGTGATCCGCCTGCCTTGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACTGTGCCTGGCCTATTTTTTTCTGTTTAAATGGAGCTCATAACAGCTTCTAATGCTGCCGTTTGTTCATCTGTTCATTGAAGGTTCTGTTGATCTCTAGATGGGTAAGAAGATACAGTCCTTAAGGACCACACAACCTTGTGGGGAGACTGATACATATATGCTATTTGTAACTCTATTACAGCACTCATTATATTGTATACTATTTGTTTATGTGTGTCTCTCCTCTACTTGATTGTGAGCCATAAGTGCTCAGTGAATGTTCATTGAATGAATGAATTATGTATGATATATGTATGTGTGTATGTATATATAATGACAAGAACTATGAAAAATGAACAGAGGCCAACCATGGTAGCTCATGCCTGTATTCCCAGTGCATTGGGAGGCCAAGGCAGGAGGATTGCTTGAGGCCAGGATTTGAAGAACAGCCTGAGCAACAGGAGAGATCTGTCTCTTAAAAAAATAAGATAAATTCGTTGGACATGGTGGCGCATGCCTGTAGTCCTAGCAACACAGGAGGCTGAGGCAGGTGGATCGCTTGAGCCCACGAGTTTCAGGCTACAGTGAGCCATGATTGCACCACTGTACTCCAGCCTGGGTGACAGAGCAAGACCGTGTCTCTAAAAAAGAAAAAAAAAGAAAAGGAAAAAATGAAGAGAATACTATAGGAAACAGAGAAGTACTAATTGAAGTACTAATTTTATTAAGAGCCAGTTTGGGAAAGATTTATAGGAGATGACATTAGAGTTTGGCTTTGAAGGGCATTCTTGCCAGATAAAATTTAGGAAAAGACAATGTGGAGATAATGAAATTTTCTGGATATTCAGAAAGTGTACTCTAGATATACCTATTCGACAGGAAAGTAGTTATACAGTTTTGCTTAAAAAAAAATAAAAAGATGAGACTGTAAATTAATGTTGTGAAGGTCTCCAAGTGTCCATCATGCTAAGGAGTTTGGAGCTCTGTCCCATAAATGATAGTCATTGAAGAGTTATAAAATCAGAATTGCTAGGTTCTGTTTTAGAACGATAACTTTGATTGCAGTATGGAGAATGGATTTAGGAGCAGGATAATTTAGGAAGCTTTTACCATAATTTAAGCTAAAGATGATGCTTTGAACAGAGGCAAGGGGAAGAAAGAAAGAAATTAGAGAAACATCTTAATTTTTTTTTTTACAAATAGCTTTATGTATGAGGCATATTTTACATATAATAAAGAAACAAAGGTATAAATGATATGATAACTAAGTATATGATGGAACGTGGGAGAGATAGGAATGTCAAAGATGAATTGGATGCATGTTAATCCAATTTATTAAGAGTACAAGAAGAAAATTTTGAAGATAGATAATGAATTTATTTGGTTTTGAACATGAAAAATTTGAGGTGCCAATGGTGCTTGGTGATAGACATCTAATAGACTGAAAACTTGCATCTAGAGCTCAGAAATAGATCTAAGGAGGAAATAGAACATTTAGGAAATAAAGCCATGAAAGTACGTGAGGTTGCCTGGGTAGAAGGAGCTAGAGAACAAGAGGCCCAGGAAAAAGCCTTGGAAATGTTACCATTGATGAATTAGACAGAAGTAGAGAACCCATCAGGAGTTAAGAGAGGTAAGAGAAAAAAACAATCAAGAGGCAAGGGAAACAGGAGAAATATGTTTTACGAATCCCAAGAGTAGTGAGAGTTTTAAGGTGAAAGTTTAAGATCCCAAATATTCATAACTTAGTGTTGTATGCCTGTTGACTGGAAAATTGTATTACTAGAATGTAGGCATATTTCCTAGGTTTTCTAGACAGATACAATTCTTTTTTAATGTTTAAATTGTTCTAAGACTTTCATAAGAATAAGCCTGTAAATCAGAAAAAAAAATAGCTTTTTAGACTATACTCCAATTTTTGGTTTAAAAAACATAGACTTTTTAGTTGAAGGACTGTTCATAACAAAGTGAAAGATAAAAGCAACATTTTATCAGTGGGCTTAAAAGAAAACATTTTCTGATGACCAGAGAAAACCCTGCACACATCTCTGAGAGGAACTGTCATTCATTTTGGGGTATTTTCTCTATTCTGTCTCTTAAAAACAGTACATTTCTTAACACTTTTTTTGTAGTATATTTGAGATTCATTTCCAGTATGCTACAAAAAATGTGTTAAGAAATATATGTCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATGGCATGAACCTGGGAGGTGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTACAGCCTGGGCGACAGTGAGACTCCATCTCAAAAAAAAAAAAAAGAAATATATGAATGTAAGAAATACAAGAAATGAAATATAGGAAATGAATTTCAAATATATTTTTCATTAATAGCTTTTTTCTCTCTAAGTACATGATTTAGAGGTTTGAAATGAATGATACTATTTTGTTGGATCCTCTTTTTCTGTGGGTTTTTATCTTTTATTTTTCCCCCGTAAGGAAACTCAAATAAAGGGAAATTTGTTCTAAGACTTCAGAGAGGCCTCAAAGAATGCAAGAGCAGAAAATGAGGTACTAGCATTATGGGAACTGGAAAATTGTCAAGAACAAGACCACTTTTTCCTTTTCTTTGGGGCCACTTTTATTTGTTTTTCTTTGATCTTCCTCTTTGTTTTTCATTGTATTTTTCAGCTTGGCTCCTTATAGCTCATGTTCCGTTGCTTATCCTTTCAATTCCAAAACTCTGCAGTCTGTTTTTCCGGGACTGGATAGGATGCTTTTGAGTCGCTTCCTTAGGCCTTATCCAGTCAGCTAGAGCTGGGGTGAGGGGCAGGTCTTTGTGGATCACAGTAGTGGATAGGGGTAATTCCCTAGAGAAAAGGACTAGGGGACTATGGGTGATTGATATCTCTAGTATTTACACAGCTTCATTTAACAAAGCTTTTTGGGAACATTTTGAAGGCATAATTTGTGATCAGTTGTTTTGTGCTTAGTGGGAAGGTGGCCTTTTTGACTTTTGAAAAATGTTAGAGTTGAAGGAAGATGAACATGGTAAATACTGAGTTCAGGCTAGAGGTTATTCTTTCTGGATGAGGGGAGAGAGATGAGTATGAGCAGAGTATACACAGAGGACTCCAGCCAAATTGGTAATAGTTTATTCCTTAAACCGAGAGGTAGCTACATGGGTGTTTATTATATTATTTCTTATATCACTTGATGTGAAGAAGAAAGCTATTTGCAGATAATTGAAGATCATTTGAATTTAGACCCTGGCAGACAGTAAGCCTCAGAAAATAAAGCAGTCCAAAAGGCTACTGAAGGCTTAGGACATTACTGTACACTATGTAGACTTTATAAACACTGTGCACTTAGGCTACACTAAGTTGATTAAAAAATCTTTCTTCAGTTATAACCTTAGTTTACTGTAACTTTTTTTACTTTATAAACTTTAATTTTTTAAAACTTTTTGACTCTTTTGTAATAATTCTTAGCTTAAAACAAATATACAGCTGTACAGCAATAGTTTATATCCTTATTCTGTAAGCTTTTGTTCAATTTAAATTTTTATGTTTTTACTTTTTATACTTTAAAAATATTTTTATCTATAAAAATATATTTTATAAAATAAAATATGTATAAAATATGTATTTTTATTTTTAGTAAATAGCAGGAGCACACTCTATAATAACAATAAAAAGTATGGTATAGTAAATCCTAGGCGATAGAAATTTTTCAGCTCCAATATAGTTTTAAGAGACCACTGTCATATATGTGGTTTGTGTTGACCAAAACATCATTGCATGGCGCATGACTGTACAGAGAACTTTGGCTATCTTCCATTATGAAAGTAGAATGTTGTGACAAGCTTATTACTATATAACAGAATGCTTAACTTTATCTTTTAAATTAAAAAAAAATTTAAATTGACACAATTGTTGCACATATTCATGGGGCAATAGTGATGTTTTGATACAATGTGTAGTGATCAGATCAGGGTGATTAGCACATCCATCATCTCAAATATTTATCACTTGTTTGTGTTGGGGACATTCAGTATCTTCTCTTCTGTCAACTTGAAAATATGCAATATGTTCTCGTTAACTGTGGCCATTCTACAGTGCTATAGAACACTAGAACTTATTCCTTCTATCTAGCTGTAATTTTGTATCCTTTAGCAAATCTCTCCTAAGCCCCCTTAACCCCCCTCCCAGCCCCTACTAACCTCTGTTCTACTTTTCTATGAGAAGAGCTCTTCTATGAGCTTCTGCATAGGAGTGAGAACATGCTGTGTTTAGCTTTCTGTTCCTGGCTTATTGACTTAACATAACGTCCTCCAGGCTCATCCATGTTGTTGTAAATGACATAATTTTTTTCTTTTTAAGGCTGAATGGTATTCCATTGTGCATATATATTTTTCTTTTCTTTTTTATGGAAATGTGAATGGTTTTGACTTAAAATGTGTTTTGCCTAATATAATTGTGACCTCCCGTGTTCTCATTTGGTTACTGTTTACATGCAGTTTCTTTTTTCAGCCTGCTAACTCCAGTCTATTTTTGTCATTAGATCTAAAGTGAGTCTCTTATAGACAGGATATAGTTCCATATTTTCTTCATCTCTTCATCTGTTGTCGGACACTTGGGTTGATTCCATTTCTTGGTTATTGTGAATAATGCTGCTATAAACATGGAAGTGCAGATATCTCTTCAGTATACCGATATTCTTTCCCTTGGATAAATATGCAGTAGTGGGATTGCTGGATCATATGACAGTTCTATTTGTAATTTTATGAGGAATCTCCATACTGTTTTCCATAGTGGATATACTAGTTTACATTCCCACCAATAGTATGTAAGAGTTCCTTTTTCTGTATATCCTCGCTAGCATTTGTTGTTGTCTTTTTGGTAATAGCCATTCTAACTGGAGTGAGATGATATCTTATTGTGGTATTGATTTGCATTTCCCTGATGCTTAGTGATGTTGAGCTTTTTAAAAACATGTTTCTTGACCATTTGTGTTTGCCCATTTTTAATTAGATTATTTGGGGTGTTTTTTTTTTTTTTTGCTCTTAAGTTGTTTGAGTTTTGTATGTTCTAGATATCAATTCCCTGTCAGGTGAGTAGTTTGCAAATATTTCCTCCAATTTTTTAGGTTGTTTTTTCATTCTGTTGATTGTTTCCCTTGCTGTGCAGAAGCTTTTTAGGTTGATATAATGCCATTTGTTAATTTTTGCTTTTGTTGCTTGTGCTTTTGGGGGTCTTGTTCATAAAATCTTTTACTAGACCAGTGCCCTGAACTGTTTCCCTTATGTTTTCTTCTAATTATTTTATAGTTTTGGGTTCTTATATTTAAGTCTTTAATCCCTTTCAAGTTGATTTTTGTGTAGTGTGAGAGATGGGGGTCTGTTTTCATTCTTCTGCATATGGAGATTCAGTTTTCCCAGCACTATTTATTGAAGAGACTATCCCCAATGAATGTTCTTGGCACCATTGTAAAAAATCAGTTGGCTATAAATACATGGATGTTCTCTGTTCTGTTCCATTGGTCTATATGTCTGATTTTATGGCAGTACTATGCTGTTTTGGTTACTGTAGCTTTGTAGTATATTTTGAAATCTGGTAGTGTGATGCCACCAGGTTTGTTTTCTTTGCTCAGGATTGCTTTGGCAATTTGGTGTCTTTTGTGTTTCCATATGGATTTCAGGATTGTTTTTTCTATTTCTGTGAAGAATACCATTGGTATTTTGATAAGGATTGTATTGAACCTGTAGATCACTTTGGGTTGTATGGTCATTTTAACAATATTAATTTTTCCAGTTCATGAACACAGGGTGTTTCTCTTTTTTGTGTGTTTTCTTTAACTTTTCTCATCAGTGTTTTACAGTTTCCTTTGTTGAGATTTTTTACCTCCTTGGTTAAACTTATTCCTAGGCAGATTTTTTTTTTCTTTTTTTTTGTATTTACTGTGGATGGGATTGCTTTCTTGATTTCTTTTTAAGGTAGTTCATTATTTATATATAGAATCACTACTGATTTTTGTATGTTGATTTTGTATCCTACAAATTTACTGAATTCATTTATTAGTTCTAATTTTTTTTTGGTGGAGTCTTTGGGATTTCCTATATGTATGATTATATTGTCTGAAAACAGGGACTATTTGACTTCCTCCTTTGCAATTTGGATGCCCTTTATTTCTTTCTCTTGCCTAATTACTGTGGCTAAGGTTTCCAGAACTATGTTGAACAAAAGTGGTGGAAGTAGGCATCCTTGTCTTGTTTCAGATCTTAGAGGAAAAGCTCATTCATGGATTACTAGGTTAAGAATACTTGTACTAGATTATGAGATTGTAGAGGTCAGAACCATTTCTCTCTCTCTCTTTCTCTCTCTCTCCCTCTCAGTTTTTTCATCTGTAAGGACCATTTCTTATTTTTAATTTTTTTAATAACATACAATGCTTTGCACAAAGTAGGTACTCAAAAAACATTGTTAGCACAAGCAGTTGATATGCAAGATTACTTTCCTCATTCTTAGTCCTAAAAGACAAATATTATTACATTGGGTTAAAAAATTCTATTTTCTATATACAATATATGTTTAAAAGAAAGTGATCCACAAAGCTGAGAATTAAAAGGAAAGGGCAGATTCAATAAAGCAATGTTTGCTTTTTATTAATAGCATGCAAAGTAGAATTAAATACTAAAACACATTTCTCTAGGGATAAACATAGTCATTTTTTATTGATGAAATATAAAATCCATTGACAGCATAACCATGAACTTTTATATACCAAATAATATCAAGTTATTTAGAAGTTAAATATTTTAGAAATATGAAATGGACAAAAATTCAGTACTTCCAGTGAATTTTTGTCCATCTCTTAATTCTTGGCAAATTAAATAGACAAAAGATAATACACAAATGTCCAACAAATAAATAAATCCAAATTTCCCATATGACTGTAGAATAAAGTTATTGTCATGATATCTCTAAAACATTTAACTGGGAAGGTGCTTTCATTATTCCGCAAGAAGTGATCTTGAATAATCTGTTTGACTTGGTAAACACTTAAAACTGGTTGGGAAAAAGGCCTTATGCAAATTGAGATTTATACTTTAAAGATAGTAAACCTTTTGATAGATTATTTTATTATAAGTAAGTTTTTTATTGCACTTTTTTCACTAATTCAATTTTCAGTTTCTACAAAAATGGTGATTTAAAAAATGTATTCTTAAAGAACAGACCAACTTTAACACATTATGTTTAAAATATTCACTTACAGGGCATGGGATAACCAAAGATTTTGAAAGATAGATTATTGGCAAATTCAGGTGCTTCATGAGAATAGATATCATAAAACTGATTTAAAACTATAATTATAGGCTGCATCTGGTGGCTCACACCCATGATCCCAGCACTTTAGGCTGAGGTGGGAGGATCACTTGAACCCACGAATTTGAGGCCAGCCTGGGCAACATAGTGAGACCTCATCTCTATAAAAATAAAAAAAAAATTAACCGGCCATGGTGGTGTGTACCTGTAGTCTCAGCTACTTGGGAGGCTGAGGTGGGAGGAACACCTGAGCCCAGGAGGTTGAGGCTGCAGTGAGCCTTGATTGCACTACTACACTCCAGCCTGGGTAACAGAGTGACACCTTGTCTCTTAAAAAAGCCTCAAAAACTAGAATTATAAACCTAAAACAGGGCATAATTCATGTTATCTGAAAATCCTTTGAAGTATAATATACATTGAAAAAGCACATAAATGTGCAAGTTGATAAAGTCTCTTAAACTGAACGTACCATTGAAATCAGCACCCAAATCCATGAATATTAGCAGTACCTCAGAAGCACTTCCTTGCTCCAATTTCAGCCACCTCTCTAGCCACTTCAAGGGTGATGACTATCCTCACTTCTAGCACTATAGATTCCTTTTGCCTATTTTATACTTTATTTAATGGAACTATATGGTATATACTCTTGTTTCTGGCATTCTTCACAACAGGTTTGTGATATTCAACTATCTTTTTGCATGTAGTTGTAATTCATTCTCATTATTAGTGGGGGATACTCTACAGTATGAAGTTTTACATAAATGTCATACTTTTAATATTTAGTAATTTATATTAAAGATGATTGCCCATTAATAGGGTATTTATTTTGAAGCAGATAAATTTTCTCTTGTTTATTCATTTATACTGTCTCAGCCTTAGTGTTAGCTAAGAAGACATAGAGATTGTGGACACATATGTAGTTTTTGATTATTTATTATATCCTTTCCATTGGAAAATTGATAATCACTATATATTTCAGCACTTATTTTAGTTTTAAAGCAAAAAGATATCTGCCCAAAGAAATTTTTTGTCTAATCATTGACAGATATTCTAGTAATATAACTGATGTTTACTTAATTTTTCAGTTGGCTAGATTAAAAATATTTTCATCAAAAGGAAAATGAGTTTTTTGTTGTTGTTGAGATGGAGTCTTGCTCTGTTGCAAGGCTGGAGTCTTGTGGCGCGATCTTGGCTCACTGCAACCTCCGACTCCCTGGTTCAAACGATTCTCCTGCCTCAGCCTCGTGAGTAGCTGGGATTACAGGCACCCGCCACCATGCCTGGCTAATTTTTGTATTTTTAGTAGAGACAGAGTTTCACTGTGTTGGCCAGGATGGTCTTGATCTCTTGACCTCATGATCCACCTCCTGCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCAGCCAAAATGAGTTAATTTCTTAAACTGTGTTGGCAGATTAAAGTTTTATGGTTGATAAAATGTTGTTTTTCTACCATGCTCCTTTATTTTCAGGAAACTTGTTTACTTTTCTTAGAGCCATATTCTATAATCTATGAAAAAATTCTAGTTTTGGAACTTTTTTTTTTGAGACGGAGTCTCGCTCTCTTGCTCAGGCTGGAGTGCAGTGGTGCGATCTTGGCTCACTGTAACCTCTGCCTCCTGGGTTCAAGCAATTCTTCTGCCTCAGCCTCCTGAGTAGCTGGGCTTACAGATGTGCGCCACCATGCCTGGCTAATTTTTGTATTTTTAGTAGAGACAGGGTTTCGCCATGTTGGCCAGGCTGGTCTTGAACTCCTGACCTCAGGTGATCCACCTGCCTCAGCCTCCCAAAGTGCTGGGATTACAGGCTTGAGCCACTGTGCCTGGCCCTATTTTTAGAACATTTTTGAGGTGTTTTCTAGGAAATTTATGGATCAGGATTATTTACCCTGATAAATCAGATTAAATTAACTTCCTTAAAAATTAACCATACCTTTTAGTTTTACATATTTGGACTTAATAGGTGAGCCCTTCTTTTTCTTAAATTTAAATTGGTGTTTTTAATATTCAACATACCGTACTGCCTCTGCATACAACTAAGTGGAGAAAATAATGTTTATCAAATATGACAGGTTTAATATATTTGTAAGGCAAATAACTCATTAACTTGATAAGAAAAACATTAAGATCCTGGTAGGTACGTGGGCAAAGGATAGAAAAGATAATTCACCAATGAAGCACTAAGGTTAAGATATAGAGCAGAACGTTTAATCCTGCTAGCATCTACCAAATTCAAATTAAAGCATAGATTTTTTTCCTGTTAAATCAGCAAAAGTAGAAAGAATAGATAATCTTAAGATATTTATGAAACAGGTAGGCCCTTCTTGGAAGACAATTTGACATTGTATTTCAAGAGCAGCATGTTCTTCACATTCTTTTGTTCAGTTTTAGCACACTTTGCAGATACTGTTCTAGGTAACCAAATAAGCAGAAATTCATAAATACAAATATGTTCATTTCAGTGTTATTTGTAATAGAGAAAATGGGACATAATATTGTATTGGTTATCTCTTGTTGTCTAACTCATAACCTCCCACCCCAATTTAAGCATTTTAAAACACCAAGCATTTGTGATCTCACATTTTCTGTATGTAAGGAATTTGGGAGCAGCTCATCTGGGTGATTTTGGCTCAATGTCTCTCATAAGGTTGCAGTCAGGATGTCAGGTAGAGCTGTAGTCCTGTGAAGATTTGACTTGGTCTGGAGGATCTACTTGCAAGGGGATAGCTCACCCACATAGCTATAGGCAGGCGGCTTCAGTTCTTTGCCACATGGACATCTATAGGATTGCTTGACTGTTCTTAGGACATGGCACCTGTATCCAAAGAGAGATGAGAATGCGTTTTATAGCCTAGTCTCAGAAGTCACATATTGTCACTTTCACCACATGTGATATTAGTTTAAAATAAGTCTCTAAGTCCCACCATACTCAATGGGAGAATAATTAATTTCCACATTTTGAAGAAAGGAGTATCAAAGAAATTTTGGATATATTTTAAAACCACCACAAGTGTCCAGCAGTTAAGTAAAGTAGGGCACAAAATTGAAGTTACAGTGATGATAGTAATAAAATTACTTATTAATGATTATATTAGGATTTAGGGAGTAATGAGTTTTTATTCTTATATTTTTCAAACTTTCTGTGTTGTTACATGATTTTTATAAAAGACTAGTTCCTCTATATTCATAGGAAAATACAGTATTACTTATCACTGAGCTAACCATGGGCAGATCATGATCTTTCCCTTTGAGTTGTGATATGATTTTTTTCCTTTCTTTTTTTATTGTGGTAACATATACATAACAAAATTTACCATTACAACCATTTCTTAAGTGTACAGTTCAGTGGTGTTAAGTACATTCATATTGTTGTACAACCATCACCACCATCCATTTCCAGATTTCTTTTCATTTTGCAAAACTGAAACCCATTAAACAGTAACTCTCCATTCTTCTCTTTTCCCAGCCCCTGACAACCACCATTCTACTTTCTGTCTCTATGTTTTTGACTATTCTGACTGTCTCATGTAAGTGGCACCATACAGTATTTGTCTTTTTGTGATTAGCTTATTTCACTTAGTATAATGTCCTCCAGATTCTTTCATGTTATAGCACATATCAGATTTCCTTCCTTTTAAGGCTGAATAATATTCTATTTTATATATAAATGTCATTTTGCTTATTCATTCATCTGTTGATAGACATTTGGGCTACTTCCACGTTAGCTACTATGAATAATGCTGCTTTGAACATGGTTGTACAAATAGTTCCTTGAGCCTCCACAAGCATTTTGGTTGTATACTGTATTAGTCCATTCTTACACTGCTATAAATTATCCAAGACTGGGTAATTTATAAAGAAAAGATGTTTAATTGATTTACAGCTCTGCATGGCTGGGGAGGCCTCAGGAAACTTACAATCATGGCAGAAGGCACCTACCTCTTTGCAGGGTGGCAGGAGAGAGAATTAGTGCAAGCAGGGGAAATGCCAGACGCTTATAAAACCGTTAGATCTCACGAGACTCACTTATTATCATGAGAACAGCATAGGGGAAACTGCCCCCCTGATCCAATTACTTCTACCTGGTCCCGCCTATGACACTTGGGGATTATGGAGATTACAATTCAGGGTGAGATTTGGGTGGGGACACAGAACCAAACCATATCATATATCTAGAAGTGGAATTGCTGGATCTTACAGTAATTCTAATTTTACTTTTTTTAGGAACTGCCACACTGTTTTCCATAGTGATTGTGTCATACATTCCCACGAACAGCACACACGGGTTCTAGTTTTGCAAATCTTTGCCAGCTGTTATTATTTTGTGTGTTTTTGATAGTAGCTTGTCCTATTTACATCCTAATGGGTGTGAAGTAGTCATAATATAATTTTTAATGATTTATTTTTCTAAGGAAAAGTAATATTATAAGATTGATGTTCGGTGAATACTGATGTGATTGCAGGTAAATTGGCAGGATGAGGAGCAGCAAATCAAAAGAGGTGCCTTTACCAAATCCAAGGAACTCTCAAAGCAAGGATACTGTTCAAGGTATGATTTTGTTTTTTTAAACAGAACTTAATACCTCATTTAGGTCATGTGTATTAAATTTTGGTCTAATAAATGAAAATTATACCTATGCTTTGATTTATAGAAAGTGGATAAATCATTTTAAAATTCTCTTATATTTTAAATATACCTGGTGCTTCATGTTAAAAGGCATCCTATGGTAAACTGTTTTGTAAAGCAGGAATCCTTTTTTAAAATGAACCAGTATATGTAAAATTTTTCAGAAATTCTGTCCTAATTTCAATAAATGATATAGCCCCAATATAGTTACAGACTTATAGGGGTACAGATAATAATAGTAACAGTCTCATTGAATTCTTTTTAAAATAATGCATTAGAGGTAATTTTAAGGCATTGATATTTTTTGAGGGAGATTGACAAACTAGAGTACATCTGGAGGAAGATTATTAGGTTTGTGAGAGATCTGGAAATCAAGTCATATGAAGAGTAGTTGGAGGAAGTGGAAATACTTATCCTGGAGAAGGATGACTTATGGAAAATACACGACCTTAAACTATTTGAAGGGTTGTTATCTGAAAGACTGAATGGATTTTACCTATGTTTTTACTGGTAGATCAAAATCTGGAATAGCTTCTAATCTTTTTCACATTGTATTATAAGATACAGTACAATATAAAATTTTATATTTGTTCAGCTCATTGGGGTAAACTGAGGATGCTGCTTTTTTTTGGAGCTGTCTGATACAGGTGTGAGTGGAAATGCCAGTCCAATGGTAGCAGTGTCTGTTTACCTTTAAACCGTATGTTGGGAATGTTGGCTTAATAATCTTATAGATAAAAACCCATTTTGACTCAATATAATAAAGCATTTGCTGAGGAATTGTTCTGCCTAATAAGGTAGTTTTCAGACTTTTGAGTGCATTTGGGAGCTTATTAAAATTCAGATTCTGAGGTCTAATTCAGAGATTCTGGTTGATTGATAAAGTAAATATCATTTTTTTTCCCAGCAGATATAACCACATCGTGGGATGCACTTTCTCAAACCAAGGCTGCTGTAAGTAGTTTTAGCTTCCATTTATTTATTTCTGGATTAAAATCAATTTCCTGTTATGATTGAATAGATTAAAATGTTTTCCATACTGTTCTTTAAATACTATATATTACTATAGAATAATACTATATATTATTATATATACTATATAAAATATGTATACACTATATATAATATATATACTATATACAGTATATACTTATATACATAGTATATAGTAATATTATAGTATTATCATATATAATATATAATATATTATAAAATGAACATAAATGACAACTATTAAATGGTAAAGTCATAAACAAGATACTACTTTAGGAATTAGGAGCACAAGCATCATGCCATATAAACTATTAAAATTGGCTGGTTGCGGTGGCTTACTCCTGTAATCCTAGCACTTTGGGAGGTTGAGGTGGGTAGATAGCTTGAGCCCAGGAGTTTGAGATTAGCCTGGGCAACATGGCAAAGCCCTGTCTCTACCAAAAAACAAAATTTAGCCAGGCCTGGTGGTGTGCACCTGTAGTCCCTGCTACTCAGGAGGCTGAGGCAGGAGGATTGCTTGAGCCCAGGGAGGTTGAGGCTTCGAGGTTGCAGTGACCCTTGATTGTGCCACTGCACTTTAGCCTGAATGACAGAGCAAGACCCTGTCTCAATAATAATAATAATAATAATTTGTACTGTTTATTCACCTTTAATAAATTTGTTCATTTGGTATGTGGAAAAGAGAATTTAGACTGAATGCAGTGGCTCATCCTGTAATCCCAGCACTTTGGAAGGCCAAGGTGAGAGGATCGTTTTGAGCCCAGGCATTCAAGACCAGCCTGTAAAACATGGCAAAACCCTGTCTCTACAAAAGTACAAAAATTAGCTGGGTGTGGTGGCATGTGCCTGTGGTCCCAGCCACTCGGGAGGCTGAGGTGGGATGATTACCTGAGCCCAGGAGCTCAAGACTGCAGTGAGCTATGATCACACCACTGCATTCTAGCCTGGGCGACAGAGTGAGATCCTGTCTCAAAAAAAAAAAAAAAAAATTTAAAACAGAGGATGAAGGTAGGGCTTTGTCTATGTAGAAAGCTCATTGTTAGGGAGGAAATGAAGGAAACCAACTTCTTGAATAGACTGTGAATTATTCTTCCTTTCCAATTTGGATGCCCTTTATTTCTTTCTGTGGTCTGATTGCTCTAGCTGGGACTTCCAGTACTATGTTGAATAACAGTGTTAAAAGAGTTACTTCTTTTTTGTATCAATAACAGGTTGTTTTTTTTCTCTTACCACTGAAGAAGAATGACTATGGGCAACCTGTTCCCATTTGTGTTCTGGCCAATCTTGTCTAAAGTTTAGGAATCCAAATTTGGCAGGTGATGCATAGTTCATTCAAAGGGAAACCCTTTTTTTTTTTTAATTTGCTCTCTGGGCAAATTACTTTCCTAATATACCCTTTTTTGTTTTCTATGCCACCTGATTTTGAATCACCTATTTCTTTACTAAGTATTTCAGTGATATTTTCTGTCCTCTAATCTTTAGAACTAGTTTGAAAATAAATCTGGAAAAGGTGTTATAACCTAGATTTGGTCCCAATGCAGAGTCATATTTTGGCCCTATAAATATAGCAATATATATTCTCTTTCTTCTACCTGAGCATTATGTTTTCTATAGTGTAGAAATAGTTTTTATCCCCCAGCTTAACTTAGTGCTAAACTAATTTAGAGTTTAGCACTAAGTTAATAACTGTTAATAGTTGTTTTTCAAAACTATTGGGTGTTTTTGAAGAACATTTAAGGATGATAAAGGAATATAAGTAGGGGAGACAATCTGATGTAGAACTATTAATATATCATAAATTATTCCAATTTTTAAAATACATTGGGATGCAAAAATCTGAAGGAAGGCTGGGCGCAGTGTCTCACGCCTGTAATCGCAGCACTTTGGGAGGCTGAGGCACGTCAATCACTTGAGCTCAGGAGTTCGAGATCAACCTGGACATTTGGGAGCTTATTAAAATTCAGATTCTGAATTTTAAACTCCATCTCTACAAAAAAAAATACAAAAATTAGCTGGGCATGGTGGTGTACATCTGTAGTCCCAGCTATTTGGGAGGCTGAGGTGGGAGGATCGCTTGAGCCCAGGAGGCAGAGGTTGCAGTTAGCTGAGATCATGTCACTGCACTCCAGCCTGGGTGACAGAGTGAGACCCTGTCTCAAAAAAAAAAAAAAAAATCTGAAGGAAGATATATCAAAATTTTATATAAAATGTCTATATGTGATGGAATTTTGATTTTAATTTCATTTTTTATATTGTTTCTGTTTTCCTAGTTTTTTGTTGTAAGCATGTAATACTTTACCATTTAGGAGAAAATGCTTTTCAGTCATTGAGAACTTTGGTGGTTATAATCAAAAGTTGACTTATATTGCAGAGTGTGTATGTGAAACCTGAATTCTATCTTGCGAATTTTGTCTTTATTGTGAAACTTATTTTATAGACATAAAACTGTAAAACACTGTGGCTGTTTGAGTACCTATAGTCCTAAATTTTGGTTCAGGTAAATAATTCTTTCCTCATTTTTGCGTGTAAAATGCACAGGATATATTGCATAATGGTATGTCACATATAAAGTCAGATAATTGGCTGGGCATAGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCTGAGGTGTGTGGATCACCTGGGGTCAGGAGATCGAGACCAGCCTGGCCAACATGTGGAAACTCCATCTCTACTAAAACTACAAAAATTAGCCCGGCGTGGTGGCGCACGCCTGTAATCCCAGCCTCTAGGGAGGCTGAGGGAGAATCTCTTGAACCCGTGAGGTGGAGGTTGCAGTGACCTGAGATCACGCCACTGCACTCCAGCCTGCGCGACAGACCGAGACCGTCTCAAAAAAAAAAAAGTCAGATAATTTACTACTTTGAAAATACTCAACAATTTAAAAAATAGATCAACATGTCAAATGGAAACTCTTTGATATTCTTCAGAATCTTACTGATTTTCTTGTCTAAATTTGGGGAGCATTGTGCCAGATGTTATTATAGCTAAAGTAAAAAATGATGAAAAACAAATATTATAAGAAAGGATTCCTGAGTTTTTTTCTACTTTGGAAAAAGGTATTGGGCAATATACTTTTAAATAACCATGTAAATTGATGGTAGCTGGATTATTTACAGAATTAACTTATGTCTTATGATGGTGTTTAAACTTTAGCTGAGACACATTGAAAATAAATTAGAAGTAGCCCCTACAAGTACAGCTGTGTGTGATTCTGTCATGGATACCAAGAAGTCTTCTACAAGTGCTACTCGAAAAATAAGTAGAAAAGGTATGTATGAAGTTACTTCCAAATATATACATATATTAAATATAAATTTTGACAAAACAGACAGCTAGATCTACAGGGGTGTGCCACCACACCTGGCTAATCTTTTTTTATCTTTTGTAGAGACGAAGTCTCTCTGTGCTGCCCCAGCTTGTCTTGAACTCCTGGCCTTAAGTAATCCTCTCATCTTGGACTCCTTTTAAAGTGCTGGGATTATAGGTGTAGCTACCACACTCTACCTTAAATATAAGTTGTTGATGTAATAGTTTTATTTCAGTTTGATTTCTTCTTACAGTCTCGTGTTGTTTGCTTTTTTGTCAGAGACACTGGCATTTAATTTTTTTCCTTTGATATTTGGTAGCTATGAAGCAGCCTCATGTTCAGTAACCGTAGATGGTGTTTTCACCATCTTTTGACCAGGGACTTTAGCACTATTGCATCACCCTTGACAATGAGTTTTCATTTGTTTTTAGTGATGAGAGTTTCATTTTGAGGTTATATATGTCTGATTTTTTTTTTTAATTCACAGATATATAGCACTTTGGAAATGAAATTTAGTTGACTTAACGTATATTTTTCTTGGGATTTAAAATCATAAATGGAATCTCTTTCATATGGTGACTGTCATAAATTGATACAATGCAAATTTTTCTTGGGTTACCTTTTGTATTTTGATTCTGTGGTTTCTTTCTTTTCTTTTTTTTTTTTTTTTTTTTGAGACAGAGTCTCACTCTGTCGCCCAGGCTGAAGTGCAGTGGCATGATCTCTACTCACTGCAACCTCTACGTCCTGGGTTTAAGCGATTCTCCTGAGTAGCCGCGCCACCGGGCCTTTCTAATTTTTGTATTTTTAGAAGAGATGGGGTTTCACCATGTTGGTCAGGCTGGTCTCAAACTCCTGACCTCGTGATCCGCCTGCCTTGGCCTCCCAAAGTGCAGGGATTACAGGCGTGAGCCACCATGCCAGGCCATTCTGTGGTTTCTTAATTGGTGGGTGAAAACCTCCATGTCAGTACTCTCTACCTAAACTGTGATATGTATATATATATATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTATATATATATATATATACATTTTTTTTTAACCTGGGAGTTTGTTAGGTTGTATATGATACTTCTCTCTTAAGTGAGTGAGAGCCATGACAAAAAATAAATTTTATGTGTAGAATTTTATTTTGCATTGTTTCCTCAGAATACAATTTCTGTATCAATGGAGTTACTTTAACAATTTAACCTGAATGTTTTTCATGTACTTATATAATATCCAACTAGAGATACTAGAGATAACTTATATAATATCCAACTATTAGAAATAATAGATACCTAATTTTTTTTTTAAGATACCAGCCTAATTCAGAGGCAGTTGTATTTTATATTATATGTTCTCACAGATTTCCTTTTCCTTCAGATGGTAGATACCTGGATGATTCTTGGGTTAATGCTCCAATCTCCAAATCCACTAAATCACGAAAAGAGAAATCTCGTAGTCCTCTCAGGGCCACCACCCTGGAGAGTAATGTGAAGAAAAATAATCGTGTGGAATTTCGTGAACCTTTGGTTTCTTATAGGTTAGTATTGAGAAAAAAAAAAGGTATCAAATAGGTTTAGCCTGTAATATTTACAGTAAGAGTATACGTTGTTTTCTGAGTTTTTTATAGCACTCTGGTT',
			'insseq'				=> 'G',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		my $expectedRnaWtSeq = 'GCCCGACUGCUCCCGCCCUGCCUCGCGCCGCCGCCACAGCCGCUGCCGCUGGCGCCGCUCCUCCUCCGUGUCAGUUGUUGGGCUGUAAUGGCGACUGGGCCGCCCCUGACGAAGUGACUCCCGGGCCGGGGAGCGGGGCCAGACCUGCGCCAGAGAGAACUGCAGGGAGCCGCAGCUCGGGGGGUGGCUUGCCCUGAGGGAGGGGAGGCAGCCUUUCCGCCUUGUCUUCCUUCCCAGCGGACCGGCGGAUCCCCGGAGCCGGUGCGAGGAGGGCACCCGGUGCGUCCCCGGAGCGGGGAGGCCAGGCCGGGCAGCCCUGGGGCCGGUCGGGGCGGCGUCACUGCACCCUCCGCCAGGCUCCGCGGGAUGCACCGUGGUAGCCGAGGGCGGAGGCGACACUCUCAGGUGAGCUCCUGUAGGACCUGGCUGGGACCGCGGAGUCUCGCUCAGCCUCAACUGUCUCUGUCCGCGGUCCCGGGUCCCCGGUGGCGGCAAGAGAGGCGGGGCCGGGCGCGGAGAGUCUUGGCCCUGCUUGGCGUCCCCCUCGGGGAACCUCCCUGGGUUCUGAGGCUUUUCGGGGACCUAGGUUUCGGGUCCCACCAGUGCCUGCUGAAGACGGAAUUUCCUUCCUUACGUGCCAGUGUCCUCCCGCCUCCUCCACCCCUCCUCUAGCCUGCAGAGGUAAGAUACCUCUAGCCUUUGCUGCCCGGAGAAGUGUCCCUUCUUCAUCUGGCGAUCAGGCGUAUUUGGGGAGAUUCCGGGGACUGGCUCCAGUGGCAGGGAUGAUCAGUUUCUGUUCCUUUUCAGUGGGGAUCACAUUCUCAAGUAUAGGGAAAACUUUGAUCUACUCUAGUCUAGCGGCUGCUGAAGGUAGUUUAUUCACAGAUAAUUUGGGGAAAGAUUUAAAAGGCUUUAUGUUUUGGAAGUUGUCAGUAUUUUCCUCUGCCUUUUAAAUGCAGCUCAACAUGCUAAGUAAUAUUUAGCUAAUGUCACAAGAUAUGAAUUGGAAUUUGCCUGGAUUUUAAAUAUUCAGAUAGAAAGUUGCACAUAGGCAUAGAUUUUAAUGUUUUAAUGCAGUAGCAAAUCCUUCAUAUGGGAAAAAACAGGUUGUCAGAAAUAGUGAUCAAUAUCCAUGUACUCCACCUGAUUGUUUUGUUUUGAUGGUGAUCCUUGUUUUGAAUUUGGGGAAGUCGUUAUUGUUGCCAUGUUAACGAUACCUCCUGCCCUGUAUUUGAAGAGUUUUAAAAGAAAGACUAUCUUUGGCAGAGCUAGUAACAAUCUGAAAAGUUCUGUAGAAACAAUGAUAAAAUACUAAAAGUUGCCCUCCAUUCUGAAAAUUUCUUCGUUUAAUAUAUUGUUUUUAAACAAGGAAAUGCACUGCUAAUGACAUGGAAUUGUGACUUUGAGGGAUUGUUAGACCCAUGCUGGAAGUUGCUUCAAAGCCAUGGGCGUAUAAUAAUUAGUUUCGGUGUGCUAAAUUGACAAAUUUUAUUUAAUGUAAUGAAAUAAAUCUGAAUCAUUUGAACUGUCAGUGCAUUAAAUUUAAGAUACUGUAAUCUGAAAAUAAACUUCAAGCUAGAAUACAUUGUAUAUGUGGAUAAAUACAUGUGUAUCUCUGUAUCUCUUCAUCCUCAGUAGCAAAAUUUUGUGUAUGUAGAGAUGGCAAGUUGUAUAUGUAAAAAAUUUCAUACGUAUGUGAACAUAUUUCUUUGAUUUUCUGAGAUAAAUUGGUGUGUGUGUCAUGUAAGAGCGUUAAUACUUAAACAGCCAACCUUCUGUGGGGCAAUCAAUUUAUGAAUUAUAAAGUUUUUGGAUUUUUUCCUUCUGCUUCUGAUUGUCAGUCUUUUGCCCAUUUUAUAACACUUCUACCAUAUAUGAGCAGAAGAAAGGAGAUGAAGUACUGCUAUAAUCAAUUUUAUUUGUUAAAUUGCUAUAAUUUAUAGUGGGAGCAAAAAAGAAGUCUAAAAUAACUGACAUCAGUUUUUGGUGUGCUUCUAUAUACAUUAUUUUAUUUAGUUUUCACAGUGAUCUUGUGAAGUAGGUAUGGCAGAUAUCCUUAUUUUGCAGUCGUGGAAACUGACACACAAGGAAGAUCAUUUAUGUGAAUGCCUCUCAUUUCUGCAGAUAAUGUGGAAUUAUUUUUAUUUUAAACUAAAUAUUAAUAUAGGUGCUCUGGCAUUUGGUGUUAUUUUAAAAUAUUAAAAAAUUAUUCUAUUUUUGUUCCUUUUUACAUGCCACUUGCCCCUUAUGCUGAUGAUCGUCAGUUCUUUAAAAACAUUAGUAUUUUUAUAACGUUUUGUUAAUUAUUGACAAGUCAUGAAUUAUGAUGGAUCACUUACUUCUUUAUGUAUAUACUUUUUGCUUGAACAAGAUUAUUAAAGUGUAGUUCAUUUAUUAUCUUUCACUUUGGAUACUCAGUUCUGUCAUCUUUAUUCUUUCAGUCAUUUAUACUCUUAGUUAUGUAAAACAAACAUUCAAGCCACAUGAUCAGUAACAUAACUGAAAAGUAGUAAUAUCUUGAUGAUUUUAUUGUUGAGACAACUGAUGAUUUAAAAAGUUAAUUCUCUGCAUAUAAUUGUGAUUGUUAUUUUUGCAACAUGAGGCAAAUUUUAGUGUUAUAGUGUAAAAAAAUCAUCAGAUUGCAUAGUAUUUUUUUCUGUAUUGAUUUGUUGUUGUUGUUGUUCAUUGUUAUCAUAAGUGUGAUUUUACAUUUACUUGUGGUUAUUUUGAGUUAUUCUUUGCCUUUCUCAUCUGGUCUGUAAGUUCCAUGAUGGCAGAGAUCCUUUAAUCACUAUUAUGCACCCAGUGCCUAGCAUAUGGUAACAUAGUCAAUAAAUACUAAUUGAAGAAAUUGACUAAAGAAUAGUCCCAGAAUAAUAUUUUAAAUGCAAAACAAAACAAAAGACUAGUUCUUAGUUCCUAAGCACACCUCAACUCUUGGGAAUAUAUUCCAGAGAUAUGCUCAGAUUUGCUGGGUAGCAAGGGAAUGACAGCCUGGGAUUAAAAAGCCACUGGGAAUAAGCCACAGUGGAAACUUCAUUUACUAUUUUAAUGAUUUGCCUAGGUGCUCUGGCAACAGGUUACGGUAAUUUAAAGAAAAUAUAAACUUUAAGACCCUAACCUGUUAUACUAGUACUAUAAUUAUAGCCUUGUAAUAACUUUAAGUCACUUCUUAAUCUCAGUAUUGGGAGUUUAAAUUCUUCAAGGCCAGCUAGUCCAAUUCCCAUAUUUUAUACAUAAGAAAACAGAAGGUUGUAUAGGGGCGAUAUCUCAUGCAUAGUUGAUUUCAUAACUGGGAGAAGAAUCCAUUGCACAGCUGCUAAGUCUACAUUUGUUUAAUAUGUCAGAUAAUUAUAAAAUCAGAACAUAAUUGUGAUCAACUUGGCAGUCUAUUCAGUAUUAAUUAAAAACUCAAAAUGGGCCUUUUUUUUCAUGAAUCCAAAGUGUUAACUUAAAAAUGAAUCUGAACCCUUUCCUCUGUAGUUGGGACUAUAAAUUGAUACCAUUGGAGCAACCUUCUGGAAGGAAAUUUGGCAUUCCCUAUCAAGAUUAAAAGUAUAUAUAUUGUUUGAUAGAAUCAUUUCCUGCUUAAAAUUUAUGUUACAUAAAUUUUAACAGCAUAGGAAUUGUAUAAGAAUGUUUAUAUUAGCAUUAUUUCAAGAUGUCAAUAACAUUUAUAAAUAAUUUAAAAAUGUCCUAAAUGUCUGUUGGUAGGCUAAUAGUUAAAAGUCUAGUAUUUUUAUGCUGUGAAUUACUGUGCAGUUGUUAACAGGUUUCAGGGUAUCUAUAUUUACUUCAUGAGCAGCAGCUGAGCAUAUCCUUCAGCUCCUGGAACCAUGCGGUGGACCUUCUGGAUUCAAAACCUAACUGCACUACUUUCUAGCUGUGACUCUUGGCGAAUUACUUAAAGUCAACUAUAAAUUGGAUAUACUAGUAGUUCUCAUGUUACAGGGGUACCAUGAGAAUUAAAUGAGUUAACAGUUAGGCCCUAAAUAAGCAUUAAUAAGCGCCACAUAAGUGUUUGCUGUUACAAUUAUGUAUUGACUUGGAAGAAUGUUCAUGAUAUAAAAUUAAGUGAAAAAAUAAGUUGCACUAUACUAUAGAAUAUGAUUUAUUUAAAAAACACAUAGACAAAUGUAGUUUUGUUUUCAUUUUGCAUGUGAAUAGAAAAAGGUCUAUAAGGAUAGACAUUAAACAUGUUAUUCCUGAGGGAGUGGAUGACUUUUCUUUUUUACUUUAGACAUUUUUGUAUGGUUUUAUUGUAAGAAUCACUUUUAUAAUUUUUGAAACUUUUUUCUAAAAAAAGGAAGUAAAUUUUUUCUUCUAUUAAAUAGUGAUAUAAAUGUAUUUUUGAAUGUGAUUUUUAGAGAUCUACACAUUUCAUCAACAGAGUGUCUCAUAAUUAGAUCCUGUGUCAAAACAUUUGUUGAUGGAAGGGCUGGGACUUGGUGGCUCACGCCUGUAAUCCCAGCUCUUAGAGAGGCCGAGGCAAGUGGGUCACUUGAGGUCAGGAGUUCGAGACCGGCCUGGCCAACAAGGCGAAACCCAUCUCUACUAAAAAUACGAAAACUAGCAAGGCGUGGUGGCAGGCCUGUAAUCCCAGCUACUGGGGAGGCCAAGGCAGGAGAAUUGCUUGAACCAGGGAGGUGGAGGUUGCAGUGAGCCAAGACUGCCACUGCACUCCAGCCUGGGCAACAGAGCAAGACUCCAUCUCAAAACAAACAAAACUUUUGUUGGCUGGGUGCGGUGUCUCAUGCCUAUAAUCCCAGCACUUUGGGAGGCUGAGGCGGGCAGAUCACCUGAGAUUAGGAGCUUCAGACCAGCCUGGCCAACAUAAGGAAACCCCAUCUCUACUAAAAGUAAAAAACUAGCUGGGCAUGGUGGCAGGCUCCUAUAAUCCCAGCUGCUGGGGAAGCUGAGGCAGGAGAAUUGCUUGAACCUGAGUGGCAGAGGUUGCAGUGAGCCAAGAUUGCACCACUGCACUCCAGCCUAGGUGACAGAGGAAGACUCUGGCUCAAAACAGACAAACAAACAAAAAAACUUUUGUUGGCAGAAGGUUGAAAUUUAGGAAAAUCUAUAAUAAAUAUAUUUGUUUUUCUCCUAGCUUUGUAGCCUUAAGUCUUAUUUCUUACUUUACAAAGUAAACCCGUAGGCUAGAUCAGUGUUUCUUAAAUCCAGUAUUCCUAAACUGUGCACUGAAGUAAAUCCUCACUUAAUGUUGUCCAUAAGUACUUGGAAACUGCAACUUCUAAGACCCCAAACACCUCUAAUAUUAAACAUUGAAAUAAAUGGGUCCUAUGCAUACAUUUAAGAAAUACUAAAAACAAGAUAAUUACCCACUUAUUUCAGUUCAGGAUGGCAGGUGGCUGGAGCCUGUCCCGGCAAGUAGGGCACAAGAUAGGCUGCCAUUCCAUUGCAGGUAGCAUUUACAGACACACACACACCCACCCACCCACCCACAUUCACUCAGACUGGGACCAUAUAGACAUGCCAAUUCAUGUCAGUUCACCUAAUGCGUGCACCUCUCUGAUGUGGGAGGAAACUGGAGCAUACGGAGAAAAGCCACACAGACAUGGGGAGAACAUGCAGGCUCCACAUAGAUGGUAGCCCUGGUGGGGAAAUCAAAAUUUUUUUCUCAUCAAAAUUGUGAAGAAAUGAUAUUGAAUGAAACUAUGUUAUGGAGGACCUGCCGUACUGUGGAGCACCAUAGUAACUCACCUGGACACUAUGCAUUAUUUGAGAAAUUUUAAUGGAAACAGCAAUAUUGGAUAGCUGUGGGAUAUUUGACAUCUAUAGGACACUACUCAAACUACCAGUUUGAGAUACCACACAAGAGAUUGUGUUCCAUUUUUUUUGAUGACACUGUGUCUUCCAAACGCUAGGUUUUGGGAGUUUGCUGUGAUAAAAAGCAAGUACUAUGCAAAAAUCAAUUUGGAACAGUAAAUGAAGGUGGUGGUGUCCAAAAUGAUUUCAAGAUUUGAGAACUUGUGUAAUAUCUAGUAGGUACACACAUCCUAUUAAGUAGUUGUGGUAAUUUAGUAAUAAAAAUAUUUAUUUCAAUUUAUACAUAUUUAUUUUUCGAUCUGCUACUAAGUUGGUAGGCUAUAAAUACUUUUUAGGUUGUGGACUGAACUACUUAGUAAAUGGAACUGUUAGGUAUUUUUUUGGCCUAGGGGUUCCAUGAAAAAAAUAACUGGCAUACUGAGGGGUCUGUGAACUGAGAAAGUUUAGGAACCUCUUGAACUUUUGACCAGGACACCUAUAUAUAUUUUCCAUGUUCUGUGCUCCCAAAGGUUGAUUUUUUACAUUUAUCUUUAUUUUUUAUUUUUUUAAAUAAAAAGAAUAGAUGUGUGUUAUACGUGUAUGUGUAUACAUAUAAUGAAUCAAAAGUUUCACUUUUCAAAGUUUUACUUUAAACAGCCUCAUGGCCAGGCGCGGUGGCUCAUGCCUAUAAUGCCAGCACUUUGGAAGGCCGAGGUAGGCAGAUCAUUUGAGGUCAGGAGUUCAAGACCAACCUGGCCAACGUGGUGAAACACCAUCUCUACUAAAAAUACAAAAAAAUUAGCUCGGCGUGGUGGCGCGCACCUGUAAUCCCAGUUACUCAGUAGGCUAAGGCAGGAGAAUCACUUGAACCUGGGAGGUGUAGGUUGCAGUGAGCCGAGAUCGCGUCACUGCACUCCAGCCAGGGCGACAGAGAGACUGUGUCUCAAAAAUAAAAUAAAAUAAAAUAAAGCGCUCACUACCUUCAAUACACUUGAUAUUUGCAAAUCAAUCCUAUUUAGUUCCAUUAAAAAAAUGCUAGCUCCAAUCCACGAUAUUGUUUUGUUUUGUUUUAAUCUACUAACGGGUCAUGACUACACACUUGAUAAUAUGAUGCUGGUUGACUCCUAACGUCCAGCUUUAAUACGGUUUGUUUCUCUGGUUAGUUAUAAAUGAAGUUAGACUACGCUUUAGCUUUUAUUUUAGAUUUGGGUUGUGGGGGAUACACGUGCAGGUUCGUUACAUGGGUAUAUUGCAUCAUGCUGAGGUUUGGGCUUCUGUUGAACCUAUCACUGAAAUAGCGAACAUAGUACCCAGUAGGUUGUUUUUCAACCUUUACCCCUUACCUGUCACUUCCCACCUUUUUUUUUGAGAUGGAGUCUCUGUCACCCAGGCUGGAGUGCAGUGGCUUGAUCUAGGCUCAUUGCAACUUCUGCCUCCCAGGUUUAAGCGAUUCUCCUGCCUCAGCCUCCCGAGUAGCUUGGACUACAGGUACGAGCCACCAUGCCCGGCUAAUUUUUUGUAUUUUUAGGAGAGAUGGGGUUUCACCAUGUUGGCCAGGCUGGUCUCAAACUUGUGAGCUCAAGUGGUCUGCCCACCUUGGCCUCCCAAGGUACUGGUAUUACAGGCGUGAGCCACCACGCCCGGCCCUUUCCCACUUUUGGAAUCCACAGUAUCUUAUUUUCAUCUUUAUGUCUAUAUUUACCUAUUGUUUAGACCCCACUUACAAGUGAGAGAAGGUGUGAUAUUUGAUUUUCUGUUUCUGCAUUAAUUUGCUUAGGAUAAUGACCUCCAGCUGCAUCCAUGUUGCUCCAAAGGACAUGAUUUCAUUCAUUUUUAUGGCUACGUAGUAUUCCACGGUGUAUAUGUACCACAUUUUCUUUAUCUAAUCUACUGUUUUUUGUUUGUUUGUUUGUUUUGUUUUGUUUUUGAGGCAGAGUCUCGCUCUGUUGCUCACGCUGGAGCGCAGUGACGUGAUCUUGGCUCACUGCAACUUCCGCCUCCCGGGUCCAAGCAAUUCUCCUGCCUCAGCCUCCUGAGUAGCUGGGAUUACAGGCACAUGCCACCAUGCCUGGCUGAUUUUUGUAUUUUUAGUAGAGAUGGAGUUUCACUAUGUUGGUCAGGCUGGUCUCAAACUCCUGACCUCAUGAUCUGCCCACCUCGGCCUCCCAAAGUGCUGGGAUUACAGGCAUGAGCCACUGCUCCCGGCCUAUCUGAUCUACUGUUGAUAGGUACUUAGUUUGAUUCAAUGACUUUGCUGUUGUGAAUAGUGCUGUGAUAAACAUACUAGUGCAGUGUCUUUUUGAUAGAAUGAUUUCUUUUCCUUUGGAUAGAUACCCAGUAGUGGGAAUGCUGGGUCUAAUGAUAGUUCUAAGUUCUUUGAGAAAUCCCUCUACUCUUUUCCACAGGGGUUGAACUAAUUUACAUUCCCACCAGCAGUGUAUGAUACUCCCUUUUCUCUGCAUCCUUGUCAUCUGUUAUUUUUUAACUCUUUAAUAAUAGCAAUUCUGACUGGUGUGAGAUGGUAUCAGUGUGGUUUUAAUUUGUAUUUUAUUGAUGAUUAGUGAUUCUGAGCAUUUUUUUCAUAUGUUUGUUGGCCACUUGUAUGUCUUCUUUUGAGAAGUAUCCGUUUAUGUCCAUUGCUCACUUUUUAAUGGAGUUAUUUUUUUUUUCUUGUUCAUUUGUUUGAGAUUCCUUGUAGAUUCUGGAUCCUAGUCCUUUGUCAGAUGUACAGUUUGCAAAUAUUUUUUUCCCAUUCUUUAGGUUACCUCUUUACUUUGUUGAUUCUUUUGCUAUGUAGAAGCUUUUUAGUUUAAUUAAAUCUUGUUUAUGAAUUUUUGUUUUGUUGCAUUUGCCUUUGAGGUUUUAGUCAUAAAUUCUUUGCCUAGGCCCAUGUGCAGAAGAGUUUUUCCUAGGUUUUUUAAAAUAGGAUCUUUAUAGUUUGAGGUUUUAUGUUUAAGUCUUUAAUACGUCUUGAGUUAAUUUUUGUAUGUGGUGAAAGGUAGGGGGUCCAGUUUCAUUCUUCUGCACAUGCUAGCCAGUUUUUUUGUGUGUUUUUUGGGUUUUUUUUUCUUUUUAGCACUGUUUAUUGAAUAGGUGUCCUUUCCCCAUUUUUGCUUUUGUCAGCUUUGUUCAAGAUCAACUGGUUAUAGGUGUAUGGCUUUAUUUCUGGGUUCUCUAUUCUGUUACAUUGAUACAGAUGUCUGUUUUUGUAUUGGUAACAUGCUGUUUUGGUUAGGAUAUCCUUGUAGUAUAGUUUGAAGUUGGGUAAUGUGAUGCCUCCUGCUGUGUUCUUUUUGCUUGGCUAUUCAGGCUCUUUUUUGGUUCCAUAUGAAUUUUAGAAUUGUUUUUUCUAAUUUGUGAAAAAUGAUGUUGCUAAUUUGAUGGGAAUAGAAUUUUUGGAUUGUUUUGGGCAGUAUGGUCAUUUAAAUUUUGAUUCUUCCAAUCCAUGAGCAUGGAAUGUUUUUCCAUUUGUUUGUGUCAUCUGUUAUUUCCUUCAUCAGUGUUUUGUAGUUCUCCUUGUAGAAAUCUUUCACCUCCCUGGUUAAAUUUAUUCCUGGGUAUUUUAUUUUUCAUGUGGCUAUUAUAGAUGGGAUUGCAUUUUUGAUUUGGUUCUCAGCUUGAACAUUAUUGGCAUAUAGAAAUGCCACUGACUGUUGUAAGUUGAUUUUGCAUCCUGAAACUUUACUGAAGUCUUUUUCUUAGGUCUAGGAGUUUUUUGAAGGACUCCUUAGGGUUGAACAGAGAUAAUUUGAGUUUCUCUUUUCCUAUAUGGAUGCCUUUUUUUAAAUUUUAAUUUUCUGCCUGAUUCCUAUGGGUAGGAUGGAAAGUAUUAUGUUGAAUAGGAGUCGUAAGAGUGCACAUCCUUGUCUUUUUUCAGUUCAUAGAGGGGAUGCUUCUAACUUUUGUCCAUUCAGUGUGAUCUUGAUUUUUUGGAGUAGUUUCAGUAAGAUUGAUACCAGCUCUUCUUUGUUUGUCUGGUAUAAUUUGGCUGUGAAUGCAUCUGAUCCUGGUCUUUUUUUGGUUUGUAGAUUUGUUUUUAUAACUGAUUCAAUUUUGUAAAUUGUUAUGAUCUGUUCAGGAUUUCAGUUUCUUCCUUGUUCAAUCUUGGGAGGUUGUGUGUUUCCAAGAAUUUAUCCAUUUCCUCUGGUUUUUCUAGUUUGUGCACAUAAAGAUGUUCAUAGUAGUCUCUGAGGCUCUUUUGUAUUUUUGUAGUAUGAGUUCUGUCACCUUUGUCAUUUCUGAUUAUGCUUCUUAGGAUCUUCUCACUUUUUUUUUUUUUAAUUUACUUAGAGGUUUAUCUAUCUUGCUUAUCCCUUCAAAGAACUAACUUUGUUCUGUUGAUCCUUUGUAUGUUUUUUUUGUUGUUGUUUGUUUGUUUGUUUGUUUUGGUCUCAAUUAAUUUAGUUCUGCUCUGAUCUUAGUUACUUCUUUUUUUUCUGCUAGCUUUGGGUUUAGUUUGUUCUUUUUCUGGUUCCUUUAAGUUCAACUUUAAGUUGUCAAUUUGAUAUCUUACUGUCUUUUUGAUGUAGGCAUUUAGUGCUGUAAACUUUCUCCUUAACAUUGCUUUUGCUGUAUCCCAGAGGUUUUGGUAUGUUUGUUUCUCUAUUUUCAUUUGUUUCAAAGAAUUUUUUUAUUUCUGCCUUAAUUUGUUGUUUAACCAAGUUAGUUAGGAGCAAGUUGUUUAGUUUCCAUGUAUUUGUGGUUUUGAGAGUUCCUUUUGGUAUUGAUUUCUAAUUUUAUUCUACUGUAGUCUGAUAAUAUGCUUGACAUGAUUUCAUUUUUUAAGAAUGUAUUGAGACUUGCUUUAUGACUGAGCUUGUGGUCCAUCUUAGAGAAUGUUCCAUGUGCAGAUGAGAGUGUAUAUUCUGUGAUUUUUGGGUACAGUAUUCUGUAGAAGUCUGUUAGGUCCGUUUGGUCAAGUGUCUAAGUCCAGAGUUUCUUUGUUAGUUUUAUGCUUUGACGUUCUGUUUCUAAUGUAUAUAGGCACAUACAUUUACAAAUGGAAAUCUCCUUUAUAGAUGUAAAUUUCUUUUACAGAAGGAUUUCAGAAUAGCCAGCUAAAUGCCAGAAAGGUUUAUUUUAAUUGGACAAGUGAUCUUUUUUCUUUUUCUUUUCUUUUUUUUUUUUUUUUUUUUUUGAGAUAGAGUCUUACUCUUGUCACCUAGGCUGGAGUGCAGUUGUGCAAUCUCGGCUUACUGCAACCUCCGCCUCCUGGGUUCAGGCGAUACUCCUGUCUCAGCCUAUUGAAUAGCUGGGAUUACAGGUGCCCGCCACCACACCUGGCUAAUUUUUUUUUUGUACUUUUAGUAGAGACGGGGUUUCACCAUGUUGGCCAGGCUGGUCUUGAACUCCUGACCUCAGGUGAUCCACUUUCCUUGGCCUUCCAAAGUGCUGGGAUUACAGGCAUGAGCCACUGUGCUCGGCCUACAGGUGAUCUUUUUAACCUAGCUACUAUUUCUUAGCUGAAAUUACUAAGUUCAGGAUGGAGCCCAUUAAGGAAUAGGGCAAAGAAAGCCUUCUCUGUACCUGGACUCAGCAAAAAUAGAUCUGAAAAAGGAGGAAGCUUACUUUACCUGGCAACCUACCUUUUAUAAACGCUAUCUAGGAUAACUUUCUUUUCACCUUCAGGGAAGAGUAGUAACGGAGCUGAAAGAUUAGCAGAUUUAAUUUUUCUUAACAAUUAGUUGCUUAAGCUUUUUAUUUGCCUUCUAUAAAGAGUCUUUUUAUAAAAAGCAAUAAAAAUAUUGAAAUCUUUUUAGAAGCUUCUGCACACCAACAGGCAUCCCUAGAUGAGACUGAUUUGGGAGCCCUCAUUUUCAAAUGUACUUCUUAAAGUGCAGUGUUCAUUUGGAAUGUACCGUUGUGAUUUUAAAUUACCUUUAGUAAGAUUUCAUCAUUUUUGUAAGCAUUUGCUGCUUCUGAGGCCUAAUAAUAAUGCAUGUAUAAGCUGGAAGGUGUAGUACUCAGGUCUUCAGAAAUUAAGGAUCCCAUUUUUACCUUGAAUAGCAGCUUUGGUCCUUGGAUCCCUUUGAUCAACUUAGCCAAUGAUUUUUCCCUCCCUAAGCAUGCAAGAAAAAUAAACAAAGGGGGCUAGAACACAAAAAUCUAUGCAAAUUUCUGAAAGCUGAAAUUGGUACCCCCUUCAGUAUUAUCAUUUAUACUAAUUUCUGUCUAAACCAGCCAGACAUAUGAGGCCUCUAACUAGAUCCAAGACAGUUAGAUACAGUUGGAUCCUGGACCCAGUCCAAUUUCUGUUGUAACUUCCAAACCUAGUUUGGAUCAGAAAUUUGCUCAAACUCUGAUAGCUCAAAUGCACACUGAUGGAGCUUCAGAAUCUGAGAGAGAACUUACCCACAAUCUCCAGUUGCUUUGAGAGAGCAGUGGACAAAAUGAGUCCAGUGGGUACCUUGCUUUGUCACUUGGUGUUCCUGGGAGUUGCUAGAAGCUUUAGUUUGGAUCCCACUUCUGACACCAUCUAUUAAAAUUAAAAACCUUAGACAAAUUAAAUUUAACAGAAUUUGAUUGAGCAAAGAAGAGUUCAUAAAUCAGCCAGCUUCCAAACCAGAAUAGGUUCAGAGAGGCUUCUUAGACUAUACUUUAUACAAUAGAUCAUCUCUGAAAAUGUGUACAGUUGCUCUAUUUUUGUAUGAUAAGUUAUUAAAUGUUAUAGAGCAAUAAUGAUGUAAAGAAGACCUUUGUACAAAUAGAAGGAUAUUUUUCCAAAGAAAAUCUGAUUUAUCUUUGAUUAGUUUGGGUUUUGGUGAUUUUUUUAAAAAGCCAUUUAUAGUUUUUGUUGUUGUUAGUGGGCUUGUAGAAGUUGUUUCUGGAAUUGUGGCACCCUCCUUUUUAAAAACUAUAUCUAAAAAUCAGUGAUAAGUAAUAGUAAUAGCACUGGCUAUAACUUAAUGUUUCAUUUCUUUUCUUUUUUUGAGACUGAGUCUCACUCUGUCGCCCAGGAGGGAGUACAGUAGCGUGAUCUUGGCUCACUGCAACCUCCGCCUUCCGGGUUCAAGCAGUUCUCCUGUCUUAGCCUCCCAAGUAGCUGGGAUUACAGACAUGCACCACCACGUCCAGCCAAUUUUUGUAUUUUUAAUAGAGAUGGGGUUUCACCAUGUUGGCCAGGCUGGUCUCCAACUCCUGACCCCAGGUGAUCGCCCACCUCGGCCUCCCAAAGUGCUGGUAUUACAGGUGUAAGCCACCGGGCCCAGCUGGCUUCAUUUCAUUCAUUAUUUCAUGUCUAAAUUUUAAAAGUGAAAAUAUUUUAUGAAAGAAUUUAUUUUUAAUUUUGAGAUUGUUUUAAAAUAUAAAAUAUGUAUGCACACGUAUACUUAUUUUUUUACGUGUUACCGUGUAAGAUAUUUUCUGCUUAGGAAUGUCUAACCUUUUCAUGUUAUAUGAAAUACCUGGCUUCAUAUGUGUAUUUGUGCCUUUUUGCAUUUAAGAAUUGGAUUCUGAGACUUAUCAAUGGUACUGUCAUUUGGCAAACAGGUCUGGAGCACUCAAUGUGGCUAUAUACUGUGCGUGAACUAGGGAUAUAAAGAGGAAUAAGAGCAGGUGCUGUGGCUCAGGCCUAUAGUCCCAGCACUCUGGGAGGCCAAGGCGGGUGCAUCACCUGAGGUCAGGAGAUCGAGACCAGCCUGGCCAACAUGGCGAAACCCCAUCUCUACUAAAAAUAUAAAAAUUAGCUGGGCAUGAUGGUAUGUGCCUGUAAUCCCAGCUACUAGGGAGGCUGAGGCAAGAGAAUCGCUUGAUCCCGGGAGGUAGAGGUUGCAGUGAGCCGAGAUUGCGCCACUGCACUCCCGCCUGGAUGACAAGAGCAAGACUCUGUCUCCAAAAAAAAAAAAAAAAAGAGGAAUAAAAUAUCUCUGCAUCUCUGCAGUUUUAUGUAAGAGACAGAUAUGUAAAUAAUAGUACAGAGGGGUAAAGAGAAUUGAUAACAUUCUGUGGAUGGAUUCAUAAAGAGGAAGCAAGUAUACUGCAUUAAGGAGUCAUUAAAGGCUUCCAAUGAAAGAAAAAUUGUUUCCAUUUCCUUAAACAGCCUCUCAGUUGACUAAGUUAUUGAACGUGGAAGGAGACUUUUGAAAAUUAUCUUGCUUUCCACGCUACCCAGAGAGGGGUGCAUACAGUGUUGUUCUGGAUUCCCAUUGUAACUUAAAGGGAAACUUUUACAAUAUCCAGAGGCCUUGAUGUCCUUAAGUUCCUUGCAGGAGGAACCCACUUAGAUGGCACCAGCCUUGACUUCCAAAUGGAACAGUACAUCUAUAAAAGGAAAAGUAAUGGCAUCUACAUCAUAAAUCUGAAGAGGACCUGAAAGAAGCUUCUGCUGGCAGCUCGUGCCAUUGUUACCUUUGAAAACCCCGCUGAUGUCAGUGUCAUAUUCUGCAGAAAUACUGGCCAGAAGGCUGUGUUGAAGUUUGCUGCUGCCACUGGAGCCACUCCAAUUGCUGACUGCUUUACUCCUGGAACCUUCACUAACCAGAUCCAGGCAGCCUUCCGGGAGCCAUGGCUUGUGGUGGUUUGACCACCAGCCUCAAACUGAAGCAUCUUAGGUUAACCUACCUACCAUUGCUCUGUGUAACACAGAGUCUCCUUUGCGCUGUGUGGACAUUGCCAUCCCGUGCAACAAUAAGGAAGCUCACUCAGUGGGUUUGAUAUGGUGGGUGCUGGCCGGGGAAGUUCUGUGUGUGUGUGCCACCAUUUCCUGUGAACACCCGUGGGAAGUCAUGCCUGAUCUCUACUUCUACAGAGAUCCUGAAGAGACUGAAAAAGCAGAGCAGGCUGCUGCUGAAAACGCUGUGACCAAGGAAGAAUUUUAGGGUGAAUGGACUGUGUCAGCUCUGGAGUUUAACACUACUCAGCCUGAGGUUGCAGACUGGUCUGAGGGCGUGCGGGUGCCCUCUGUGCCUAUUCAGCAGUUUCCUACUAAAGACUGGAGUGCUCAGAAGACUGGUCUACAGCUCCCACUGCUCAGGCCACUGAAUGGGUAGGAAUAAGCACUGAAUUGGCCUUAAGCUGUUCUUGCAUGGACUCUUAAGCAACAUGGAAACAGGGUUGAUGGAAAAUAAGCAGCAGUUUCUAAAAAAGACAAAAACAAAUUAUUUCUUUUUUCCUUGAAGUUUUCCUCAAAAAGUGAUCAGUGUGUCUGCUUUAAAAAAAGAAUUCAGAAACAGAAGCUGGGUGUGGUGGUGUGUGCUGUAAUCCCAGCUACUUGGGAGGCUGAAGUGGGAGGGUUGCUUGAGCCCAGGAAUUCAGCAAGACUAGCCUGGGCAACAAUGAGACCCUGUCUUCUAAAAAAUAAAAAAAAUAUACACUGAAAAAUCUACUUCCUAUCCUGGACCCUCACUUCCUCUCCUUGGAGUUAACUAAUAGUAACACUUUCUUAUGUAUAAGGGAAAAUUUUAUAGUGAAACUUUUAUACUUUUAGGCAAAAAAUAAUCCUGACUUGCAAACUGUAAGACACCUUUUUGUAGAUAUAUUCCUUGAAAUUAAUUUAGAAAAGUACUCAGGAUAACUAGAUAGUGAGUAAGCAAAAGUGUAGGCCUAUUUUUAGGGCUCCCCCCUUCUUUUGGUUAUUUCUGGUCCUGUUUGCUGCCCUCAGCUCUCCUAAAAUAAAAAACUAUGGGUGUCUUGUUUACGAAGUUUAGCUCUUUAUUUCAUCAGUGUGGGACAUUGUGUUCUCUUUUAUCCUUACUCUACUUCCUAUGCUGCAAAGGCCAAAAUAGUACUAAUGUUACAUAGUAUCUUCUGGAUUUCCUGGAUUACAGUAAGUCACUCAGUUUUUGUAAAAGUAGUUUUUGAGGUUUAUAUUUUAUAAGUAGAAACAGGUGAGACUAAAUUCUCAGAUUAGUGGUAUAUUAGAAUAUGUUAUAGGACUGUUGGUUAUUCUGGAAUUUACUUGGAACAGUGUAAUAGAUCUGAGUUAUAACUGAGUAUUUUCCACUUUUUUGUCUGCUGGUCCUUACUAAAUUUACCAUUUGUCCAUUUUAUGUUGAGGAUAGAUAGUAGGAAUAUUGCAGAGUUGUCUUUUCUUAACCUUUCAGUUAACUUUUUCUUAACCUAUCAUAGAUCUGAUAUAGAUCAGAGAACUGAAAUGUAGCAAAGGUGUGGAAAAUAUUCUAAGAGCCAGAAUGAAGUCAUUUUUUUCAUCUUUUGUCCCUUAAUGCCUGGCACUUAGUAGGAGGCCAAUAAAUAUUGAGCAAAAUGAAACUAUGUGAUUGGAUGAAAACACUGUCAAUACAUGUGAAAGAAAAAUUUGCUGUUAGAAGAUAAUUGAGAAUUGGCUAUUUCACAAUGUCAUCUUAUAUCUUACUUGUGUAUACACACCACAGGUAUAGAUCAUUGUUUCACAGAUGUUUUGGUAUAAGACUGCCUCAUGCAGUUGACUUAUGGAAGACUCAAAAAGAUUUUGUUCAUGUGGGAUUUAUCUGUUGAUGUUUACUGUAUUUGAAAUUUUCAGUGUACUUGUUAGUUUAUUUAAAAUAAUGAUAAACCUAUUACAUGUUAAAAUAAUAACUUUUUGUUGUUGUUGUUUGUUUUUUUUUUUUUUGAGACAGAGUCUUGCUCUGUCGCUCAGGCUGGAGUGCAGCAGUGCAAUCUCGGCUCAUUGCAACCUCUGCCUCCCGGGUUUAAGCAAUUCUCCUGCUUUAGCCUCCCAAGUAGCUGGGACUACAGGCACAUGCCACCACACCCAGCUAAUUUUUGUAUUUUUAGUAGAGACUGGGUUUCACCAUGUUGGCCAGGCUGGUCUCGAACUCCUGACUGCAGGUGAUCUGCCUGCCUCGACCUCCCAAAGUGCUGAGAUUAUAGGCGUGAGCCACUGCGCCCUGCCUUGUUUUGUUUUUAAGACAGGGUCUCACCCCAUUGCCCAGGCUGGAGUGUAGUAGUGCAAUUAUAGCUCACUGCAGUCUUGACGUUUAAGGCUCAAGUGAUCCUCCAACCUCAAUCUCCUGAGUAGUGGGGACUACAGGUGCAUGCCACCACACCUGGCUAAUUAAAAAUUUUUUUUUGUGGACACUGGGUCUCACUUUGUUGUCCAGGCUGGUCUUGAACUCUUGGACUCAAGCAAUCUUCCUAUCUUGGCCUCCCAAAGCGCUGGGAUUACAGGUGUGACCCACUGCACCUGGGCUUAACAUUUUAAAAUGGAAAAUAACUAACUGUAUCUCAACAAUAAAAAGCAGUGAGUGGCAUUGUUUGACAUUUUUCAACUCUUUUUAGUGUCUGGUUUAAUAGAAGACAUCUGGAUUCUCAUAUCUGCUUCUGUUUGCAUUCAGUCUGCUGUGAUAUCAUAUGUCAUGUGGCCUUUGGAAAACUCCAUUGUACAUUCAUAAGAGAAUGAGAAUAAAAAAAGACAAAGUCUUAAUAUUAUUAUGGAAAUAGUUUUGAUUUCAAAGAUGCUUCCUGAAAGAGUGUGAGAGACCCUCAAAGAUCCCUGAACCAUACUUUGAGAACUGGUAGUAUAGAUUAUAACUUCAGAGUAGAAAAGAACCUUACAGGUCAUAUGGUCAUCUCUUCUUGAGGUCACCUUUUCUUGUAGUUGUUUCUGGUUUAUGCAGCUGGAUGACUGAGAUAAGGAAUACCUCGUUUGUUUUUUGUGCGAUAAUGAGCUGGCUCUUGGAGAUGGUGAGUUUAAUGUGCCUUUCAGAUUUCUAAAUAAGCAGUUUGAUAUGUGGGUUAGGAUAUUAGAGGUAAAAUCUGAUGUUUUGGUAGUCAUUGGCAUGUAGGUGGUUUUGAAAGCCUUUUGCAAGAGUGAGGUCGCCCAGCGAGAGAGUAUAGUGUUAGAAGCAAAAAAGCUGAAGGUCAAAUCCUCCUGAACUUCAAGGUUUAAAAGGCAUAUGGAGAAAGACUAGCUAGUAAAGGAGACUGAGAAAGAGAACUAGAGAAGUGGGAAGAUGUAUGUAUUGUAUGUAUGGAGUGUUGUGUAUAGAAGUCAAGGGAAAAGAGUAUGUAUGGAAAAUGUGUGUAUGGAGUGUAGUGUAUAGAAGUCAAGGGAAAAGAGUGCUGUUUAUUUUAUUUUUUAUUUUUUGAGAGGGAGUCUUGCUCUGUCACCCAGGCUGGAGUGCAGUGGUGCAAUAUCAGCUCACUGCAACCUCCGACUCCCUGGUUCAAGCAAUUCUCCUGCCUCAGCCUCCUGAGUAGCUGGGAUUACAGGCACACGCCACCACACCCUGCCGAAAAGAGUGUUUUUAGAGAGUGGUCAAUUGUAAUGAAGCUGCCAGGAGGUCAAGUAGGUUGAGGAUUGAAAAUACCCAUUAGAUUUGGAGAUGAAGAGGUUUUUGGUGAUUUUAAUGAGAGUUAUUCAGAAAUGGGGAUGGCAACCACAUUGGAGUUGCUUGAAGAGUCAGAGUCUAAAUGGACACUGAGUGUAUUUAACUUUUUUUUUUUUUUUUGAAACAGAGUCUCAUUCUGUCGCCAGGCUGGAGUGCAGUGGCAUGAUCUCGGCUCACUGCAAGCUCCGCCUCCGGGGUUCACGCCAUUCUCCUGCCUCAGCCUCCCGAGUAGCUGGGACUACAGGCGCCCACCACCACACCCCUUUAAUUUUUUGUAUUUUUAGUAGAGACGGGGUUUCACCAUGUUAGCCAGGAUGGUCUUGAUCUCCUGACCUCGUGAUCCAUCUGCCUCAGCCUCCCAAAGUGCUGGGAUUACAGGCGUGAGCCACCUCUCCCGGCCGCAAUAGGUACAACAUUUUUAACAAGGCAAAUAUGUUUUCUAAGUUGUAGCCUAAUAAUCAGCAGUGGUACCCCUCUAGUAAUUCCUGCAAAUAAUCUAUGUACAUCAUUCCCUUUUAAGAAACACUAACAUUUUUUGUGAGGUUGAAGCCAUAGACCUAGUAGUCAAUCUUUGCUUUGCAUGGGUAGCUGAGAGAAAUGUCUUCCCUUCUGAGUAGUUUCUGCUGUUCAGUGAAUCUAGAUGGCAGUGUUUUCAUUGUUCCAGUUUAACUGUUAAUCACCCUGUUCCUGUCAGUGAGUGUAUAAAUGAGUAUCCAGUAACGUUUAUAUGUAUUUCUUAUACUUUUUUUUAGUGUAAUCAUAAUUAGGCUUUUACUUGGUGAUAGAGCAGUAGUAUUAUCUUUUACUCCAGUUUUCUGAUUUAUUGCACAUUAAUGAUGACAUGCCUUUAUCAAACGGCUCCCAUUAUUUCCAUGACAUCAUAUAUUUAGCUUAUUAUGUGACUACUUGUAGCGCUUUUCUUUGCCCUGUUAUCAGCGGGUCAAAUUUUGUGUGUGUGUUCCUUGUAGUUUAUUUUACUAUGUAUUCUUACAGGCUACUUGCAAAAUAGUAUGUGAAUUUCAUUUCCUUAUUUUUCAUUUCUUUUCAUUUUCUAUUUUUUCAAUUUCUUUUUCUCCCCCUCAUUUACUGCCGUAGUGAUAAAUAUGUUAAUUGUUCCUCUUUCUUUCUUAUGAGUGGUUUUGCAGCCUGGCCCUUUAACUACCUCUUCUGCUUCUGAGAUCAUGUUUCUUAACGCUUCUGUUGACCCAUCAUACCUCCUUCUUCACCAUUCAUCCUUCCACUGGCAGCUGUGCCUUUGAUUCUCCUUUAUUUUUUUCUGUUUUUUUUUUGUUGUUGUUGCUGGUGUUUGUUUGUUUGUUUUGUUUUGUUUUGAGAUGGAGUCUCGCUCUGUUGCCCAGACUGGAGUGCAGUGGCGUGAUCUCGGCUCACUGCAAGCUCUGCCUCUCGGGUUCAUGCCAUUCUCCUGCCUCAGCCUCCCACGUAGCUGGGACCACAGGUGCCUGCCACCACGUCUGGCUAAUUUUUUGGUAUUUUUAGUAGAGACAGAGUUUCACCAUGUUAGCCAGGAUGGUCUCGAUCUCCUGACCUUGUGAUCCGCCUGCCUUGGCCUCCCAAAGUGCUGGGAUUACAGGCGUGAGCCACUGUGCCUGGCCUAUUUUUUUCUGUUUAAAUGGAGCUCAUAACAGCUUCUAAUGCUGCCGUUUGUUCAUCUGUUCAUUGAAGGUUCUGUUGAUCUCUAGAUGGGUAAGAAGAUACAGUCCUUAAGGACCACACAACCUUGUGGGGAGACUGAUACAUAUAUGCUAUUUGUAACUCUAUUACAGCACUCAUUAUAUUGUAUACUAUUUGUUUAUGUGUGUCUCUCCUCUACUUGAUUGUGAGCCAUAAGUGCUCAGUGAAUGUUCAUUGAAUGAAUGAAUUAUGUAUGAUAUAUGUAUGUGUGUAUGUAUAUAUAAUGACAAGAACUAUGAAAAAUGAACAGAGGCCAACCAUGGUAGCUCAUGCCUGUAUUCCCAGUGCAUUGGGAGGCCAAGGCAGGAGGAUUGCUUGAGGCCAGGAUUUGAAGAACAGCCUGAGCAACAGGAGAGAUCUGUCUCUUAAAAAAAUAAGAUAAAUUCGUUGGACAUGGUGGCGCAUGCCUGUAGUCCUAGCAACACAGGAGGCUGAGGCAGGUGGAUCGCUUGAGCCCACGAGUUUCAGGCUACAGUGAGCCAUGAUUGCACCACUGUACUCCAGCCUGGGUGACAGAGCAAGACCGUGUCUCUAAAAAAGAAAAAAAAAGAAAAGGAAAAAAUGAAGAGAAUACUAUAGGAAACAGAGAAGUACUAAUUGAAGUACUAAUUUUAUUAAGAGCCAGUUUGGGAAAGAUUUAUAGGAGAUGACAUUAGAGUUUGGCUUUGAAGGGCAUUCUUGCCAGAUAAAAUUUAGGAAAAGACAAUGUGGAGAUAAUGAAAUUUUCUGGAUAUUCAGAAAGUGUACUCUAGAUAUACCUAUUCGACAGGAAAGUAGUUAUACAGUUUUGCUUAAAAAAAAAUAAAAAGAUGAGACUGUAAAUUAAUGUUGUGAAGGUCUCCAAGUGUCCAUCAUGCUAAGGAGUUUGGAGCUCUGUCCCAUAAAUGAUAGUCAUUGAAGAGUUAUAAAAUCAGAAUUGCUAGGUUCUGUUUUAGAACGAUAACUUUGAUUGCAGUAUGGAGAAUGGAUUUAGGAGCAGGAUAAUUUAGGAAGCUUUUACCAUAAUUUAAGCUAAAGAUGAUGCUUUGAACAGAGGCAAGGGGAAGAAAGAAAGAAAUUAGAGAAACAUCUUAAUUUUUUUUUUUACAAAUAGCUUUAUGUAUGAGGCAUAUUUUACAUAUAAUAAAGAAACAAAGGUAUAAAUGAUAUGAUAACUAAGUAUAUGAUGGAACGUGGGAGAGAUAGGAAUGUCAAAGAUGAAUUGGAUGCAUGUUAAUCCAAUUUAUUAAGAGUACAAGAAGAAAAUUUUGAAGAUAGAUAAUGAAUUUAUUUGGUUUUGAACAUGAAAAAUUUGAGGUGCCAAUGGUGCUUGGUGAUAGACAUCUAAUAGACUGAAAACUUGCAUCUAGAGCUCAGAAAUAGAUCUAAGGAGGAAAUAGAACAUUUAGGAAAUAAAGCCAUGAAAGUACGUGAGGUUGCCUGGGUAGAAGGAGCUAGAGAACAAGAGGCCCAGGAAAAAGCCUUGGAAAUGUUACCAUUGAUGAAUUAGACAGAAGUAGAGAACCCAUCAGGAGUUAAGAGAGGUAAGAGAAAAAAACAAUCAAGAGGCAAGGGAAACAGGAGAAAUAUGUUUUACGAAUCCCAAGAGUAGUGAGAGUUUUAAGGUGAAAGUUUAAGAUCCCAAAUAUUCAUAACUUAGUGUUGUAUGCCUGUUGACUGGAAAAUUGUAUUACUAGAAUGUAGGCAUAUUUCCUAGGUUUUCUAGACAGAUACAAUUCUUUUUUAAUGUUUAAAUUGUUCUAAGACUUUCAUAAGAAUAAGCCUGUAAAUCAGAAAAAAAAAUAGCUUUUUAGACUAUACUCCAAUUUUUGGUUUAAAAAACAUAGACUUUUUAGUUGAAGGACUGUUCAUAACAAAGUGAAAGAUAAAAGCAACAUUUUAUCAGUGGGCUUAAAAGAAAACAUUUUCUGAUGACCAGAGAAAACCCUGCACACAUCUCUGAGAGGAACUGUCAUUCAUUUUGGGGUAUUUUCUCUAUUCUGUCUCUUAAAAACAGUACAUUUCUUAACACUUUUUUUGUAGUAUAUUUGAGAUUCAUUUCCAGUAUGCUACAAAAAAUGUGUUAAGAAAUAUAUGUCCCAGCUACUUGGGAGGCUGAGGCAGGAGAAUGGCAUGAACCUGGGAGGUGGAGCUUGCAGUGAGCCGAGAUCGCGCCACUGCACUACAGCCUGGGCGACAGUGAGACUCCAUCUCAAAAAAAAAAAAAAGAAAUAUAUGAAUGUAAGAAAUACAAGAAAUGAAAUAUAGGAAAUGAAUUUCAAAUAUAUUUUUCAUUAAUAGCUUUUUUCUCUCUAAGUACAUGAUUUAGAGGUUUGAAAUGAAUGAUACUAUUUUGUUGGAUCCUCUUUUUCUGUGGGUUUUUAUCUUUUAUUUUUCCCCCGUAAGGAAACUCAAAUAAAGGGAAAUUUGUUCUAAGACUUCAGAGAGGCCUCAAAGAAUGCAAGAGCAGAAAAUGAGGUACUAGCAUUAUGGGAACUGGAAAAUUGUCAAGAACAAGACCACUUUUUCCUUUUCUUUGGGGCCACUUUUAUUUGUUUUUCUUUGAUCUUCCUCUUUGUUUUUCAUUGUAUUUUUCAGCUUGGCUCCUUAUAGCUCAUGUUCCGUUGCUUAUCCUUUCAAUUCCAAAACUCUGCAGUCUGUUUUUCCGGGACUGGAUAGGAUGCUUUUGAGUCGCUUCCUUAGGCCUUAUCCAGUCAGCUAGAGCUGGGGUGAGGGGCAGGUCUUUGUGGAUCACAGUAGUGGAUAGGGGUAAUUCCCUAGAGAAAAGGACUAGGGGACUAUGGGUGAUUGAUAUCUCUAGUAUUUACACAGCUUCAUUUAACAAAGCUUUUUGGGAACAUUUUGAAGGCAUAAUUUGUGAUCAGUUGUUUUGUGCUUAGUGGGAAGGUGGCCUUUUUGACUUUUGAAAAAUGUUAGAGUUGAAGGAAGAUGAACAUGGUAAAUACUGAGUUCAGGCUAGAGGUUAUUCUUUCUGGAUGAGGGGAGAGAGAUGAGUAUGAGCAGAGUAUACACAGAGGACUCCAGCCAAAUUGGUAAUAGUUUAUUCCUUAAACCGAGAGGUAGCUACAUGGGUGUUUAUUAUAUUAUUUCUUAUAUCACUUGAUGUGAAGAAGAAAGCUAUUUGCAGAUAAUUGAAGAUCAUUUGAAUUUAGACCCUGGCAGACAGUAAGCCUCAGAAAAUAAAGCAGUCCAAAAGGCUACUGAAGGCUUAGGACAUUACUGUACACUAUGUAGACUUUAUAAACACUGUGCACUUAGGCUACACUAAGUUGAUUAAAAAAUCUUUCUUCAGUUAUAACCUUAGUUUACUGUAACUUUUUUUACUUUAUAAACUUUAAUUUUUUAAAACUUUUUGACUCUUUUGUAAUAAUUCUUAGCUUAAAACAAAUAUACAGCUGUACAGCAAUAGUUUAUAUCCUUAUUCUGUAAGCUUUUGUUCAAUUUAAAUUUUUAUGUUUUUACUUUUUAUACUUUAAAAAUAUUUUUAUCUAUAAAAAUAUAUUUUAUAAAAUAAAAUAUGUAUAAAAUAUGUAUUUUUAUUUUUAGUAAAUAGCAGGAGCACACUCUAUAAUAACAAUAAAAAGUAUGGUAUAGUAAAUCCUAGGCGAUAGAAAUUUUUCAGCUCCAAUAUAGUUUUAAGAGACCACUGUCAUAUAUGUGGUUUGUGUUGACCAAAACAUCAUUGCAUGGCGCAUGACUGUACAGAGAACUUUGGCUAUCUUCCAUUAUGAAAGUAGAAUGUUGUGACAAGCUUAUUACUAUAUAACAGAAUGCUUAACUUUAUCUUUUAAAUUAAAAAAAAAUUUAAAUUGACACAAUUGUUGCACAUAUUCAUGGGGCAAUAGUGAUGUUUUGAUACAAUGUGUAGUGAUCAGAUCAGGGUGAUUAGCACAUCCAUCAUCUCAAAUAUUUAUCACUUGUUUGUGUUGGGGACAUUCAGUAUCUUCUCUUCUGUCAACUUGAAAAUAUGCAAUAUGUUCUCGUUAACUGUGGCCAUUCUACAGUGCUAUAGAACACUAGAACUUAUUCCUUCUAUCUAGCUGUAAUUUUGUAUCCUUUAGCAAAUCUCUCCUAAGCCCCCUUAACCCCCCUCCCAGCCCCUACUAACCUCUGUUCUACUUUUCUAUGAGAAGAGCUCUUCUAUGAGCUUCUGCAUAGGAGUGAGAACAUGCUGUGUUUAGCUUUCUGUUCCUGGCUUAUUGACUUAACAUAACGUCCUCCAGGCUCAUCCAUGUUGUUGUAAAUGACAUAAUUUUUUUCUUUUUAAGGCUGAAUGGUAUUCCAUUGUGCAUAUAUAUUUUUCUUUUCUUUUUUAUGGAAAUGUGAAUGGUUUUGACUUAAAAUGUGUUUUGCCUAAUAUAAUUGUGACCUCCCGUGUUCUCAUUUGGUUACUGUUUACAUGCAGUUUCUUUUUUCAGCCUGCUAACUCCAGUCUAUUUUUGUCAUUAGAUCUAAAGUGAGUCUCUUAUAGACAGGAUAUAGUUCCAUAUUUUCUUCAUCUCUUCAUCUGUUGUCGGACACUUGGGUUGAUUCCAUUUCUUGGUUAUUGUGAAUAAUGCUGCUAUAAACAUGGAAGUGCAGAUAUCUCUUCAGUAUACCGAUAUUCUUUCCCUUGGAUAAAUAUGCAGUAGUGGGAUUGCUGGAUCAUAUGACAGUUCUAUUUGUAAUUUUAUGAGGAAUCUCCAUACUGUUUUCCAUAGUGGAUAUACUAGUUUACAUUCCCACCAAUAGUAUGUAAGAGUUCCUUUUUCUGUAUAUCCUCGCUAGCAUUUGUUGUUGUCUUUUUGGUAAUAGCCAUUCUAACUGGAGUGAGAUGAUAUCUUAUUGUGGUAUUGAUUUGCAUUUCCCUGAUGCUUAGUGAUGUUGAGCUUUUUAAAAACAUGUUUCUUGACCAUUUGUGUUUGCCCAUUUUUAAUUAGAUUAUUUGGGGUGUUUUUUUUUUUUUUUGCUCUUAAGUUGUUUGAGUUUUGUAUGUUCUAGAUAUCAAUUCCCUGUCAGGUGAGUAGUUUGCAAAUAUUUCCUCCAAUUUUUUAGGUUGUUUUUUCAUUCUGUUGAUUGUUUCCCUUGCUGUGCAGAAGCUUUUUAGGUUGAUAUAAUGCCAUUUGUUAAUUUUUGCUUUUGUUGCUUGUGCUUUUGGGGGUCUUGUUCAUAAAAUCUUUUACUAGACCAGUGCCCUGAACUGUUUCCCUUAUGUUUUCUUCUAAUUAUUUUAUAGUUUUGGGUUCUUAUAUUUAAGUCUUUAAUCCCUUUCAAGUUGAUUUUUGUGUAGUGUGAGAGAUGGGGGUCUGUUUUCAUUCUUCUGCAUAUGGAGAUUCAGUUUUCCCAGCACUAUUUAUUGAAGAGACUAUCCCCAAUGAAUGUUCUUGGCACCAUUGUAAAAAAUCAGUUGGCUAUAAAUACAUGGAUGUUCUCUGUUCUGUUCCAUUGGUCUAUAUGUCUGAUUUUAUGGCAGUACUAUGCUGUUUUGGUUACUGUAGCUUUGUAGUAUAUUUUGAAAUCUGGUAGUGUGAUGCCACCAGGUUUGUUUUCUUUGCUCAGGAUUGCUUUGGCAAUUUGGUGUCUUUUGUGUUUCCAUAUGGAUUUCAGGAUUGUUUUUUCUAUUUCUGUGAAGAAUACCAUUGGUAUUUUGAUAAGGAUUGUAUUGAACCUGUAGAUCACUUUGGGUUGUAUGGUCAUUUUAACAAUAUUAAUUUUUCCAGUUCAUGAACACAGGGUGUUUCUCUUUUUUGUGUGUUUUCUUUAACUUUUCUCAUCAGUGUUUUACAGUUUCCUUUGUUGAGAUUUUUUACCUCCUUGGUUAAACUUAUUCCUAGGCAGAUUUUUUUUUUCUUUUUUUUUGUAUUUACUGUGGAUGGGAUUGCUUUCUUGAUUUCUUUUUAAGGUAGUUCAUUAUUUAUAUAUAGAAUCACUACUGAUUUUUGUAUGUUGAUUUUGUAUCCUACAAAUUUACUGAAUUCAUUUAUUAGUUCUAAUUUUUUUUUGGUGGAGUCUUUGGGAUUUCCUAUAUGUAUGAUUAUAUUGUCUGAAAACAGGGACUAUUUGACUUCCUCCUUUGCAAUUUGGAUGCCCUUUAUUUCUUUCUCUUGCCUAAUUACUGUGGCUAAGGUUUCCAGAACUAUGUUGAACAAAAGUGGUGGAAGUAGGCAUCCUUGUCUUGUUUCAGAUCUUAGAGGAAAAGCUCAUUCAUGGAUUACUAGGUUAAGAAUACUUGUACUAGAUUAUGAGAUUGUAGAGGUCAGAACCAUUUCUCUCUCUCUCUUUCUCUCUCUCUCCCUCUCAGUUUUUUCAUCUGUAAGGACCAUUUCUUAUUUUUAAUUUUUUUAAUAACAUACAAUGCUUUGCACAAAGUAGGUACUCAAAAAACAUUGUUAGCACAAGCAGUUGAUAUGCAAGAUUACUUUCCUCAUUCUUAGUCCUAAAAGACAAAUAUUAUUACAUUGGGUUAAAAAAUUCUAUUUUCUAUAUACAAUAUAUGUUUAAAAGAAAGUGAUCCACAAAGCUGAGAAUUAAAAGGAAAGGGCAGAUUCAAUAAAGCAAUGUUUGCUUUUUAUUAAUAGCAUGCAAAGUAGAAUUAAAUACUAAAACACAUUUCUCUAGGGAUAAACAUAGUCAUUUUUUAUUGAUGAAAUAUAAAAUCCAUUGACAGCAUAACCAUGAACUUUUAUAUACCAAAUAAUAUCAAGUUAUUUAGAAGUUAAAUAUUUUAGAAAUAUGAAAUGGACAAAAAUUCAGUACUUCCAGUGAAUUUUUGUCCAUCUCUUAAUUCUUGGCAAAUUAAAUAGACAAAAGAUAAUACACAAAUGUCCAACAAAUAAAUAAAUCCAAAUUUCCCAUAUGACUGUAGAAUAAAGUUAUUGUCAUGAUAUCUCUAAAACAUUUAACUGGGAAGGUGCUUUCAUUAUUCCGCAAGAAGUGAUCUUGAAUAAUCUGUUUGACUUGGUAAACACUUAAAACUGGUUGGGAAAAAGGCCUUAUGCAAAUUGAGAUUUAUACUUUAAAGAUAGUAAACCUUUUGAUAGAUUAUUUUAUUAUAAGUAAGUUUUUUAUUGCACUUUUUUCACUAAUUCAAUUUUCAGUUUCUACAAAAAUGGUGAUUUAAAAAAUGUAUUCUUAAAGAACAGACCAACUUUAACACAUUAUGUUUAAAAUAUUCACUUACAGGGCAUGGGAUAACCAAAGAUUUUGAAAGAUAGAUUAUUGGCAAAUUCAGGUGCUUCAUGAGAAUAGAUAUCAUAAAACUGAUUUAAAACUAUAAUUAUAGGCUGCAUCUGGUGGCUCACACCCAUGAUCCCAGCACUUUAGGCUGAGGUGGGAGGAUCACUUGAACCCACGAAUUUGAGGCCAGCCUGGGCAACAUAGUGAGACCUCAUCUCUAUAAAAAUAAAAAAAAAAUUAACCGGCCAUGGUGGUGUGUACCUGUAGUCUCAGCUACUUGGGAGGCUGAGGUGGGAGGAACACCUGAGCCCAGGAGGUUGAGGCUGCAGUGAGCCUUGAUUGCACUACUACACUCCAGCCUGGGUAACAGAGUGACACCUUGUCUCUUAAAAAAGCCUCAAAAACUAGAAUUAUAAACCUAAAACAGGGCAUAAUUCAUGUUAUCUGAAAAUCCUUUGAAGUAUAAUAUACAUUGAAAAAGCACAUAAAUGUGCAAGUUGAUAAAGUCUCUUAAACUGAACGUACCAUUGAAAUCAGCACCCAAAUCCAUGAAUAUUAGCAGUACCUCAGAAGCACUUCCUUGCUCCAAUUUCAGCCACCUCUCUAGCCACUUCAAGGGUGAUGACUAUCCUCACUUCUAGCACUAUAGAUUCCUUUUGCCUAUUUUAUACUUUAUUUAAUGGAACUAUAUGGUAUAUACUCUUGUUUCUGGCAUUCUUCACAACAGGUUUGUGAUAUUCAACUAUCUUUUUGCAUGUAGUUGUAAUUCAUUCUCAUUAUUAGUGGGGGAUACUCUACAGUAUGAAGUUUUACAUAAAUGUCAUACUUUUAAUAUUUAGUAAUUUAUAUUAAAGAUGAUUGCCCAUUAAUAGGGUAUUUAUUUUGAAGCAGAUAAAUUUUCUCUUGUUUAUUCAUUUAUACUGUCUCAGCCUUAGUGUUAGCUAAGAAGACAUAGAGAUUGUGGACACAUAUGUAGUUUUUGAUUAUUUAUUAUAUCCUUUCCAUUGGAAAAUUGAUAAUCACUAUAUAUUUCAGCACUUAUUUUAGUUUUAAAGCAAAAAGAUAUCUGCCCAAAGAAAUUUUUUGUCUAAUCAUUGACAGAUAUUCUAGUAAUAUAACUGAUGUUUACUUAAUUUUUCAGUUGGCUAGAUUAAAAAUAUUUUCAUCAAAAGGAAAAUGAGUUUUUUGUUGUUGUUGAGAUGGAGUCUUGCUCUGUUGCAAGGCUGGAGUCUUGUGGCGCGAUCUUGGCUCACUGCAACCUCCGACUCCCUGGUUCAAACGAUUCUCCUGCCUCAGCCUCGUGAGUAGCUGGGAUUACAGGCACCCGCCACCAUGCCUGGCUAAUUUUUGUAUUUUUAGUAGAGACAGAGUUUCACUGUGUUGGCCAGGAUGGUCUUGAUCUCUUGACCUCAUGAUCCACCUCCUGCGGCCUCCCAAAGUGCUGGGAUUACAGGCGUGAGCCACCGCGCCCAGCCAAAAUGAGUUAAUUUCUUAAACUGUGUUGGCAGAUUAAAGUUUUAUGGUUGAUAAAAUGUUGUUUUUCUACCAUGCUCCUUUAUUUUCAGGAAACUUGUUUACUUUUCUUAGAGCCAUAUUCUAUAAUCUAUGAAAAAAUUCUAGUUUUGGAACUUUUUUUUUUGAGACGGAGUCUCGCUCUCUUGCUCAGGCUGGAGUGCAGUGGUGCGAUCUUGGCUCACUGUAACCUCUGCCUCCUGGGUUCAAGCAAUUCUUCUGCCUCAGCCUCCUGAGUAGCUGGGCUUACAGAUGUGCGCCACCAUGCCUGGCUAAUUUUUGUAUUUUUAGUAGAGACAGGGUUUCGCCAUGUUGGCCAGGCUGGUCUUGAACUCCUGACCUCAGGUGAUCCACCUGCCUCAGCCUCCCAAAGUGCUGGGAUUACAGGCUUGAGCCACUGUGCCUGGCCCUAUUUUUAGAACAUUUUUGAGGUGUUUUCUAGGAAAUUUAUGGAUCAGGAUUAUUUACCCUGAUAAAUCAGAUUAAAUUAACUUCCUUAAAAAUUAACCAUACCUUUUAGUUUUACAUAUUUGGACUUAAUAGGUGAGCCCUUCUUUUUCUUAAAUUUAAAUUGGUGUUUUUAAUAUUCAACAUACCGUACUGCCUCUGCAUACAACUAAGUGGAGAAAAUAAUGUUUAUCAAAUAUGACAGGUUUAAUAUAUUUGUAAGGCAAAUAACUCAUUAACUUGAUAAGAAAAACAUUAAGAUCCUGGUAGGUACGUGGGCAAAGGAUAGAAAAGAUAAUUCACCAAUGAAGCACUAAGGUUAAGAUAUAGAGCAGAACGUUUAAUCCUGCUAGCAUCUACCAAAUUCAAAUUAAAGCAUAGAUUUUUUUCCUGUUAAAUCAGCAAAAGUAGAAAGAAUAGAUAAUCUUAAGAUAUUUAUGAAACAGGUAGGCCCUUCUUGGAAGACAAUUUGACAUUGUAUUUCAAGAGCAGCAUGUUCUUCACAUUCUUUUGUUCAGUUUUAGCACACUUUGCAGAUACUGUUCUAGGUAACCAAAUAAGCAGAAAUUCAUAAAUACAAAUAUGUUCAUUUCAGUGUUAUUUGUAAUAGAGAAAAUGGGACAUAAUAUUGUAUUGGUUAUCUCUUGUUGUCUAACUCAUAACCUCCCACCCCAAUUUAAGCAUUUUAAAACACCAAGCAUUUGUGAUCUCACAUUUUCUGUAUGUAAGGAAUUUGGGAGCAGCUCAUCUGGGUGAUUUUGGCUCAAUGUCUCUCAUAAGGUUGCAGUCAGGAUGUCAGGUAGAGCUGUAGUCCUGUGAAGAUUUGACUUGGUCUGGAGGAUCUACUUGCAAGGGGAUAGCUCACCCACAUAGCUAUAGGCAGGCGGCUUCAGUUCUUUGCCACAUGGACAUCUAUAGGAUUGCUUGACUGUUCUUAGGACAUGGCACCUGUAUCCAAAGAGAGAUGAGAAUGCGUUUUAUAGCCUAGUCUCAGAAGUCACAUAUUGUCACUUUCACCACAUGUGAUAUUAGUUUAAAAUAAGUCUCUAAGUCCCACCAUACUCAAUGGGAGAAUAAUUAAUUUCCACAUUUUGAAGAAAGGAGUAUCAAAGAAAUUUUGGAUAUAUUUUAAAACCACCACAAGUGUCCAGCAGUUAAGUAAAGUAGGGCACAAAAUUGAAGUUACAGUGAUGAUAGUAAUAAAAUUACUUAUUAAUGAUUAUAUUAGGAUUUAGGGAGUAAUGAGUUUUUAUUCUUAUAUUUUUCAAACUUUCUGUGUUGUUACAUGAUUUUUAUAAAAGACUAGUUCCUCUAUAUUCAUAGGAAAAUACAGUAUUACUUAUCACUGAGCUAACCAUGGGCAGAUCAUGAUCUUUCCCUUUGAGUUGUGAUAUGAUUUUUUUCCUUUCUUUUUUUAUUGUGGUAACAUAUACAUAACAAAAUUUACCAUUACAACCAUUUCUUAAGUGUACAGUUCAGUGGUGUUAAGUACAUUCAUAUUGUUGUACAACCAUCACCACCAUCCAUUUCCAGAUUUCUUUUCAUUUUGCAAAACUGAAACCCAUUAAACAGUAACUCUCCAUUCUUCUCUUUUCCCAGCCCCUGACAACCACCAUUCUACUUUCUGUCUCUAUGUUUUUGACUAUUCUGACUGUCUCAUGUAAGUGGCACCAUACAGUAUUUGUCUUUUUGUGAUUAGCUUAUUUCACUUAGUAUAAUGUCCUCCAGAUUCUUUCAUGUUAUAGCACAUAUCAGAUUUCCUUCCUUUUAAGGCUGAAUAAUAUUCUAUUUUAUAUAUAAAUGUCAUUUUGCUUAUUCAUUCAUCUGUUGAUAGACAUUUGGGCUACUUCCACGUUAGCUACUAUGAAUAAUGCUGCUUUGAACAUGGUUGUACAAAUAGUUCCUUGAGCCUCCACAAGCAUUUUGGUUGUAUACUGUAUUAGUCCAUUCUUACACUGCUAUAAAUUAUCCAAGACUGGGUAAUUUAUAAAGAAAAGAUGUUUAAUUGAUUUACAGCUCUGCAUGGCUGGGGAGGCCUCAGGAAACUUACAAUCAUGGCAGAAGGCACCUACCUCUUUGCAGGGUGGCAGGAGAGAGAAUUAGUGCAAGCAGGGGAAAUGCCAGACGCUUAUAAAACCGUUAGAUCUCACGAGACUCACUUAUUAUCAUGAGAACAGCAUAGGGGAAACUGCCCCCCUGAUCCAAUUACUUCUACCUGGUCCCGCCUAUGACACUUGGGGAUUAUGGAGAUUACAAUUCAGGGUGAGAUUUGGGUGGGGACACAGAACCAAACCAUAUCAUAUAUCUAGAAGUGGAAUUGCUGGAUCUUACAGUAAUUCUAAUUUUACUUUUUUUAGGAACUGCCACACUGUUUUCCAUAGUGAUUGUGUCAUACAUUCCCACGAACAGCACACACGGGUUCUAGUUUUGCAAAUCUUUGCCAGCUGUUAUUAUUUUGUGUGUUUUUGAUAGUAGCUUGUCCUAUUUACAUCCUAAUGGGUGUGAAGUAGUCAUAAUAUAAUUUUUAAUGAUUUAUUUUUCUAAGGAAAAGUAAUAUUAUAAGAUUGAUGUUCGGUGAAUACUGAUGUGAUUGCAGGUAAAUUGGCAGGAUGAGGAGCAGCAAAUCAAAAGAGGUGCCUUUACCAAAUCCAAGGAACUCUCAAAGCAAGGAUACUGUUCAAGGUAUGAUUUUGUUUUUUUAAACAGAACUUAAUACCUCAUUUAGGUCAUGUGUAUUAAAUUUUGGUCUAAUAAAUGAAAAUUAUACCUAUGCUUUGAUUUAUAGAAAGUGGAUAAAUCAUUUUAAAAUUCUCUUAUAUUUUAAAUAUACCUGGUGCUUCAUGUUAAAAGGCAUCCUAUGGUAAACUGUUUUGUAAAGCAGGAAUCCUUUUUUAAAAUGAACCAGUAUAUGUAAAAUUUUUCAGAAAUUCUGUCCUAAUUUCAAUAAAUGAUAUAGCCCCAAUAUAGUUACAGACUUAUAGGGGUACAGAUAAUAAUAGUAACAGUCUCAUUGAAUUCUUUUUAAAAUAAUGCAUUAGAGGUAAUUUUAAGGCAUUGAUAUUUUUUGAGGGAGAUUGACAAACUAGAGUACAUCUGGAGGAAGAUUAUUAGGUUUGUGAGAGAUCUGGAAAUCAAGUCAUAUGAAGAGUAGUUGGAGGAAGUGGAAAUACUUAUCCUGGAGAAGGAUGACUUAUGGAAAAUACACGACCUUAAACUAUUUGAAGGGUUGUUAUCUGAAAGACUGAAUGGAUUUUACCUAUGUUUUUACUGGUAGAUCAAAAUCUGGAAUAGCUUCUAAUCUUUUUCACAUUGUAUUAUAAGAUACAGUACAAUAUAAAAUUUUAUAUUUGUUCAGCUCAUUGGGGUAAACUGAGGAUGCUGCUUUUUUUUGGAGCUGUCUGAUACAGGUGUGAGUGGAAAUGCCAGUCCAAUGGUAGCAGUGUCUGUUUACCUUUAAACCGUAUGUUGGGAAUGUUGGCUUAAUAAUCUUAUAGAUAAAAACCCAUUUUGACUCAAUAUAAUAAAGCAUUUGCUGAGGAAUUGUUCUGCCUAAUAAGGUAGUUUUCAGACUUUUGAGUGCAUUUGGGAGCUUAUUAAAAUUCAGAUUCUGAGGUCUAAUUCAGAGAUUCUGGUUGAUUGAUAAAGUAAAUAUCAUUUUUUUUCCCAGCAGAUAUAACCACAUCGUGGGAUGCACUUUCUCAAACCAAGGCUGCUGUAAGUAGUUUUAGCUUCCAUUUAUUUAUUUCUGGAUUAAAAUCAAUUUCCUGUUAUGAUUGAAUAGAUUAAAAUGUUUUCCAUACUGUUCUUUAAAUACUAUAUAUUACUAUAGAAUAAUACUAUAUAUUAUUAUAUAUACUAUAUAAAAUAUGUAUACACUAUAUAUAAUAUAUAUACUAUAUACAGUAUAUACUUAUAUACAUAGUAUAUAGUAAUAUUAUAGUAUUAUCAUAUAUAAUAUAUAAUAUAUUAUAAAAUGAACAUAAAUGACAACUAUUAAAUGGUAAAGUCAUAAACAAGAUACUACUUUAGGAAUUAGGAGCACAAGCAUCAUGCCAUAUAAACUAUUAAAAUUGGCUGGUUGCGGUGGCUUACUCCUGUAAUCCUAGCACUUUGGGAGGUUGAGGUGGGUAGAUAGCUUGAGCCCAGGAGUUUGAGAUUAGCCUGGGCAACAUGGCAAAGCCCUGUCUCUACCAAAAAACAAAAUUUAGCCAGGCCUGGUGGUGUGCACCUGUAGUCCCUGCUACUCAGGAGGCUGAGGCAGGAGGAUUGCUUGAGCCCAGGGAGGUUGAGGCUUCGAGGUUGCAGUGACCCUUGAUUGUGCCACUGCACUUUAGCCUGAAUGACAGAGCAAGACCCUGUCUCAAUAAUAAUAAUAAUAAUAAUUUGUACUGUUUAUUCACCUUUAAUAAAUUUGUUCAUUUGGUAUGUGGAAAAGAGAAUUUAGACUGAAUGCAGUGGCUCAUCCUGUAAUCCCAGCACUUUGGAAGGCCAAGGUGAGAGGAUCGUUUUGAGCCCAGGCAUUCAAGACCAGCCUGUAAAACAUGGCAAAACCCUGUCUCUACAAAAGUACAAAAAUUAGCUGGGUGUGGUGGCAUGUGCCUGUGGUCCCAGCCACUCGGGAGGCUGAGGUGGGAUGAUUACCUGAGCCCAGGAGCUCAAGACUGCAGUGAGCUAUGAUCACACCACUGCAUUCUAGCCUGGGCGACAGAGUGAGAUCCUGUCUCAAAAAAAAAAAAAAAAAAUUUAAAACAGAGGAUGAAGGUAGGGCUUUGUCUAUGUAGAAAGCUCAUUGUUAGGGAGGAAAUGAAGGAAACCAACUUCUUGAAUAGACUGUGAAUUAUUCUUCCUUUCCAAUUUGGAUGCCCUUUAUUUCUUUCUGUGGUCUGAUUGCUCUAGCUGGGACUUCCAGUACUAUGUUGAAUAACAGUGUUAAAAGAGUUACUUCUUUUUUGUAUCAAUAACAGGUUGUUUUUUUUCUCUUACCACUGAAGAAGAAUGACUAUGGGCAACCUGUUCCCAUUUGUGUUCUGGCCAAUCUUGUCUAAAGUUUAGGAAUCCAAAUUUGGCAGGUGAUGCAUAGUUCAUUCAAAGGGAAACCCUUUUUUUUUUUUAAUUUGCUCUCUGGGCAAAUUACUUUCCUAAUAUACCCUUUUUUGUUUUCUAUGCCACCUGAUUUUGAAUCACCUAUUUCUUUACUAAGUAUUUCAGUGAUAUUUUCUGUCCUCUAAUCUUUAGAACUAGUUUGAAAAUAAAUCUGGAAAAGGUGUUAUAACCUAGAUUUGGUCCCAAUGCAGAGUCAUAUUUUGGCCCUAUAAAUAUAGCAAUAUAUAUUCUCUUUCUUCUACCUGAGCAUUAUGUUUUCUAUAGUGUAGAAAUAGUUUUUAUCCCCCAGCUUAACUUAGUGCUAAACUAAUUUAGAGUUUAGCACUAAGUUAAUAACUGUUAAUAGUUGUUUUUCAAAACUAUUGGGUGUUUUUGAAGAACAUUUAAGGAUGAUAAAGGAAUAUAAGUAGGGGAGACAAUCUGAUGUAGAACUAUUAAUAUAUCAUAAAUUAUUCCAAUUUUUAAAAUACAUUGGGAUGCAAAAAUCUGAAGGAAGGCUGGGCGCAGUGUCUCACGCCUGUAAUCGCAGCACUUUGGGAGGCUGAGGCACGUCAAUCACUUGAGCUCAGGAGUUCGAGAUCAACCUGGACAUUUGGGAGCUUAUUAAAAUUCAGAUUCUGAAUUUUAAACUCCAUCUCUACAAAAAAAAAUACAAAAAUUAGCUGGGCAUGGUGGUGUACAUCUGUAGUCCCAGCUAUUUGGGAGGCUGAGGUGGGAGGAUCGCUUGAGCCCAGGAGGCAGAGGUUGCAGUUAGCUGAGAUCAUGUCACUGCACUCCAGCCUGGGUGACAGAGUGAGACCCUGUCUCAAAAAAAAAAAAAAAAAUCUGAAGGAAGAUAUAUCAAAAUUUUAUAUAAAAUGUCUAUAUGUGAUGGAAUUUUGAUUUUAAUUUCAUUUUUUAUAUUGUUUCUGUUUUCCUAGUUUUUUGUUGUAAGCAUGUAAUACUUUACCAUUUAGGAGAAAAUGCUUUUCAGUCAUUGAGAACUUUGGUGGUUAUAAUCAAAAGUUGACUUAUAUUGCAGAGUGUGUAUGUGAAACCUGAAUUCUAUCUUGCGAAUUUUGUCUUUAUUGUGAAACUUAUUUUAUAGACAUAAAACUGUAAAACACUGUGGCUGUUUGAGUACCUAUAGUCCUAAAUUUUGGUUCAGGUAAAUAAUUCUUUCCUCAUUUUUGCGUGUAAAAUGCACAGGAUAUAUUGCAUAAUGGUAUGUCACAUAUAAAGUCAGAUAAUUGGCUGGGCAUAGUGGCUCACGCCUGUAAUCCCAGCACUUUGGGAGGCUGAGGUGUGUGGAUCACCUGGGGUCAGGAGAUCGAGACCAGCCUGGCCAACAUGUGGAAACUCCAUCUCUACUAAAACUACAAAAAUUAGCCCGGCGUGGUGGCGCACGCCUGUAAUCCCAGCCUCUAGGGAGGCUGAGGGAGAAUCUCUUGAACCCGUGAGGUGGAGGUUGCAGUGACCUGAGAUCACGCCACUGCACUCCAGCCUGCGCGACAGACCGAGACCGUCUCAAAAAAAAAAAAGUCAGAUAAUUUACUACUUUGAAAAUACUCAACAAUUUAAAAAAUAGAUCAACAUGUCAAAUGGAAACUCUUUGAUAUUCUUCAGAAUCUUACUGAUUUUCUUGUCUAAAUUUGGGGAGCAUUGUGCCAGAUGUUAUUAUAGCUAAAGUAAAAAAUGAUGAAAAACAAAUAUUAUAAGAAAGGAUUCCUGAGUUUUUUUCUACUUUGGAAAAAGGUAUUGGGCAAUAUACUUUUAAAUAACCAUGUAAAUUGAUGGUAGCUGGAUUAUUUACAGAAUUAACUUAUGUCUUAUGAUGGUGUUUAAACUUUAGCUGAGACACAUUGAAAAUAAAUUAGAAGUAGCCCCUACAAGUACAGCUGUGUGUGAUUCUGUCAUGGAUACCAAGAAGUCUUCUACAAGUGCUACUCGAAAAAUAAGUAGAAAAGGUAUGUAUGAAGUUACUUCCAAAUAUAUACAUAUAUUAAAUAUAAAUUUUGACAAAACAGACAGCUAGAUCUACAGGGGUGUGCCACCACACCUGGCUAAUCUUUUUUUAUCUUUUGUAGAGACGAAGUCUCUCUGUGCUGCCCCAGCUUGUCUUGAACUCCUGGCCUUAAGUAAUCCUCUCAUCUUGGACUCCUUUUAAAGUGCUGGGAUUAUAGGUGUAGCUACCACACUCUACCUUAAAUAUAAGUUGUUGAUGUAAUAGUUUUAUUUCAGUUUGAUUUCUUCUUACAGUCUCGUGUUGUUUGCUUUUUUGUCAGAGACACUGGCAUUUAAUUUUUUUCCUUUGAUAUUUGGUAGCUAUGAAGCAGCCUCAUGUUCAGUAACCGUAGAUGGUGUUUUCACCAUCUUUUGACCAGGGACUUUAGCACUAUUGCAUCACCCUUGACAAUGAGUUUUCAUUUGUUUUUAGUGAUGAGAGUUUCAUUUUGAGGUUAUAUAUGUCUGAUUUUUUUUUUUAAUUCACAGAUAUAUAGCACUUUGGAAAUGAAAUUUAGUUGACUUAACGUAUAUUUUUCUUGGGAUUUAAAAUCAUAAAUGGAAUCUCUUUCAUAUGGUGACUGUCAUAAAUUGAUACAAUGCAAAUUUUUCUUGGGUUACCUUUUGUAUUUUGAUUCUGUGGUUUCUUUCUUUUCUUUUUUUUUUUUUUUUUUUUGAGACAGAGUCUCACUCUGUCGCCCAGGCUGAAGUGCAGUGGCAUGAUCUCUACUCACUGCAACCUCUACGUCCUGGGUUUAAGCGAUUCUCCUGAGUAGCCGCGCCACCGGGCCUUUCUAAUUUUUGUAUUUUUAGAAGAGAUGGGGUUUCACCAUGUUGGUCAGGCUGGUCUCAAACUCCUGACCUCGUGAUCCGCCUGCCUUGGCCUCCCAAAGUGCAGGGAUUACAGGCGUGAGCCACCAUGCCAGGCCAUUCUGUGGUUUCUUAAUUGGUGGGUGAAAACCUCCAUGUCAGUACUCUCUACCUAAACUGUGAUAUGUAUAUAUAUAUAUGUGUGUGUGUGUGUGUGUGUGUGUGUGUGUGUGUGUGUGUAUAUAUAUAUAUAUAUACAUUUUUUUUUAACCUGGGAGUUUGUUAGGUUGUAUAUGAUACUUCUCUCUUAAGUGAGUGAGAGCCAUGACAAAAAAUAAAUUUUAUGUGUAGAAUUUUAUUUUGCAUUGUUUCCUCAGAAUACAAUUUCUGUAUCAAUGGAGUUACUUUAACAAUUUAACCUGAAUGUUUUUCAUGUACUUAUAUAAUAUCCAACUAGAGAUACUAGAGAUAACUUAUAUAAUAUCCAACUAUUAGAAAUAAUAGAUACCUAAUUUUUUUUUUAAGAUACCAGCCUAAUUCAGAGGCAGUUGUAUUUUAUAUUAUAUGUUCUCACAGAUUUCCUUUUCCUUCAGAUGGUAGAUACCUGGAUGAUUCUUGGGUUAAUGCUCCAAUCUCCAAAUCCACUAAAUCACGAAAAGAGAAAUCUCGUAGUCCUCUCAGGGCCACCACCCUGGAGAGUAAUGUGAAGAAAAAUAAUCGUGUGGAAUUUCGUGAACCUUUGGUUUCUUAUAGGUUAGUAUUGAGAAAAAAAAAAGGUAUCAAAUAGGUUUAGCCUGUAAUAUUUACAGUAAGAGUAUACGUUGUUUUCUGAGUUUUUUAUAGCACUCUGGUU';
		my $expectedCdsWtSeq = 'ATGAGGAGCAGCAAATCAAAAGAGGTGCCTTTACCAAATCCAAGGAACTCTCAAAGCAAGGATACTGTTCAAGGTATGATTTTGTTTTTTTAAACAGAACTTAATACCTCATTTAGGTCATGTGTATTAAATTTTGGTCTAATAAATGAAAATTATACCTATGCTTTGATTTATAGAAAGTGGATAAATCATTTTAAAATTCTCTTATATTTTAAATATACCTGGTGCTTCATGTTAAAAGGCATCCTATGGTAAACTGTTTTGTAAAGCAGGAATCCTTTTTTAAAATGAACCAGTATATGTAAAATTTTTCAGAAATTCTGTCCTAATTTCAATAAATGATATAGCCCCAATATAGTTACAGACTTATAGGGGTACAGATAATAATAGTAACAGTCTCATTGAATTCTTTTTAAAATAATGCATTAGAGGTAATTTTAAGGCATTGATATTTTTTGAGGGAGATTGACAAACTAGAGTACATCTGGAGGAAGATTATTAGGTTTGTGAGAGATCTGGAAATCAAGTCATATGAAGAGTAGTTGGAGGAAGTGGAAATACTTATCCTGGAGAAGGATGACTTATGGAAAATACACGACCTTAAACTATTTGAAGGGTTGTTATCTGAAAGACTGAATGGATTTTACCTATGTTTTTACTGGTAGATCAAAATCTGGAATAGCTTCTAATCTTTTTCACATTGTATTATAAGATACAGTACAATATAAAATTTTATATTTGTTCAGCTCATTGGGGTAAACTGAGGATGCTGCTTTTTTTTGGAGCTGTCTGATACAGGTGTGAGTGGAAATGCCAGTCCAATGGTAGCAGTGTCTGTTTACCTTTAAACCGTATGTTGGGAATGTTGGCTTAATAATCTTATAGATAAAAACCCATTTTGACTCAATATAATAAAGCATTTGCTGAGGAATTGTTCTGCCTAATAAGGTAGTTTTCAGACTTTTGAGTGCATTTGGGAGCTTATTAAAATTCAGATTCTGAGGTCTAATTCAGAGATTCTGGTTGATTGATAAAGTAAATATCATTTTTTTTCCCAGCAGATATAACCACATCGTGGGATGCACTTTCTCAAACCAAGGCTGCTGTAAGTAGTTTTAGCTTCCATTTATTTATTTCTGGATTAAAATCAATTTCCTGTTATGATTGAATAGATTAAAATGTTTTCCATACTGTTCTTTAAATACTATATATTACTATAGAATAATACTATATATTATTATATATACTATATAAAATATGTATACACTATATATAATATATATACTATATACAGTATATACTTATATACATAGTATATAGTAATATTATAGTATTATCATATATAATATATAATATATTATAAAATGAACATAAATGACAACTATTAAATGGTAAAGTCATAAACAAGATACTACTTTAGGAATTAGGAGCACAAGCATCATGCCATATAAACTATTAAAATTGGCTGGTTGCGGTGGCTTACTCCTGTAATCCTAGCACTTTGGGAGGTTGAGGTGGGTAGATAGCTTGAGCCCAGGAGTTTGAGATTAGCCTGGGCAACATGGCAAAGCCCTGTCTCTACCAAAAAACAAAATTTAGCCAGGCCTGGTGGTGTGCACCTGTAGTCCCTGCTACTCAGGAGGCTGAGGCAGGAGGATTGCTTGAGCCCAGGGAGGTTGAGGCTTCGAGGTTGCAGTGACCCTTGATTGTGCCACTGCACTTTAGCCTGAATGACAGAGCAAGACCCTGTCTCAATAATAATAATAATAATAATTTGTACTGTTTATTCACCTTTAATAAATTTGTTCATTTGGTATGTGGAAAAGAGAATTTAGACTGAATGCAGTGGCTCATCCTGTAATCCCAGCACTTTGGAAGGCCAAGGTGAGAGGATCGTTTTGAGCCCAGGCATTCAAGACCAGCCTGTAAAACATGGCAAAACCCTGTCTCTACAAAAGTACAAAAATTAGCTGGGTGTGGTGGCATGTGCCTGTGGTCCCAGCCACTCGGGAGGCTGAGGTGGGATGATTACCTGAGCCCAGGAGCTCAAGACTGCAGTGAGCTATGATCACACCACTGCATTCTAGCCTGGGCGACAGAGTGAGATCCTGTCTCAAAAAAAAAAAAAAAAAATTTAAAACAGAGGATGAAGGTAGGGCTTTGTCTATGTAGAAAGCTCATTGTTAGGGAGGAAATGAAGGAAACCAACTTCTTGAATAGACTGTGAATTATTCTTCCTTTCCAATTTGGATGCCCTTTATTTCTTTCTGTGGTCTGATTGCTCTAGCTGGGACTTCCAGTACTATGTTGAATAACAGTGTTAAAAGAGTTACTTCTTTTTTGTATCAATAACAGGTTGTTTTTTTTCTCTTACCACTGAAGAAGAATGACTATGGGCAACCTGTTCCCATTTGTGTTCTGGCCAATCTTGTCTAAAGTTTAGGAATCCAAATTTGGCAGGTGATGCATAGTTCATTCAAAGGGAAACCCTTTTTTTTTTTTAATTTGCTCTCTGGGCAAATTACTTTCCTAATATACCCTTTTTTGTTTTCTATGCCACCTGATTTTGAATCACCTATTTCTTTACTAAGTATTTCAGTGATATTTTCTGTCCTCTAATCTTTAGAACTAGTTTGAAAATAAATCTGGAAAAGGTGTTATAACCTAGATTTGGTCCCAATGCAGAGTCATATTTTGGCCCTATAAATATAGCAATATATATTCTCTTTCTTCTACCTGAGCATTATGTTTTCTATAGTGTAGAAATAGTTTTTATCCCCCAGCTTAACTTAGTGCTAAACTAATTTAGAGTTTAGCACTAAGTTAATAACTGTTAATAGTTGTTTTTCAAAACTATTGGGTGTTTTTGAAGAACATTTAAGGATGATAAAGGAATATAAGTAGGGGAGACAATCTGATGTAGAACTATTAATATATCATAAATTATTCCAATTTTTAAAATACATTGGGATGCAAAAATCTGAAGGAAGGCTGGGCGCAGTGTCTCACGCCTGTAATCGCAGCACTTTGGGAGGCTGAGGCACGTCAATCACTTGAGCTCAGGAGTTCGAGATCAACCTGGACATTTGGGAGCTTATTAAAATTCAGATTCTGAATTTTAAACTCCATCTCTACAAAAAAAAATACAAAAATTAGCTGGGCATGGTGGTGTACATCTGTAGTCCCAGCTATTTGGGAGGCTGAGGTGGGAGGATCGCTTGAGCCCAGGAGGCAGAGGTTGCAGTTAGCTGAGATCATGTCACTGCACTCCAGCCTGGGTGACAGAGTGAGACCCTGTCTCAAAAAAAAAAAAAAAAATCTGAAGGAAGATATATCAAAATTTTATATAAAATGTCTATATGTGATGGAATTTTGATTTTAATTTCATTTTTTATATTGTTTCTGTTTTCCTAGTTTTTTGTTGTAAGCATGTAATACTTTACCATTTAGGAGAAAATGCTTTTCAGTCATTGAGAACTTTGGTGGTTATAATCAAAAGTTGACTTATATTGCAGAGTGTGTATGTGAAACCTGAATTCTATCTTGCGAATTTTGTCTTTATTGTGAAACTTATTTTATAGACATAAAACTGTAAAACACTGTGGCTGTTTGAGTACCTATAGTCCTAAATTTTGGTTCAGGTAAATAATTCTTTCCTCATTTTTGCGTGTAAAATGCACAGGATATATTGCATAATGGTATGTCACATATAAAGTCAGATAATTGGCTGGGCATAGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCTGAGGTGTGTGGATCACCTGGGGTCAGGAGATCGAGACCAGCCTGGCCAACATGTGGAAACTCCATCTCTACTAAAACTACAAAAATTAGCCCGGCGTGGTGGCGCACGCCTGTAATCCCAGCCTCTAGGGAGGCTGAGGGAGAATCTCTTGAACCCGTGAGGTGGAGGTTGCAGTGACCTGAGATCACGCCACTGCACTCCAGCCTGCGCGACAGACCGAGACCGTCTCAAAAAAAAAAAAGTCAGATAATTTACTACTTTGAAAATACTCAACAATTTAAAAAATAGATCAACATGTCAAATGGAAACTCTTTGATATTCTTCAGAATCTTACTGATTTTCTTGTCTAAATTTGGGGAGCATTGTGCCAGATGTTATTATAGCTAAAGTAAAAAATGATGAAAAACAAATATTATAAGAAAGGATTCCTGAGTTTTTTTCTACTTTGGAAAAAGGTATTGGGCAATATACTTTTAAATAACCATGTAAATTGATGGTAGCTGGATTATTTACAGAATTAACTTATGTCTTATGATGGTGTTTAAACTTTAGCTGAGACACATTGAAAATAAATTAGAAGTAGCCCCTACAAGTACAGCTGTGTGTGATTCTGTCATGGATACCAAGAAGTCTTCTACAAGTGCTACTCGAAAAATAAGTAGAAAAGGTATGTATGAAGTTACTTCCAAATATATACATATATTAAATATAAATTTTGACAAAACAGACAGCTAGATCTACAGGGGTGTGCCACCACACCTGGCTAATCTTTTTTTATCTTTTGTAGAGACGAAGTCTCTCTGTGCTGCCCCAGCTTGTCTTGAACTCCTGGCCTTAAGTAATCCTCTCATCTTGGACTCCTTTTAAAGTGCTGGGATTATAGGTGTAGCTACCACACTCTACCTTAAATATAAGTTGTTGATGTAATAGTTTTATTTCAGTTTGATTTCTTCTTACAGTCTCGTGTTGTTTGCTTTTTTGTCAGAGACACTGGCATTTAATTTTTTTCCTTTGATATTTGGTAGCTATGAAGCAGCCTCATGTTCAGTAACCGTAGATGGTGTTTTCACCATCTTTTGACCAGGGACTTTAGCACTATTGCATCACCCTTGACAATGAGTTTTCATTTGTTTTTAGTGATGAGAGTTTCATTTTGAGGTTATATATGTCTGATTTTTTTTTTTAATTCACAGATATATAGCACTTTGGAAATGAAATTTAGTTGACTTAACGTATATTTTTCTTGGGATTTAAAATCATAAATGGAATCTCTTTCATATGGTGACTGTCATAAATTGATACAATGCAAATTTTTCTTGGGTTACCTTTTGTATTTTGATTCTGTGGTTTCTTTCTTTTCTTTTTTTTTTTTTTTTTTTTGAGACAGAGTCTCACTCTGTCGCCCAGGCTGAAGTGCAGTGGCATGATCTCTACTCACTGCAACCTCTACGTCCTGGGTTTAAGCGATTCTCCTGAGTAGCCGCGCCACCGGGCCTTTCTAATTTTTGTATTTTTAGAAGAGATGGGGTTTCACCATGTTGGTCAGGCTGGTCTCAAACTCCTGACCTCGTGATCCGCCTGCCTTGGCCTCCCAAAGTGCAGGGATTACAGGCGTGAGCCACCATGCCAGGCCATTCTGTGGTTTCTTAATTGGTGGGTGAAAACCTCCATGTCAGTACTCTCTACCTAAACTGTGATATGTATATATATATATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTATATATATATATATATACATTTTTTTTTAACCTGGGAGTTTGTTAGGTTGTATATGATACTTCTCTCTTAAGTGAGTGAGAGCCATGACAAAAAATAAATTTTATGTGTAGAATTTTATTTTGCATTGTTTCCTCAGAATACAATTTCTGTATCAATGGAGTTACTTTAACAATTTAACCTGAATGTTTTTCATGTACTTATATAATATCCAACTAGAGATACTAGAGATAACTTATATAATATCCAACTATTAGAAATAATAGATACCTAATTTTTTTTTTAAGATACCAGCCTAATTCAGAGGCAGTTGTATTTTATATTATATGTTCTCACAGATTTCCTTTTCCTTCAGATGGTAGATACCTGGATGATTCTTGGGTTAATGCTCCAATCTCCAAATCCACTAAATCACGAAAAGAGAAATCTCGTAGTCCTCTCAGGGCCACCACCCTGGAGAGTAATGTGAAGAAAAATAATCGTGTGGAATTTCGTGAACCTTTGGTTTCTTATAGGTTAGTATTGAGAAAAAAAAAAGGTATCAAATAGGTTTAGCCTGTAATATTTACAGTAAGAGTATACGTTGTTTTCTGAGTTTTTTATAGCACTCTGGTT';

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getCDSClass,$a->getExonClass,$a->getEssentialSpliceSiteClass,$a->getSpliceRegionClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1,0,813,100,$expectedRnaWtSeq,'G','r.1_813+100del37584insg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get2KBUpStreamVariantClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1,0,395,100,$expectedCdsWtSeq,'G','c.1_395+100del6140insG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass,$a->getStartLostVariantClass);
		done_testing();
	};

}
sub testRearHalfTranscriptRemoval_1_CEP350 {
	my $file = shift;

	subtest 'Testing CEP350 Rear Half Transcript Removal + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 180061649,
			'maxpos'				=> 180085015,
			'delseq' 				=> 'TTTTTTTCACTAAAATTTGGCTTTTCTCTAAACTTTGGTTCACAAAATTGAAAGATAGTGATTAAATAATTTTGATGTGGTGATTTGTTTTTTCTTTAAGATGCTTCACTGTCTGGTTCTGAGAGATCAGTATCAGAAAGGTCTTTATCTGCATATGCAAAGAGAGTAAATGAATGGGACAGTCGAACAGAAGATTTTCAGACCCCATCTCCAGTTCTCAGATCATCAAGGAAAATCAGAGAAGAATCTGGAGATTCTCTAGAAAATGTACCTGCATTACATCTTCTCAAAGAATTAAATGCCACTAGTAGAATTCTTGATATGTCAGATGGCAAGGTTGGAGAATCTAGTAAAAAATCAGAAATAAAAGAAATAGAGTATACAAAATTGAAGAAGAGTAAGATTGAAGATGCCTTTTCTAAAGAAGGTAAATCTGATGTCTTACTGAAATTAGTCCTAGAACAGGGAGATTCATCTGAAATTCTTTCAAAGAAAGATCTTCCTTTAGATTCTGAAAATGTTCAGAAAGACCTAGTTGGATTAGCTATTGAAAATCTCCATAAAAGTGAGGAAATGTTGAAAGAGAGACAGTCAGATCAAGATATGAATCATAGTCCAAACATCCAATCAGGAAAAGACATTCACGAACAAAAGAACACAAAGGAAAAAGATTTGTCTTGGTCAGAACATCTTTTTGCTCCTAAAGAGATACCATACTCTGAAGATTTTGAAGTGTCTTCTTTCAAGAAAGAAATTTCAGCTGAATTGTACAAAGATGATTTTGAGGTGTCATCTTTGCTGTCACTCAGGAAAGACTCTCAGTCTTGCAGAGATAAGCCACAGCCAATGAGGAGCTCTACAAGTGGAGCCACTAGCTTTGGTAGTAATGAGGAAATCAGTGAGTGCCTAAGTGAGAAAAGCCTTTCTATCCATAGCAATGTTCATTCTGACAGGCTGTTGGAACTCAAGTCCCCTACTGAGCTGATGAAAAGTAAGGAGCGCAGTGATGTGGAGCATGAACAGCAAGTTACTGAATCCCCTTCCTTGGCTTCAGTTCCTACTGCAGACGAGTTATTTGATTTCCACATTGGTGATAGGGTGTTGATTGGAAATGTTCAGCCAGGAATTCTTCGATTCAAAGGTGAGACTAGTTTTGCTAAAGGATTTTGGGCCGGAGTGGAGTTAGATAAACCTGAAGGAAATAACAATGGAACATATGATGGTATTGCATATTTTGAGTGCAAAGAAAAGCATGGTATTTTTGCTCCTCCTCAAAAAATATCTCACATTCCAGAAAACTTTGATGACTATGTAGACATTAATGAAGATGAAGACTGTTACTCAGATGAACGATATCAGTGCTATAATCAAGAGCAAAATGATACAGAGGGTCCAAAAGACAGAGAAAAGGATGTCAGTGAATATTTTTATGAGAAATCCCTACCTAGTGTGAATGATATAGAAGCCTCAGTTAATAGAAGTAGAAGCCTTAAAATAGAAACAGACAATGTACAGGACATTTCTGGGGTACTTGAAGCCCATGTTCACCAGCAGTCTTCAGTGGATTCACAGATTTCTTCAAAGGAAAACAAAGACCTCATTTCTGATGCCACAGAAAAGGTTTCCATCGCTGCAGAAGATGACACTTTAGACAATACCTTTTCCGAAGAATTGGAGAAGCAACAGCAGTTTACAGAAGAGGAAGACAACCTATATGCTGAAGCTTCAGAAAAGCTTTGTACACCACTTCTGGATCTTTTAACAAGAGAAAAAAACCAACTGGAAGCCCAGCTGAAGTCATCACTAAATGAGGAAAAAAAGTCAAAACAACAACTGGAAAAAATCAGCTTACTGACAGACAGTTTACTAAAAGTCTTTGTAAAGGACACAGTCAATCAACTACAACAAATCAAAAAAACCAGGGATGAGAAAATCCAGCTTAGCAATCAGGAGCTTCTTGGTGATGACCAAAAGAAAGTAACACCCCAAGACCTATCCCAAAATGTTGAGGAACAGTCGCCAAGTATTTCAGGTTGCTTCTTAAGTTCTGAATTGGAAGATGAAAAAGAAGAGATTTCCTCTCCAGATATGTGTCCCAGACCGGTGAGTATTTCGTCTAGAAAAAGTATCAGTATACCTTTTGACACTTTTAACAAGGCAACTTACCTTGCCTATATAGTACATTTTAATTCTGTGACTCCTGATAGGTTTTTTTTCCCTTGACGTGTTGTTTATTAGAGAAGGATCACTGTACTGGTTTTTAAAGGCCTTGTAGCTGTTTTTCCTGGCTGAAGTGTAATTGTCCTCTATTATGTCATGAGTAGTCAATATAGGGTTGTCACAAGCACTGGCAGGTTTTGGGGAACTAAAACATCCATTTAAACTGCTGAGTTGACAACCCCTGTACTGTCCATAGATTACATCTCCTTCCAGATTCACTGAACAAATTCGTTCCATGGAGCTTTCCAGTGGCTCACCACCCAATTCTTTTTTGCCTGTGAAATTGTTTACATATCTACTCAGTCATGGTTGGCCTGCAGAGGTTAGCAGTTGTGATGTTGATTCTGGTAAAATGAAATAGATCCACTATGACATGCACTCAGACCATCAGCTTGTTTATGTTGCCCATGGACAGTTGAGCAAAAGTAATGACTTTGAGTAACTTGAGCCAGATGGTCACCAGAGGTGAGAAATTTTATATCCCATTTCCATATCCCATAGACAGTTCTTTAGAAACAGAATAATGATAATACAAAAAGCGGGGAGAGCCCTAATCAATAAACATTTTATAAAAATTCAGTTATATGCAAGGTTTATATAAGTAGAATTGAGGGAACAGGATAACTTTTTCTTAGATTGTTTATTAGTAGCCACCATATTGTTTTTTATTCTGCTTACTTATAGAGAAATATTAGATATATAAATAATGCATTTACTCTGTCTTTTTTTAAAAAAAAATGATATGAGGTCCCTAAATGGTGCTCTGTTTGAATGGTTCCATTTCTATAGGAGAGCCCAGTATTTGGTGCCAGTGGGCAGGAAGAACTTGCTAAGAGACTTGCTGAACTTGAACTCAGCCGGGAGTTCCTGAGCGCGTTAGGAGATGATCAAGACTGGTTTGATGAAGACTTTGGTTTGAGCTCTTCTCACAAGATCCAAAAAAATAAGGCAGAAGAAACCATTGTACCTCTAATGGCAGAACCTAAAAGAGTAACCCAACAACCATGTGAAACATTATTGGCAGTCCCCCATACTGCAGAAGAAGTAGAGATTCTTGTACATAATGCAGCAGAAGAACTTTGGAAATGGAAAGAATTAGGCCACGATCTTCATAGCATCAGTATTCCTACAAAACTGCTTGGCTGTGCCAGTAAAGGTCTAGATATAGAAAGCACTAGTAAAAGGGTCTACAAACAGGTAGGTGAAATAAAAGGATAATTTTAGTTTTTAAACTTTTCTCCCATTGTCTGAAATTTGGATTCTTAGCCATTTTGTTTTGTTTTGTTTTGTTTTTCTCTATCAAGGCGGTTTTTGATTTAACAAAAGAGATTTTTGAGGAAATATTTGCTGAGGATCCCAACTTAAATCAACCTGTCTGGATGAAGCCATGTAGAATCAACTCTAGTTATTTCCGACGAGTGAAAAATCCAAATAACCTTGATGAAATCAAGGTAAACTGCAAACTATAAAGTGTCTTCTTTTTTGACTTGCTGTTCATTTACACATATCTGTAAATGTTGTAAAAGCAGGTGCTCTTCTGTATGATTAACCTTTGTGTCTTCAGGACCTGTTCCATAGTGGGTATTCAATAAAGATTTGTTGACTTACTGAATCCATGAACAGTGGAGGGATTGGAATGAGTGAAAGGGCTGAATCCTCAAGTTTGTGCAAAAAAAAAAAAAAATTCTTTTACCAGTCATTTTCTCTATATACATAATTATAATATTTAAGAAAAATTAGGAAATACCTAGGGAGATAAATTATTGCAACTAATAGGAGTTTTAGGTAAGTTTTGGGGGCTTTAGAGATAGGTATAATGTTCTAAAAGTACATATGAATACGTGTCATATATATGTCTTTTCAATATGCCAGCCCTGGACACATAGGTGTTAACTTCTTTATCTTTGACTTTGGAGGTTTCAAATTGGTAATGTGGTTTTGAACTATGTAGGATAGTAATGTGATACGTAAAAATATTTATTTTGTTATCTAAGAAACCTTTGTAAACCTTTGTCCAGAAGCTACAAGGACTATGAACTTTAAAGTCACATTGCTGGCAATCTGCTCTGCCTAGTCACATGTTTCTTTACGGCAAAGGTTATATTTCCAGAGCAGTGATTCTCAACTTTTAATGTGCATGTGAATGATCTGAGGATCTTGTTGAAATGGAGATTCTGATTTAGTAGGTCTGGAGTAGGACCTGAAATTCTGCATTTCCAACAAGCTCTCAAGTCATGCTGGTCCTGCTGGTCCTGCTGGTCCTGCTGGTCCCAGGACCACGCTGAGTAGCAAAGCTTTAAACTTGTATTTCCCAAACTTATTTGATCACAGAAAGCCTTTAAAAATATTTTATACTTATAGAACACCAACAGATATTTATGAGATGGGATGTATGCTTACCTCAGGGGTTGGCAAACTATGGCCCATAAGCCAAATCCAGGCCCCTGTATGTTTTTATAAATAAAGTTTTGTTGGACCACAGCCATACCTAGTTGTTTACATGTTACAACTACATGTAACTACAACAGTGGAATTGAATAGTTGTGACAAAGACCACATGGCCCGCAAAGCCTAAAATATTTACAGTCTGATTTTTTACAGAAAAAATTGGCCATCCCCTGGCTTAAAGTAACTAGCTTTTTAAAAACACATTGTATACTCCATTCCTTCTACCCTTTCACAATTGAAAATGCTTGAATGCTATGTGTGTCTGTGAACAGGTAGGGAGTGATGGAATTTCCGGGCACGTTGGTACTTTGAGAATCTTTCCAGGTGATTCTGATTCATAAATGTGTCCCTTCCCCTGTTAAGAACATTGTTTTTTTAACAGAACAACCATAAAAGTAACAAAATAAATTATTAATAACATTAAGTTGAGGTAAGACAAGATATATAATATATACCAACAATCAAACTACATGCATCAATATAGCATACAATAAATGCCCACTAATTTTAATATTGACATAACATTTTGTAAAGTTTTTAAGGGTGATAGTTCAAATATTGTTAGTAGTAAAATTATTATGCTTCTCATACTGAAGTGTTGATATTTTTTAATTTAACATATTTGAGTTTTGGTGCTTCAAAATAAACTATCCCAAAGTCTTTTGCTTTTCTCACATTACATTGTCCCCAAAAATTTCTTTGAGATTTTTCCAATTCAGAGATTTATTTTATGCTTTCTTTAGTAGTTTTGCCAGCTTGAAGTTAGCTAATGGTAGTTTTGTTGTTGTTTTTTCTTTGTTTTTCTTGAGGTGGAGTTTCACTCTTGTTGCCCAGGCTGGAGTACAATGGCGTGATCTCAGCCCACTGCAACCTCCGCCTCCTGGATTCAAGTGATTCTCCTGCCTCAGCCTCCTGAGTAGCTGGGATTACAGGCATGCGCTACCTCACCCAGCTAATTTTGTATTTTTAGTAGAGACAGGGTTTCTCCATGTTGGTCAGGCTGGTCTCGAACTCCCGACCTCAGGTGATCCGCCGCCTTGGCCTTTCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCACCCAGCATGGTAGTTCTTAAGTAAACTCTAAAGTGATTTTTCTTTGCAGATTAATCAGTGGGTTTTACTTTGTTGTTGTTGTTGTTGTTGTTATTTTTAGATAAGATCTCGCTTTGTCACCTAGGCTGGAGTGCAGTGGTGCAATCTCAGCTCACTGCAGTCTGGACCTCCTGGGCTCAATCAGTTCTCCCACCTCACCCTCCTGAGTAGCCAGAACTACAGTCATGTACCACCATGCCCAGCTAATTTTATATATTTTTTGTAGAGATGAGGTCTCTCACTATGTTGCTCAGGCTGGTCTCAAACTCCTGGGCTCAAGTGATTCTCCTGCCTTGGCTTCAGAGAGTGCTGGGATTACAGACGTGAGTCACCGTGCCCAGCTAGGTTTTACTTTTATATGTTGCGTAAGTGAATAAAGTCTCAGCATTTGGAATGCTTTTCCCTTCAAATTGCTATGTGGAACTGTAAGTGGAAGAAATTTGTATGATAATTTTATTGTAAATATTATATTCTTATCGTGAATCTGTTTACTGCATTGTAAACACATGAAACTCAGGAATTTCTGTTGTTTGTGTTTCGTGTATTTGTCTTGTTTCCACAGAGCTTCATAGCAAGTGAAGTACTCAAGTTGTTCAGTCTTAAAAAGGAGCCAAACCACAAAACAGATTGGCAGAAAATGATGAAATTTGGAAGAAAGAAAAGAGACCGAGTGGATCATATCCTGGTCAGTGTATACAACCAAACTGTTTTTATTTTGACCATATCTTTTAAACTCAGAATGCAGTTAGATTCTAAAAGGAAAAGGGATAGACTTATTCTTAAATCTATGGTATTAAAATCTTAAAGCCTTTACTTGCTGCATTTCTACCCCTCAGTTTTGTTAACTATTTTCTTTTGTACTCAGTAGATGAACATTTAACTTGATCAATAGTATTTTATGTAGAACTTGAGAAATTTACTCAAATGGAAGCAAAGTTTAAGGTAGTGAATTGAAAAGAGACCAGTGTGGACTCTGGAGCCTGACTGCCCATGTTTAAATCTCTTGATGTGTAATCTTGGACAAATTTGTCTAACATTTTTGTGCAAAGTGAAGGTAACAATAATAGTACTTACCTGATAAAATTATTAAGAAAATTAAATAAGTATATGTGTAAAGCACTTAGAATAGTTCCTCATATACAACTCAGTAAATAATAGCTGTTACCTCATCATCATCACTATTATCATTATCCAGAACATTCCCTGGATAGGATTTATTTGGTTTTGCAATATATTTCTGTAACATTAGGAAAGTGTTTATCAAAAAAGAAAAAAAGAAAAATTCAGGTTTCACTCTAAACATTTTAATTCTGGACCTTTTATCTATATTGACATTTATATAATGAAAAAAAGTGTAATACTTAAATTTATCACCGTTCTCATGTATGTACCTCTCTAGTGTTTCATTTACTCAGTAATAACTGTACTAGATAACATTTGCTTAGTATTTTATAGTTTATCATTCTTATGCCTTATTTGATTTTTTTTTTTTTTTTTTTTTTTGAGACACAGTCTCGCTCTGTCACCCTGGCTGGAGTGTGGTGGCACAATCTCTGCTCACTGCAGCCTCTGCCTCCCAGGTTCAAGTGATTCTCATGCCTCAACCTCCCGAGTAGCTGGGATTACAGGCACCTGCAACCATGTCCAGCTAATTTTTGTATTTTTAGTAGAGATGGGATTTCGCCATGTTGGCCAGGCTGGTCTCAAACTCCTGACCACAAGTGATCCGCCTGCCTTGGCCTCCCAAAGTGCTGGGATTACAGGCATGAGCCACCGCGCCCAGTGTCTTATGCCTTATTTGATCCTCACGATAATTACATGAGGTAGGAAGAACTAATGTTATCCTTTGTTTTACAAATGAGGTGTCTTACGCCTTATTTGATCCTCACAATAATTACATGAGGTAGGAAGGACTAATGTTATTCTTTGTTTTACAAATGAGGAAATTTAGATTCAAAGATGTTAAGTGACTTGCCCAGTAAAATAAGAGACGAGATTTTAATCAACTAATTCTTCATAGATCTTCACATCTCAAGAAAGTGAATTTGTTCTCCTATTTCTATCATGTCTTAAAGTCTGAATTGTTATTAATTTCCAAGTTGTTACTTCCTTTACTTCTTAGGAATTAAAGAAAAAAGCATATTAAGAGGCATTCATCTAACCTATCAGCTAAATTGATGATATCCCCTATAAATCACCTCTGTATATACATTGAGTCAAAAATTGTGCTAATTAAGGACTGATGTTAACAATAGTAAGCAAAACAGACACAGACCCCACCTTCATAGATTTCCAGTTTGAAGGCAGAAACCAGCCATTAATATATCGTTTACACACAAGTGTGTAAGTTATGCATAAAAATGGAAGAATAGGGAGCTGTTAAAACCACAGGAAGCTGATATAATGAGGATGGTCAAGGCAGGCTTCTCTGAGGAAGTGACATTTGAGCTCAGATCTTAAAGATGAATACATGGCAGCAGTAAGAGAAAATGAGGAGGAAGCAAAAGTGGAAACCCCTGATAAACCCATCAGATCTTGTGAGACTTAATTGCTATCACGAGAATAGCATGGGAAAGACTCCCATGATACAATTACCTCCCCCTGGGTCCCTCCCACACACGTGGGAATTCTGGGAGATACAATTCAAGTTGAGATTTAGGTGGGGACACAGCCAAACCATATCAGTACCCTGACAGCAGAGGCTTCCAGTAGGTGCCAAACAGGGTCAGTTAAGCAGACATCTCTCTAACAGTTATCAGAAATTCTCAAATACTAAACAAAATATTTCTTTTGCAGGGGTAAGGTTATTATAGCTTGGGTGAGGGGAGTACTGTGCATAAGAATTGAAGCAATTTAGCAACTGGTAAACAATAATAATAAAGTATATATTATGTGACTATCCCCTTTGCTGGTTGTAAAATATTTTTTAGCTTCCTACAGAGAGACCTGACCTATATAACTTTTTGCATGAGCTCGTTCTTCTTAATTTTTTTTTTTAACAAAATGATACTGTCAGCTTCTTGAGCTAATTATTAACAAGAAATGCTAAATGAACTAGTCTCCTCTCCCATAATTAGATTTGGTAAACCTTGTATTTGCAGATTCTAATGTTTTTACATTTTAGGGGGATTGAGATTTAAAACAATTTTTTAATAACACGTAAAAGTTTTCTCAAACCACTTCCTCTACTCTCACATTACAACAATCATCACCACACAAGAAAACTTCTGTGACCATATGTGTGTGGACTTTTCCCCCACACACCAAGCAGTGGACACCAGCTGGGTATCCTCTGAATCAGTTCTGATACTCCCTACCCAGATATAGTGTCACATCCCACAGAGTGAGTGCTCAGTCCCCAAGACTGCCCCTGCCCCCTCCACACACACATACCAGTTGCAAGTCCAGGCCTCCAGAACTTCTGATGGACTAGCTTCAAGTTGGGGCTCCCAAGATCCACCTTTGGGTTTAATTAATTTGATGGAGTAGCTCACAGAACTCAGTAAACACTTACTTAGGTTTACTGGTTTATTAGGAAGGATATTGCAAAGGATACAGAGACGCATAGGGTGGAGTATGGGAAAGGGGCACAGAACTTCCATGCCTTCCCTGGGCTGCCACCCTCCAGGAACCTCCAGGTTTAGCTATCCAGAAGCTACCTGAACTCTTTCCTCTTGGGTTTTTCTGAAAGCTTCATGACATCAGCATTCCTTCCCCCAAGGTATTGGGTAGGACCCTCTCATGGGAGGGCCTTAAGACCCACAGTCAGAAAGGTAGGGGAACATTAGAGTGAAAGGAGGGCGGGCCTGCCACTGAGGCCTAACACCCTTGACATTTTTTTTTTTTTTTTTTTGAGACAGAGTCTCGCTTTGTTACCCAGGCTGGAGTGCAGTGGCACCACCTCGGCTCACTACAACCTCCGCCTCCCGGGTTCAAGCAATTCTCCTGTCTCAGTCTCCGGAGTAGCTGGGATTAGAGGCGCATGCCACCATGCCTGGCTAATTTTTGTATTTTTAGTAGAGATGGGGTTTCACCATGTTGGCCAGGCTGGTCTCAAGCTCCTGACCTCAGGTGATCCACCTGTCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGGGAGCCACCGCACCTGGCCAACCCCCCACATTTTAACAAGACTATAACAAGGGCTATGGGAATTATTAGCCATGTTTTAAGTTTTTACTAAAAGAATGAGCTTATTTTGCTCTATGATTTTGTTACACTTTTCACCCATTATTGGTTATTTTTCTAAGTCTTTTACATGTTTAATTATGAAATATTTTTAAAGTATAAAAAAGAATAGAAAATAATATAGGTGCATTGTTTCTGCATCCTTTTCCTTCATTCTGAGCTCGGTTGCATCTACCTAAAGTACACACTTAGTATTTCATTAACTTCAGGCCTTTAATGTAAACTCTCAGTCTTTAAAATGTCTTTATTTCAGCTCTTACTTTTGAGTTATATTTTAGCTGTATATAGAATTTTAAGTTGATAGCTCTTTTTACTTAAGCAGTTTGAAGATATCTCATTATAATTTTGAAGCAGCATCATTGTCTGCGGTGAATACCCAGGGTTCGTCGTCTCGAGCTGAGAAAATTAACGACACAGACACACATACATGAAGTGGGTTAAGGAGTGGAAAGTTTAATAGGCAGAAGAAATGAGAGAGGAGAGCAGCTCTCTCTTTTATGAAAGAGAGGTGTTCCGAAAAAGGAAAAGCGGTGGACCGCAGCAGATTTTATAGGCAGGCTTTAAGAGGTGGTATCTGATTTACATAGGGCCCACAGATTGGTTCCACCAGGTGTGATGTTTACTGGTGTGTGGGGAAGGCTGGTCACCCCACCCTAATCTTATTATGCAAATAGATTTTCTACTTGGCCAGTGCCATCTTGCCTGCTCCTTACTGTGCATGTGGCAAAGAGAAGAGAAGATGGAGCTGCCATTTTGAACAGGCCTATTCCACAGGTAGTTTTTTCTATTGGCACAATTGCCGGCATTCGCTTGTGCAAGCTTCCAGCTTGTTTGTCTGTGTCTGCAGCTTGATTTTACAGGCTGCTCTTCGTTAGAAAAAAAAAAATCATTTGGGGACTGCTTTTCATTAAAAGGAAAACCTTACCAAGGACTTCCTTACCCTATCTGCTTAAGTAATTTCTTTTTTAACTTTTATATCATTCTGACCTCTCGTTGCTTTTGAGAAATCTGACAGTTCTTTATGGGGAATTGGTCTTTTCTGTTTTCCATTTAAGATTGTGTGTGTGTGTGCGTGTATTCTGTAATTTTTCTGGATATGTCTAAATGTGGATTTGCTTTTATTCTGCTTTAGACTCCTTGTGCTTATTGAATCTATAGATTCATATCTTTACTCAGTTCTGGAAAATTCTCAGCTATTGTCTTTTAAAAGATGCTTCTTCCAGCTTTCCCAATTCTCTATGGAACTTCTTTAGACATCTTTGGCCTTAATTCTCTTTCCATGGCTCGTAACTCCTCCCTTAGTTAGTCATAGCGCTTGCCCCTCTTTTAAATCACTATTTGTTACTTAGGATTTTCCCTCCTTTCTTAGAATTTCTTTTGTTATATTTTTAAATATACAAATTTATATTTTATATACAAATATATGTTTTAAAATATATACAAAATATATACAAATATATATTTTAAAATATATACAAATATATATTTTAAAATATATACAAATATATATTTTAAAATATATACAAAATATATACAAATATATATTTTAAATATACACAAATATATATATTTTAAATATACAAATATTTTTGTTATATTTTATCTAATACATCTACATATTTATATGAGGACTTTTTTGGCATTATCTTTACCATTTCTGCTGGAACTGATTCTTGACAGCTGTAGTTTTCTTGTAGAGTTGTTTTTTGTTGTTTTTTGATATATTCCAAGAGTTTTTCAAATTTTATCTATAATTAATATTTGGATATTGTAGTTCAAAAGAGCAGGAGAAGACAGTGGTTAATATTTTGGACAGCCCAATACTATAGATAAAATGAAGTCACCTTTTCCTATAAAATGGATTTTTTTTATTGTTTTTCTAGAGCTTTGACTATTGGGAATTTTGCTTAGGATCTCTAGTTTTAAAAGTTTCTGAATCTTCCTGAGCCTCCTTCAGTATTTTTGCATGTAATAAACCCCAAAGTTCTCTCATGGGGATTCTTTTATAAATGGAGAGAATGGCATAAGCAGTTAAGAGCAAGGTCTCTGGAATCAAACAGCTTGTAAGTGTCCAGTCTTGTCTTTACCACTTATTAACATTGTGATTTTGGGTGAGTTCTCTCACCTCTCTTTCTCAGTCTCCCGCTCTTTAAAGTAGCAACCTCGTAGAGTTGGGATGATTAACCAAGATCTAATATGTATTAAGTGGCATCTAAGGAGCTCTCAAAAATATTAGATATTATTAATATAATTATTTTCATTGATTCACAGTAACTTTTTCCTTGCACTTAAAGCTATTTCAGGTGTAGCAGTTTTAATCCGCCACCTGATACAGCATTCTTGAAGGCAAAAGCACTTGGACTAATAATAAGCATAGGCCGTCACTTCCCCCGCATTATACATTTTCTCTACTGAATCTTTTCCATTAGCATACCAACTCATATTTCTTTAAAGCTAAAATAACTCTTTTGTCATCACTTCCCCCTCCAATTTTCCCTTTACAATAAAACTCCTTTAATATGTCCATTTTTACACTCCCTAATCTGTCTTTCATGAAACCAAACTCCAGTCAAGGTTTTACCTCATCCATGCCACCAAATCTGCTTGGATTAATAAACCCAGTGGTAAATCTCAGTTCCTGTGTTACTTGACCCATCAGCAGCATTCAACACAGTTGATCCCTCCCTTTTCCTTGAAACCTGGCTTTCTGCTCCCTCACTGGCTGCTGCTCCTTGGTCTCCTTTGTTGGTTCATCCTCACCTTGAAATGTTGGCATGCTCCTGGAATCCGCCCACGGACCTCTTCTTTATTTGTCTACACTCCTTCTCAGTCTCATGGCTTTGACTACGATCTGTATGTTCACCACTCCCAAATTCATATCTCCAGCCTAGGACTCTTGCTTGACTCCAGATGCATATAGCCATTTGCCTACTCAGTATCTCCATTTGGATGTCTATTAGACATCTCAAATTTCACATGTTCAAAACAACTTCTGAGAATTATGTTAGCTCTGGGTTTGTCATACATGGCCTTTCTACTGAGGTATTATTTCCTCTATACCAATCTGTTGAGAGTTTTTATCATGAAAGGGTATCAAAGGCTTTTTCTGCATTTATTTAAATGATCCTCTGATTTTTAGCCTTCATTTTGTTAATGTGGTGTATCATATTTATTGATTGTGTATGCAATAAATCCCATTTGATTGTGGTGATTCTTTTATCAAAATCTTAATTACCCCCTTTTAAAACATTTGTAACATGATAAACTATCTTACTAAACATTTCATAATAAACAGTTGCTTCTTCCAATTTAAGAAAATGTATATTTTAGTTGTTGCTTGGGCTGAGGATGTCCCTAATCAATCTAACTATATTTGTGTCTTCAGGCATTATTTGAGGTGAATAAAATGCAAATAATTTTTTAGACATAAAGGCCTCTGTTAGCCAAATGAGAAATCTTTCTTTCACCTCTAGTTCTCATTCACTCCCCCTTTTTCTTTCACAAAGGCATTGATTAATTTATATCTTAGGTTAGCTAAGAATTTGAGCTATTTTTAGCAGTACATTCTCAAAGTATAGATGCTGTGGAATGAGTAGCAAGTAATTAGAGCACCGTTTGTATAACCTTTATGTACAGTAGCCTTGAAAAATCTCCGTACTGTAGGTAGAATTGCACCTAAGATGACAAATCTTGTAGTAGGACCACCTAGCTAATGATTGAGTCTCCTTTATAGCCTGGGAACTGATGAGAAGGACACTCAGCACTCCTGAGTTCCTCCTGTCTCTACGTCTTTCTCATATTAAGCCTAAGTTGGCTTTCACACAATTAATTACAATTTGCTTCTATAGCTTTGGGGCTTTACGTAACAACTTTAATTTTTTCTACACCCTTTTGAAATAGCTATCATGTTTGACCCGGGTCTTTTCTCCTTATAAAACATTAGTATTGTCTCCTGTTGTCTGATATAGTTTTGAGTCCCTTCAACATCAACTGTAGCTCTGTCAGTGTCCCTCTTTAAATACAGTGCTCAGAGTTAGCCTTATAATGTTTCATGGTTTTGTTTTTGTTTTTAGATACAAGGTCTTGTTCTGTCACCCAGGCTGGAGTGCAGTGGCATGATGGTAGCTCACTGCAGCCTTGGACTCCTGGGCTCAAGGGATCCTCCTGCCTCGGTCTTCCAAATCACTAGGATTATAGGGATGCATCACCGCGACCAGCCTCGTGTATTATTTAATTAGTAGTCATTCAACCTATTTTTTTCTAGATAGTACATTTTGATTAGTACAGTTTATGATTTTTTTTTTTTTTACCTGTATCTTACTTATAGTACTTACAGGAGCTGAGGAAAACCCCAAGACCAACCTGTGCTCAACCTGTGCTTTTGTAGTGATTTTTTTGGACCCTATTACAGAACCTTATTTATTCCTTCATCTTGTTAAATGTGGCCCATTGTTCCAGGCTACCACCACATTTCTTGATTCTTTCTTTGAACATAATTATTATATTTCCCAGCTTCGTGTCATTTCTAAAATTGATGTACCTGCCTTCATCTAGGTCATTGCTAAAAACAATGAAGTAGAGCTAAAGACAGCCTTTCAGCATTCACAAAACTTATCTATTCTATATGCCATATATTTATATAGCACCTTCTATTTGCACTGTCCTAGTCATGAGGGAATCATAAATAAAAAGAAATGGCCCTTTAAGAGCTCACAGCCTAGTGGGAAAACAGATACACAGTATATAATTTAGCATGATAATTGCTGTAATGGAATATGAAGAAGATAGAATGGGAGCACAGAGGAAAAACTCCGAGTCTTTGAGGCAGAAACTTCCTTGCAGGTTGAGAGCTGTCCTTTAATCAAAAAAGTGTAAATTAGCCAGGCACGGTGGCTCACCTAATCCTAGCACTTTGGGAGGCCGGGGCAGATGGATCACCTGAGGTCAGGAGTTTGAGACCAGCCTGGCCAACATGGTGAAACCCCATCTCTCTGAAAATACAAAAATTAGCCGGGCATAGTGGCGCATGCCTGTAATCCCAGCTACTCGGGAGACTGAGGCAGGAGGATCACTTGAACCCGGGAGGTGGAGGTTGCAGTGAGCCGAGATGGCGCCACTGCACTCCAGCCTGGGTGACAGAGCGAGACTCCCGTCTCAAACAAGTGTAAATCCACCTAACCTGTCTCTTTTTTCCCAACTTGTTCACAAGGATACAGTGAAAGACCCTGGCATCTGCTTGCTGAAATTGGGCCACATTCTATTTACTACATTTCTCTATTCCACCAGCCTAATACCCTGTCCCCTCCTCTGCTCCCAAAAAAACAAAACAAAACAAAACAAAACAAAAAACTGTTAGTCTTTTAGATGTAATTAATTAAGCTGTTTATATTTTAAATGTATTAATTTTTACATATACATACATATATCTTCAGGAGCAGGATCGGGGAAATGAGTCATAGAATATAATAGACTCTGAAGATCAGAAGGAAAAAAAATAGAGGCAAATTTGCAGGCTATGAGATCCCCCAGTTGTTATAGTTGAAGATCCAAATTGGGCTCTGAGCTTTCTAAAGGCCAAAGAGAAATGGAAAAGTTCATGTTCTTACTGACTGACAGAGAAATAGTTCCTGTTCTTAAGTGGGAACAAAGCTTTTGCCAACACCGAAATCTTACATTTCAGTGGGGGCTGCATTGAGGTTACACCGAGAGAGAGAATAAAAAATGTCATGAGCATATATTTCTACTGCAAGTGCACAAAAGATTTCCTTAGTGTTTCCTAAATTTTGTTATTGTGGAGGAAATTGCTCAAAGAACAGATTTAATATTACAGGCCAAGAGAGGCTGGGCAGCATGGCTCATGCCTGTAATCCCAGCACTTTGGAAGGCCAAGGCAGGCTAGATAGCTTGAGCCCGGGAGGTGGAGGTTACAGTGAGCTGAGATCACGCCACTGCACTCCAGCCTGGGTGACAGAGTGAGACTCTGTCTCAAAAAAAAAGTGGGGCGGGGTGGGGGACTAAAAGAACTGTTATAGAAAGGAAAGCTGCGCAACTTTGAATAACAAATTCATTGGGGATGTGAAATTTTTTGTGCACACACCTAGACCATCCCTCAGCACAAGTTCAAGAAGTCTTATATTAGGCTGTTTATTATAGCATCCTTCATAAGAGAGGAACACCAATTCCATGGTAGGTTGACAAGTATTTCGTAGATACCCATGATTCTGTGGTCAAATAAATATAGGAAACGCTGAGTTAAACAAGTTTCTTTACTGCTGGTGGGCTTCTGTCTCATAGGCTTAAGTACAATGTAACTCTCCCAGAGGGAGTTCTCCAAAACTTTATTTGGCCATAGAACCCTTTTTGGTGGTACATTAGGGAACAGAGGAATAAATGTTTACATATCCTCCAAACAACCTTTCTAAAAGAGAAAGAATGTGGCAGCCACCCTGTTACATATTTGTCAACTCCAGCCAAGGATATATGTTGGACAATACTGTTCTTCCCCTCTCTTACTCTCTCAAGCTTCAGGATGATGTAGTCTCTTTTTTGCTAGGGTCTTCTCGTTTCCTTATAAGGCAAACCAGAGTACAAGTTCCCCTCTTCACTTTCTCTTTTTTTTTTTTTTTGAGATGGAGTCACTCTGTCGCCTAGGCTGGAGTGCAGTGGCGCAATCTTGGCTCACTGCAACCTCCGCCTCCCGGGTTCAAGCGATTCTCCTGCCTCAGCCTCTCTGAGTAGCTGCGATTACAGGTGCACACAACCACACCCGGCTAATTTTTGTATTTTTAGTAGAAACGGGGTTTCACCACATTGGTCAGGGTGGTCTTGAACTCCTGACCTCGTGATCCACCCGTCTCCGCCTCCCAAAGTGCTAGGATTACAGGCTTGAGCCACCGCGCTCGGCCCACTTTCTCTGCTTTTTTAAAGATGAAACTATGATCAGGGCAAGCCAGGAGCATGTTGTCACTTGGCTTTTTTGCCAAGTAAGACTTCTAGTTAGAGACCCCTGTCACTGCAATTAACTTTGCTGAATTTGTAACTTGAGATCAGGGTTCACTTCCTTATCTGGCCTCCATCCCATTCACAATGTCATTTTGTATAGCAAATCTAATGTGTTTTTTTTGTTTGTTTGTTTGTTTTTTAGACGGAGTTTCACTCTGTCGCCCAGGCTACAGTGAAGTGGCGCAATCTTGGCCCACTGCAACCTCCGCCCCCTGGGTTCAAGCGATTCTCCTGCCTCAGCCTCCCTAGTAGCTGGGACTACAGGCGCGTGCCACCATGCCTGGCAGATTTTGTATTTTTAGTAGAGATGGGGTTTCGCCATGTTGGCCAGGCTGATCTCCAACTTCTGACCTCAGGTGATCCACCCGCCTAGGCCTCCTGAAGTGCTGGGATTACAGGCGTGAGCCACCGCACCCAGCCAATCTAATGCGTTTCTTTTAGAATAAAGTAGGTTTTTTTCCCATGTATTATAAACAAGCAATAAATCTGGGGAGGTTATTTCTACTTGATTAGTTTTTATTTACCTAATAAGAATCCTAAAATGCAAAAATAGTATTGTATTAATATAACAGTGACTTTAAAATAAACTTCAGCAATATTACATTTGAAAAATACCCTAAGGTTAAATTTATTTATAACAAGTGTTATTCATTTATTATAATACCCCAGAGTCACATTCAGTATTGTCTTTCGATATATGCATTTTATAATTGTAAATCCATATACAGATGTCCCCTACTTAATGAGGGTTCGACTTACTGTTTTTCAACATCACGATGGTATGAAAGCAACACACATTCAGTAGAAACAGTACTTCGATTTTTGAATTTTGATCTTTTCCTGGGCCAGTGATAAGTGGTATGGGGGCAGCAGCAGCGAGCCACAGCTTGCAGTCAGCCACACCATCATGAGGGTAAACACCGATACTCTCCTGTGGACTGTGATGCAGATGATTCTGCCCGACTGTAGACCAATATACATGCTCTGAGCATGTTTCAGATGGACTAGGCTAAGCTATATGTTCAGTAAGTTAGGTGTATTCACTGCATTTTTGACTTAATGATATTTTCAATTTACTGCGAATTTATTGGACATAACTCCATTGTAAGTAGAGGAACATTTCTGTCTAATATTTTGTGTTTTCCTCCCATTAAATTCCTTATTCTCTTATGTTATGGAATTCCTTATGGAATACATATTCCCTTACGTTGACATTGACAACATACTGGTAATCTTTAATATATATGAAATACTTGAATATTGTCAAATGTTATTTTTGAAGTCATTTTGATATTTGGCAATTGAGGTTTTTTCAACTTTTTGCAATTTTAAAAGGGGTATTGATATCATATGAATTATCTTTTTTCTCTGGTTTATTTCCTTGGGCTGTATTCCTATGGATAGGACTACAAAGTCAGCAATTATATTTTCTAGCTTTTGTATGACAGGAATTTTCTAACCTTATTGTCTTCTAAATCCTAGGTTCAGGAGCTCCATGAGGAGGAGGCACAGTGGGTGAACTATGATGAGGATGAGTTGTGTGTGAAAATGCAGCTAGCCGACGGGATCTTTGAGACCCTGATCAAAGATACTATTGATGTTCTGAATCAGATCAGTGAAAAGCAGGGGAGAATGCTACTTGTGTGACATCTTGCAAATAAATCGAACGCTGAGTGCTAATGTGAGTCCTGGGCCTTTCTGCCTCCTGATGTACACCCATCGCCATCATAGCAAGAGTGCTTCTGGACCTTGTACTTATTCTTAAAGACTACCAGTATGGAGTTCATAGGACAATGTGGTACACCTGGTATTACAGCCTTTGCCTTTCGAGACTATCCACTGGATTAATGGGTTATTTTCAGTGGGCAGGGTTGCACAGTGTAATCCTACACCTTTTGCTAACACCCCTACTAGGTCCCAGAGGGCCAGAAACACCTGACTTACCTCTGAGTTTAGACTAGGGTATCACTTCTTTGAGTCTGAAGTCAAGTGAGAGAGGTATAGATAAATGCATCACATCACTTTTGAAATGTAATTCTGGTCTATACCATGGAAGTCATAAATGGACATTATAGTCTCTAAACAGTATTAAACCCTTAACCACTTCTAAAATAGGCAAGCTCAATAATGTCTGCCAACTTCACATTCTGGAGTTTATTTCATTTCTTTTTGAAGACCATTTTCTTCCATTATTGTAGTTGAGCAGCACCAAGTGGACTGTCAGGCTAACAGGAATAAGTGGTAGCCTTGCTTTCTGAGCACCATCTAAAGAATTTTAAACCTCTGCATTATGTTTAGTGTTCTCTGTGTGGGCATGAAAACAAAGAAATGCCCCTACTGAAGACTGGGCTCAAAGGACCAATGCAGGGCTGGTTCTTTTACCCCTTGGTTTACTCCTGCACTTGTCTTACTCATTCTGAATCTCTACCAGCTGCTCCCAGAATCACAGATACTCAGGACCATCTCAGGCATCCAGGCATGGCAAAGTGGAAATAATTCATTTTGGCCTGCAAACCTTCAACTCCCATCTTTCCCCAAGAAGTCAGAGATGCTGTTACTTGAATGATTTAGGAAATGAGTGTGTGCACCAAAGACAAAAAGATATTGTCTATTGTTTGTGTGCTTGTTTTGCGCTATGGAACATTTTTTAATTTATTTTAAGATAAATTATTAAGTTGAAAATGTGTGTCCCTATTCAGAAGTGAAAGATTCATCTTGTAATAGTTAAACCTCCATCTTGAAGCTTCTATGGTTCATAGTCTTTGCACAGGAACCTGTGGTTTTAACAAACCAATACACATATTGAAGAAGTCATTTTAATTCAGTGAAACGAAGATGGGCTTTTCCAGATCACCTGCAATAGCAGCAGTGGGATAAAATGATTTAAAAACACTGTACAATTTAACTCTGCCTCTCTTGCAGCATTGCTTCTCACAACTATTACCTGCATCTGAAAAAAAATCTATAGACTCCAGCTGCTACATTAGAGCATAAGAGATGCTCTCCTGGGACCTCAGTACCCTGCCTTCTTGACTGGTTTCCGTTCATCAGTCCTGTCCCTCTTCAAGTAATCTAGAAGAATGTGGATACTCTTAGGCGTGAATGTAAATGCCTTAATATTGAAGGTCCTGGTTAGAAGCATGATACAAGACATCTACTGGATTCATATTTACAAATATCCTGGAATGTTATAGCTTCAAAGTATATTAGAAAAACCCCAAAGATGGTATAATCTTTAAGTGTGCACGTTCGTTTATTTCTGCATCTTCCCTCCAAACTTGCCTTTGCATCTTAAATATTTCACTATGCACACTCCCATTCCTCTTGGGTTTCATCTTGTCGTTTAAGAAATGTACTGAAATAATCATTGGAATATTTGCATTTTGCACAATGACTGGTATGATAGCTCTTGACAAATAAGGAAAGCACTGAAATGTTGTGATTGGGTCTCGGGAAATGCTCAGATTGATGTCTTACCAGCATTTCTTCTGGGCTTGTGATGTTGAGCTGTAGTCTTGTAGCCATAATGAGCAAATTGACTAAGAGAAGCAAAGGTTTCTTGGGGTTATTAACCAGTAGTGTGGAAATACTAGTTTTATGTGGCCAAGGAAAAGCAAAGGCTTTTCTTTTCAGTTTGTGTTATTTGGAAGACAGAAAAACATCTTGTCTACATCCTTTGGCTGTTTGTAGGATCACGTTGTCCTTACGATACTGAAACTTTACAGCTGCTGTAAATTTTTTATAAATGAATTTCAAAATGTTATAATGGGACTGTAGGTTGTTTTTCTACATCTTCATTATTTGGACCTAAAACCAGTTTTTAATAAGAAAGTTTATCTTTACTCTTTCTGAAATTATGACTCCAGAAAAAGAAAAAAAAAATACAAGTCATGGAATCAGCAATCTGGTAAGAAATGCTGCCAAGAATGTGGCAGTAGCTGTCCTGACAGACTCCAACTGTCTTTACTATCTGAAGAATCCTAGGCTCCACATGAGAGGCAGAAATGGATCAGTCTTATTCTTTTCTAGAAATGGTTATCTGTAGTTTGGTAGCAAAAAAAAAGAAAAAAGAATCCATAATTAGCAGATTTCTTATTAACTATTTGGATCTAATTGAAATGGCTTTATTCTTAGGATTAAGAAAGATAGATGTGGATACCCAGCCACTCGTTCCATATTGGTATCTTTTTAAATCAGCTCTGCCTCTTAATCAAGAACCTAAATATTCCCTCTTTCTAATCTTTGTTCCTTCTCCCTACACCCTCATCCTCTTTCACTCTTCCTTCATAATTCCTCTAAGAAAAATATCTTTGCATCAGCAGTAATATCTTTTAGAATAGCACTATCAGAATTTAGCAGTAAACCAACATACAGGCTTCAGATTTACTTCTGAGTCCAAAACAATTTGTGCTATCCAGGGTAGTTAACTCTGGGTTAAACAAGTACAGGGTATAGATTCCCTCTTCAGGTCTACACAGGAATTTTTACCATAGGGAAAAGTGGGGAGAGCTCAAACGTAGTTAATAAGGAAGGTAATTTGTTTTTCTTTTACCTAAAAGAAAAGAAAATTCCTTCTGTGACTACAGGTCTCTGAGAAATTATCTTTCAAAAGAGATTTCATTGCTCATAAGAGTGTTGTGGCCTATTGATAAAAACAATTTTGTTCAGTTTCTTGTCTTGAAAAAAAAGTGGCCTTAGCTTTTTGCAATACTTGAATAAAGTGTGTACTCGCAAAAGAATTTCTGTAGCACAGCATTAGAGACTCATAACTTTTCTGCAAGAAATACAAACTTACATCTTCCTTTTACTACCTTAAGAATACTAGTGAATAAAACATTAATTCAAAGAGCAAATTATAGAAACTACAATGACATTTAATGCAAATTGTAGGAATTTACATGTTTACAAATCATCTTCAACTGGTTGTGCAGCAATTCAATAAAATATCTTTGTATTATAAAAATGTGAAGAAAAAATGTAAACTGATGTAAAGGAGGTACTGTCATTTTAATTAACCTATGTTTAATAGCTTTTCCTTCTGGACTTTGCAAAGCCTTCTTGGCAAACACATTGCAAAGCATTCTCTGGGAGGTTCAGCCTCCTTGTGTGTACTGTACTGTGCAGACATGAAAAAATAAACCCGTTTACTGTGTGCGTGTAAATAGCCTGGTCATCAGGCCATTTTCAGCCAATAGTCACATCCAGTGCAATTTTGCACCGAACACTTAAGGGTGTGGTTTGTAAGTACGATCTGTAAAATAACTGGGATGAATTCCCATGTATACCTGTGTAAATAGATTTGTTAACTGAAATATACTTTAAGAAAGATAAAATCTGTAAATAAACTGATTTATAAATTAATTTCATGTCAGGTCTCCATTTTTTCCCTCCTGAGCCAACTTCAGAAAATTAGTCCCCATAAAATACATCTGATTGCTTATTTAGTCATAGGGTGCTGATGCTTGGACAGTACAATCAGGTAAGCTGTGGATAAGAGCTATTTACCTGCACATAAATGCCATCAGCACATACATTTGAGCAGCCTGGGCTTTGCTCCCACACTCATTGATTCAGCAGATAGCCACTGAGCAACTTTGTGTTCTGAGTACTGTGCTGAGTGCTGGAAACATAGCAAAAAAGACAGACATGGTCCCTGATCTCAAGAAACTTAATGTCTAATGAAAGCCACCAAAAGATATTTTGAAAACCTGTGTGACAAGTTCTAACAAGTGCTGCACAGGAAGCATACATAAGATTTGAGAATATCTTAATAAAAGGACCCAACAAGGGTAGAGTGAGTCAGGGAAGAGTTCTCTGATAAATTATAACTGTTTTATGAGCTTCTTCCTTCTAATTCATGGCAAATGCCATTTCTCCCAAATGACCTAAATTTTGATAACCCTTGCAAAATAACACATGGAAGAGTCACACTCAAGTGGACCACAAAGAAGAAGTGCAACCAGCTGATGCCCAGGTGTGGAGGCTGTAACTTATTTCTGGAAATATGAGCTTGGAGGAGAATGAAGCTGAAGTGATAGAAGAATAAAAACAATCCAAGGATGGAAGTGTACAGCACTCAAGATAATGAATTGTCTCTGGTCTGACAAATCGGAAGTGGTAACCACCATCATCTTACCTGAATTCCTCTAGTTTCCCCCAGATGGCCATCCCTGCTATATCTTTAATCCCTACCACAACCATGAACATTTCGGAAATTCTCAGAAAAAGATTGGCCCCAAATTACAGTAGGTCCTAAGAAATTCTTCCTCTCATAAAACAGATAGTTGAGACTGGGCATGGTGGCTCACACCTGTAATCCCAGCACTTTGGGATGCCAAGGCAGGAAGATTCCTTGAGCCCAGAAGT',
			'insseq'				=> 'G',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		my $expectedRnaWtSeq = 'UUUUUUUCACUAAAAUUUGGCUUUUCUCUAAACUUUGGUUCACAAAAUUGAAAGAUAGUGAUUAAAUAAUUUUGAUGUGGUGAUUUGUUUUUUCUUUAAGAUGCUUCACUGUCUGGUUCUGAGAGAUCAGUAUCAGAAAGGUCUUUAUCUGCAUAUGCAAAGAGAGUAAAUGAAUGGGACAGUCGAACAGAAGAUUUUCAGACCCCAUCUCCAGUUCUCAGAUCAUCAAGGAAAAUCAGAGAAGAAUCUGGAGAUUCUCUAGAAAAUGUACCUGCAUUACAUCUUCUCAAAGAAUUAAAUGCCACUAGUAGAAUUCUUGAUAUGUCAGAUGGCAAGGUUGGAGAAUCUAGUAAAAAAUCAGAAAUAAAAGAAAUAGAGUAUACAAAAUUGAAGAAGAGUAAGAUUGAAGAUGCCUUUUCUAAAGAAGGUAAAUCUGAUGUCUUACUGAAAUUAGUCCUAGAACAGGGAGAUUCAUCUGAAAUUCUUUCAAAGAAAGAUCUUCCUUUAGAUUCUGAAAAUGUUCAGAAAGACCUAGUUGGAUUAGCUAUUGAAAAUCUCCAUAAAAGUGAGGAAAUGUUGAAAGAGAGACAGUCAGAUCAAGAUAUGAAUCAUAGUCCAAACAUCCAAUCAGGAAAAGACAUUCACGAACAAAAGAACACAAAGGAAAAAGAUUUGUCUUGGUCAGAACAUCUUUUUGCUCCUAAAGAGAUACCAUACUCUGAAGAUUUUGAAGUGUCUUCUUUCAAGAAAGAAAUUUCAGCUGAAUUGUACAAAGAUGAUUUUGAGGUGUCAUCUUUGCUGUCACUCAGGAAAGACUCUCAGUCUUGCAGAGAUAAGCCACAGCCAAUGAGGAGCUCUACAAGUGGAGCCACUAGCUUUGGUAGUAAUGAGGAAAUCAGUGAGUGCCUAAGUGAGAAAAGCCUUUCUAUCCAUAGCAAUGUUCAUUCUGACAGGCUGUUGGAACUCAAGUCCCCUACUGAGCUGAUGAAAAGUAAGGAGCGCAGUGAUGUGGAGCAUGAACAGCAAGUUACUGAAUCCCCUUCCUUGGCUUCAGUUCCUACUGCAGACGAGUUAUUUGAUUUCCACAUUGGUGAUAGGGUGUUGAUUGGAAAUGUUCAGCCAGGAAUUCUUCGAUUCAAAGGUGAGACUAGUUUUGCUAAAGGAUUUUGGGCCGGAGUGGAGUUAGAUAAACCUGAAGGAAAUAACAAUGGAACAUAUGAUGGUAUUGCAUAUUUUGAGUGCAAAGAAAAGCAUGGUAUUUUUGCUCCUCCUCAAAAAAUAUCUCACAUUCCAGAAAACUUUGAUGACUAUGUAGACAUUAAUGAAGAUGAAGACUGUUACUCAGAUGAACGAUAUCAGUGCUAUAAUCAAGAGCAAAAUGAUACAGAGGGUCCAAAAGACAGAGAAAAGGAUGUCAGUGAAUAUUUUUAUGAGAAAUCCCUACCUAGUGUGAAUGAUAUAGAAGCCUCAGUUAAUAGAAGUAGAAGCCUUAAAAUAGAAACAGACAAUGUACAGGACAUUUCUGGGGUACUUGAAGCCCAUGUUCACCAGCAGUCUUCAGUGGAUUCACAGAUUUCUUCAAAGGAAAACAAAGACCUCAUUUCUGAUGCCACAGAAAAGGUUUCCAUCGCUGCAGAAGAUGACACUUUAGACAAUACCUUUUCCGAAGAAUUGGAGAAGCAACAGCAGUUUACAGAAGAGGAAGACAACCUAUAUGCUGAAGCUUCAGAAAAGCUUUGUACACCACUUCUGGAUCUUUUAACAAGAGAAAAAAACCAACUGGAAGCCCAGCUGAAGUCAUCACUAAAUGAGGAAAAAAAGUCAAAACAACAACUGGAAAAAAUCAGCUUACUGACAGACAGUUUACUAAAAGUCUUUGUAAAGGACACAGUCAAUCAACUACAACAAAUCAAAAAAACCAGGGAUGAGAAAAUCCAGCUUAGCAAUCAGGAGCUUCUUGGUGAUGACCAAAAGAAAGUAACACCCCAAGACCUAUCCCAAAAUGUUGAGGAACAGUCGCCAAGUAUUUCAGGUUGCUUCUUAAGUUCUGAAUUGGAAGAUGAAAAAGAAGAGAUUUCCUCUCCAGAUAUGUGUCCCAGACCGGUGAGUAUUUCGUCUAGAAAAAGUAUCAGUAUACCUUUUGACACUUUUAACAAGGCAACUUACCUUGCCUAUAUAGUACAUUUUAAUUCUGUGACUCCUGAUAGGUUUUUUUUCCCUUGACGUGUUGUUUAUUAGAGAAGGAUCACUGUACUGGUUUUUAAAGGCCUUGUAGCUGUUUUUCCUGGCUGAAGUGUAAUUGUCCUCUAUUAUGUCAUGAGUAGUCAAUAUAGGGUUGUCACAAGCACUGGCAGGUUUUGGGGAACUAAAACAUCCAUUUAAACUGCUGAGUUGACAACCCCUGUACUGUCCAUAGAUUACAUCUCCUUCCAGAUUCACUGAACAAAUUCGUUCCAUGGAGCUUUCCAGUGGCUCACCACCCAAUUCUUUUUUGCCUGUGAAAUUGUUUACAUAUCUACUCAGUCAUGGUUGGCCUGCAGAGGUUAGCAGUUGUGAUGUUGAUUCUGGUAAAAUGAAAUAGAUCCACUAUGACAUGCACUCAGACCAUCAGCUUGUUUAUGUUGCCCAUGGACAGUUGAGCAAAAGUAAUGACUUUGAGUAACUUGAGCCAGAUGGUCACCAGAGGUGAGAAAUUUUAUAUCCCAUUUCCAUAUCCCAUAGACAGUUCUUUAGAAACAGAAUAAUGAUAAUACAAAAAGCGGGGAGAGCCCUAAUCAAUAAACAUUUUAUAAAAAUUCAGUUAUAUGCAAGGUUUAUAUAAGUAGAAUUGAGGGAACAGGAUAACUUUUUCUUAGAUUGUUUAUUAGUAGCCACCAUAUUGUUUUUUAUUCUGCUUACUUAUAGAGAAAUAUUAGAUAUAUAAAUAAUGCAUUUACUCUGUCUUUUUUUAAAAAAAAAUGAUAUGAGGUCCCUAAAUGGUGCUCUGUUUGAAUGGUUCCAUUUCUAUAGGAGAGCCCAGUAUUUGGUGCCAGUGGGCAGGAAGAACUUGCUAAGAGACUUGCUGAACUUGAACUCAGCCGGGAGUUCCUGAGCGCGUUAGGAGAUGAUCAAGACUGGUUUGAUGAAGACUUUGGUUUGAGCUCUUCUCACAAGAUCCAAAAAAAUAAGGCAGAAGAAACCAUUGUACCUCUAAUGGCAGAACCUAAAAGAGUAACCCAACAACCAUGUGAAACAUUAUUGGCAGUCCCCCAUACUGCAGAAGAAGUAGAGAUUCUUGUACAUAAUGCAGCAGAAGAACUUUGGAAAUGGAAAGAAUUAGGCCACGAUCUUCAUAGCAUCAGUAUUCCUACAAAACUGCUUGGCUGUGCCAGUAAAGGUCUAGAUAUAGAAAGCACUAGUAAAAGGGUCUACAAACAGGUAGGUGAAAUAAAAGGAUAAUUUUAGUUUUUAAACUUUUCUCCCAUUGUCUGAAAUUUGGAUUCUUAGCCAUUUUGUUUUGUUUUGUUUUGUUUUUCUCUAUCAAGGCGGUUUUUGAUUUAACAAAAGAGAUUUUUGAGGAAAUAUUUGCUGAGGAUCCCAACUUAAAUCAACCUGUCUGGAUGAAGCCAUGUAGAAUCAACUCUAGUUAUUUCCGACGAGUGAAAAAUCCAAAUAACCUUGAUGAAAUCAAGGUAAACUGCAAACUAUAAAGUGUCUUCUUUUUUGACUUGCUGUUCAUUUACACAUAUCUGUAAAUGUUGUAAAAGCAGGUGCUCUUCUGUAUGAUUAACCUUUGUGUCUUCAGGACCUGUUCCAUAGUGGGUAUUCAAUAAAGAUUUGUUGACUUACUGAAUCCAUGAACAGUGGAGGGAUUGGAAUGAGUGAAAGGGCUGAAUCCUCAAGUUUGUGCAAAAAAAAAAAAAAAUUCUUUUACCAGUCAUUUUCUCUAUAUACAUAAUUAUAAUAUUUAAGAAAAAUUAGGAAAUACCUAGGGAGAUAAAUUAUUGCAACUAAUAGGAGUUUUAGGUAAGUUUUGGGGGCUUUAGAGAUAGGUAUAAUGUUCUAAAAGUACAUAUGAAUACGUGUCAUAUAUAUGUCUUUUCAAUAUGCCAGCCCUGGACACAUAGGUGUUAACUUCUUUAUCUUUGACUUUGGAGGUUUCAAAUUGGUAAUGUGGUUUUGAACUAUGUAGGAUAGUAAUGUGAUACGUAAAAAUAUUUAUUUUGUUAUCUAAGAAACCUUUGUAAACCUUUGUCCAGAAGCUACAAGGACUAUGAACUUUAAAGUCACAUUGCUGGCAAUCUGCUCUGCCUAGUCACAUGUUUCUUUACGGCAAAGGUUAUAUUUCCAGAGCAGUGAUUCUCAACUUUUAAUGUGCAUGUGAAUGAUCUGAGGAUCUUGUUGAAAUGGAGAUUCUGAUUUAGUAGGUCUGGAGUAGGACCUGAAAUUCUGCAUUUCCAACAAGCUCUCAAGUCAUGCUGGUCCUGCUGGUCCUGCUGGUCCUGCUGGUCCCAGGACCACGCUGAGUAGCAAAGCUUUAAACUUGUAUUUCCCAAACUUAUUUGAUCACAGAAAGCCUUUAAAAAUAUUUUAUACUUAUAGAACACCAACAGAUAUUUAUGAGAUGGGAUGUAUGCUUACCUCAGGGGUUGGCAAACUAUGGCCCAUAAGCCAAAUCCAGGCCCCUGUAUGUUUUUAUAAAUAAAGUUUUGUUGGACCACAGCCAUACCUAGUUGUUUACAUGUUACAACUACAUGUAACUACAACAGUGGAAUUGAAUAGUUGUGACAAAGACCACAUGGCCCGCAAAGCCUAAAAUAUUUACAGUCUGAUUUUUUACAGAAAAAAUUGGCCAUCCCCUGGCUUAAAGUAACUAGCUUUUUAAAAACACAUUGUAUACUCCAUUCCUUCUACCCUUUCACAAUUGAAAAUGCUUGAAUGCUAUGUGUGUCUGUGAACAGGUAGGGAGUGAUGGAAUUUCCGGGCACGUUGGUACUUUGAGAAUCUUUCCAGGUGAUUCUGAUUCAUAAAUGUGUCCCUUCCCCUGUUAAGAACAUUGUUUUUUUAACAGAACAACCAUAAAAGUAACAAAAUAAAUUAUUAAUAACAUUAAGUUGAGGUAAGACAAGAUAUAUAAUAUAUACCAACAAUCAAACUACAUGCAUCAAUAUAGCAUACAAUAAAUGCCCACUAAUUUUAAUAUUGACAUAACAUUUUGUAAAGUUUUUAAGGGUGAUAGUUCAAAUAUUGUUAGUAGUAAAAUUAUUAUGCUUCUCAUACUGAAGUGUUGAUAUUUUUUAAUUUAACAUAUUUGAGUUUUGGUGCUUCAAAAUAAACUAUCCCAAAGUCUUUUGCUUUUCUCACAUUACAUUGUCCCCAAAAAUUUCUUUGAGAUUUUUCCAAUUCAGAGAUUUAUUUUAUGCUUUCUUUAGUAGUUUUGCCAGCUUGAAGUUAGCUAAUGGUAGUUUUGUUGUUGUUUUUUCUUUGUUUUUCUUGAGGUGGAGUUUCACUCUUGUUGCCCAGGCUGGAGUACAAUGGCGUGAUCUCAGCCCACUGCAACCUCCGCCUCCUGGAUUCAAGUGAUUCUCCUGCCUCAGCCUCCUGAGUAGCUGGGAUUACAGGCAUGCGCUACCUCACCCAGCUAAUUUUGUAUUUUUAGUAGAGACAGGGUUUCUCCAUGUUGGUCAGGCUGGUCUCGAACUCCCGACCUCAGGUGAUCCGCCGCCUUGGCCUUUCAAAGUGCUGGGAUUACAGGCGUGAGCCACCGCACCCAGCAUGGUAGUUCUUAAGUAAACUCUAAAGUGAUUUUUCUUUGCAGAUUAAUCAGUGGGUUUUACUUUGUUGUUGUUGUUGUUGUUGUUAUUUUUAGAUAAGAUCUCGCUUUGUCACCUAGGCUGGAGUGCAGUGGUGCAAUCUCAGCUCACUGCAGUCUGGACCUCCUGGGCUCAAUCAGUUCUCCCACCUCACCCUCCUGAGUAGCCAGAACUACAGUCAUGUACCACCAUGCCCAGCUAAUUUUAUAUAUUUUUUGUAGAGAUGAGGUCUCUCACUAUGUUGCUCAGGCUGGUCUCAAACUCCUGGGCUCAAGUGAUUCUCCUGCCUUGGCUUCAGAGAGUGCUGGGAUUACAGACGUGAGUCACCGUGCCCAGCUAGGUUUUACUUUUAUAUGUUGCGUAAGUGAAUAAAGUCUCAGCAUUUGGAAUGCUUUUCCCUUCAAAUUGCUAUGUGGAACUGUAAGUGGAAGAAAUUUGUAUGAUAAUUUUAUUGUAAAUAUUAUAUUCUUAUCGUGAAUCUGUUUACUGCAUUGUAAACACAUGAAACUCAGGAAUUUCUGUUGUUUGUGUUUCGUGUAUUUGUCUUGUUUCCACAGAGCUUCAUAGCAAGUGAAGUACUCAAGUUGUUCAGUCUUAAAAAGGAGCCAAACCACAAAACAGAUUGGCAGAAAAUGAUGAAAUUUGGAAGAAAGAAAAGAGACCGAGUGGAUCAUAUCCUGGUCAGUGUAUACAACCAAACUGUUUUUAUUUUGACCAUAUCUUUUAAACUCAGAAUGCAGUUAGAUUCUAAAAGGAAAAGGGAUAGACUUAUUCUUAAAUCUAUGGUAUUAAAAUCUUAAAGCCUUUACUUGCUGCAUUUCUACCCCUCAGUUUUGUUAACUAUUUUCUUUUGUACUCAGUAGAUGAACAUUUAACUUGAUCAAUAGUAUUUUAUGUAGAACUUGAGAAAUUUACUCAAAUGGAAGCAAAGUUUAAGGUAGUGAAUUGAAAAGAGACCAGUGUGGACUCUGGAGCCUGACUGCCCAUGUUUAAAUCUCUUGAUGUGUAAUCUUGGACAAAUUUGUCUAACAUUUUUGUGCAAAGUGAAGGUAACAAUAAUAGUACUUACCUGAUAAAAUUAUUAAGAAAAUUAAAUAAGUAUAUGUGUAAAGCACUUAGAAUAGUUCCUCAUAUACAACUCAGUAAAUAAUAGCUGUUACCUCAUCAUCAUCACUAUUAUCAUUAUCCAGAACAUUCCCUGGAUAGGAUUUAUUUGGUUUUGCAAUAUAUUUCUGUAACAUUAGGAAAGUGUUUAUCAAAAAAGAAAAAAAGAAAAAUUCAGGUUUCACUCUAAACAUUUUAAUUCUGGACCUUUUAUCUAUAUUGACAUUUAUAUAAUGAAAAAAAGUGUAAUACUUAAAUUUAUCACCGUUCUCAUGUAUGUACCUCUCUAGUGUUUCAUUUACUCAGUAAUAACUGUACUAGAUAACAUUUGCUUAGUAUUUUAUAGUUUAUCAUUCUUAUGCCUUAUUUGAUUUUUUUUUUUUUUUUUUUUUUUGAGACACAGUCUCGCUCUGUCACCCUGGCUGGAGUGUGGUGGCACAAUCUCUGCUCACUGCAGCCUCUGCCUCCCAGGUUCAAGUGAUUCUCAUGCCUCAACCUCCCGAGUAGCUGGGAUUACAGGCACCUGCAACCAUGUCCAGCUAAUUUUUGUAUUUUUAGUAGAGAUGGGAUUUCGCCAUGUUGGCCAGGCUGGUCUCAAACUCCUGACCACAAGUGAUCCGCCUGCCUUGGCCUCCCAAAGUGCUGGGAUUACAGGCAUGAGCCACCGCGCCCAGUGUCUUAUGCCUUAUUUGAUCCUCACGAUAAUUACAUGAGGUAGGAAGAACUAAUGUUAUCCUUUGUUUUACAAAUGAGGUGUCUUACGCCUUAUUUGAUCCUCACAAUAAUUACAUGAGGUAGGAAGGACUAAUGUUAUUCUUUGUUUUACAAAUGAGGAAAUUUAGAUUCAAAGAUGUUAAGUGACUUGCCCAGUAAAAUAAGAGACGAGAUUUUAAUCAACUAAUUCUUCAUAGAUCUUCACAUCUCAAGAAAGUGAAUUUGUUCUCCUAUUUCUAUCAUGUCUUAAAGUCUGAAUUGUUAUUAAUUUCCAAGUUGUUACUUCCUUUACUUCUUAGGAAUUAAAGAAAAAAGCAUAUUAAGAGGCAUUCAUCUAACCUAUCAGCUAAAUUGAUGAUAUCCCCUAUAAAUCACCUCUGUAUAUACAUUGAGUCAAAAAUUGUGCUAAUUAAGGACUGAUGUUAACAAUAGUAAGCAAAACAGACACAGACCCCACCUUCAUAGAUUUCCAGUUUGAAGGCAGAAACCAGCCAUUAAUAUAUCGUUUACACACAAGUGUGUAAGUUAUGCAUAAAAAUGGAAGAAUAGGGAGCUGUUAAAACCACAGGAAGCUGAUAUAAUGAGGAUGGUCAAGGCAGGCUUCUCUGAGGAAGUGACAUUUGAGCUCAGAUCUUAAAGAUGAAUACAUGGCAGCAGUAAGAGAAAAUGAGGAGGAAGCAAAAGUGGAAACCCCUGAUAAACCCAUCAGAUCUUGUGAGACUUAAUUGCUAUCACGAGAAUAGCAUGGGAAAGACUCCCAUGAUACAAUUACCUCCCCCUGGGUCCCUCCCACACACGUGGGAAUUCUGGGAGAUACAAUUCAAGUUGAGAUUUAGGUGGGGACACAGCCAAACCAUAUCAGUACCCUGACAGCAGAGGCUUCCAGUAGGUGCCAAACAGGGUCAGUUAAGCAGACAUCUCUCUAACAGUUAUCAGAAAUUCUCAAAUACUAAACAAAAUAUUUCUUUUGCAGGGGUAAGGUUAUUAUAGCUUGGGUGAGGGGAGUACUGUGCAUAAGAAUUGAAGCAAUUUAGCAACUGGUAAACAAUAAUAAUAAAGUAUAUAUUAUGUGACUAUCCCCUUUGCUGGUUGUAAAAUAUUUUUUAGCUUCCUACAGAGAGACCUGACCUAUAUAACUUUUUGCAUGAGCUCGUUCUUCUUAAUUUUUUUUUUUAACAAAAUGAUACUGUCAGCUUCUUGAGCUAAUUAUUAACAAGAAAUGCUAAAUGAACUAGUCUCCUCUCCCAUAAUUAGAUUUGGUAAACCUUGUAUUUGCAGAUUCUAAUGUUUUUACAUUUUAGGGGGAUUGAGAUUUAAAACAAUUUUUUAAUAACACGUAAAAGUUUUCUCAAACCACUUCCUCUACUCUCACAUUACAACAAUCAUCACCACACAAGAAAACUUCUGUGACCAUAUGUGUGUGGACUUUUCCCCCACACACCAAGCAGUGGACACCAGCUGGGUAUCCUCUGAAUCAGUUCUGAUACUCCCUACCCAGAUAUAGUGUCACAUCCCACAGAGUGAGUGCUCAGUCCCCAAGACUGCCCCUGCCCCCUCCACACACACAUACCAGUUGCAAGUCCAGGCCUCCAGAACUUCUGAUGGACUAGCUUCAAGUUGGGGCUCCCAAGAUCCACCUUUGGGUUUAAUUAAUUUGAUGGAGUAGCUCACAGAACUCAGUAAACACUUACUUAGGUUUACUGGUUUAUUAGGAAGGAUAUUGCAAAGGAUACAGAGACGCAUAGGGUGGAGUAUGGGAAAGGGGCACAGAACUUCCAUGCCUUCCCUGGGCUGCCACCCUCCAGGAACCUCCAGGUUUAGCUAUCCAGAAGCUACCUGAACUCUUUCCUCUUGGGUUUUUCUGAAAGCUUCAUGACAUCAGCAUUCCUUCCCCCAAGGUAUUGGGUAGGACCCUCUCAUGGGAGGGCCUUAAGACCCACAGUCAGAAAGGUAGGGGAACAUUAGAGUGAAAGGAGGGCGGGCCUGCCACUGAGGCCUAACACCCUUGACAUUUUUUUUUUUUUUUUUUUGAGACAGAGUCUCGCUUUGUUACCCAGGCUGGAGUGCAGUGGCACCACCUCGGCUCACUACAACCUCCGCCUCCCGGGUUCAAGCAAUUCUCCUGUCUCAGUCUCCGGAGUAGCUGGGAUUAGAGGCGCAUGCCACCAUGCCUGGCUAAUUUUUGUAUUUUUAGUAGAGAUGGGGUUUCACCAUGUUGGCCAGGCUGGUCUCAAGCUCCUGACCUCAGGUGAUCCACCUGUCUCGGCCUCCCAAAGUGCUGGGAUUACAGGCGGGAGCCACCGCACCUGGCCAACCCCCCACAUUUUAACAAGACUAUAACAAGGGCUAUGGGAAUUAUUAGCCAUGUUUUAAGUUUUUACUAAAAGAAUGAGCUUAUUUUGCUCUAUGAUUUUGUUACACUUUUCACCCAUUAUUGGUUAUUUUUCUAAGUCUUUUACAUGUUUAAUUAUGAAAUAUUUUUAAAGUAUAAAAAAGAAUAGAAAAUAAUAUAGGUGCAUUGUUUCUGCAUCCUUUUCCUUCAUUCUGAGCUCGGUUGCAUCUACCUAAAGUACACACUUAGUAUUUCAUUAACUUCAGGCCUUUAAUGUAAACUCUCAGUCUUUAAAAUGUCUUUAUUUCAGCUCUUACUUUUGAGUUAUAUUUUAGCUGUAUAUAGAAUUUUAAGUUGAUAGCUCUUUUUACUUAAGCAGUUUGAAGAUAUCUCAUUAUAAUUUUGAAGCAGCAUCAUUGUCUGCGGUGAAUACCCAGGGUUCGUCGUCUCGAGCUGAGAAAAUUAACGACACAGACACACAUACAUGAAGUGGGUUAAGGAGUGGAAAGUUUAAUAGGCAGAAGAAAUGAGAGAGGAGAGCAGCUCUCUCUUUUAUGAAAGAGAGGUGUUCCGAAAAAGGAAAAGCGGUGGACCGCAGCAGAUUUUAUAGGCAGGCUUUAAGAGGUGGUAUCUGAUUUACAUAGGGCCCACAGAUUGGUUCCACCAGGUGUGAUGUUUACUGGUGUGUGGGGAAGGCUGGUCACCCCACCCUAAUCUUAUUAUGCAAAUAGAUUUUCUACUUGGCCAGUGCCAUCUUGCCUGCUCCUUACUGUGCAUGUGGCAAAGAGAAGAGAAGAUGGAGCUGCCAUUUUGAACAGGCCUAUUCCACAGGUAGUUUUUUCUAUUGGCACAAUUGCCGGCAUUCGCUUGUGCAAGCUUCCAGCUUGUUUGUCUGUGUCUGCAGCUUGAUUUUACAGGCUGCUCUUCGUUAGAAAAAAAAAAAUCAUUUGGGGACUGCUUUUCAUUAAAAGGAAAACCUUACCAAGGACUUCCUUACCCUAUCUGCUUAAGUAAUUUCUUUUUUAACUUUUAUAUCAUUCUGACCUCUCGUUGCUUUUGAGAAAUCUGACAGUUCUUUAUGGGGAAUUGGUCUUUUCUGUUUUCCAUUUAAGAUUGUGUGUGUGUGUGCGUGUAUUCUGUAAUUUUUCUGGAUAUGUCUAAAUGUGGAUUUGCUUUUAUUCUGCUUUAGACUCCUUGUGCUUAUUGAAUCUAUAGAUUCAUAUCUUUACUCAGUUCUGGAAAAUUCUCAGCUAUUGUCUUUUAAAAGAUGCUUCUUCCAGCUUUCCCAAUUCUCUAUGGAACUUCUUUAGACAUCUUUGGCCUUAAUUCUCUUUCCAUGGCUCGUAACUCCUCCCUUAGUUAGUCAUAGCGCUUGCCCCUCUUUUAAAUCACUAUUUGUUACUUAGGAUUUUCCCUCCUUUCUUAGAAUUUCUUUUGUUAUAUUUUUAAAUAUACAAAUUUAUAUUUUAUAUACAAAUAUAUGUUUUAAAAUAUAUACAAAAUAUAUACAAAUAUAUAUUUUAAAAUAUAUACAAAUAUAUAUUUUAAAAUAUAUACAAAUAUAUAUUUUAAAAUAUAUACAAAAUAUAUACAAAUAUAUAUUUUAAAUAUACACAAAUAUAUAUAUUUUAAAUAUACAAAUAUUUUUGUUAUAUUUUAUCUAAUACAUCUACAUAUUUAUAUGAGGACUUUUUUGGCAUUAUCUUUACCAUUUCUGCUGGAACUGAUUCUUGACAGCUGUAGUUUUCUUGUAGAGUUGUUUUUUGUUGUUUUUUGAUAUAUUCCAAGAGUUUUUCAAAUUUUAUCUAUAAUUAAUAUUUGGAUAUUGUAGUUCAAAAGAGCAGGAGAAGACAGUGGUUAAUAUUUUGGACAGCCCAAUACUAUAGAUAAAAUGAAGUCACCUUUUCCUAUAAAAUGGAUUUUUUUUAUUGUUUUUCUAGAGCUUUGACUAUUGGGAAUUUUGCUUAGGAUCUCUAGUUUUAAAAGUUUCUGAAUCUUCCUGAGCCUCCUUCAGUAUUUUUGCAUGUAAUAAACCCCAAAGUUCUCUCAUGGGGAUUCUUUUAUAAAUGGAGAGAAUGGCAUAAGCAGUUAAGAGCAAGGUCUCUGGAAUCAAACAGCUUGUAAGUGUCCAGUCUUGUCUUUACCACUUAUUAACAUUGUGAUUUUGGGUGAGUUCUCUCACCUCUCUUUCUCAGUCUCCCGCUCUUUAAAGUAGCAACCUCGUAGAGUUGGGAUGAUUAACCAAGAUCUAAUAUGUAUUAAGUGGCAUCUAAGGAGCUCUCAAAAAUAUUAGAUAUUAUUAAUAUAAUUAUUUUCAUUGAUUCACAGUAACUUUUUCCUUGCACUUAAAGCUAUUUCAGGUGUAGCAGUUUUAAUCCGCCACCUGAUACAGCAUUCUUGAAGGCAAAAGCACUUGGACUAAUAAUAAGCAUAGGCCGUCACUUCCCCCGCAUUAUACAUUUUCUCUACUGAAUCUUUUCCAUUAGCAUACCAACUCAUAUUUCUUUAAAGCUAAAAUAACUCUUUUGUCAUCACUUCCCCCUCCAAUUUUCCCUUUACAAUAAAACUCCUUUAAUAUGUCCAUUUUUACACUCCCUAAUCUGUCUUUCAUGAAACCAAACUCCAGUCAAGGUUUUACCUCAUCCAUGCCACCAAAUCUGCUUGGAUUAAUAAACCCAGUGGUAAAUCUCAGUUCCUGUGUUACUUGACCCAUCAGCAGCAUUCAACACAGUUGAUCCCUCCCUUUUCCUUGAAACCUGGCUUUCUGCUCCCUCACUGGCUGCUGCUCCUUGGUCUCCUUUGUUGGUUCAUCCUCACCUUGAAAUGUUGGCAUGCUCCUGGAAUCCGCCCACGGACCUCUUCUUUAUUUGUCUACACUCCUUCUCAGUCUCAUGGCUUUGACUACGAUCUGUAUGUUCACCACUCCCAAAUUCAUAUCUCCAGCCUAGGACUCUUGCUUGACUCCAGAUGCAUAUAGCCAUUUGCCUACUCAGUAUCUCCAUUUGGAUGUCUAUUAGACAUCUCAAAUUUCACAUGUUCAAAACAACUUCUGAGAAUUAUGUUAGCUCUGGGUUUGUCAUACAUGGCCUUUCUACUGAGGUAUUAUUUCCUCUAUACCAAUCUGUUGAGAGUUUUUAUCAUGAAAGGGUAUCAAAGGCUUUUUCUGCAUUUAUUUAAAUGAUCCUCUGAUUUUUAGCCUUCAUUUUGUUAAUGUGGUGUAUCAUAUUUAUUGAUUGUGUAUGCAAUAAAUCCCAUUUGAUUGUGGUGAUUCUUUUAUCAAAAUCUUAAUUACCCCCUUUUAAAACAUUUGUAACAUGAUAAACUAUCUUACUAAACAUUUCAUAAUAAACAGUUGCUUCUUCCAAUUUAAGAAAAUGUAUAUUUUAGUUGUUGCUUGGGCUGAGGAUGUCCCUAAUCAAUCUAACUAUAUUUGUGUCUUCAGGCAUUAUUUGAGGUGAAUAAAAUGCAAAUAAUUUUUUAGACAUAAAGGCCUCUGUUAGCCAAAUGAGAAAUCUUUCUUUCACCUCUAGUUCUCAUUCACUCCCCCUUUUUCUUUCACAAAGGCAUUGAUUAAUUUAUAUCUUAGGUUAGCUAAGAAUUUGAGCUAUUUUUAGCAGUACAUUCUCAAAGUAUAGAUGCUGUGGAAUGAGUAGCAAGUAAUUAGAGCACCGUUUGUAUAACCUUUAUGUACAGUAGCCUUGAAAAAUCUCCGUACUGUAGGUAGAAUUGCACCUAAGAUGACAAAUCUUGUAGUAGGACCACCUAGCUAAUGAUUGAGUCUCCUUUAUAGCCUGGGAACUGAUGAGAAGGACACUCAGCACUCCUGAGUUCCUCCUGUCUCUACGUCUUUCUCAUAUUAAGCCUAAGUUGGCUUUCACACAAUUAAUUACAAUUUGCUUCUAUAGCUUUGGGGCUUUACGUAACAACUUUAAUUUUUUCUACACCCUUUUGAAAUAGCUAUCAUGUUUGACCCGGGUCUUUUCUCCUUAUAAAACAUUAGUAUUGUCUCCUGUUGUCUGAUAUAGUUUUGAGUCCCUUCAACAUCAACUGUAGCUCUGUCAGUGUCCCUCUUUAAAUACAGUGCUCAGAGUUAGCCUUAUAAUGUUUCAUGGUUUUGUUUUUGUUUUUAGAUACAAGGUCUUGUUCUGUCACCCAGGCUGGAGUGCAGUGGCAUGAUGGUAGCUCACUGCAGCCUUGGACUCCUGGGCUCAAGGGAUCCUCCUGCCUCGGUCUUCCAAAUCACUAGGAUUAUAGGGAUGCAUCACCGCGACCAGCCUCGUGUAUUAUUUAAUUAGUAGUCAUUCAACCUAUUUUUUUCUAGAUAGUACAUUUUGAUUAGUACAGUUUAUGAUUUUUUUUUUUUUUACCUGUAUCUUACUUAUAGUACUUACAGGAGCUGAGGAAAACCCCAAGACCAACCUGUGCUCAACCUGUGCUUUUGUAGUGAUUUUUUUGGACCCUAUUACAGAACCUUAUUUAUUCCUUCAUCUUGUUAAAUGUGGCCCAUUGUUCCAGGCUACCACCACAUUUCUUGAUUCUUUCUUUGAACAUAAUUAUUAUAUUUCCCAGCUUCGUGUCAUUUCUAAAAUUGAUGUACCUGCCUUCAUCUAGGUCAUUGCUAAAAACAAUGAAGUAGAGCUAAAGACAGCCUUUCAGCAUUCACAAAACUUAUCUAUUCUAUAUGCCAUAUAUUUAUAUAGCACCUUCUAUUUGCACUGUCCUAGUCAUGAGGGAAUCAUAAAUAAAAAGAAAUGGCCCUUUAAGAGCUCACAGCCUAGUGGGAAAACAGAUACACAGUAUAUAAUUUAGCAUGAUAAUUGCUGUAAUGGAAUAUGAAGAAGAUAGAAUGGGAGCACAGAGGAAAAACUCCGAGUCUUUGAGGCAGAAACUUCCUUGCAGGUUGAGAGCUGUCCUUUAAUCAAAAAAGUGUAAAUUAGCCAGGCACGGUGGCUCACCUAAUCCUAGCACUUUGGGAGGCCGGGGCAGAUGGAUCACCUGAGGUCAGGAGUUUGAGACCAGCCUGGCCAACAUGGUGAAACCCCAUCUCUCUGAAAAUACAAAAAUUAGCCGGGCAUAGUGGCGCAUGCCUGUAAUCCCAGCUACUCGGGAGACUGAGGCAGGAGGAUCACUUGAACCCGGGAGGUGGAGGUUGCAGUGAGCCGAGAUGGCGCCACUGCACUCCAGCCUGGGUGACAGAGCGAGACUCCCGUCUCAAACAAGUGUAAAUCCACCUAACCUGUCUCUUUUUUCCCAACUUGUUCACAAGGAUACAGUGAAAGACCCUGGCAUCUGCUUGCUGAAAUUGGGCCACAUUCUAUUUACUACAUUUCUCUAUUCCACCAGCCUAAUACCCUGUCCCCUCCUCUGCUCCCAAAAAAACAAAACAAAACAAAACAAAACAAAAAACUGUUAGUCUUUUAGAUGUAAUUAAUUAAGCUGUUUAUAUUUUAAAUGUAUUAAUUUUUACAUAUACAUACAUAUAUCUUCAGGAGCAGGAUCGGGGAAAUGAGUCAUAGAAUAUAAUAGACUCUGAAGAUCAGAAGGAAAAAAAAUAGAGGCAAAUUUGCAGGCUAUGAGAUCCCCCAGUUGUUAUAGUUGAAGAUCCAAAUUGGGCUCUGAGCUUUCUAAAGGCCAAAGAGAAAUGGAAAAGUUCAUGUUCUUACUGACUGACAGAGAAAUAGUUCCUGUUCUUAAGUGGGAACAAAGCUUUUGCCAACACCGAAAUCUUACAUUUCAGUGGGGGCUGCAUUGAGGUUACACCGAGAGAGAGAAUAAAAAAUGUCAUGAGCAUAUAUUUCUACUGCAAGUGCACAAAAGAUUUCCUUAGUGUUUCCUAAAUUUUGUUAUUGUGGAGGAAAUUGCUCAAAGAACAGAUUUAAUAUUACAGGCCAAGAGAGGCUGGGCAGCAUGGCUCAUGCCUGUAAUCCCAGCACUUUGGAAGGCCAAGGCAGGCUAGAUAGCUUGAGCCCGGGAGGUGGAGGUUACAGUGAGCUGAGAUCACGCCACUGCACUCCAGCCUGGGUGACAGAGUGAGACUCUGUCUCAAAAAAAAAGUGGGGCGGGGUGGGGGACUAAAAGAACUGUUAUAGAAAGGAAAGCUGCGCAACUUUGAAUAACAAAUUCAUUGGGGAUGUGAAAUUUUUUGUGCACACACCUAGACCAUCCCUCAGCACAAGUUCAAGAAGUCUUAUAUUAGGCUGUUUAUUAUAGCAUCCUUCAUAAGAGAGGAACACCAAUUCCAUGGUAGGUUGACAAGUAUUUCGUAGAUACCCAUGAUUCUGUGGUCAAAUAAAUAUAGGAAACGCUGAGUUAAACAAGUUUCUUUACUGCUGGUGGGCUUCUGUCUCAUAGGCUUAAGUACAAUGUAACUCUCCCAGAGGGAGUUCUCCAAAACUUUAUUUGGCCAUAGAACCCUUUUUGGUGGUACAUUAGGGAACAGAGGAAUAAAUGUUUACAUAUCCUCCAAACAACCUUUCUAAAAGAGAAAGAAUGUGGCAGCCACCCUGUUACAUAUUUGUCAACUCCAGCCAAGGAUAUAUGUUGGACAAUACUGUUCUUCCCCUCUCUUACUCUCUCAAGCUUCAGGAUGAUGUAGUCUCUUUUUUGCUAGGGUCUUCUCGUUUCCUUAUAAGGCAAACCAGAGUACAAGUUCCCCUCUUCACUUUCUCUUUUUUUUUUUUUUUGAGAUGGAGUCACUCUGUCGCCUAGGCUGGAGUGCAGUGGCGCAAUCUUGGCUCACUGCAACCUCCGCCUCCCGGGUUCAAGCGAUUCUCCUGCCUCAGCCUCUCUGAGUAGCUGCGAUUACAGGUGCACACAACCACACCCGGCUAAUUUUUGUAUUUUUAGUAGAAACGGGGUUUCACCACAUUGGUCAGGGUGGUCUUGAACUCCUGACCUCGUGAUCCACCCGUCUCCGCCUCCCAAAGUGCUAGGAUUACAGGCUUGAGCCACCGCGCUCGGCCCACUUUCUCUGCUUUUUUAAAGAUGAAACUAUGAUCAGGGCAAGCCAGGAGCAUGUUGUCACUUGGCUUUUUUGCCAAGUAAGACUUCUAGUUAGAGACCCCUGUCACUGCAAUUAACUUUGCUGAAUUUGUAACUUGAGAUCAGGGUUCACUUCCUUAUCUGGCCUCCAUCCCAUUCACAAUGUCAUUUUGUAUAGCAAAUCUAAUGUGUUUUUUUUGUUUGUUUGUUUGUUUUUUAGACGGAGUUUCACUCUGUCGCCCAGGCUACAGUGAAGUGGCGCAAUCUUGGCCCACUGCAACCUCCGCCCCCUGGGUUCAAGCGAUUCUCCUGCCUCAGCCUCCCUAGUAGCUGGGACUACAGGCGCGUGCCACCAUGCCUGGCAGAUUUUGUAUUUUUAGUAGAGAUGGGGUUUCGCCAUGUUGGCCAGGCUGAUCUCCAACUUCUGACCUCAGGUGAUCCACCCGCCUAGGCCUCCUGAAGUGCUGGGAUUACAGGCGUGAGCCACCGCACCCAGCCAAUCUAAUGCGUUUCUUUUAGAAUAAAGUAGGUUUUUUUCCCAUGUAUUAUAAACAAGCAAUAAAUCUGGGGAGGUUAUUUCUACUUGAUUAGUUUUUAUUUACCUAAUAAGAAUCCUAAAAUGCAAAAAUAGUAUUGUAUUAAUAUAACAGUGACUUUAAAAUAAACUUCAGCAAUAUUACAUUUGAAAAAUACCCUAAGGUUAAAUUUAUUUAUAACAAGUGUUAUUCAUUUAUUAUAAUACCCCAGAGUCACAUUCAGUAUUGUCUUUCGAUAUAUGCAUUUUAUAAUUGUAAAUCCAUAUACAGAUGUCCCCUACUUAAUGAGGGUUCGACUUACUGUUUUUCAACAUCACGAUGGUAUGAAAGCAACACACAUUCAGUAGAAACAGUACUUCGAUUUUUGAAUUUUGAUCUUUUCCUGGGCCAGUGAUAAGUGGUAUGGGGGCAGCAGCAGCGAGCCACAGCUUGCAGUCAGCCACACCAUCAUGAGGGUAAACACCGAUACUCUCCUGUGGACUGUGAUGCAGAUGAUUCUGCCCGACUGUAGACCAAUAUACAUGCUCUGAGCAUGUUUCAGAUGGACUAGGCUAAGCUAUAUGUUCAGUAAGUUAGGUGUAUUCACUGCAUUUUUGACUUAAUGAUAUUUUCAAUUUACUGCGAAUUUAUUGGACAUAACUCCAUUGUAAGUAGAGGAACAUUUCUGUCUAAUAUUUUGUGUUUUCCUCCCAUUAAAUUCCUUAUUCUCUUAUGUUAUGGAAUUCCUUAUGGAAUACAUAUUCCCUUACGUUGACAUUGACAACAUACUGGUAAUCUUUAAUAUAUAUGAAAUACUUGAAUAUUGUCAAAUGUUAUUUUUGAAGUCAUUUUGAUAUUUGGCAAUUGAGGUUUUUUCAACUUUUUGCAAUUUUAAAAGGGGUAUUGAUAUCAUAUGAAUUAUCUUUUUUCUCUGGUUUAUUUCCUUGGGCUGUAUUCCUAUGGAUAGGACUACAAAGUCAGCAAUUAUAUUUUCUAGCUUUUGUAUGACAGGAAUUUUCUAACCUUAUUGUCUUCUAAAUCCUAGGUUCAGGAGCUCCAUGAGGAGGAGGCACAGUGGGUGAACUAUGAUGAGGAUGAGUUGUGUGUGAAAAUGCAGCUAGCCGACGGGAUCUUUGAGACCCUGAUCAAAGAUACUAUUGAUGUUCUGAAUCAGAUCAGUGAAAAGCAGGGGAGAAUGCUACUUGUGUGACAUCUUGCAAAUAAAUCGAACGCUGAGUGCUAAUGUGAGUCCUGGGCCUUUCUGCCUCCUGAUGUACACCCAUCGCCAUCAUAGCAAGAGUGCUUCUGGACCUUGUACUUAUUCUUAAAGACUACCAGUAUGGAGUUCAUAGGACAAUGUGGUACACCUGGUAUUACAGCCUUUGCCUUUCGAGACUAUCCACUGGAUUAAUGGGUUAUUUUCAGUGGGCAGGGUUGCACAGUGUAAUCCUACACCUUUUGCUAACACCCCUACUAGGUCCCAGAGGGCCAGAAACACCUGACUUACCUCUGAGUUUAGACUAGGGUAUCACUUCUUUGAGUCUGAAGUCAAGUGAGAGAGGUAUAGAUAAAUGCAUCACAUCACUUUUGAAAUGUAAUUCUGGUCUAUACCAUGGAAGUCAUAAAUGGACAUUAUAGUCUCUAAACAGUAUUAAACCCUUAACCACUUCUAAAAUAGGCAAGCUCAAUAAUGUCUGCCAACUUCACAUUCUGGAGUUUAUUUCAUUUCUUUUUGAAGACCAUUUUCUUCCAUUAUUGUAGUUGAGCAGCACCAAGUGGACUGUCAGGCUAACAGGAAUAAGUGGUAGCCUUGCUUUCUGAGCACCAUCUAAAGAAUUUUAAACCUCUGCAUUAUGUUUAGUGUUCUCUGUGUGGGCAUGAAAACAAAGAAAUGCCCCUACUGAAGACUGGGCUCAAAGGACCAAUGCAGGGCUGGUUCUUUUACCCCUUGGUUUACUCCUGCACUUGUCUUACUCAUUCUGAAUCUCUACCAGCUGCUCCCAGAAUCACAGAUACUCAGGACCAUCUCAGGCAUCCAGGCAUGGCAAAGUGGAAAUAAUUCAUUUUGGCCUGCAAACCUUCAACUCCCAUCUUUCCCCAAGAAGUCAGAGAUGCUGUUACUUGAAUGAUUUAGGAAAUGAGUGUGUGCACCAAAGACAAAAAGAUAUUGUCUAUUGUUUGUGUGCUUGUUUUGCGCUAUGGAACAUUUUUUAAUUUAUUUUAAGAUAAAUUAUUAAGUUGAAAAUGUGUGUCCCUAUUCAGAAGUGAAAGAUUCAUCUUGUAAUAGUUAAACCUCCAUCUUGAAGCUUCUAUGGUUCAUAGUCUUUGCACAGGAACCUGUGGUUUUAACAAACCAAUACACAUAUUGAAGAAGUCAUUUUAAUUCAGUGAAACGAAGAUGGGCUUUUCCAGAUCACCUGCAAUAGCAGCAGUGGGAUAAAAUGAUUUAAAAACACUGUACAAUUUAACUCUGCCUCUCUUGCAGCAUUGCUUCUCACAACUAUUACCUGCAUCUGAAAAAAAAUCUAUAGACUCCAGCUGCUACAUUAGAGCAUAAGAGAUGCUCUCCUGGGACCUCAGUACCCUGCCUUCUUGACUGGUUUCCGUUCAUCAGUCCUGUCCCUCUUCAAGUAAUCUAGAAGAAUGUGGAUACUCUUAGGCGUGAAUGUAAAUGCCUUAAUAUUGAAGGUCCUGGUUAGAAGCAUGAUACAAGACAUCUACUGGAUUCAUAUUUACAAAUAUCCUGGAAUGUUAUAGCUUCAAAGUAUAUUAGAAAAACCCCAAAGAUGGUAUAAUCUUUAAGUGUGCACGUUCGUUUAUUUCUGCAUCUUCCCUCCAAACUUGCCUUUGCAUCUUAAAUAUUUCACUAUGCACACUCCCAUUCCUCUUGGGUUUCAUCUUGUCGUUUAAGAAAUGUACUGAAAUAAUCAUUGGAAUAUUUGCAUUUUGCACAAUGACUGGUAUGAUAGCUCUUGACAAAUAAGGAAAGCACUGAAAUGUUGUGAUUGGGUCUCGGGAAAUGCUCAGAUUGAUGUCUUACCAGCAUUUCUUCUGGGCUUGUGAUGUUGAGCUGUAGUCUUGUAGCCAUAAUGAGCAAAUUGACUAAGAGAAGCAAAGGUUUCUUGGGGUUAUUAACCAGUAGUGUGGAAAUACUAGUUUUAUGUGGCCAAGGAAAAGCAAAGGCUUUUCUUUUCAGUUUGUGUUAUUUGGAAGACAGAAAAACAUCUUGUCUACAUCCUUUGGCUGUUUGUAGGAUCACGUUGUCCUUACGAUACUGAAACUUUACAGCUGCUGUAAAUUUUUUAUAAAUGAAUUUCAAAAUGUUAUAAUGGGACUGUAGGUUGUUUUUCUACAUCUUCAUUAUUUGGACCUAAAACCAGUUUUUAAUAAGAAAGUUUAUCUUUACUCUUUCUGAAAUUAUGACUCCAGAAAAAGAAAAAAAAAAUACAAGUCAUGGAAUCAGCAAUCUGGUAAGAAAUGCUGCCAAGAAUGUGGCAGUAGCUGUCCUGACAGACUCCAACUGUCUUUACUAUCUGAAGAAUCCUAGGCUCCACAUGAGAGGCAGAAAUGGAUCAGUCUUAUUCUUUUCUAGAAAUGGUUAUCUGUAGUUUGGUAGCAAAAAAAAAGAAAAAAGAAUCCAUAAUUAGCAGAUUUCUUAUUAACUAUUUGGAUCUAAUUGAAAUGGCUUUAUUCUUAGGAUUAAGAAAGAUAGAUGUGGAUACCCAGCCACUCGUUCCAUAUUGGUAUCUUUUUAAAUCAGCUCUGCCUCUUAAUCAAGAACCUAAAUAUUCCCUCUUUCUAAUCUUUGUUCCUUCUCCCUACACCCUCAUCCUCUUUCACUCUUCCUUCAUAAUUCCUCUAAGAAAAAUAUCUUUGCAUCAGCAGUAAUAUCUUUUAGAAUAGCACUAUCAGAAUUUAGCAGUAAACCAACAUACAGGCUUCAGAUUUACUUCUGAGUCCAAAACAAUUUGUGCUAUCCAGGGUAGUUAACUCUGGGUUAAACAAGUACAGGGUAUAGAUUCCCUCUUCAGGUCUACACAGGAAUUUUUACCAUAGGGAAAAGUGGGGAGAGCUCAAACGUAGUUAAUAAGGAAGGUAAUUUGUUUUUCUUUUACCUAAAAGAAAAGAAAAUUCCUUCUGUGACUACAGGUCUCUGAGAAAUUAUCUUUCAAAAGAGAUUUCAUUGCUCAUAAGAGUGUUGUGGCCUAUUGAUAAAAACAAUUUUGUUCAGUUUCUUGUCUUGAAAAAAAAGUGGCCUUAGCUUUUUGCAAUACUUGAAUAAAGUGUGUACUCGCAAAAGAAUUUCUGUAGCACAGCAUUAGAGACUCAUAACUUUUCUGCAAGAAAUACAAACUUACAUCUUCCUUUUACUACCUUAAGAAUACUAGUGAAUAAAACAUUAAUUCAAAGAGCAAAUUAUAGAAACUACAAUGACAUUUAAUGCAAAUUGUAGGAAUUUACAUGUUUACAAAUCAUCUUCAACUGGUUGUGCAGCAAUUCAAUAAAAUAUCUUUGUAUUAUAAAAAUGUGAAGAAAAAAUGUAAACUGAUGUAAAGGAGGUACUGUCAUUUUAAUUAACCUAUGUUUAAUAGCUUUUCCUUCUGGACUUUGCAAAGCCUUCUUGGCAAACACAUUGCAAAGCAUUCUCUGGGAGGUUCAGCCUCCUUGUGUGUACUGUACUGUGCAGACAUGAAAAAAUAAACCCGUUUACUGUGUGCGUGUAAAUAGCCUGGUCAUCAGGCCAUUUUCAGCCAAUAGUCACAUCCAGUGCAAUUUUGCACCGAACACUUAAGGGUGUGGUUUGUAAGUACGAUCUGUAAAAUAACUGGGAUGAAUUCCCAUGUAUACCUGUGUAAAUAGAUUUGUUAACUGAAAUAUACUUUAAGAAAGAUAAAAUCUGUAAAUAAACUGAUUUAUAAAUUAAUUUCA';
		my $expectedCdsWtSeq = 'TTTTTTTCACTAAAATTTGGCTTTTCTCTAAACTTTGGTTCACAAAATTGAAAGATAGTGATTAAATAATTTTGATGTGGTGATTTGTTTTTTCTTTAAGATGCTTCACTGTCTGGTTCTGAGAGATCAGTATCAGAAAGGTCTTTATCTGCATATGCAAAGAGAGTAAATGAATGGGACAGTCGAACAGAAGATTTTCAGACCCCATCTCCAGTTCTCAGATCATCAAGGAAAATCAGAGAAGAATCTGGAGATTCTCTAGAAAATGTACCTGCATTACATCTTCTCAAAGAATTAAATGCCACTAGTAGAATTCTTGATATGTCAGATGGCAAGGTTGGAGAATCTAGTAAAAAATCAGAAATAAAAGAAATAGAGTATACAAAATTGAAGAAGAGTAAGATTGAAGATGCCTTTTCTAAAGAAGGTAAATCTGATGTCTTACTGAAATTAGTCCTAGAACAGGGAGATTCATCTGAAATTCTTTCAAAGAAAGATCTTCCTTTAGATTCTGAAAATGTTCAGAAAGACCTAGTTGGATTAGCTATTGAAAATCTCCATAAAAGTGAGGAAATGTTGAAAGAGAGACAGTCAGATCAAGATATGAATCATAGTCCAAACATCCAATCAGGAAAAGACATTCACGAACAAAAGAACACAAAGGAAAAAGATTTGTCTTGGTCAGAACATCTTTTTGCTCCTAAAGAGATACCATACTCTGAAGATTTTGAAGTGTCTTCTTTCAAGAAAGAAATTTCAGCTGAATTGTACAAAGATGATTTTGAGGTGTCATCTTTGCTGTCACTCAGGAAAGACTCTCAGTCTTGCAGAGATAAGCCACAGCCAATGAGGAGCTCTACAAGTGGAGCCACTAGCTTTGGTAGTAATGAGGAAATCAGTGAGTGCCTAAGTGAGAAAAGCCTTTCTATCCATAGCAATGTTCATTCTGACAGGCTGTTGGAACTCAAGTCCCCTACTGAGCTGATGAAAAGTAAGGAGCGCAGTGATGTGGAGCATGAACAGCAAGTTACTGAATCCCCTTCCTTGGCTTCAGTTCCTACTGCAGACGAGTTATTTGATTTCCACATTGGTGATAGGGTGTTGATTGGAAATGTTCAGCCAGGAATTCTTCGATTCAAAGGTGAGACTAGTTTTGCTAAAGGATTTTGGGCCGGAGTGGAGTTAGATAAACCTGAAGGAAATAACAATGGAACATATGATGGTATTGCATATTTTGAGTGCAAAGAAAAGCATGGTATTTTTGCTCCTCCTCAAAAAATATCTCACATTCCAGAAAACTTTGATGACTATGTAGACATTAATGAAGATGAAGACTGTTACTCAGATGAACGATATCAGTGCTATAATCAAGAGCAAAATGATACAGAGGGTCCAAAAGACAGAGAAAAGGATGTCAGTGAATATTTTTATGAGAAATCCCTACCTAGTGTGAATGATATAGAAGCCTCAGTTAATAGAAGTAGAAGCCTTAAAATAGAAACAGACAATGTACAGGACATTTCTGGGGTACTTGAAGCCCATGTTCACCAGCAGTCTTCAGTGGATTCACAGATTTCTTCAAAGGAAAACAAAGACCTCATTTCTGATGCCACAGAAAAGGTTTCCATCGCTGCAGAAGATGACACTTTAGACAATACCTTTTCCGAAGAATTGGAGAAGCAACAGCAGTTTACAGAAGAGGAAGACAACCTATATGCTGAAGCTTCAGAAAAGCTTTGTACACCACTTCTGGATCTTTTAACAAGAGAAAAAAACCAACTGGAAGCCCAGCTGAAGTCATCACTAAATGAGGAAAAAAAGTCAAAACAACAACTGGAAAAAATCAGCTTACTGACAGACAGTTTACTAAAAGTCTTTGTAAAGGACACAGTCAATCAACTACAACAAATCAAAAAAACCAGGGATGAGAAAATCCAGCTTAGCAATCAGGAGCTTCTTGGTGATGACCAAAAGAAAGTAACACCCCAAGACCTATCCCAAAATGTTGAGGAACAGTCGCCAAGTATTTCAGGTTGCTTCTTAAGTTCTGAATTGGAAGATGAAAAAGAAGAGATTTCCTCTCCAGATATGTGTCCCAGACCGGTGAGTATTTCGTCTAGAAAAAGTATCAGTATACCTTTTGACACTTTTAACAAGGCAACTTACCTTGCCTATATAGTACATTTTAATTCTGTGACTCCTGATAGGTTTTTTTTCCCTTGACGTGTTGTTTATTAGAGAAGGATCACTGTACTGGTTTTTAAAGGCCTTGTAGCTGTTTTTCCTGGCTGAAGTGTAATTGTCCTCTATTATGTCATGAGTAGTCAATATAGGGTTGTCACAAGCACTGGCAGGTTTTGGGGAACTAAAACATCCATTTAAACTGCTGAGTTGACAACCCCTGTACTGTCCATAGATTACATCTCCTTCCAGATTCACTGAACAAATTCGTTCCATGGAGCTTTCCAGTGGCTCACCACCCAATTCTTTTTTGCCTGTGAAATTGTTTACATATCTACTCAGTCATGGTTGGCCTGCAGAGGTTAGCAGTTGTGATGTTGATTCTGGTAAAATGAAATAGATCCACTATGACATGCACTCAGACCATCAGCTTGTTTATGTTGCCCATGGACAGTTGAGCAAAAGTAATGACTTTGAGTAACTTGAGCCAGATGGTCACCAGAGGTGAGAAATTTTATATCCCATTTCCATATCCCATAGACAGTTCTTTAGAAACAGAATAATGATAATACAAAAAGCGGGGAGAGCCCTAATCAATAAACATTTTATAAAAATTCAGTTATATGCAAGGTTTATATAAGTAGAATTGAGGGAACAGGATAACTTTTTCTTAGATTGTTTATTAGTAGCCACCATATTGTTTTTTATTCTGCTTACTTATAGAGAAATATTAGATATATAAATAATGCATTTACTCTGTCTTTTTTTAAAAAAAAATGATATGAGGTCCCTAAATGGTGCTCTGTTTGAATGGTTCCATTTCTATAGGAGAGCCCAGTATTTGGTGCCAGTGGGCAGGAAGAACTTGCTAAGAGACTTGCTGAACTTGAACTCAGCCGGGAGTTCCTGAGCGCGTTAGGAGATGATCAAGACTGGTTTGATGAAGACTTTGGTTTGAGCTCTTCTCACAAGATCCAAAAAAATAAGGCAGAAGAAACCATTGTACCTCTAATGGCAGAACCTAAAAGAGTAACCCAACAACCATGTGAAACATTATTGGCAGTCCCCCATACTGCAGAAGAAGTAGAGATTCTTGTACATAATGCAGCAGAAGAACTTTGGAAATGGAAAGAATTAGGCCACGATCTTCATAGCATCAGTATTCCTACAAAACTGCTTGGCTGTGCCAGTAAAGGTCTAGATATAGAAAGCACTAGTAAAAGGGTCTACAAACAGGTAGGTGAAATAAAAGGATAATTTTAGTTTTTAAACTTTTCTCCCATTGTCTGAAATTTGGATTCTTAGCCATTTTGTTTTGTTTTGTTTTGTTTTTCTCTATCAAGGCGGTTTTTGATTTAACAAAAGAGATTTTTGAGGAAATATTTGCTGAGGATCCCAACTTAAATCAACCTGTCTGGATGAAGCCATGTAGAATCAACTCTAGTTATTTCCGACGAGTGAAAAATCCAAATAACCTTGATGAAATCAAGGTAAACTGCAAACTATAAAGTGTCTTCTTTTTTGACTTGCTGTTCATTTACACATATCTGTAAATGTTGTAAAAGCAGGTGCTCTTCTGTATGATTAACCTTTGTGTCTTCAGGACCTGTTCCATAGTGGGTATTCAATAAAGATTTGTTGACTTACTGAATCCATGAACAGTGGAGGGATTGGAATGAGTGAAAGGGCTGAATCCTCAAGTTTGTGCAAAAAAAAAAAAAAATTCTTTTACCAGTCATTTTCTCTATATACATAATTATAATATTTAAGAAAAATTAGGAAATACCTAGGGAGATAAATTATTGCAACTAATAGGAGTTTTAGGTAAGTTTTGGGGGCTTTAGAGATAGGTATAATGTTCTAAAAGTACATATGAATACGTGTCATATATATGTCTTTTCAATATGCCAGCCCTGGACACATAGGTGTTAACTTCTTTATCTTTGACTTTGGAGGTTTCAAATTGGTAATGTGGTTTTGAACTATGTAGGATAGTAATGTGATACGTAAAAATATTTATTTTGTTATCTAAGAAACCTTTGTAAACCTTTGTCCAGAAGCTACAAGGACTATGAACTTTAAAGTCACATTGCTGGCAATCTGCTCTGCCTAGTCACATGTTTCTTTACGGCAAAGGTTATATTTCCAGAGCAGTGATTCTCAACTTTTAATGTGCATGTGAATGATCTGAGGATCTTGTTGAAATGGAGATTCTGATTTAGTAGGTCTGGAGTAGGACCTGAAATTCTGCATTTCCAACAAGCTCTCAAGTCATGCTGGTCCTGCTGGTCCTGCTGGTCCTGCTGGTCCCAGGACCACGCTGAGTAGCAAAGCTTTAAACTTGTATTTCCCAAACTTATTTGATCACAGAAAGCCTTTAAAAATATTTTATACTTATAGAACACCAACAGATATTTATGAGATGGGATGTATGCTTACCTCAGGGGTTGGCAAACTATGGCCCATAAGCCAAATCCAGGCCCCTGTATGTTTTTATAAATAAAGTTTTGTTGGACCACAGCCATACCTAGTTGTTTACATGTTACAACTACATGTAACTACAACAGTGGAATTGAATAGTTGTGACAAAGACCACATGGCCCGCAAAGCCTAAAATATTTACAGTCTGATTTTTTACAGAAAAAATTGGCCATCCCCTGGCTTAAAGTAACTAGCTTTTTAAAAACACATTGTATACTCCATTCCTTCTACCCTTTCACAATTGAAAATGCTTGAATGCTATGTGTGTCTGTGAACAGGTAGGGAGTGATGGAATTTCCGGGCACGTTGGTACTTTGAGAATCTTTCCAGGTGATTCTGATTCATAAATGTGTCCCTTCCCCTGTTAAGAACATTGTTTTTTTAACAGAACAACCATAAAAGTAACAAAATAAATTATTAATAACATTAAGTTGAGGTAAGACAAGATATATAATATATACCAACAATCAAACTACATGCATCAATATAGCATACAATAAATGCCCACTAATTTTAATATTGACATAACATTTTGTAAAGTTTTTAAGGGTGATAGTTCAAATATTGTTAGTAGTAAAATTATTATGCTTCTCATACTGAAGTGTTGATATTTTTTAATTTAACATATTTGAGTTTTGGTGCTTCAAAATAAACTATCCCAAAGTCTTTTGCTTTTCTCACATTACATTGTCCCCAAAAATTTCTTTGAGATTTTTCCAATTCAGAGATTTATTTTATGCTTTCTTTAGTAGTTTTGCCAGCTTGAAGTTAGCTAATGGTAGTTTTGTTGTTGTTTTTTCTTTGTTTTTCTTGAGGTGGAGTTTCACTCTTGTTGCCCAGGCTGGAGTACAATGGCGTGATCTCAGCCCACTGCAACCTCCGCCTCCTGGATTCAAGTGATTCTCCTGCCTCAGCCTCCTGAGTAGCTGGGATTACAGGCATGCGCTACCTCACCCAGCTAATTTTGTATTTTTAGTAGAGACAGGGTTTCTCCATGTTGGTCAGGCTGGTCTCGAACTCCCGACCTCAGGTGATCCGCCGCCTTGGCCTTTCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCACCCAGCATGGTAGTTCTTAAGTAAACTCTAAAGTGATTTTTCTTTGCAGATTAATCAGTGGGTTTTACTTTGTTGTTGTTGTTGTTGTTGTTATTTTTAGATAAGATCTCGCTTTGTCACCTAGGCTGGAGTGCAGTGGTGCAATCTCAGCTCACTGCAGTCTGGACCTCCTGGGCTCAATCAGTTCTCCCACCTCACCCTCCTGAGTAGCCAGAACTACAGTCATGTACCACCATGCCCAGCTAATTTTATATATTTTTTGTAGAGATGAGGTCTCTCACTATGTTGCTCAGGCTGGTCTCAAACTCCTGGGCTCAAGTGATTCTCCTGCCTTGGCTTCAGAGAGTGCTGGGATTACAGACGTGAGTCACCGTGCCCAGCTAGGTTTTACTTTTATATGTTGCGTAAGTGAATAAAGTCTCAGCATTTGGAATGCTTTTCCCTTCAAATTGCTATGTGGAACTGTAAGTGGAAGAAATTTGTATGATAATTTTATTGTAAATATTATATTCTTATCGTGAATCTGTTTACTGCATTGTAAACACATGAAACTCAGGAATTTCTGTTGTTTGTGTTTCGTGTATTTGTCTTGTTTCCACAGAGCTTCATAGCAAGTGAAGTACTCAAGTTGTTCAGTCTTAAAAAGGAGCCAAACCACAAAACAGATTGGCAGAAAATGATGAAATTTGGAAGAAAGAAAAGAGACCGAGTGGATCATATCCTGGTCAGTGTATACAACCAAACTGTTTTTATTTTGACCATATCTTTTAAACTCAGAATGCAGTTAGATTCTAAAAGGAAAAGGGATAGACTTATTCTTAAATCTATGGTATTAAAATCTTAAAGCCTTTACTTGCTGCATTTCTACCCCTCAGTTTTGTTAACTATTTTCTTTTGTACTCAGTAGATGAACATTTAACTTGATCAATAGTATTTTATGTAGAACTTGAGAAATTTACTCAAATGGAAGCAAAGTTTAAGGTAGTGAATTGAAAAGAGACCAGTGTGGACTCTGGAGCCTGACTGCCCATGTTTAAATCTCTTGATGTGTAATCTTGGACAAATTTGTCTAACATTTTTGTGCAAAGTGAAGGTAACAATAATAGTACTTACCTGATAAAATTATTAAGAAAATTAAATAAGTATATGTGTAAAGCACTTAGAATAGTTCCTCATATACAACTCAGTAAATAATAGCTGTTACCTCATCATCATCACTATTATCATTATCCAGAACATTCCCTGGATAGGATTTATTTGGTTTTGCAATATATTTCTGTAACATTAGGAAAGTGTTTATCAAAAAAGAAAAAAAGAAAAATTCAGGTTTCACTCTAAACATTTTAATTCTGGACCTTTTATCTATATTGACATTTATATAATGAAAAAAAGTGTAATACTTAAATTTATCACCGTTCTCATGTATGTACCTCTCTAGTGTTTCATTTACTCAGTAATAACTGTACTAGATAACATTTGCTTAGTATTTTATAGTTTATCATTCTTATGCCTTATTTGATTTTTTTTTTTTTTTTTTTTTTTGAGACACAGTCTCGCTCTGTCACCCTGGCTGGAGTGTGGTGGCACAATCTCTGCTCACTGCAGCCTCTGCCTCCCAGGTTCAAGTGATTCTCATGCCTCAACCTCCCGAGTAGCTGGGATTACAGGCACCTGCAACCATGTCCAGCTAATTTTTGTATTTTTAGTAGAGATGGGATTTCGCCATGTTGGCCAGGCTGGTCTCAAACTCCTGACCACAAGTGATCCGCCTGCCTTGGCCTCCCAAAGTGCTGGGATTACAGGCATGAGCCACCGCGCCCAGTGTCTTATGCCTTATTTGATCCTCACGATAATTACATGAGGTAGGAAGAACTAATGTTATCCTTTGTTTTACAAATGAGGTGTCTTACGCCTTATTTGATCCTCACAATAATTACATGAGGTAGGAAGGACTAATGTTATTCTTTGTTTTACAAATGAGGAAATTTAGATTCAAAGATGTTAAGTGACTTGCCCAGTAAAATAAGAGACGAGATTTTAATCAACTAATTCTTCATAGATCTTCACATCTCAAGAAAGTGAATTTGTTCTCCTATTTCTATCATGTCTTAAAGTCTGAATTGTTATTAATTTCCAAGTTGTTACTTCCTTTACTTCTTAGGAATTAAAGAAAAAAGCATATTAAGAGGCATTCATCTAACCTATCAGCTAAATTGATGATATCCCCTATAAATCACCTCTGTATATACATTGAGTCAAAAATTGTGCTAATTAAGGACTGATGTTAACAATAGTAAGCAAAACAGACACAGACCCCACCTTCATAGATTTCCAGTTTGAAGGCAGAAACCAGCCATTAATATATCGTTTACACACAAGTGTGTAAGTTATGCATAAAAATGGAAGAATAGGGAGCTGTTAAAACCACAGGAAGCTGATATAATGAGGATGGTCAAGGCAGGCTTCTCTGAGGAAGTGACATTTGAGCTCAGATCTTAAAGATGAATACATGGCAGCAGTAAGAGAAAATGAGGAGGAAGCAAAAGTGGAAACCCCTGATAAACCCATCAGATCTTGTGAGACTTAATTGCTATCACGAGAATAGCATGGGAAAGACTCCCATGATACAATTACCTCCCCCTGGGTCCCTCCCACACACGTGGGAATTCTGGGAGATACAATTCAAGTTGAGATTTAGGTGGGGACACAGCCAAACCATATCAGTACCCTGACAGCAGAGGCTTCCAGTAGGTGCCAAACAGGGTCAGTTAAGCAGACATCTCTCTAACAGTTATCAGAAATTCTCAAATACTAAACAAAATATTTCTTTTGCAGGGGTAAGGTTATTATAGCTTGGGTGAGGGGAGTACTGTGCATAAGAATTGAAGCAATTTAGCAACTGGTAAACAATAATAATAAAGTATATATTATGTGACTATCCCCTTTGCTGGTTGTAAAATATTTTTTAGCTTCCTACAGAGAGACCTGACCTATATAACTTTTTGCATGAGCTCGTTCTTCTTAATTTTTTTTTTTAACAAAATGATACTGTCAGCTTCTTGAGCTAATTATTAACAAGAAATGCTAAATGAACTAGTCTCCTCTCCCATAATTAGATTTGGTAAACCTTGTATTTGCAGATTCTAATGTTTTTACATTTTAGGGGGATTGAGATTTAAAACAATTTTTTAATAACACGTAAAAGTTTTCTCAAACCACTTCCTCTACTCTCACATTACAACAATCATCACCACACAAGAAAACTTCTGTGACCATATGTGTGTGGACTTTTCCCCCACACACCAAGCAGTGGACACCAGCTGGGTATCCTCTGAATCAGTTCTGATACTCCCTACCCAGATATAGTGTCACATCCCACAGAGTGAGTGCTCAGTCCCCAAGACTGCCCCTGCCCCCTCCACACACACATACCAGTTGCAAGTCCAGGCCTCCAGAACTTCTGATGGACTAGCTTCAAGTTGGGGCTCCCAAGATCCACCTTTGGGTTTAATTAATTTGATGGAGTAGCTCACAGAACTCAGTAAACACTTACTTAGGTTTACTGGTTTATTAGGAAGGATATTGCAAAGGATACAGAGACGCATAGGGTGGAGTATGGGAAAGGGGCACAGAACTTCCATGCCTTCCCTGGGCTGCCACCCTCCAGGAACCTCCAGGTTTAGCTATCCAGAAGCTACCTGAACTCTTTCCTCTTGGGTTTTTCTGAAAGCTTCATGACATCAGCATTCCTTCCCCCAAGGTATTGGGTAGGACCCTCTCATGGGAGGGCCTTAAGACCCACAGTCAGAAAGGTAGGGGAACATTAGAGTGAAAGGAGGGCGGGCCTGCCACTGAGGCCTAACACCCTTGACATTTTTTTTTTTTTTTTTTTGAGACAGAGTCTCGCTTTGTTACCCAGGCTGGAGTGCAGTGGCACCACCTCGGCTCACTACAACCTCCGCCTCCCGGGTTCAAGCAATTCTCCTGTCTCAGTCTCCGGAGTAGCTGGGATTAGAGGCGCATGCCACCATGCCTGGCTAATTTTTGTATTTTTAGTAGAGATGGGGTTTCACCATGTTGGCCAGGCTGGTCTCAAGCTCCTGACCTCAGGTGATCCACCTGTCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGGGAGCCACCGCACCTGGCCAACCCCCCACATTTTAACAAGACTATAACAAGGGCTATGGGAATTATTAGCCATGTTTTAAGTTTTTACTAAAAGAATGAGCTTATTTTGCTCTATGATTTTGTTACACTTTTCACCCATTATTGGTTATTTTTCTAAGTCTTTTACATGTTTAATTATGAAATATTTTTAAAGTATAAAAAAGAATAGAAAATAATATAGGTGCATTGTTTCTGCATCCTTTTCCTTCATTCTGAGCTCGGTTGCATCTACCTAAAGTACACACTTAGTATTTCATTAACTTCAGGCCTTTAATGTAAACTCTCAGTCTTTAAAATGTCTTTATTTCAGCTCTTACTTTTGAGTTATATTTTAGCTGTATATAGAATTTTAAGTTGATAGCTCTTTTTACTTAAGCAGTTTGAAGATATCTCATTATAATTTTGAAGCAGCATCATTGTCTGCGGTGAATACCCAGGGTTCGTCGTCTCGAGCTGAGAAAATTAACGACACAGACACACATACATGAAGTGGGTTAAGGAGTGGAAAGTTTAATAGGCAGAAGAAATGAGAGAGGAGAGCAGCTCTCTCTTTTATGAAAGAGAGGTGTTCCGAAAAAGGAAAAGCGGTGGACCGCAGCAGATTTTATAGGCAGGCTTTAAGAGGTGGTATCTGATTTACATAGGGCCCACAGATTGGTTCCACCAGGTGTGATGTTTACTGGTGTGTGGGGAAGGCTGGTCACCCCACCCTAATCTTATTATGCAAATAGATTTTCTACTTGGCCAGTGCCATCTTGCCTGCTCCTTACTGTGCATGTGGCAAAGAGAAGAGAAGATGGAGCTGCCATTTTGAACAGGCCTATTCCACAGGTAGTTTTTTCTATTGGCACAATTGCCGGCATTCGCTTGTGCAAGCTTCCAGCTTGTTTGTCTGTGTCTGCAGCTTGATTTTACAGGCTGCTCTTCGTTAGAAAAAAAAAAATCATTTGGGGACTGCTTTTCATTAAAAGGAAAACCTTACCAAGGACTTCCTTACCCTATCTGCTTAAGTAATTTCTTTTTTAACTTTTATATCATTCTGACCTCTCGTTGCTTTTGAGAAATCTGACAGTTCTTTATGGGGAATTGGTCTTTTCTGTTTTCCATTTAAGATTGTGTGTGTGTGTGCGTGTATTCTGTAATTTTTCTGGATATGTCTAAATGTGGATTTGCTTTTATTCTGCTTTAGACTCCTTGTGCTTATTGAATCTATAGATTCATATCTTTACTCAGTTCTGGAAAATTCTCAGCTATTGTCTTTTAAAAGATGCTTCTTCCAGCTTTCCCAATTCTCTATGGAACTTCTTTAGACATCTTTGGCCTTAATTCTCTTTCCATGGCTCGTAACTCCTCCCTTAGTTAGTCATAGCGCTTGCCCCTCTTTTAAATCACTATTTGTTACTTAGGATTTTCCCTCCTTTCTTAGAATTTCTTTTGTTATATTTTTAAATATACAAATTTATATTTTATATACAAATATATGTTTTAAAATATATACAAAATATATACAAATATATATTTTAAAATATATACAAATATATATTTTAAAATATATACAAATATATATTTTAAAATATATACAAAATATATACAAATATATATTTTAAATATACACAAATATATATATTTTAAATATACAAATATTTTTGTTATATTTTATCTAATACATCTACATATTTATATGAGGACTTTTTTGGCATTATCTTTACCATTTCTGCTGGAACTGATTCTTGACAGCTGTAGTTTTCTTGTAGAGTTGTTTTTTGTTGTTTTTTGATATATTCCAAGAGTTTTTCAAATTTTATCTATAATTAATATTTGGATATTGTAGTTCAAAAGAGCAGGAGAAGACAGTGGTTAATATTTTGGACAGCCCAATACTATAGATAAAATGAAGTCACCTTTTCCTATAAAATGGATTTTTTTTATTGTTTTTCTAGAGCTTTGACTATTGGGAATTTTGCTTAGGATCTCTAGTTTTAAAAGTTTCTGAATCTTCCTGAGCCTCCTTCAGTATTTTTGCATGTAATAAACCCCAAAGTTCTCTCATGGGGATTCTTTTATAAATGGAGAGAATGGCATAAGCAGTTAAGAGCAAGGTCTCTGGAATCAAACAGCTTGTAAGTGTCCAGTCTTGTCTTTACCACTTATTAACATTGTGATTTTGGGTGAGTTCTCTCACCTCTCTTTCTCAGTCTCCCGCTCTTTAAAGTAGCAACCTCGTAGAGTTGGGATGATTAACCAAGATCTAATATGTATTAAGTGGCATCTAAGGAGCTCTCAAAAATATTAGATATTATTAATATAATTATTTTCATTGATTCACAGTAACTTTTTCCTTGCACTTAAAGCTATTTCAGGTGTAGCAGTTTTAATCCGCCACCTGATACAGCATTCTTGAAGGCAAAAGCACTTGGACTAATAATAAGCATAGGCCGTCACTTCCCCCGCATTATACATTTTCTCTACTGAATCTTTTCCATTAGCATACCAACTCATATTTCTTTAAAGCTAAAATAACTCTTTTGTCATCACTTCCCCCTCCAATTTTCCCTTTACAATAAAACTCCTTTAATATGTCCATTTTTACACTCCCTAATCTGTCTTTCATGAAACCAAACTCCAGTCAAGGTTTTACCTCATCCATGCCACCAAATCTGCTTGGATTAATAAACCCAGTGGTAAATCTCAGTTCCTGTGTTACTTGACCCATCAGCAGCATTCAACACAGTTGATCCCTCCCTTTTCCTTGAAACCTGGCTTTCTGCTCCCTCACTGGCTGCTGCTCCTTGGTCTCCTTTGTTGGTTCATCCTCACCTTGAAATGTTGGCATGCTCCTGGAATCCGCCCACGGACCTCTTCTTTATTTGTCTACACTCCTTCTCAGTCTCATGGCTTTGACTACGATCTGTATGTTCACCACTCCCAAATTCATATCTCCAGCCTAGGACTCTTGCTTGACTCCAGATGCATATAGCCATTTGCCTACTCAGTATCTCCATTTGGATGTCTATTAGACATCTCAAATTTCACATGTTCAAAACAACTTCTGAGAATTATGTTAGCTCTGGGTTTGTCATACATGGCCTTTCTACTGAGGTATTATTTCCTCTATACCAATCTGTTGAGAGTTTTTATCATGAAAGGGTATCAAAGGCTTTTTCTGCATTTATTTAAATGATCCTCTGATTTTTAGCCTTCATTTTGTTAATGTGGTGTATCATATTTATTGATTGTGTATGCAATAAATCCCATTTGATTGTGGTGATTCTTTTATCAAAATCTTAATTACCCCCTTTTAAAACATTTGTAACATGATAAACTATCTTACTAAACATTTCATAATAAACAGTTGCTTCTTCCAATTTAAGAAAATGTATATTTTAGTTGTTGCTTGGGCTGAGGATGTCCCTAATCAATCTAACTATATTTGTGTCTTCAGGCATTATTTGAGGTGAATAAAATGCAAATAATTTTTTAGACATAAAGGCCTCTGTTAGCCAAATGAGAAATCTTTCTTTCACCTCTAGTTCTCATTCACTCCCCCTTTTTCTTTCACAAAGGCATTGATTAATTTATATCTTAGGTTAGCTAAGAATTTGAGCTATTTTTAGCAGTACATTCTCAAAGTATAGATGCTGTGGAATGAGTAGCAAGTAATTAGAGCACCGTTTGTATAACCTTTATGTACAGTAGCCTTGAAAAATCTCCGTACTGTAGGTAGAATTGCACCTAAGATGACAAATCTTGTAGTAGGACCACCTAGCTAATGATTGAGTCTCCTTTATAGCCTGGGAACTGATGAGAAGGACACTCAGCACTCCTGAGTTCCTCCTGTCTCTACGTCTTTCTCATATTAAGCCTAAGTTGGCTTTCACACAATTAATTACAATTTGCTTCTATAGCTTTGGGGCTTTACGTAACAACTTTAATTTTTTCTACACCCTTTTGAAATAGCTATCATGTTTGACCCGGGTCTTTTCTCCTTATAAAACATTAGTATTGTCTCCTGTTGTCTGATATAGTTTTGAGTCCCTTCAACATCAACTGTAGCTCTGTCAGTGTCCCTCTTTAAATACAGTGCTCAGAGTTAGCCTTATAATGTTTCATGGTTTTGTTTTTGTTTTTAGATACAAGGTCTTGTTCTGTCACCCAGGCTGGAGTGCAGTGGCATGATGGTAGCTCACTGCAGCCTTGGACTCCTGGGCTCAAGGGATCCTCCTGCCTCGGTCTTCCAAATCACTAGGATTATAGGGATGCATCACCGCGACCAGCCTCGTGTATTATTTAATTAGTAGTCATTCAACCTATTTTTTTCTAGATAGTACATTTTGATTAGTACAGTTTATGATTTTTTTTTTTTTTACCTGTATCTTACTTATAGTACTTACAGGAGCTGAGGAAAACCCCAAGACCAACCTGTGCTCAACCTGTGCTTTTGTAGTGATTTTTTTGGACCCTATTACAGAACCTTATTTATTCCTTCATCTTGTTAAATGTGGCCCATTGTTCCAGGCTACCACCACATTTCTTGATTCTTTCTTTGAACATAATTATTATATTTCCCAGCTTCGTGTCATTTCTAAAATTGATGTACCTGCCTTCATCTAGGTCATTGCTAAAAACAATGAAGTAGAGCTAAAGACAGCCTTTCAGCATTCACAAAACTTATCTATTCTATATGCCATATATTTATATAGCACCTTCTATTTGCACTGTCCTAGTCATGAGGGAATCATAAATAAAAAGAAATGGCCCTTTAAGAGCTCACAGCCTAGTGGGAAAACAGATACACAGTATATAATTTAGCATGATAATTGCTGTAATGGAATATGAAGAAGATAGAATGGGAGCACAGAGGAAAAACTCCGAGTCTTTGAGGCAGAAACTTCCTTGCAGGTTGAGAGCTGTCCTTTAATCAAAAAAGTGTAAATTAGCCAGGCACGGTGGCTCACCTAATCCTAGCACTTTGGGAGGCCGGGGCAGATGGATCACCTGAGGTCAGGAGTTTGAGACCAGCCTGGCCAACATGGTGAAACCCCATCTCTCTGAAAATACAAAAATTAGCCGGGCATAGTGGCGCATGCCTGTAATCCCAGCTACTCGGGAGACTGAGGCAGGAGGATCACTTGAACCCGGGAGGTGGAGGTTGCAGTGAGCCGAGATGGCGCCACTGCACTCCAGCCTGGGTGACAGAGCGAGACTCCCGTCTCAAACAAGTGTAAATCCACCTAACCTGTCTCTTTTTTCCCAACTTGTTCACAAGGATACAGTGAAAGACCCTGGCATCTGCTTGCTGAAATTGGGCCACATTCTATTTACTACATTTCTCTATTCCACCAGCCTAATACCCTGTCCCCTCCTCTGCTCCCAAAAAAACAAAACAAAACAAAACAAAACAAAAAACTGTTAGTCTTTTAGATGTAATTAATTAAGCTGTTTATATTTTAAATGTATTAATTTTTACATATACATACATATATCTTCAGGAGCAGGATCGGGGAAATGAGTCATAGAATATAATAGACTCTGAAGATCAGAAGGAAAAAAAATAGAGGCAAATTTGCAGGCTATGAGATCCCCCAGTTGTTATAGTTGAAGATCCAAATTGGGCTCTGAGCTTTCTAAAGGCCAAAGAGAAATGGAAAAGTTCATGTTCTTACTGACTGACAGAGAAATAGTTCCTGTTCTTAAGTGGGAACAAAGCTTTTGCCAACACCGAAATCTTACATTTCAGTGGGGGCTGCATTGAGGTTACACCGAGAGAGAGAATAAAAAATGTCATGAGCATATATTTCTACTGCAAGTGCACAAAAGATTTCCTTAGTGTTTCCTAAATTTTGTTATTGTGGAGGAAATTGCTCAAAGAACAGATTTAATATTACAGGCCAAGAGAGGCTGGGCAGCATGGCTCATGCCTGTAATCCCAGCACTTTGGAAGGCCAAGGCAGGCTAGATAGCTTGAGCCCGGGAGGTGGAGGTTACAGTGAGCTGAGATCACGCCACTGCACTCCAGCCTGGGTGACAGAGTGAGACTCTGTCTCAAAAAAAAAGTGGGGCGGGGTGGGGGACTAAAAGAACTGTTATAGAAAGGAAAGCTGCGCAACTTTGAATAACAAATTCATTGGGGATGTGAAATTTTTTGTGCACACACCTAGACCATCCCTCAGCACAAGTTCAAGAAGTCTTATATTAGGCTGTTTATTATAGCATCCTTCATAAGAGAGGAACACCAATTCCATGGTAGGTTGACAAGTATTTCGTAGATACCCATGATTCTGTGGTCAAATAAATATAGGAAACGCTGAGTTAAACAAGTTTCTTTACTGCTGGTGGGCTTCTGTCTCATAGGCTTAAGTACAATGTAACTCTCCCAGAGGGAGTTCTCCAAAACTTTATTTGGCCATAGAACCCTTTTTGGTGGTACATTAGGGAACAGAGGAATAAATGTTTACATATCCTCCAAACAACCTTTCTAAAAGAGAAAGAATGTGGCAGCCACCCTGTTACATATTTGTCAACTCCAGCCAAGGATATATGTTGGACAATACTGTTCTTCCCCTCTCTTACTCTCTCAAGCTTCAGGATGATGTAGTCTCTTTTTTGCTAGGGTCTTCTCGTTTCCTTATAAGGCAAACCAGAGTACAAGTTCCCCTCTTCACTTTCTCTTTTTTTTTTTTTTTGAGATGGAGTCACTCTGTCGCCTAGGCTGGAGTGCAGTGGCGCAATCTTGGCTCACTGCAACCTCCGCCTCCCGGGTTCAAGCGATTCTCCTGCCTCAGCCTCTCTGAGTAGCTGCGATTACAGGTGCACACAACCACACCCGGCTAATTTTTGTATTTTTAGTAGAAACGGGGTTTCACCACATTGGTCAGGGTGGTCTTGAACTCCTGACCTCGTGATCCACCCGTCTCCGCCTCCCAAAGTGCTAGGATTACAGGCTTGAGCCACCGCGCTCGGCCCACTTTCTCTGCTTTTTTAAAGATGAAACTATGATCAGGGCAAGCCAGGAGCATGTTGTCACTTGGCTTTTTTGCCAAGTAAGACTTCTAGTTAGAGACCCCTGTCACTGCAATTAACTTTGCTGAATTTGTAACTTGAGATCAGGGTTCACTTCCTTATCTGGCCTCCATCCCATTCACAATGTCATTTTGTATAGCAAATCTAATGTGTTTTTTTTGTTTGTTTGTTTGTTTTTTAGACGGAGTTTCACTCTGTCGCCCAGGCTACAGTGAAGTGGCGCAATCTTGGCCCACTGCAACCTCCGCCCCCTGGGTTCAAGCGATTCTCCTGCCTCAGCCTCCCTAGTAGCTGGGACTACAGGCGCGTGCCACCATGCCTGGCAGATTTTGTATTTTTAGTAGAGATGGGGTTTCGCCATGTTGGCCAGGCTGATCTCCAACTTCTGACCTCAGGTGATCCACCCGCCTAGGCCTCCTGAAGTGCTGGGATTACAGGCGTGAGCCACCGCACCCAGCCAATCTAATGCGTTTCTTTTAGAATAAAGTAGGTTTTTTTCCCATGTATTATAAACAAGCAATAAATCTGGGGAGGTTATTTCTACTTGATTAGTTTTTATTTACCTAATAAGAATCCTAAAATGCAAAAATAGTATTGTATTAATATAACAGTGACTTTAAAATAAACTTCAGCAATATTACATTTGAAAAATACCCTAAGGTTAAATTTATTTATAACAAGTGTTATTCATTTATTATAATACCCCAGAGTCACATTCAGTATTGTCTTTCGATATATGCATTTTATAATTGTAAATCCATATACAGATGTCCCCTACTTAATGAGGGTTCGACTTACTGTTTTTCAACATCACGATGGTATGAAAGCAACACACATTCAGTAGAAACAGTACTTCGATTTTTGAATTTTGATCTTTTCCTGGGCCAGTGATAAGTGGTATGGGGGCAGCAGCAGCGAGCCACAGCTTGCAGTCAGCCACACCATCATGAGGGTAAACACCGATACTCTCCTGTGGACTGTGATGCAGATGATTCTGCCCGACTGTAGACCAATATACATGCTCTGAGCATGTTTCAGATGGACTAGGCTAAGCTATATGTTCAGTAAGTTAGGTGTATTCACTGCATTTTTGACTTAATGATATTTTCAATTTACTGCGAATTTATTGGACATAACTCCATTGTAAGTAGAGGAACATTTCTGTCTAATATTTTGTGTTTTCCTCCCATTAAATTCCTTATTCTCTTATGTTATGGAATTCCTTATGGAATACATATTCCCTTACGTTGACATTGACAACATACTGGTAATCTTTAATATATATGAAATACTTGAATATTGTCAAATGTTATTTTTGAAGTCATTTTGATATTTGGCAATTGAGGTTTTTTCAACTTTTTGCAATTTTAAAAGGGGTATTGATATCATATGAATTATCTTTTTTCTCTGGTTTATTTCCTTGGGCTGTATTCCTATGGATAGGACTACAAAGTCAGCAATTATATTTTCTAGCTTTTGTATGACAGGAATTTTCTAACCTTATTGTCTTCTAAATCCTAGGTTCAGGAGCTCCATGAGGAGGAGGCACAGTGGGTGAACTATGATGAGGATGAGTTGTGTGTGAAAATGCAGCTAGCCGACGGGATCTTTGAGACCCTGATCAAAGATACTATTGATGTTCTGAATCAGATCAGTGAAAAGCAGGGGAGAATGCTACTTGTGTGA';

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get3PrimeUtrClass,$a->getCDSClass,$a->getExonClass,$a->getEssentialSpliceSiteClass,$a->getSpliceRegionClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									6927,-100,13491,0,$expectedRnaWtSeq,'G','r.6927-100_13491del22367insg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get500BPDownStreamVariantClass,$a->get5KBDownStreamVariantClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									6509,-100,9354,0,$expectedCdsWtSeq,'G','c.6509-100_9354del18648insG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','p.?',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getUnknownVariantClass,$a->getStopLostVariantClass);
		done_testing();
	};

}



# TOR1AIP2 protein coding gene with 5 prime utr exons on - strand of genome

sub testUpsteamMilesAway_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Upstream Miles Away 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179856934,
			'maxpos'				=> 179856938,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');
		done_testing();
	};
}
sub testEndsUpsteam5001bp_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Upstream Miles Away 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179851935,
			'maxpos'				=> 179851939,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testEndsUpsteam5000bp_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Ends Upstream 5000 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179851934,
			'maxpos'				=> 179851938,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);


		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get5KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam2001bp_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Ends Upstream 2001 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179848935,
			'maxpos'				=> 179848939,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);


		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get5KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam2000bp_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Ends Upstream 2000 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179848934,
			'maxpos'				=> 179848938,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);


		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get2KBUpStreamVariantClass,$a->get5KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam1997bp_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Ends Upstream 1997 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179848931,
			'maxpos'				=> 179848935,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);


		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get2KBUpStreamVariantClass,$a->get5KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam1996bp_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Ends Upstream 1996 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179848930,
			'maxpos'				=> 179848934,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);


		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get2KBUpStreamVariantClass);

		done_testing();
	};
}
sub testEndsUpsteam1bp_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Ends Upstream 1 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846935,
			'maxpos'				=> 179846939,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);


		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get2KBUpStreamVariantClass);

		done_testing();
	};
}

sub test5PrimeUTRIntronic_1bp_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Intronic 1bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846817,
			'maxpos'				=> 179846817,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic_1bp_2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Intronic 1bp - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846650,
			'maxpos'				=> 179846650,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic_1bp_3_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Intronic 1bp - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846519,
			'maxpos'				=> 179846519,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic_10bp_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Intronic 10bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846808,
			'maxpos'				=> 179846817,
			'delseq' 				=> 'ATATATATAT',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic_10bp_2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Intronic 10bp - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846519,
			'maxpos'				=> 179846528,
			'delseq' 				=> 'ATATATATAT',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}

sub testCDSIntronic_1bp_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 CDS Intronic 1bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821756,
			'maxpos'				=> 179821756,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testCDSIntronic_10bp_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 CDS Intronic 10bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821747,
			'maxpos'				=> 179821756,
			'delseq' 				=> 'ATATATATAT',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testCDSIntronic_1bp_2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 CDS Intronic 1bp - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820946,
			'maxpos'				=> 179820946,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testCDSIntronic_10bp_2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 CDS Intronic 10bp - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820946,
			'maxpos'				=> 179820955,
			'delseq' 				=> 'ATATATATAT',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testCDSIntronic_1bp_3_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 CDS Intronic 1bp - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820509,
			'maxpos'				=> 179820509,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testCDSIntronic_10bp_3_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 CDS Intronic 10bp - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820509,
			'maxpos'				=> 179820518,
			'delseq' 				=> 'ATATATATAT',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}

sub testStartsDownsteam1bp_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream 1 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179809097,
			'maxpos'				=> 179809101,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);


		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get500BPDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam496bp_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream 496 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179808602,
			'maxpos'				=> 179808606,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);


		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get500BPDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam497bp_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream 497 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179808601,
			'maxpos'				=> 179808605,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);


		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get500BPDownStreamVariantClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam498bp_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream 498 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179808600,
			'maxpos'				=> 179808604,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);


		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get500BPDownStreamVariantClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam499bp_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream 499 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179808599,
			'maxpos'				=> 179808603,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);


		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get500BPDownStreamVariantClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam500bp_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream 500 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179808598,
			'maxpos'				=> 179808602,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);


		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get500BPDownStreamVariantClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam501bp_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream 501 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179808599,
			'maxpos'				=> 179808603,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);


		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get5KBDownStreamVariantClass,$a->get500BPDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam5000bp_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream 5000 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179804098,
			'maxpos'				=> 179804102,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);


		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->get5KBDownStreamVariantClass);

		done_testing();
	};
}
sub testStartsDownsteam5001bp_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream 5001 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179804097,
			'maxpos'				=> 179804101,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);


		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}
sub testStartsDownsteamMilesAway_5bp_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Starts Downstream Miles Away 5bp long - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 169804097,
			'maxpos'				=> 169804101,
			'delseq' 				=> 'CCCCC',
			'insseq'				=> 'TTTT',);


		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		ok(!defined($res[0]),'annotation not defined');

		done_testing();
	};
}

sub test5PrimeUtrExon_1bp_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Exon 1bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846829,
			'maxpos'				=> 179846829,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									106,0,106,0,'U','AAAA','r.106delUinsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUtrExon_1bp_2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Exon 1bp - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846828,
			'maxpos'				=> 179846828,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									107,0,107,0,'U','AAAA','r.107delUinsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUtrEssentialSpliceSite_1bp_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Essential Splice Site 1bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846827,
			'maxpos'				=> 179846827,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									107,1,107,1,'G','AAAA','r.107+1delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUtrEssentialSpliceSite_1bp_3_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Essential Splice Site 1bp - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846825,
			'maxpos'				=> 179846825,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									107,3,107,3,'G','AAAA','r.107+3delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUtrEssentialSpliceSite_1bp_5_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Essential Splice Site 1bp - strand 5' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846823,
			'maxpos'				=> 179846823,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									107,5,107,5,'G','AAAA','r.107+5delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUtrSpliceRegion_1bp_6_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Splice Region 1bp - strand 6' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846822,
			'maxpos'				=> 179846822,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									107,6,107,6,'G','AAAA','r.107+6delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUtrSpliceRegion_1bp_10_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Splice Region 1bp - strand 10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846818,
			'maxpos'				=> 179846818,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									107,10,107,10,'G','AAAA','r.107+10delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUTRIntronic_1bp_11_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Intronic 1bp - strand 11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846817,
			'maxpos'				=> 179846817,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUTRIntronic_1bp_n11_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Intronic 1bp - strand -11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846519,
			'maxpos'				=> 179846519,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test5PrimeUtrSpliceRegion_1bp_n10_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Splice Region 1bp - strand -10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846518,
			'maxpos'				=> 179846518,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									108,-10,108,-10,'G','AAAA','r.108-10delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUtrSpliceRegion_1bp_n3_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Splice Region 1bp - strand -3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846511,
			'maxpos'				=> 179846511,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									108,-3,108,-3,'G','AAAA','r.108-3delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUtrEssentialSpliceSite_1bp_n2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Essential Splice Site 1bp - strand -2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846510,
			'maxpos'				=> 179846510,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									108,-2,108,-2,'G','AAAA','r.108-2delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUtrEssentialSpliceSite_1bp_n1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Essential Splice Site 1bp - strand -1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846509,
			'maxpos'				=> 179846509,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get5PrimeUtrClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									108,-1,108,-1,'G','AAAA','r.108-1delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUtrExon_1bp_3_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Exon 1bp - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846508,
			'maxpos'				=> 179846508,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									108,0,108,0,'C','AAAA','r.108delCinsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUtrExon_1bp_4_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Exon 1bp - strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846507,
			'maxpos'				=> 179846507,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									109,0,109,0,'G','AAAA','r.109delGinsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get5PrimeUtrVariantClass);

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

sub test5PrimeUtrExon_2bp_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Exon 2bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846828,
			'maxpos'				=> 179846829,
			'delseq' 				=> 'AA',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									106,0,107,0,'UU','AAAA','r.106_107delUUinsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get5PrimeUtrVariantClass);

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
sub test5PrimeUtrExon_2bp_2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 5 Prime UTR Exon 1bp - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179846507,
			'maxpos'				=> 179846508,
			'delseq' 				=> 'CG',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									108,0,109,0,'CG','AAAA','r.108_109delCGinsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get5PrimeUtrVariantClass);

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

sub testEssentialSpliceSite_1bp_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice Site 1bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821766,
			'maxpos'				=> 179821766,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									422,1,422,1,'G','AAAA','r.422+1delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									34,1,34,1,'G','AAAA','c.34+1delGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

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
sub testEssentialSpliceSite_1bp_3_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice Site 1bp - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821764,
			'maxpos'				=> 179821764,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									422,3,422,3,'G','AAAA','r.422+3delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									34,3,34,3,'G','AAAA','c.34+3delGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

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
sub testEssentialSpliceSite_1bp_5_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice Site 1bp - strand 5' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821762,
			'maxpos'				=> 179821762,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									422,5,422,5,'G','AAAA','r.422+5delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									34,5,34,5,'G','AAAA','c.34+5delGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

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
sub testSpliceRegion_1bp_6_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region 1bp - strand 6' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821761,
			'maxpos'				=> 179821761,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									422,6,422,6,'G','AAAA','r.422+6delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									34,6,34,6,'G','AAAA','c.34+6delGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

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
sub testSpliceRegion_1bp_10_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region 1bp - strand 10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821757,
			'maxpos'				=> 179821757,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									422,10,422,10,'G','AAAA','r.422+10delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									34,10,34,10,'G','AAAA','c.34+10delGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

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
sub testIntronic_1bp_11_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Intronic 1bp - strand 11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821756,
			'maxpos'				=> 179821756,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic_1bp_n11_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Intronic 1bp - strand -11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820509,
			'maxpos'				=> 179820509,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testSpliceRegion_1bp_n10_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region 1bp - strand -10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820508,
			'maxpos'				=> 179820508,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-10,423,-10,'G','AAAA','r.423-10delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-10,35,-10,'G','AAAA','c.35-10delGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

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
sub testSpliceRegion_1bp_n3_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region 1bp - strand -3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820501,
			'maxpos'				=> 179820501,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-3,423,-3,'G','AAAA','r.423-3delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-3,35,-3,'G','AAAA','c.35-3delGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

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
sub testEssentialSpliceSite_1bp_n2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice Site 1bp - strand -2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820500,
			'maxpos'				=> 179820500,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-2,423,-2,'G','AAAA','r.423-2delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-2,35,-2,'G','AAAA','c.35-2delGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

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
sub testEssentialSpliceSite_1bp_n1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice Site 1bp - strand -1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820499,
			'maxpos'				=> 179820499,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-1,423,-1,'G','AAAA','r.423-1delginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-1,35,-1,'G','AAAA','c.35-1delGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

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

sub testEssentialSpliceSite_3bp_1t3_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice Site 3bp - strand 1 to 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821764,
			'maxpos'				=> 179821766,
			'delseq' 				=> 'CCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									422,1,422,3,'GGG','AAAA','r.422+1_422+3delggginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									34,1,34,3,'GGG','AAAA','c.34+1_34+3delGGGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

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
sub testEssentialSpliceSite_3bp_3t5_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice Site 3bp - strand 3 to 5' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821762,
			'maxpos'				=> 179821764,
			'delseq' 				=> 'CCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									422,3,422,5,'GGG','AAAA','r.422+3_422+5delggginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									34,3,34,5,'GGG','AAAA','c.34+3_34+5delGGGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

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
sub testEssentialSpliceSite_3bp_4t6_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice Site 3bp - strand 4 to 6' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821761,
			'maxpos'				=> 179821763,
			'delseq' 				=> 'CCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getEssentialSpliceSiteClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									422,4,422,6,'GGG','AAAA','r.422+4_422+6delggginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									34,4,34,6,'GGG','AAAA','c.34+4_34+6delGGGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->getSpliceRegionVariantClass);

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
sub testEssentialSpliceSite_3bp_5t7_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice Site 3bp - strand 5 to 7' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821760,
			'maxpos'				=> 179821762,
			'delseq' 				=> 'CCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getEssentialSpliceSiteClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									422,5,422,7,'GGG','AAAA','r.422+5_422+7delggginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									34,5,34,7,'GGG','AAAA','c.34+5_34+7delGGGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->getSpliceRegionVariantClass);

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
sub testSpliceRegion_3bp_6t8_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region 3bp - strand 6 to 8' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821759,
			'maxpos'				=> 179821761,
			'delseq' 				=> 'CCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									422,6,422,8,'GGG','AAAA','r.422+6_422+8delggginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									34,6,34,8,'GGG','AAAA','c.34+6_34+8delGGGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

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
sub testSpliceRegion_3bp_8t10_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region 3bp - strand 6 to 8' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821757,
			'maxpos'				=> 179821759,
			'delseq' 				=> 'CCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									422,8,422,10,'GGG','AAAA','r.422+8_422+10delggginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									34,8,34,10,'GGG','AAAA','c.34+8_34+10delGGGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

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
sub testSpliceRegion_3bp_9t11_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region 3bp - strand 9 to 11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821756,
			'maxpos'				=> 179821758,
			'delseq' 				=> 'CCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									422,9,422,11,'GGG','AAAA','r.422+9_422+11delggginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getIntronVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									34,9,34,11,'GGG','AAAA','c.34+9_34+11delGGGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getIntronVariantClass);

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
sub testSpliceRegion_3bp_10t12_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region 3bp - strand 10 to 12' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821755,
			'maxpos'				=> 179821757,
			'delseq' 				=> 'CCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									422,10,422,12,'GGG','AAAA','r.422+10_422+12delggginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getIntronVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									34,10,34,12,'GGG','AAAA','c.34+10_34+12delGGGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getIntronVariantClass);

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
sub testIntronic_3bp_11t13_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Intronic 3bp - strand 11 to 13' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179821754,
			'maxpos'				=> 179821756,
			'delseq' 				=> 'CCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic_3bp_n13tn11_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Intronic 3bp - strand -13 to -11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820509,
			'maxpos'				=> 179820511,
			'delseq' 				=> 'CCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),1,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(!defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getUnknownAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffSequenceSubtype(),
									0,0,0,0,'?','?','r.?',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testSpliceRegion_3bp_n12tn10_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region 3bp - strand -12 to -10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820508,
			'maxpos'				=> 179820510,
			'delseq' 				=> 'CCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-12,423,-10,'GGG','AAAA','r.423-12_423-10delggginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getIntronVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-12,35,-10,'GGG','AAAA','c.35-12_35-10delGGGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getIntronVariantClass);

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
sub testSpliceRegion_3bp_n11tn9_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region 3bp - strand -11 to -9' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820507,
			'maxpos'				=> 179820509,
			'delseq' 				=> 'CCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-11,423,-9,'GGG','AAAA','r.423-11_423-9delggginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getIntronVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-11,35,-9,'GGG','AAAA','c.35-11_35-9delGGGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getIntronVariantClass);

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
sub testSpliceRegion_3bp_n10tn8_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region 3bp - strand -10 to -8' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820506,
			'maxpos'				=> 179820508,
			'delseq' 				=> 'CCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-10,423,-8,'GGG','AAAA','r.423-10_423-8delggginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-10,35,-8,'GGG','AAAA','c.35-10_35-8delGGGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

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
sub testSpliceRegion_3bp_n5tn3_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region 3bp - strand -5 to -3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820501,
			'maxpos'				=> 179820503,
			'delseq' 				=> 'CCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-5,423,-3,'GGG','AAAA','r.423-5_423-3delggginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-5,35,-3,'GGG','AAAA','c.35-5_35-3delGGGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

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
sub testEssentialSpliceSite_3bp_n4tn2_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice Site 3bp - strand -4 to -2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820500,
			'maxpos'				=> 179820502,
			'delseq' 				=> 'CCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-4,423,-2,'GGG','AAAA','r.423-4_423-2delggginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-4,35,-2,'GGG','AAAA','c.35-4_35-2delGGGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getEssentialSpliceSiteVariantClass);

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
sub testEssentialSpliceSite_3bp_n3tn1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Essential Splice Site 3bp - strand -3 to -1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820499,
			'maxpos'				=> 179820501,
			'delseq' 				=> 'CCC',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getSpliceRegionClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									423,-3,423,-1,'GGG','AAAA','r.423-3_423-1delggginsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getEssentialSpliceSiteVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									35,-3,35,-1,'GGG','AAAA','c.35-3_35-1delGGGinsAAAA',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->getEssentialSpliceSiteVariantClass);

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

sub testInFrameInPhasePreserveLength_3bp_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 In Frame In Phase Preserve Length - strand 3bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820494,
			'maxpos'				=> 179820496,
			'delseq' 				=> 'AGA',
			'insseq'				=> 'CCC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									425,0,427,0,'UCU','GGG','r.425_427delUCUinsggg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									37,0,39,0,'TCT','GGG','c.37_39delTCTinsGGG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									13,0,13,0,'S','G','p.S13G',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getNonSynonymousVariantClass,);
		done_testing();
	};
}
sub testInFrameOutOfPhasePreserveLength_2bp_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 In Frame Out Of Phase Preserve Length - strand 2bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820494,
			'maxpos'				=> 179820495,
			'delseq' 				=> 'AG',
			'insseq'				=> 'CC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									426,0,427,0,'CU','GG','r.426_427delCUinsgg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									38,0,39,0,'CT','GG','c.38_39delCTinsGG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getSubstitutionAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									13,0,13,0,'S','W','p.S13W',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getSubstitutionClass,$a->getNonSynonymousVariantClass,);
		done_testing();
	};
}
sub testInFrameOutOfPhasePreserveLength_4bp_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 In Frame Out Of Phase Preserve Length - strand 4bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820494,
			'maxpos'				=> 179820497,
			'delseq' 				=> 'AGAG',
			'insseq'				=> 'CCCC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									424,0,427,0,'CUCU','GGGG','r.424_427delCUCUinsgggg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									36,0,39,0,'CTCT','GGGG','c.36_39delCTCTinsGGGG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									12,0,13,0,'DS','EG','p.D12_S13delinsEG',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,);
		done_testing();
	};
}

sub testInFrameInPhaseShorten_9to3_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 In Frame In Phase Shorten - strand 9bp to 3bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820488,
			'maxpos'				=> 179820496,
			'delseq' 				=> 'CTTTTGAGA',
			'insseq'				=> 'CCC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									425,0,433,0,'UCUCAAAAG','GGG','r.425_433delUCUCAAAAGinsggg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									37,0,45,0,'TCTCAAAAG','GGG','c.37_45delTCTCAAAAGinsGGG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									13,0,15,0,'SQK','G','p.S13_K15delinsG',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,);
		done_testing();
	};
}
sub testInFrameInPhaseLengthen_3to9_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 In Frame In Phase Lengthen - strand 3bp to 9bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820494,
			'maxpos'				=> 179820496,
			'delseq' 				=> 'AGA',
			'insseq'				=> 'CCCCCCCCC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									425,0,427,0,'UCU','GGGGGGGGG','r.425_427delUCUinsggggggggg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									37,0,39,0,'TCT','GGGGGGGGG','c.37_39delTCTinsGGGGGGGGG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									13,0,13,0,'S','GGG','p.S13delinsGGG',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,);
		done_testing();
	};
}

sub testInFrameOutOfPhaseShorten_8to2_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 In Frame Out Of Phase Shorten - strand 8bp to 2bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820489,
			'maxpos'				=> 179820496,
			'delseq' 				=> 'TTTTGAGA',
			'insseq'				=> 'CC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									425,0,432,0,'UCUCAAAA','GG','r.425_432delUCUCAAAAinsgg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									37,0,44,0,'TCTCAAAA','GG','c.37_44delTCTCAAAAinsGG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									13,0,15,0,'SQK','G','p.S13_K15delinsG',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,);
		done_testing();
	};
}
sub testInFrameOutOfPhaseLengthen_2to8_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 In Frame Out Of Phase Lengthen - strand 2bp to 8bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820495,
			'maxpos'				=> 179820496,
			'delseq' 				=> 'GA',
			'insseq'				=> 'CCCCCCCC');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									425,0,426,0,'UC','GGGGGGGG','r.425_426delUCinsgggggggg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									37,0,38,0,'TC','GGGGGGGG','c.37_38delTCinsGGGGGGGG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getInFrameVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									13,0,13,0,'S','GGG','p.S13delinsGGG',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,);
		done_testing();
	};
}

sub testFrameShiftShorten_9to1_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 In Frame Out Of Phase Shorten - strand 8bp to 2bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820488,
			'maxpos'				=> 179820496,
			'delseq' 				=> 'CTTTTGAGA',
			'insseq'				=> 'G');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									425,0,433,0,'UCUCAAAAG','C','r.425_433delUCUCAAAAGinsc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									37,0,45,0,'TCTCAAAAG','C','c.37_45delTCTCAAAAGinsC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									13,0,13,0,'S','RFGK*','p.S13fs*5',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getFrameShiftVariantClass,);
		done_testing();
	};
}
sub testFrameShiftLengthen_1to9_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 In Frame Out Of Phase Shorten - strand 8bp to 2bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179820496,
			'maxpos'				=> 179820496,
			'delseq' 				=> 'A',
			'insseq'				=> 'GGGGGGGGG');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									425,0,425,0,'U','CCCCCCCCC','r.425delUinsccccccccc',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									37,0,37,0,'T','CCCCCCCCC','c.37delTinsCCCCCCCCC',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getFrameShiftVariantClass);

		AnnotationTestUtils::checkAnnotation('examine Protein annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getFrameShiftAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									13,0,13,0,'S','PPPLKRIWKMIHQ*','p.S13fs*14',$t[0]->getProteinAccession,$t[0]->getProteinAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getFrameShiftVariantClass,);
		done_testing();
	};
}

sub testIntronToExonBoundry_15to1_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Intron To Exon Boundry - strand 15bp to 1bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179816770,
			'maxpos'				=> 179816783,
			'delseq' 				=> 'CTCTGGAGAGGAAA',
			'insseq'				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass,$a->getSpliceRegionClass,$a->getEssentialSpliceSiteClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									942,-12,943,0,'UUUCCUCUCCAGAG','G','r.942-12_943deluuuccucuccagAGinsg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									554,-12,555,0,'TTTCCTCTCCAGAG','G','c.554-12_555delTTTCCTCTCCAGAGinsG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

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
sub testSpliceRegionToExonBoundry_10to1_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Splice Region To Exon Boundry - strand 10bp to 1bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179816770,
			'maxpos'				=> 179816779,
			'delseq' 				=> 'CTCTGGAGAG',
			'insseq'				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass,$a->getSpliceRegionClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									942,-8,943,0,'CUCUCCAGAG','G','r.942-8_943delcucuccagAGinsg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									554,-8,555,0,'CTCTCCAGAG','G','c.554-8_555delCTCTCCAGAGinsG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

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
sub testEssSpliceToExonBoundry_4to1_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Ess Splice To Exon Boundry - strand 4bp to 1bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179816770,
			'maxpos'				=> 179816773,
			'delseq' 				=> 'CTCT',
			'insseq'				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									942,-2,943,0,'AGAG','G','r.942-2_943delagAGinsg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									554,-2,555,0,'AGAG','G','c.554-2_555delAGAGinsG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

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

sub testExonToEssSpliceBoundry_4to1_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Exon To Ess Splice Boundry - strand 4bp to 1bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179816668,
			'maxpos'				=> 179816671,
			'delseq' 				=> 'ACCA',
			'insseq'				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1042,0,1043,2,'UGGU','G','r.1042_1043+2delUGguinsg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,0,655,2,'TGGT','G','c.654_655+2delTGGTinsG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

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
sub testExonToSpliceRegionBoundry_10to1_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Exon To Splice Region Boundry - strand 10bp to 1bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179816662,
			'maxpos'				=> 179816671,
			'delseq' 				=> 'GTTCTCACCA',
			'insseq'				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass,$a->getEssentialSpliceSiteClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1042,0,1043,8,'UGGUGAGAAC','G','r.1042_1043+8delUGgugagaacinsg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,0,655,8,'TGGTGAGAAC','G','c.654_655+8delTGGTGAGAACinsG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

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
sub testExonToIntronBoundry_16to1_1_TOR1AIP2 {
	my $file = shift;

	subtest 'Testing TOR1AIP2 Exon To Intron Boundry - strand 16bp to 1bp 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 1,
			'minpos'				=> 179816656,
			'maxpos'				=> 179816671,
			'delseq' 				=> 'AGGTTTGTTCTCACCA',
			'insseq'				=> 'C');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->getCDSClass,$a->getExonClass,$a->getEssentialSpliceSiteClass,$a->getSpliceRegionClass,$a->getIntronClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1042,0,1043,14,'UGGUGAGAACAAACCU','G','r.1042_1043+14del16insg',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

		AnnotationTestUtils::checkAnnotation('examine CDS annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									654,0,655,14,'TGGTGAGAACAAACCT','G','c.654_655+14del16insG',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getComplexChangeVariantClass);

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

# MYB protein coding gene with 3 prime utr exons on + strand of genome

sub test3PrimeUTRIntronic_1bp_96_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Intronic 1bp + strand 96' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135524473,
			'maxpos'				=> 135524473,
			'delseq' 				=> 'A',
			'insseq'                => 'TTT');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic_10bp_96_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Intronic 1bp + strand 96' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135524473,
			'maxpos'				=> 135524482,
			'delseq' 				=> 'TATATATATA',
			'insseq'                => 'TTT');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic_1bp_99_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Intronic 1bp + strand 99' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135538991,
			'maxpos'				=> 135538991,
			'delseq' 				=> 'T',
			'insseq'                => 'TTT');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic_10bp_99_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Intronic 10bp + strand 99' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135538982,
			'maxpos'				=> 135538991,
			'delseq' 				=> 'TATATATATA',
			'insseq'                => 'TTT');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}

sub test3PrimeUTRExon_1bp_1_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 Prime UTR Exon 1bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135521336,
			'maxpos'				=> 135521336,
			'delseq' 				=> 'T',
			'insseq'				=> 'AAAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									1587,0,1587,0,'U','AAAA','r.1587delUinsaaaa',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUTRExon_1bp_2_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 Prime UTR Exon 1bp + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135521337,
			'maxpos'				=> 135521337,
			'delseq' 				=> 'A',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									1588,0,1588,0,'A','UUUU','r.1588delAinsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUTREssentialSpliceSite_1bp_1_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 Prime UTR Essential Splice Site 1bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135521338,
			'maxpos'				=> 135521338,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get3PrimeUtrClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1588,1,1588,1,'G','UUUU','r.1588+1delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUTREssentialSpliceSite_1bp_5_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 Prime UTR Essential Splice Site 1bp + strand 5' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135521342,
			'maxpos'				=> 135521342,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get3PrimeUtrClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1588,5,1588,5,'G','UUUU','r.1588+5delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUTRSpliceRegion_1bp_6_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 Prime UTR Splice Region 1bp + strand 6' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135521343,
			'maxpos'				=> 135521343,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get3PrimeUtrClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1588,6,1588,6,'G','UUUU','r.1588+6delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUTRSpliceRegion_1bp_10_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 Prime UTR Splice Region 1bp + strand 10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135521347,
			'maxpos'				=> 135521347,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get3PrimeUtrClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1588,10,1588,10,'G','UUUU','r.1588+10delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUTRIntronic_1bp_11_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Intronic 1bp + strand 96' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135521348,
			'maxpos'				=> 135521348,
			'delseq' 				=> 'A',
			'insseq'                => 'TTT');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRIntronic_1bp_n11_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 prime UTR Intronic 1bp + strand -11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135521417,
			'maxpos'				=> 135521417,
			'delseq' 				=> 'A',
			'insseq'                => 'TTT');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub test3PrimeUTRSpliceRegion_1bp_n10_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 Prime UTR Splice Region 1bp + strand -10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135521418,
			'maxpos'				=> 135521418,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get3PrimeUtrClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1589,-10,1589,-10,'G','UUUU','r.1589-10delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUTRSpliceRegion_1bp_n3_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 Prime UTR Splice Region 1bp + strand -3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135521425,
			'maxpos'				=> 135521425,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get3PrimeUtrClass,$a->getSpliceRegionClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1589,-3,1589,-3,'G','UUUU','r.1589-3delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUTREssentialSpliceSite_1bp_n2_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 Prime UTR Essential Splice Site 1bp + strand -2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135521426,
			'maxpos'				=> 135521426,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get3PrimeUtrClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1589,-2,1589,-2,'G','UUUU','r.1589-2delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUTREssentialSpliceSite_1bp_n1_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 Prime UTR Essential Splice Site 1bp + strand -1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135521427,
			'maxpos'				=> 135521427,
			'delseq' 				=> 'G',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getProteinCodingType(),'annotation group type - proteincoding');
		is(scalar(@{$res[0]->getAllAnnotations}),3,'annotation count for group');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext())),'has mRNA context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getCDSAnnotationContext())),'has have CDS context annotation');
		ok(defined($res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getProteinAnnotationContext())),'has have protein context annotation');

		AnnotationTestUtils::checkAnnotationGroup('examine annotation group in detail',$res[0],
									$t[0]->getGeneName,$t[0]->getCCDS,$t[0]->getAccession,$t[0]->getGeneType,
									$a->getProteinCodingClass,$a->get3PrimeUtrClass,$a->getEssentialSpliceSiteClass);

		AnnotationTestUtils::checkAnnotation('examine mRNA annotation in detail',
									$res[0]->getAnnotationByContext(Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext()),
									Sanger::CGP::Vagrent::Data::Annotation::getmRNAAnnotationContext(),
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									1589,-1,1589,-1,'G','UUUU','r.1589-1delginsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUTRExon_1bp_3_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 Prime UTR Exon 1bp + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135521428,
			'maxpos'				=> 135521428,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									1589,0,1589,0,'C','UUUU','r.1589delCinsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get3PrimeUtrVariantClass);

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
sub test3PrimeUTRExon_1bp_4_MYB {
	my $file = shift;

	subtest 'Testing MYB 3 Prime UTR Exon 1bp + strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 6,
			'minpos'				=> 135521429,
			'maxpos'				=> 135521429,
			'delseq' 				=> 'C',
			'insseq'				=> 'TTTT',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									1590,0,1590,0,'C','UUUU','r.1590delCinsuuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->get3PrimeUtrVariantClass);

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

sub testIntronic_1bp_1_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Intronic 1bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573310,
			'maxpos'				=> 91573310,
			'delseq' 				=> 'C',
			'insseq'                => 'TTT');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic_1bp_2_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Intronic 1bp + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573342,
			'maxpos'				=> 91573342,
			'delseq' 				=> 'C',
			'insseq'                => 'TTT');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic_1bp_3_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Intronic 1bp + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573482,
			'maxpos'				=> 91573482,
			'delseq' 				=> 'G',
			'insseq'                => 'TTT');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testIntronic_1bp_4_AC068831 {
	my $file = shift;

	subtest 'Testing AC068831.2 Intronic 1bp + strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 15,
			'minpos'				=> 91573489,
			'maxpos'				=> 91573489,
			'delseq' 				=> 'T',
			'insseq'                => 'AAA');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}

# AC011503.1 lincRNA gene on - strand of genome

sub testIntronic_1bp_n11_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Intronic 1bp - strand -11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24135249,
			'maxpos'				=> 24135249,
			'delseq' 				=> 'T',
			'insseq'                => 'AAA');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
sub testSpliceRegion_1bp_n10_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Splice Region 1bp - strand -10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24135248,
			'maxpos'				=> 24135248,
			'delseq' 				=> 'T',
			'insseq'				=> 'AAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - proteincoding');
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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									201,-10,201,-10,'A','UUU','r.201-10delainsuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		done_testing();
	};
}
sub testSpliceRegion_1bp_n3_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Splice Region 1bp - strand -3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24135241,
			'maxpos'				=> 24135241,
			'delseq' 				=> 'T',
			'insseq'				=> 'AAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - proteincoding');
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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									201,-3,201,-3,'A','UUU','r.201-3delainsuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		done_testing();
	};
}
sub testEssentialSpliceSite_1bp_n2_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Essential Splice Site 1bp - strand -2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24135240,
			'maxpos'				=> 24135240,
			'delseq' 				=> 'T',
			'insseq'				=> 'AAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - proteincoding');
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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									201,-2,201,-2,'A','UUU','r.201-2delainsuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		done_testing();
	};
}
sub testEssentialSpliceSite_1bp_n1_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Essential Splice Site 1bp - strand -1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24135239,
			'maxpos'				=> 24135239,
			'delseq' 				=> 'T',
			'insseq'				=> 'AAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - proteincoding');
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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									201,-1,201,-1,'A','UUU','r.201-1delainsuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		done_testing();
	};
}
sub testExon_1bp_1_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Exon 1bp - strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24135238,
			'maxpos'				=> 24135238,
			'delseq' 				=> 'C',
			'insseq'				=> 'AAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - proteincoding');
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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									201,0,201,0,'G','UUU','r.201delGinsuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getNonCodingTranscriptVariantClass);

		done_testing();
	};
}
sub testExon_1bp_2_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Exon 1bp - strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24135237,
			'maxpos'				=> 24135237,
			'delseq' 				=> 'T',
			'insseq'				=> 'AAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - proteincoding');
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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									202,0,202,0,'A','UUU','r.202delAinsuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getNonCodingTranscriptVariantClass);

		done_testing();
	};
}
sub testExon_1bp_3_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Exon 1bp - strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134973,
			'maxpos'				=> 24134973,
			'delseq' 				=> 'T',
			'insseq'				=> 'AAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - proteincoding');
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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									466,0,466,0,'A','UUU','r.466delAinsuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getNonCodingTranscriptVariantClass);

		done_testing();
	};
}
sub testExon_1bp_4_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Exon 1bp - strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134972,
			'maxpos'				=> 24134972,
			'delseq' 				=> 'C',
			'insseq'				=> 'AAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - proteincoding');
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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionKnownSubtype(),
									467,0,467,0,'G','UUU','r.467delGinsuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getNonCodingTranscriptVariantClass);

		done_testing();
	};
}
sub testEssentialSpliceSite_1bp_1_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Essential Splice Site 1bp + strand 1' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134971,
			'maxpos'				=> 24134971,
			'delseq' 				=> 'T',
			'insseq'				=> 'AAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - proteincoding');
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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									467,1,467,1,'A','UUU','r.467+1delainsuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		done_testing();
	};
}
sub testEssentialSpliceSite_1bp_2_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Essential Splice Site 1bp + strand 2' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134970,
			'maxpos'				=> 24134970,
			'delseq' 				=> 'T',
			'insseq'				=> 'AAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - proteincoding');
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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									467,2,467,2,'A','UUU','r.467+2delainsuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		done_testing();
	};
}
sub testEssentialSpliceSite_1bp_3_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Essential Splice Site 1bp + strand 3' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134969,
			'maxpos'				=> 24134969,
			'delseq' 				=> 'T',
			'insseq'				=> 'AAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - proteincoding');
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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									467,3,467,3,'A','UUU','r.467+3delainsuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		done_testing();
	};
}
sub testEssentialSpliceSite_1bp_4_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Essential Splice Site 1bp + strand 4' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134968,
			'maxpos'				=> 24134968,
			'delseq' 				=> 'T',
			'insseq'				=> 'AAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - proteincoding');
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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									467,4,467,4,'A','UUU','r.467+4delainsuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		done_testing();
	};
}
sub testEssentialSpliceSite_1bp_5_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Essential Splice Site 1bp + strand 5' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134967,
			'maxpos'				=> 24134967,
			'delseq' 				=> 'T',
			'insseq'				=> 'AAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - proteincoding');
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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									467,5,467,5,'A','UUU','r.467+5delainsuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getEssentialSpliceSiteVariantClass);

		done_testing();
	};
}
sub testSpliceRegion_1bp_6_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Splice Region 1bp + strand 6' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134966,
			'maxpos'				=> 24134966,
			'delseq' 				=> 'T',
			'insseq'				=> 'AAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - proteincoding');
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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									467,6,467,6,'A','UUU','r.467+6delainsuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		done_testing();
	};
}
sub testSpliceRegion_1bp_10_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Splice Region 1bp + strand 10' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134962,
			'maxpos'				=> 24134962,
			'delseq' 				=> 'T',
			'insseq'				=> 'AAA',);

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

		my @res = $a->getAnnotation($sub);

		is(scalar(@res),1,'annotation group count');
		is($res[0]->getType,Sanger::CGP::Vagrent::Data::Transcript::getLincRnaType(),'annotation group type - proteincoding');
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
									Sanger::CGP::Vagrent::Data::Annotation::getComplexAnnotationType(),
									Sanger::CGP::Vagrent::Data::Annotation::getPositionOffsetSubtype(),
									467,10,467,10,'A','UUU','r.467+10delainsuuu',$t[0]->getAccession,$t[0]->getAccessionVersion,$t[0]->getDatabase,$t[0]->getDatabaseVersion,
									$a->getComplexIndelClass,$a->getSpliceRegionVariantClass);

		done_testing();
	};
}
sub testIntronic_1bp_11_AC011503 {
	my $file = shift;

	subtest 'Testing AC011503.1 Intronic 1bp - strand 11' => sub {
		 my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => AnnotationTestUtils::TRANSCRIPT_CACHE);   


		my $sub = Sanger::CGP::Vagrent::Data::ComplexIndel->new(
			'species'				=> 'human',
			'genomeVersion' => 'GRCh37',
			'chr' 					=> 19,
			'minpos'				=> 24134961,
			'maxpos'				=> 24134961,
			'delseq' 				=> 'T',
			'insseq'                => 'AAA');

		my @t = $ts->getTranscripts($sub);

		my $a = Sanger::CGP::Vagrent::Annotators::ComplexIndelAnnotator->new(transcriptSource => $ts);

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
									$a->getComplexIndelClass,$a->getIntronVariantClass);
		done_testing();
	};
}
