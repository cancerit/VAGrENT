=head1 NAME

AnnotationTestUtils - utility class for testing

=head1 DESCRIPTION

Handy utility class for testing

=cut

package AnnotationTestUtils;

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

use Test::More;

1;

use constant PTEN_TRANSCRIPT_DATA_FILE => 'testData/PTEN.transcript.datadumper';
use constant BRAF_TRANSCRIPT_DATA_FILE => 'testData/BRAF.transcript.datadumper';

use constant MYB_TRANSCRIPT => 'testData/MYB.transcript.datadumper';
use constant FAM54A_TRANSCRIPT => 'testData/FAM54A.transcript.datadumper';

use constant AC068831_TRANSCRIPT => 'testData/AC068831.2.transcript.datadumper';
use constant AC011503_TRANSCRIPT => 'testData/AC011503.1.transcript.datadumper';

use constant CEP350_TRANSCRIPT => 'testData/CEP350.transcript.datadumper';
use constant TOR1AIP2_TRANSCRIPT => 'testData/TOR1AIP2.transcript.datadumper';

use constant TRANSCRIPT_CACHE => 'testData/test_transcript.cache.gz';

sub getReferenceTranscripts {
	my $ref = do shift;
	return $ref;
}

sub checkAnnotation {
	my ($name,$anno,$context,$type,$subtype,$min,$minOff,$max,$maxOff,$wt,$mt,$desc,$acc,$accVers,$db,$dbVers,@ontologies) = @_;
	subtest $name => sub {
		is($anno->getContext,$context,"check context: $context");
		is($anno->getType,$type,"check type: $type");
		is($anno->getSubtype,$subtype,"check subtype: $subtype");
		is($anno->getMinPos,$min,"check minpos: $min");
		is($anno->getMinOffset,$minOff,"check min offset: $minOff");
		is($anno->getMaxPos,$max,"check maxpos: $max");
		is($anno->getMaxOffset,$maxOff,"check max offset: $maxOff");
		if(length($wt) > 100){
			is($anno->getWt,$wt,"check wt seq length: ".length($wt));
		} else {
			is($anno->getWt,$wt,"check wt seq: $wt");
		}
		if(length($mt) > 100){
			is($anno->getMt,$mt,"check mt seq length: ".length($mt));
		} else {
			is($anno->getMt,$mt,"check mt seq: $mt");
		}

		is($anno->getDescription,$desc,"check descripton: $desc");
		is($anno->getAccession,$acc,"check accession: $acc");
		is($anno->getSequenceVersion,$accVers,"check accession version: $accVers");
		is($anno->getDatabase,$db,"check database: $db");
		is($anno->getDatabaseVersion,$dbVers,"check database version: $dbVers");
		ok(scalar($anno->getClassifications) == scalar(@ontologies),"check ontology list lengths, tested (".scalar($anno->getClassifications).") against reference (".scalar(@ontologies).")");
		my @refOnt = sort {$a cmp $b} @ontologies;
		my @testOnt = sort {$a cmp $b} $anno->getClassifications;
		for(my $i = 0; $i < scalar(@ontologies) ; $i++){
			is($testOnt[$i],$refOnt[$i],"check ontology: $i $refOnt[$i]");
		}
		done_testing();
	};
}

sub checkAnnotationGroup {
	my ($name,$group,$label,$ccds,$acc,$type,@ontologies) = @_;
	subtest $name => sub {
		is($group->getLabel,$label,"check label: $label");
		is($group->getCCDS,$ccds,"check ccds: $ccds");
		is($group->getAccession,$acc,"check accession: $acc");
		is($group->getType,$type,"check type: $type");
		ok(scalar($group->getClassifications) == scalar(@ontologies),"check ontology list lengths, tested (".scalar($group->getClassifications).") against reference (".scalar(@ontologies).")");
		my @refOnt = sort {$a cmp $b} @ontologies;
		my @testOnt = sort {$a cmp $b} $group->getClassifications;
		for(my $i = 0; $i < scalar(@ontologies) ; $i++){
			is($testOnt[$i],$refOnt[$i],"check ontology: $i $refOnt[$i]");
		}
		done_testing();
	};

}


