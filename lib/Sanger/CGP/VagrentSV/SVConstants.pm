package Sanger::CGP::VagrentSV::SVConstants;

use strict;

use warnings FATAL => 'all';
use autodie qw(:all);
use Const::Fast qw(const);

use Sanger::CGP::Vagrent; # exports VERSION
use FindBin qw($Bin);

####
#   Base constants
####
const our $PADDING_BUFFER_TR => 0;
const our $PADDING_BUFFER_SEQ => 1000;
const our $PROMOTER_SEARCH_PADDING => 500000;
const our $ENHANCER_SEARCH_PADDING => 500000;


1;
