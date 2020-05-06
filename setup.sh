#!/bin/bash

##########LICENCE##########
# Copyright (c) 2014-2020 Genome Research Ltd.
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

# ALL tool versions used by opt-build.sh
# ensure updated in Dockerfile too
export VER_BEDTOOLS="2.28.0"
export VER_VCFTOOLS="0.1.16"
export VER_BIODBHTS="2.10"
export VER_HTSLIB="1.9"
export VER_SAMTOOLS="1.9"

if [[ ($# -ne 1) ]] ; then
  echo "Please provide an installation path dependencies expected in PATH/PERL5LIB, e.g."
  echo "  ./setup.sh /opt/myBundle"
  exit 0
fi

INST_PATH=$1

# get current directory
INIT_DIR=`pwd`

set -e

# build tools from other repos
bash build/opt-build.sh $INST_PATH
bash build/opt-build-local.sh $INST_PATH

echo
echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of PERL5LIB:"
echo "  $PERLROOT"
echo
