#!/bin/bash

##########LICENCE##########
# Copyright (c) 2014-2017 Genome Research Ltd.
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

SOURCE_SAMTOOLS="https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2"
SOURCE_HTSLIB="https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2"
SOURCE_BIOBDHTS="https://github.com/Ensembl/Bio-HTS/archive/2.3.tar.gz"
SOURCE_VCFTOOLS="https://github.com/vcftools/vcftools/releases/download/v0.1.14/vcftools-0.1.14.tar.gz"
# Warning bedtools 2.24.0 and 2.25.0 have a swapped usage in coverageBed
# No upgrades until [this ticket](https://github.com/arq5x/bedtools2/issues/319) is resolved
SOURCE_BEDTOOLS="https://github.com/arq5x/bedtools2/releases/download/v2.21.0/bedtools-2.21.0.tar.gz"

get_distro () {
  EXT=""
  DECOMP=""
  if [[ $2 == *.tar.bz2* ]] ; then
    EXT="tar.bz2"
    DECOMP="-j"
  elif [[ $2 == *.tar.gz* ]] ; then
    EXT="tar.gz"
    DECOMP="-z"
  else
    echo "I don't understand the file type for $1"
    exit 1
  fi

  if hash curl 2>/dev/null; then
    curl -sS -o $1.$EXT -L $2
  else
    wget -nv -O $1.$EXT $2
  fi
  mkdir -p $1
  tar --strip-components 1 -C $1 $DECOMP -xf $1.$EXT
}

get_file () {
# output, source
  if hash curl 2>/dev/null; then
    curl -sS -o $1 -L $2
  else
    wget -nv -O $1 $2
  fi
}

if [[ ($# -ne 1 && $# -ne 2) ]] ; then
  echo "Please provide an installation path and optionally perl lib paths to allow, e.g."
  echo "  ./setup.sh /opt/myBundle"
  echo "OR all elements versioned:"
  echo "  ./setup.sh /opt/myBundle /opt/cgpVcf-X.X.X/lib/perl5:/opt/PCAP-X.X.X/lib/perl5"
  exit 0
fi

INST_PATH=$1

if [[ $# -eq 2 ]] ; then
  CGP_PERLLIBS=$2
fi

CPU=`grep -c ^processor /proc/cpuinfo`
if [[ $? -eq 0 ]]; then
  if [[ $CPU -gt 6 ]]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

INST_PATH=$1

# get current directory
INIT_DIR=`pwd`

set -e

# cleanup inst_path
mkdir -p $INST_PATH/bin
cd $INST_PATH
INST_PATH=`pwd`
cd $INIT_DIR

# make sure that build is self contained
unset PERL5LIB
PERLROOT=$INST_PATH/lib/perl5

# allows user to knowingly specify other PERL5LIB areas.
if [ -z ${CGP_PERLLIBS+x} ]; then
  export PERL5LIB="$PERLROOT"
else
  export PERL5LIB="$PERLROOT:$CGP_PERLLIBS"
fi

export PATH=$INST_PATH/bin:$PATH

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

cd $SETUP_DIR

## grab cpanm and stick in workspace, then do a self upgrade into bin:
get_file $SETUP_DIR/cpanm https://cpanmin.us/
perl $SETUP_DIR/cpanm -l $INST_PATH App::cpanminus
CPANM=`which cpanm`
echo $CPANM

perlmods=( "File::ShareDir" "File::ShareDir::Install" "Bio::Root::Version@1.006924")

for i in "${perlmods[@]}" ; do
  echo -n "Installing build prerequisite $i..."
  $CPANM --notest --mirror http://cpan.metacpan.org -l $INST_PATH $i
done

echo -n "Building bedtools2 ..."
if [ -e $SETUP_DIR/bedtools.success ]; then
  echo -n " previously installed (resumed)...";
elif [ -e $INST_PATH/bin/bedtools ]; then
  echo -n " previously installed ...";
else
  cd $SETUP_DIR
  get_distro "bedtools2" $SOURCE_BEDTOOLS
  mkdir -p bedtools2
  tar --strip-components 1 -C bedtools2 -zxf bedtools2.tar.gz
  make -C bedtools2 -j$CPU
  cp bedtools2/bin/* $INST_PATH/bin/.
  touch $SETUP_DIR/bedtools.success
fi

CURR_TOOL="vcftools"
CURR_SOURCE=$SOURCE_VCFTOOLS
echo -n "Building $CURR_TOOL ..."
if [ -e $SETUP_DIR/$CURR_TOOL.success ]; then
  echo -n " previously installed (resumed) ..."
elif [ -e $INST_PATH/bin/$CURR_TOOL ]; then
  echo -n " previously installed ...";
else
  get_distro $CURR_TOOL $CURR_SOURCE
  cd $SETUP_DIR/$CURR_TOOL
  patch src/perl/Vcf.pm < $INIT_DIR/patches/vcfToolsProcessLog.diff
  ./configure --prefix=$INST_PATH --with-pmdir=$INST_PATH/lib/perl5
  make -j$CPU
  make install
  touch $SETUP_DIR/$CURR_TOOL.success
fi

if [ -e $SETUP_DIR/htslibGet.success ]; then
  echo " already staged ...";
else
  echo
  cd $SETUP_DIR
  get_distro "htslib" $SOURCE_HTSLIB
  touch $SETUP_DIR/htslibGet.success
fi

echo -n "Building htslib ..."
if [ -e $SETUP_DIR/htslib.success ]; then
  echo " previously installed ...";
else
  echo
  mkdir -p htslib
  tar --strip-components 1 -C htslib -jxf htslib.tar.bz2
  cd htslib
  ./configure --enable-plugins --enable-libcurl --prefix=$INST_PATH
  make -j$CPU
  make install
  cd $SETUP_DIR
  touch $SETUP_DIR/htslib.success
fi

export HTSLIB=$INST_PATH

CHK=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Bio::DB::HTS`
if [[ "x$CHK" == "x" ]] ; then
  echo -n "Building Bio::DB::HTS ..."
  if [ -e $SETUP_DIR/biohts.success ]; then
    echo " previously installed ...";
  else
    echo
    cd $SETUP_DIR
    cpanm --no-interactive --notest --mirror http://cpan.metacpan.org -l $INST_PATH Module::Build
    cpanm --no-interactive --notest --mirror http://cpan.metacpan.org -l $INST_PATH Bio::Root::Version
    rm -rf bioDbHts
    get_distro "bioDbHts" $SOURCE_BIOBDHTS
    tar --strip-components 1 -C bioDbHts -zxf bioDbHts.tar.gz
    cd bioDbHts
    perl Build.PL --install_base=$INST_PATH --htslib=$INST_PATH
    ./Build test
    ./Build install
    cd $SETUP_DIR
    rm -f bioDbHts.tar.gz
    touch $SETUP_DIR/biohts.success
  fi
else
  echo "Bio::DB::HTS already installed ..."
fi

if [ -e $SETUP_DIR/samtools.success ]; then
  echo " previously installed ...";
else
echo
  cd $SETUP_DIR
  rm -rf samtools
  get_distro "samtools" $SOURCE_SAMTOOLS
  mkdir -p samtools
  tar --strip-components 1 -C samtools -xjf samtools.tar.bz2
  cd samtools
  ./configure --enable-plugins --enable-libcurl --with-htslib=$HTSLIB --prefix=$INST_PATH
  make -j$CPU all
  make install
  cd $SETUP_DIR
  rm -f samtools.tar.bz2
  touch $SETUP_DIR/samtools.success
fi

# echo -n "Building samtools ..."
# if [ -e "$SETUP_DIR/samtools.success" ]; then
#   echo -n " previously installed (resumed) ...";
# elif [ -e $INST_PATH/bin/samtools ]; then
#   echo -n " previously installed ...";
# else
#   cd $SETUP_DIR
#   get_distro "samtools" $SOURCE_SAMTOOLS
#   cd samtools
#   ./configure --enable-plugins --enable-libcurl --prefix=$INST_PATH
#   make all all-htslib
#   make install install-htslib
#   touch $SETUP_DIR/samtools.success
# fi
#
# CHK=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Bio::DB::HTS`
# if [[ "x$CHK" == "x" ]] ; then
#   echo -n "Building Bio::DB::HTS ..."
#   cd $SETUP_DIR
#   # now Bio::DB::HTS
#   get_file "INSTALL.pl" $BIODBHTS_INSTALL
#   perl -I $PERL5LIB INSTALL.pl --prefix $INST_PATH --static
#   rm -f BioDbHTS_INSTALL.pl
# else
#   echo "Bio::DB::HTS already installed"
# fi

cd $INIT_DIR

echo -n "Installing Perl prerequisites ..."
if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
    echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi

$CPANM --mirror http://cpan.metacpan.org --notest -l $INST_PATH/ --installdeps . < /dev/null

echo -n "Installing vagrent ..."
cd $INIT_DIR
perl Makefile.PL INSTALL_BASE=$INST_PATH
make
make test
make install

# cleanup all junk
rm -rf $SETUP_DIR

echo
echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of PERL5LIB:"
echo "  $PERLROOT"
echo
