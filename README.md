VAGrENT
======

Variation Annotation GENeraTor

A suite of perl modules that compares genomic variations with reference genome annotations
and generates the possible effects each variant may have on the transcripts it overlaps. It evaluates
each variation/transcript combination and describes the effects in the mRNA, CDS and protein sequence
contexts. It provides details of the sequence and position of the change within the transcript/protein
as well as Sequence Ontology terms to classify its consequences.

| Master | Dev |
|---|---|
|  [![Build Status](https://travis-ci.org/cancerit/VAGrENT.svg?branch=master)](https://travis-ci.org/cancerit/VAGrENT) | [![Build Status](https://travis-ci.org/cancerit/VAGrENT.svg?branch=dev)](https://travis-ci.org/cancerit/VAGrENT) |

---

### Dependencies/Install
Some of the code included in this package has dependencies on several packages:

 * [Samtools v1.3+](https://github.com/samtools/samtools)
 * [vcftools](https://vcftools.github.io/)
 * [Bio::DB::HTS](http://search.cpan.org/~rishidev/Bio-DB-HTS/)
 * [bedtools2](http://bedtools.readthedocs.io/en/latest/index.html)
   * Not >=2.24.0, no upgrades until [this ticket](https://github.com/arq5x/bedtools2/issues/319) is resolved (which may involve code changes)

And various perl modules.

Please use `setup.sh` to install the dependencies.  Setting the environment variable `CGP_PERLLIBS` allows you to to append to `PERL5LIB` during install.  Without this all dependancies are installed into the target area.

Please be aware that this expects basic C compilation libraries and tools to be available.

---

## Creating a release
#### Preparation
* Commit/push all relevant changes.
* Pull a clean version of the repo and use this for the following steps.

#### Cutting the release
1. Update `lib/Sanger/CGP/Vagrent.pm` to the correct version.
2. Update `CHANGES.md` to show major items.
3. Run `./prerelease.sh`
4. Check all tests and coverage reports are acceptable.
5. Commit the updated docs tree and updated module/version.
6. Push commits.
7. Use the GitHub tools to draft a release.

LICENCE
=======

Copyright (c) 2014-2018 Genome Research Ltd.

Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>

This file is part of VAGrENT.

VAGrENT is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
2009, 2011-2012’ should be interpreted as being identical to a statement that
reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012’."
