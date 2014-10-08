LICENCE
=======

Copyright (c) 2014 Genome Research Ltd.

Author: Cancer Genome Project <cgpit@sanger.ac.uk>

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


VAGrENT
======

Variation Annotation GENeraTor

A suite of perl modules that compares genomic variations with reference genome annotations
and generates the possible effects each variant may have on the transcripts it overlaps. It evaluates
each variation/transcript combination and describes the effects in the mRNA, CDS and protein sequence
contexts. It provides details of the sequence and position of the change within the transcript/protein
as well as Sequence Ontology terms to classify its consequences.

---

###Dependencies/Install
Some of the code included in this package has dependencies on several C packages:

 * [Samtools](https://github.com/samtools/samtools) - max 0.1.20 until perl bindings are updated
 * [tabix](https://github.com/samtools/samtools/tabix)
 * [vcftools](http://vcftools.sourceforge.net/)

And various perl modules.

Please use `setup.sh` to install the dependencies.  Please be aware that this expects basic C
compilation libraries and tools to be available, most are listed in `INSTALL`.

---

##Creating a release
####Preparation
* Commit/push all relevant changes.
* Pull a clean version of the repo and use this for the following steps.

####Cutting the release
1. Update `lib/Sanger/CGP/Vagrent.pm` to the correct version.
2. Update `Changes` to show major items.
3. Run `./prerelease.sh`
4. Check all tests and coverage reports are acceptable.
5. Commit the updated docs tree and updated module/version.
6. Push commits.
7. Use the GitHub tools to draft a release.
