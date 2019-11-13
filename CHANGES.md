# CHANGES

## 3.5.1

 * Handle chr not present in itree hash when running `Admin_GeneRegionBedDumper.pl`

## 3.5.0

* Loads vagrent cache into IntervalTree to speed up processing by:
  * reducing redundant/random disk access
  * reducing repeated decompression of same data when events are local to each other

## 3.4.0

* Add Dockerfile and supporting scripts
* Switch travis-ci to build and test under docker

## 3.3.5

* fixing pipefail issue for non-bash environments, splitting command into 2 separate executions.

## 3.3.4

* Don't assume working directory is writable (fixes #33)

## 3.3.3

* fixing issue in Admin_CacheFileBuilder.pl where it would fail without the optional fai file

## 3.3.2

* Actually get the new `htslib` and `samtools` links correct in `setup.sh`
* Move `samtools` install in `setup.sh` to happen before `Bio::DB::HTS`

## 3.3.1

* `samtools` and `htslib` updated to 1.7
* `Bio::DB::HTS` updated to 2.10, fixes error with GRCh38 contig names in Tabix.
* Update tabix query to use `query_full`

## 3.3.0

* Complete re-write of the reference generation code

## 3.2.2

* Add bedtools2 to `setup.sh`
* Added bedtools2 to `README.md`
* Changes `Bio::DB::HTS`, `samtools` and `HTSlib` install methods.
* Corrected condition indicating sort is required.
* Fixes #23 Changed from vcf-sort to normal linux sort to ensure multiple indels with
same start coord are sorted in a stable way.

## 3.2.0

* Allows use of ensemblgenomes.org as a datasource
* Handle genes without names, and give more useful error message
* Ensure an error code is emmitted on failure

## 3.1.0

* Adds travis testing
* Cleans up install script and adds multi versioned paths to options

## 3.0.0

* Removed use of legacy Tabix codebase and switches to samtools 1.2+
* Switches to gihub version of vcftools and makes several patches unnecessary
* Several bugfixes
