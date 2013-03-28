# Genotypes in SciDB

Adds support to [SciDB](http://www.scidb.org/) for loading genotype
data from VCF (Variant Call Format) files. It includes custom data
types and loaders.

## Compiling

To compile, add the genotypes directory to the examples
subdirectory in the SciDB source tree, and include an
'add_subdirectory("genotypes")' directive in the
examples/CMakeLists.txt file.

## Loading Data

To load a typical set of VCF files into a set of related SciDB arrays
with the common basename "foobar" you would use the following
commands on the SciDB coordinator node:

        $ loadgt.sh -d foobar_samples.csv foobar *.vcf.gz
        $ samples_load.sh foobar foobar_samples.csv

The 'foobar_samples.csv' is comma separated list of samples
containing the list of samples which matches those in the VCF header
files, the populations to which they belong, and their sex.

        #col,sample,population,founder,sex
        1,NA19625,ASW,true,F
        2,NA19700,ASW,true,M
        3,NA19701,ASW,true,F
        4,NA19703,ASW,true,M
        5,NA19704,ASW,true,F

For optimal performance of queries involving an entire population,
this file should be sorted by population, then by sample. The sample
columns in the VCF file will be reordered according to the column
numbers listed in this file. An convience script is included which may
be used to create a sample file which can then be editted by hand.

        $ samples_collect_from_vcf.sh *.vcf.gz > foobar_samples.csv

## Analysis

To create an array containing allele counts for each population, for
each variation site:

        $ allele_counts_calculate.py foobar

This will create an array with the three dimensions that denote a
variation (chrom, pos, and var) and two additional dimensions,
population and allele. In the allele dimension, an index of -1 is the
total number of alleles, 0 is the number of alleles that are the same
as the reference allele, 1 is the first alternate, 2 is the second,
and so forth.
