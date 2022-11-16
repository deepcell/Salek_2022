Extended Figure 7
=================

The data is available from s3://deepcell_public/Salek_2022_NatBiotechnol/.
This directory contains the following files::

    * README.rst
    * src/figure7ext_primary_analysis.R
    * src/figure7ext_targetedsnp_2.sh
    * src/figure7ext_genotype_and_mixture_estimation.R
    * data/fastq_files_list.txt

Set the directory containing the fastq files listed under `fastq_files_list.txt` as working directory and process the data using the commands described in `figure7ext_targetedsnp_2.sh`, followed by `figure7ext_primary_analysis.R`.
Finally use the output, `snp_allele_fractions.csv`, to run `figure7ext_genotype_and_mixture_estimation.R` to create the final plot used for extended figure 7.

..Note::
    fastq files listed under fastq_files_list.txt are available upon request only.

Code Description
----------------
The scripts listed here are also available to download from https://github.com/deepcell/Salek_2022_NatBiotechnol/.
The following lists the scripts in `src` directory used to process the sample::

    * figure7ext_targetedsnp_2.sh - This script processes the data::

        1. alignment to the reference genome using bwa mem
        2. indexes the bam files with samtools index
        3. Performs joint variant calling on all the alignment files on the SWIFT panel (97 SNPs + 8 Y)
        4. Annotates variants with their dbSNP rsIDs
        5. Extracts reference and alternate read counts into a table
        6. Finds all variants within +/- 20 b of each variant and their read counts, as a proxy indicator of PCR + sequencing error
        7. Computes summary/aggregate metrics

    * figure7ext_primary_analysis.R - This R script takes the results of e2etargetedsnp_2.sh, which are the summary measures + allele read counts, and performs the following::

        1. Sample QC
        2. Assay QC
        3. Allele fraction calculation
        4. Write out of data in long form with QCs joined

    * figure7ext_genotype_and_mixture_estimation.R - This R script plots the figure used in ExtendedFigure7

Inputs
------
Fastq files listed in `fastq_files_list.txt` are the required inputs to `figure7ext_targetedsnp_2.sh`.
These raw read files are available upon request.
