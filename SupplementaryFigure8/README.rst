Supplementary Figure 8
======================

This directory contains the following files::

    * README.rst
    * src/primary_analysis.R
    * src/targetedsnp.sh
    * src/genotype_and_mixture_estimation.R

Set the directory containing the fastq files as working directory and process the data using the commands described in `targetedsnp.sh`, followed by `primary_analysis.R`.
Finally use the output, `snp_allele_fractions.csv`, to run `genotype_and_mixture_estimation.R` to create the final plot used for Supplementary Figure 8.

Code Description
----------------
The following lists the scripts in `src` directory used to process the sample::

    * targetedsnp.sh - This script processes the data::

        1. alignment to the reference genome using bwa mem
        2. indexes the bam files with samtools index
        3. Performs joint variant calling on all the alignment files on the SWIFT panel (97 SNPs + 8 Y)
        4. Annotates variants with their dbSNP rsIDs
        5. Extracts reference and alternate read counts into a table
        6. Finds all variants within +/- 20 b of each variant and their read counts, as a proxy indicator of PCR + sequencing error
        7. Computes summary/aggregate metrics

    * primary_analysis.R - This R script takes the results of targetedsnp.sh, which are the summary measures + allele read counts, and performs the following::

        1. Sample QC
        2. Assay QC
        3. Allele fraction calculation
        4. Write out of data in long form with QCs joined

    * genotype_and_mixture_estimation.R - This R script plots the figure used in Supplementary Figure 8.
