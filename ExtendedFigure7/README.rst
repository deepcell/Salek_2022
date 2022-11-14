Extended Figure 7
=================
Description
-----------
The following lists the scripts used to process the sample::

    * e2etargetedsnp_2.sh - This script processes all the fastq.gz files in a directory through to::

        1. alignment to the reference genome using bwa mem
        2. indexes the bam files with samtools index
        3. Performs joint variant calling on all the alignment files on the SWIFT panel (97 SNPs + 8 Y)
        4. Annotates variants with their dbSNP rsIDs
        5. Extracts reference and alternate read counts into a table
        6. Finds all variants within +/- 20 b of each variant and their read counts, as a proxy indicator of PCR + sequencing error
        7. Computes summary/aggregate metrics

    * e2e_primary_analysis.R - This R script takes the results of e2etargetedsnp_2.sh, which are the summary measures + allele read counts, and performs the following::

        1. Sample QC
        2. Assay QC
        3. Allele fraction calculation
        4. Write out of data in long form with QCs joined

    * paper_custom_genotype_and_mixture_estimation.R - This R script plots the figure used in ExtendedFigure7
