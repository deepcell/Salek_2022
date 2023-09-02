Figure 3
========

This directory contains the following files::

    * README.rst
    * src/plot_cnv.R
    * src/process_snp.sh
    * data/targets.bed
    * data/target_strings.txt

Description
-----------

The following lists the scripts in `src` directory used to process and analyze the molecular data depicted in Figure 3::

    * process_snp.sh for figure 3c - This script processes the data for mutation analysis::
        
        1. alignment to the reference genome using bwa mem
        2. indexes the bam files with samtools index
        3. add readgroups
        4. generate allele counts
    
    * plot_cnv.R for figures 3d thru f - This R script takes the counts data and normalize to reference data from GM12878 and plots

        1. read the counts data (number of observed counts per 1 mb-sized bin) from whole genome sequencing pipeline
        2. normalize the data to normal reference sample with the ploidy of 2
        3. calculate counts-based coverage
