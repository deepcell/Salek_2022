Figure 3
========

The data is available from s3://deepcell-public/Salek_2022/.
This directory contains the following files::

    * README.rst
    * src/figure3_cnv.R
    * src/figure3_process_snp.sh
    * data/_1_DTC_RSEC_MolsPerCell.csv
    * data/3_DTC1.mat
    * data/3_DTC_1mb.csv
    * data/3_GM12878_1mb.csv
    * data/DTC1_RSEC_MolsPerCell.csv
    * data/fig3d_all_mutation_strings.txt
    * data/fig3d_mutation_data.csv
    * data/fig3d_mutations_found_all_unique.bed
    * allele_counts/*_allele_counts.csv
    * bams/*bam
    * bams/*bam.bai
    * bams/samples.tsv

..Note::
    fastq files and aligned reads files used for CNV as well as mutation analysis are available upon request only.

Description
-----------
The scripts listed here are also available to download from https://github.com/deepcell/Salek_2022/.
The following lists the scripts in `src` directory used to process and analyze the molecular data depicted in Figure 3::

    * figure3_process_snp.sh for figure 3d - This script processes the data for mutation analysis::
        
        1. alignment to the reference genome using bwa mem
        2. indexes the bam files with samtools index
        3. add readgroups
        4. generate allele counts, note that the data is filtered to only report the mutations used for Figure 3d
    
    * figure3_cnv.R for figures 3 e thru g - This R script takes the counts data and normalize to reference data from GM12878 and plots across the genomes, for Figure 3, e thru g

        1. process raw sequencing data using the whole genome sequencing pipeline (available from https://github.com/deepcell/Salek_2022_NatBiotechnol)
        2. read the counts data (number of observed counts per 1 mb-sized bin)
        3. normalize the data to normal reference sample with the ploidy of 2
        4. calculate counts-based coverage

The following lists the BD dataview session files (.mat) saved used in Figure 3h and 3i

    * 3_DTC1.mat - This file loads a data table file, named `_1_DTC_RSEC_MolsPerCell.csv` and `DTC1_RSEC_MolsPerCell.csv` for scRNA-seq analysis

Open DB DataView, under Load data on upper left corner, click "Load session (.mat)" button to load the session file used for 3h and 3i.
The tabular formatted counts data used for the analysis is also provided in the `data` directory.
The raw sequencing data will be available upon request.

Inputs
------
figure3_process_snp.sh::

    * fastqs are available upon request, input to figure3_process_snp.sh
    * fig3d_mutations_found_all_unique.bed with the regions used in figure 3d
    * fig3d_mutation_data.csv with the data used to plot figure 3d

figure3_cnv.r::

    * 3_DTC_1mb.csv w/ counts from DTC-sorted and pre-sort samples, where each row represents counts observed in 1mb bin
    * 3_GM12878_1mb.csv w/ counts from GM12878 sample, where each row represents counts observed in 1mb bin

Outputs
-------
allele_counts directory contains files with mutation allele counts used to plot figure 3d.
bams directory contains aligned reads overlaps to the following regions, used for figure 3d. See samples.tsv for sample type, Pre-sorted vs. Sorted::

    * chr12   25245342        25245359
    * chr17   7673774 7673791
    * chr17   7675199 7675216

.. note::
    chr17,7673774,7673791,7673783,CGGTCTCTCCCAGGACA,CGGTCTCTACCAGGACA: todo
    chr12,25245342,25245359,25245351,TACGCCACCAGCTCCAA,TACGCCACAAGCTCCAA: verified
    chr17,7675199,7675216,7675208,CCAGTTGGCAAAACATC,CCAGTTGGTAAAACATC: verified
