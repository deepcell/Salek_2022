Figure 3
========
Description
-----------
The following lists the scripts used to process and analyze the molecular data depicted in Figure 3::

    * figure3_process_snp.sh for figure 3d- This bash scripts takes the fastq file and process to get allele counts. Please note that we added a filtering step to only report the mutations used for Figure 3d, where the original analysis was performed did not include this.
    * figure3_cnv.R for figures 3 e thru g - This R script takes the counts data and normalize to reference data from GM12878 and plots across the genomes, for Figure 3, e thru g

Inputs and Outputs
------------------
figure3_process_snp.sh::

    * fastqs are available upon request, input to figure3_process_snp.sh
    * fig3d_mutations_found_all_unique.bed with the regions used in figure 3d
    * fig3d_mutation_data.csv with the data used to plot figure 3d
    * allele_counts, containing files containing mutation allele counts

figure3_cnv.r::

    * 3_DTC_1mb.csv w/ counts from DTC-sorted and pre-sort samples, where each row represents count observed in 1mb bin
    * 3_GM12878_1mb.csv w/ counts from GM12878 sample, where each row represents count observed in 1mb bin


.. note::
    TODO: chr17,7673774,7673791,7673783,CGGTCTCTCCCAGGACA,CGGTCTCTACCAGGACA
    chr12,25245342,25245359,25245351,TACGCCACCAGCTCCAA,TACGCCACAAGCTCCAA: verified
    chr17,7675199,7675216,7675208,CCAGTTGGCAAAACATC,CCAGTTGGTAAAACATC: verified

