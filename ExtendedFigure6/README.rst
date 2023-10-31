Extended Figure 6
=================

The data is available from s3://deepcell-public/Salek_2022/.
This directory contains the following files::

    * README.rst
    * src/figure6.R

Set the directory containing `snp_allele_fractions.csv` as the working directory.
Then run the following commands::

    Rscript src/figure6.R

Code Description
----------------
    * figure6.R::

        1. Maximum likelihood estimation of sample mixture based on SNP allele frequencies
        2. Plot estimated sample mixture against true sample mixture

Input
------
SNP allele fraction file in s3 named `snp_allele_fractions.csv`