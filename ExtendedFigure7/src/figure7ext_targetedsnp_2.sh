# This script processes all the fastq.gz files in a directory through to
# 1. alignment to the reference genome using bwa mem
# 2. indexes the bam files with samtools index
# 3. Performs joint variant calling on all the alignment files on the SWIFT panel (97 SNPs + 8 Y)
# 4. Annotates variants with their dbSNP rsIDs
# 5. Extracts reference and alternate read counts into a table
# 6. Finds all variants within +/- 20 b of each variant and their read counts, as a proxy indicator of PCR + sequencing error
# 7. Computes summary/aggregate metrics

# map targeted fastq files to GRCh38
fqzfiles=(*_R1_*.fastq.gz)
for fqz in ${fqzfiles[@]}; do
    prefix="${fqz/_R1_001.fastq.gz/}"
    fqz2=$prefix"_R2_001.fastq.gz"
    bamfile=$prefix".bam"
    if [[ ! -f $bamfile ]];
    then
	bwa mem -M -t 10 GRCh38/hg38.fullAnalysisSet.chroms/grch38chromsbwaidx $fqz $fqz2 \
	    | samtools view -S -b -q 60 \
	    | samtools sort \
	    > $bamfile
    fi
done

# index bam files
bamfiles=(*.bam)
for  bam in ${bamfiles[@]}; do
    samtools index $bam
done

# use bcftools to post process the alignments in the bam files, call variants and extract read counts for each allele
# call variants only at SNP loci

ls *.bam > all_bam_files
bcftools mpileup -Ou -f GRCh38/hg38.fullAnalysisSet.chroms/wg.fa \
    --min-MQ 60 --min-BQ 25 --max-depth 20000 \
    --threads 8 --skip-indels -b all_bam_files \
    --regions-file ~/code/scripts/snp.bed \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    | bcftools call -m -Ov -g 25 -o "variants_joint_deep_testgvcf.vcf" 

# compress and index variant information; annotate with dbSNP rsIDs on GRCh38
bgzip -c variants_joint_deep_testgvcf.vcf > variants_joint.vcf.gz && tabix variants_joint.vcf.gz
bcftools annotate -a ~/Documents/snpdata/dbsnp_common_snps/common_all_20180418_chr.vcf.gz -c ID variants_joint.vcf.gz  > variants_joint_annotated.vcf

# post-process to get a table with the counts
grep CHROM variants_joint_annotated.vcf | sed "s/#//" | \
    awk 'BEGIN{OFS = "\t"; FS = "\t"}{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5, $6, $7); for (i = 10; i <= NF; i++){name=$i; sub(/_.*/,"",name); printf("\tS%s_1\tS%s_2", name, name)}; printf("\n") }' > high_qual_snp_ad_deep.txt
cat variants_joint_annotated.vcf | awk '$4 ~ /[ACGT]/ && $5 ~ /[ACGT\.]/' | awk 'BEGIN{OFS = "\t"}{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5, $6, $7); \
    for (i = 10; i <= NF; i++){ad = $i; sub(/^.*:/, "", ad); n = split(ad, adj, ","); printf("\t%s", adj[1]); if(n > 1) {printf("\t%s", adj[2])} else{ printf("\t0")};}; printf("\n");}' \
      >> high_qual_snp_ad_deep.txt

# snp calling within TP53 amplicons
#bcftools mpileup -Ou -f GRCh38/hg38.fullAnalysisSet.chroms/wg.fa \
#    --min-MQ 60 --min-BQ 25 --max-depth 20000 \
#    --threads 8 --skip-indels -b all_bam_files \
#    --regions-file ~/code/scripts/hglft_genome_380c3_349310.bed \
#    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
#    | bcftools call -m -Ov -o "variants_joint_tp53.vcf" 

# post-process TP53 snps ...
#grep CHROM variants_joint_tp53.vcf | sed "s/#//" | awk 'BEGIN{OFS = "\t"; FS = "\t"}{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5, $6, $7); for (i = 10; i <= NF; i++){name=$i; sub(/_.*/,"",name); printf("\t%s_1\t%s_2", name#, name)}; printf("\n") }' > snp_ad_tp53.txt

#grep -v "^#" variants_joint_tp53.vcf | awk '$5 == "A" || $5 == "C" || $5 == "G" || $5 == "T" || length($4) != length($5)' | awk 'BEGIN{OFS = "\t"}{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5, $6, $7); \
#    for (i = 10; i <= NF; i++){ad = $i; sub(/^.*:/, "", ad); split(ad, adj, ","); printf("\t%s\t%s", adj[1], adj[2]);}; printf("\n");}' \
#      >> snp_ad_tp53.txt

# custom search for H522 indel in TP53
#wt="AAGATGCTGAGGAGGGGCCAG"
#wtrev=(`echo $wt | rev | tr "ACGT" "TGCA"`)
#del="AAGATGCTGAGAGGGGCCAG"
#delrev=(`echo $del | rev | tr "ACGT" "TGCA"`)

#[[ -f "tp53_indel.csv" ]] && rm tp53_indel.csv;

#fqzfiles=(*_R1_*.fastq.gz)
#for fqz in ${fqzfiles[@]}; do
#    prefix="${fqz/_R1_001.fastq.gz/}"
#    fqz2=$prefix"_R2_001.fastq.gz"
#    numwtfwd=(`gzcat $fqz  | grep -e $wt -e $wtrev | wc -l`)
#    numwtrev=(`gzcat $fqz2 | grep -e $wt -e $wtrev | wc -l`)
#    numdelfwd=(`gzcat $fqz  | grep -e $del -e $delrev | wc -l`)
#    numdelrev=(`gzcat $fqz2 | grep -e $del -e $delrev | wc -l`)
#    echo "$prefix,$numwtfwd,$numwtrev,$numdelfwd,$numdelrev" >> tp53_indel.csv
#done


# call variants only at SNP loci +/- 20 bases
bcftools mpileup -Ou -f ~/apps/reference_genomes/GRCh38/hg38.fullAnalysisSet.chroms/wg.fa \
    --min-MQ 60 --min-BQ 25 --max-depth 20000 \
    --threads 8 --skip-indels -b all_bam_files \
    --regions-file ~/code/scripts/snp40.bed \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    | bcftools call -m -Ov -v -o "variants_joint_neighborhood.vcf" 

# post-process neighborhood snps ...
grep CHROM variants_joint_neighborhood.vcf | sed "s/#//" | awk 'BEGIN{OFS = "\t"; FS = "\t"}{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5, $6, $7); for (i = 10; i <= NF; i++){name=$i; sub(/_.*/,"",name); printf("\t%s_1\t%s_2", name, name)}; printf("\n") }' > snp_ad_neighborhood.txt

grep -v "^#" variants_joint_neighborhood.vcf | awk '$5 == "A" || $5 == "C" || $5 == "G" || $5 == "T" ' | awk 'BEGIN{OFS = "\t"}{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5, $6, $7); \
    for (i = 10; i <= NF; i++){ad = $i; sub(/^.*:/, "", ad); split(ad, adj, ","); printf("\t%s\t%s", adj[1], adj[2]);}; printf("\n");}' \
      >> snp_ad_neighborhood.txt


# summaries
[[ -f "fastq_readcounts.csv" ]] && rm fastq_readcounts.csv;
[[ -f "all_good_alignments.csv" ]] && rm all_good_alignments.csv;
[[ -f "snp_coverage.csv" ]] && rm snp_coverage.csv; 
[[ -f "y_coverage.csv" ]] && rm y_coverage.csv; 
[[ -f "snp_reads_summary.csv" ]] && rm snp_reads_summary.csv; 
[[ -f "y_reads_summary.csv" ]] && rm y_reads_summary.csv; 
   
for fqz in ${fqzfiles[@]}; do
    prefix="${fqz/_L001_R1_001.fastq.gz/}"
    bamfile=$prefix"_L001.bam"
    gzcat $fqz | wc -l | awk -v fn=$prefix 'BEGIN{OFS = ","}{print (fn, $1/2)}' >> fastq_readcounts.csv
    samtools view $bamfile | grep "[1-9][0-9][0-9]M" | cut -f3,4 | uniq -c | sed "s/^ *//" | awk -v fn=$prefix 'BEGIN{OFS = ","; sum=0}{sum += $1}END{print (fn, sum)}' >> all_good_alignments.csv
    bedtools coverage -a ~/code/scripts/snp.bed -b $bamfile -counts | awk -v prefix=$prefix 'BEGIN{OFS = ","}{print prefix, $1, $3, $4}' | tee -a snp_coverage.csv | \
	awk 'BEGIN{FS=","; OFS=","; sum = 0; prevBam = "none"; counter = 0;}{if (prevBam != $1 && counter > 0){print (prevBam, sum); sum = 0;}; counter++; sum = sum + $4; prevBam = $1;}END {print ($1, sum)}' \
	    >> snp_reads_summary.csv
    bedtools coverage -a ~/code/scripts/hglft_swift_y_GRCH38.bed -b $bamfile -counts | awk -v prefix=$prefix 'BEGIN{OFS = ","}{print prefix, $1, $3, $4}'  | tee -a y_coverage.csv | \
	awk 'BEGIN{FS=","; OFS=","; sum = 0; prevBam = "none"; counter = 0;}{if (prevBam != $1 && counter > 0){print (prevBam, sum/8); sum = 0;}; counter++; sum = sum + $4; prevBam = $1;}END {print ($1, sum)}' \
	    >> y_reads_summary.csv
done

paste -d',' fastq_readcounts.csv <(cut -d, -f2 all_good_alignments.csv) <(cut -d, -f2 snp_reads_summary.csv) <(cut -d, -f2 y_reads_summary.csv) | \
    awk 'BEGIN{FS = ","; OFS = ","}{print $0, $4 + $5, $4/97, $5/8}' > summary_metrics.csv

# postprocess with an R script (combining here)
/usr/local/bin/Rscript ~/src/figure7ext_primary_analysis.R

# copy files
analysis_dir="analysis_files/"
[ -f $analysis_dir ] || mkdir $analysis_dir
cp summary_metrics_with_assay_noise.csv $analysis_dir
cp summary_metrics_with_assay_noise.csv $analysis_dir
cp *qc.csv $analysis_dir
cp snp_allele_fractions*.csv $analysis_dir
cp snp_ad_neighborhood.txt $analysis_dir
cp y_reads_summary.csv $analysis_dir
