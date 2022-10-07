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
	bwa mem -M -t 10 /bioinformatics/reference_data/grch38chromsbwaidx $fqz $fqz2 \
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

# compile number of good alignments per Mb in each of the bam files

[[ -f "genome_mb_coverage.csv" ]] && rm genome_mb_coverage.csv;

for bam_file in ${bamfiles[@]}; do
    prefix="${bam_file/_L001.bam/}"
    samtools view $bam_file | awk '$5 > 50 {print $3, $4}' | uniq | awk  '{print $1, int($2/1000000)}' | uniq -c | sed "s/^ *//; s/ /,/g" | awk -v bf=$prefix 'BEGIN{FS = ","; OFS = ","}{print $1, $2, $3, bf}'  >> genome_mb_coverage.csv
done;

# count up GC per file too

[[ -f "gc_histogram.csv" ]] && rm gc_histogram.csv;
for bam_file in ${bamfiles[@]}; do
    samtools view $bam_file | awk '$5 > 50' | awk '{ n=length($10); print int(gsub(/[GCgc]/,"",$10)/n * 1000)/1000;}' | sort | uniq -c | sed "s/^ *//; s/ /,/" | awk -v bf=$bam_file 'BEGIN{FS=","; OFS=","}{print $1,$2, bf}' >> gc_histogram.csv
done;

# get GC distribution of 150 b segments in the human reference genome
#grep -v ">"  ~/apps/reference_genomes/GRCh38/hg38.fullAnalysisSet.chroms/wg.fa | grep -v "N" | grep -v "[acgt]" | awk '{ n=length($1); print int(gsub(/[GCgc]/,"",$1)/n * 1000)/1000;}' #| awk 'BEGIN{i=0; sum = 0; OFS = ","}{sum = sum + $1; i++; if (i % 3 == 0) {print int(sum/3 * 1000)/1000; sum=0} }' | sort |  uniq -c  | sed "s/^ *//; s/ /,/" > genome_gc_histogram.csv

# GC and non-N bases per Mb of genome to join with cov data
#bedtools getfasta -fi /bioinformatics/reference_data/wg.fa -bed /bioinformatics/reference_data/grch38_mb.bed  -tab | tr "acgt" "ACGT" | sed 's/Chr/chr/' | awk 'BEGIN{FS = "\t"; OFS = ","}{split($1,a,":"); split(a[2],b,"-"); print (a[1],b[1], b[2], $2)}' | awk  'BEGIN{FS = ","; OFS = "\t"}{a = $4; gsub(/N/,"",a); b = a; gsub(/[AT]/,"",b); print $1, $2, $3, length(a), length(b)}' > grch38_mb_gc.bed 

# Mark duplicates and get metrics on them
#for bam_file in ${bamfiles[@]}; do
#    ofile="dup_marked_"$bam_file
#    prefix="${bam_file/.bam/}"
#    mfile="dup_metrics_"$prefix".txt"
#    java -jar ~/apps/picard/picard.jar MarkDuplicates \
#      I=$bam_file \
#      O=$ofile \
#      M=$mfile
#done;

[[ -f "duplicate_stats.csv" ]] && rm duplicate_stats.csv;
for bam_file in ${bamfiles[@]}; do
    prefix="${bam_file/_L001.bam/}"
    n=(`samtools view $bam_file | wc -l`)
    n1=(`samtools view $bam_file | cut -f3,4 | uniq | wc -l`)
    echo $prefix","$n","$n1 >> duplicate_stats.csv
done;

for bam_file in ${bamfiles[@]}; do
    prefix="${bam_file/.bam/}"
    samtools stats $bam_file > $prefix"_samtools_stats.txt"
done;

# extract mismatch rate from samtools stats
[[ -f "mismatch_stats.csv" ]] && rm mismatch_stats.csv;
[[ -f "gc_stats.csv" ]] && rm gc_stats.csv;

for bam_file in ${bamfiles[@]}; do
    prefix="${bam_file/_L001.bam/}"
    statfile=$prefix"_L001_samtools_stats.txt"
    grep "# mismatch" $statfile | cut -f3 | awk -v bf=$prefix 'BEGIN{OFS = ","}{print bf, $1}'  >> mismatch_stats.csv
    grep ^GCD $statfile | awk -v bf=$prefix 'BEGIN{OFS = ","}{print bf, $2,$3, $4, $5,$6,$7,$8}'  >> gc_stats.csv
done;

# copy files to NAS
wgs_analysis_dir="/Volumes/Biology/NGS/WGS/WGS29/analysis_files"
#mkdir $wgs_analysis_dir
#cp duplicate_stats.csv $wgs_analysis_dir
#cp gc_stats.csv $wgs_analysis_dir
#cp genome_mb_coverage.csv $wgs_analysis_dir
#cp gc_histogram.csv $wgs_analysis_dir
#cp mismatch_stats.csv $wgs_analysis_dir


# custom work for WGS15 where there was a mixup in the initial setup on the sequencer
#cp 0819*.bam ../WGS15_2
#cp 2120-S-20210304-Aries-ResolveDNA_S39_L001.bam ../WGS15_2/2120-S-20210304-Taurus-ResolveDNA_S38_L001.bam
#cp 2120-B-20210304-ResolveDNA_S40_L001.bam ../WGS15_2/2120-S-20210304-Aries-ResolveDNA_S39_L001.bam
#cp 2120-Epcam-20210304-ResolveDNA_S41_L001.bam ../WGS15_2/2120-B-20210304-ResolveDNA_S40_L001.bam
#cp 2120-CD45-20210304-ResolveDNA_S42_L001.bam ../WGS15_2/2120-Epcam-20210304-ResolveDNA_S41_L001.bam
#cp Live12878-SOP-ResolveDNA_S43_L001.bam ../WGS15_2/2120-CD45-20210304-ResolveDNA_S42_L001.bam
#cp  ../WGS15_2/Live12878-SOP-ResolveDNA_S43_L001.bam
