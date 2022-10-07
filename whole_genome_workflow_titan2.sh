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
