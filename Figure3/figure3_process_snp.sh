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


# add readgroups and process the two panels separately
swift_lung_bam_files=(*_L001.bam)
for bamfile in ${swift_lung_bam_files[@]}; do
    prefix="${bamfile/.bam/}"
    bam2=$prefix"_RG.bam"
    samp="${prefix/_L001/}"
    vcf=$samp".vcf"
    gatk AddOrReplaceReadGroups -I $bamfile --RGLB SwiftLung --RGPL Illumina --RGPU 3 --RGSM $samp -O $bam2
    samtools index $bam2
done



# generate allele counts
all_bamfiles=(*_RG.bam)
for bam in ${all_bamfiles[@]}; do
    prefix="${bam/_RG.bam/}"
    allele_counts_file=$prefix"_allele_counts.csv"
    samtools view $bam -L fig3d_mutations_found_all_unique.bed \
    | grep -f fig3d_all_mutation_strings.txt -o | sort | uniq -c | sed 's/^ *//; s/ /,/g' | awk 'BEGIN{FS = ","; OFS = ","}{print $2, $1}' > $allele_counts_file
done

ls *allele_counts.csv > allele_counts_files
