##################
# This R script takes the results of targetedsnp.sh, which are the summary measures + allele read counts, and performs the following:
# 1. Sample QC
# 2. Assay QC
# 3. Allele fraction calculation
# 4. Write out of data in long form with QCs joined
#if (!require('scales')) install.packages('scales'); library('scales')
eps <- 1e-5

# read SNP data
d <- read.table("high_qual_snp_ad_deep.txt", header = T, sep = "\t", as.is = T, check.names = F)
names(d)[8:ncol(d)] <-  substr(names(d)[8:ncol(d)], 2, nchar(names(d)[8:ncol(d)]))
sample_names <- unique(sub("_[12]","",names(d)[8:ncol(d)]))
sample_names <- sample_names[sample_names != "Undetermined"]

# read summary metrics
metrics <- read.csv("summary_metrics.csv", as.is = T, header = F)
names(metrics) <- c("sample_name","reads","alignments","snp_reads","y_reads","target_reads", "snp_ave","y_ave")
metrics$y_ratio <- with(metrics, y_ave/snp_ave)
metrics$sample_name <- sub("_S[0-9]*$","", metrics$sample_name)

# long form of SNP data
dlong <- c()
for (samp in sample_names) {
  d1s <- d[,c("CHROM","POS","REF",paste0(samp, "_1"),paste0(samp, "_2"))]
  names(d1s)[4:5] <- c("nr","na")
  d1s$sample_name <- samp
  dlong <- rbind(dlong, d1s)
}

dlong$total <- with(dlong, nr + na)
dlong$af <- with(dlong, na/total)

# sample QC on SNP data (can't use Y because not always present)
sample_medians <- with (dlong, aggregate(total, list(sample_name), median))
names(sample_medians) <- c("sample_name", "sample_median")
sample_medians$sample_qc <- ifelse(sample_medians$sample_median  < 200, "FAIL","PASS")
write.csv(sample_medians, "sample_qc.csv", row.names = F, quote = F)
dlong <- merge(dlong, sample_medians)

# Assay QC on SNP data 
dlong$scaled_total <- with(dlong, total/sample_median)
assay_medians <- with(dlong[dlong$sample_qc == "PASS",], aggregate(scaled_total, list(CHROM,POS), median))
names(assay_medians) <- c("CHROM","POS","assay_median")
assay_medians$assay_qc <- ifelse(assay_medians$assay_median < 0.1, "FAIL","PASS")
dlong <- merge(dlong, assay_medians)

# read Y assay data
y_coverage <- read.csv("y_coverage.csv", header = F, as.is = T)
names(y_coverage) <- c("sample_name","chromosome","position","nreads")
y_coverage$sample_name <- sub("_S[0-9]*$","", y_coverage$sample_name)
y_coverage <- merge(y_coverage, sample_medians)

# Y assay QC and append to SNP assay QC
male_samples <- metrics$sample_name[metrics$y_ratio > 0.2]
male_samples <- male_samples[y_coverage$sample_name %in% male_samples & y_coverage$sample_qc == "PASS"]
if (length(male_samples) > 0) {
  y_qc <- with(y_coverage[y_coverage$sample_name %in% male_samples & y_coverage$sample_qc == "PASS",],
               aggregate(nreads/(sample_median + eps), list(chromosome, position), median))
  names(y_qc) <- c("CHROM","POS", "assay_median")
  y_qc$assay_qc <- with(y_qc, ifelse(assay_median < 0.1, "FAIL","PASS"))
  assay_medians <- rbind(assay_medians, y_qc)
  
}
write.csv(assay_medians, "locus_qc.csv", row.names = F, quote = F)

# write out the long-form data
write.csv(dlong, "snp_allele_fractions.csv", row.names = F, quote = F)

# also turn the snp_allele_fractions data into a pivot for Bio folk
saf <- read.csv("snp_allele_fractions.csv",  as.is = T, check.names = F)
saf <- saf[saf$sample_qc == "PASS" & saf$assay_qc == "PASS",]
saf$af[saf$total < 50] <- NA

saf_pivot <- c()
for (samp in unique(saf$sample_name)) {
  di <- saf[saf$sample_name == samp,c("CHROM","POS","af")]
  names(di)[3] <- samp
  if (is.null(saf_pivot)) {
    saf_pivot <- di
  } else {
    saf_pivot <- merge(saf_pivot, di)
  }
}
write.csv(saf_pivot, "snp_allele_fractions_pivot.csv", row.names = F, quote = F)


# look at neighborhood and estimate assay error, join to metrics
sad <- read.table("snp_ad_neighborhood.txt", header = T, sep = "\t", as.is = T, check.names = F)

# compute the fractions for each neighborhood SNP
sample_names <- unique(sub("_[12]$", "", names(sad)[grep("_[12]$", names(sad))]))
sadp <- sad[,1:7]

for (samp in sample_names) {
  s1 <- sad[,paste0(samp, "_1")]
  s2 <- sad[,paste0(samp, "_2")]
  tot <- s1 + s2
  sadp[,samp] <- ifelse(tot < 50, NA, s2/(s1+s2))
}
sad <- sadp

purity_estimates_neigh <- apply(sad[,8:ncol(sad)], 2, function(x) 4 * sum(x, na.rm = T) / (40*89) )
pen <- data.frame(purity_estimates_neigh)
pen$sample_name <- names(purity_estimates_neigh)
info_snp_neigh_1pct <- apply(sad[,8:ncol(sad)], 2, function(x) length(x[!is.na(x) & x > 0.01]))
isn1 <- data.frame(info_snp_neigh_1pct)
isn1$sample_name <- names(info_snp_neigh_1pct)

info_snp_neigh_10pct <- apply(sad[,8:ncol(sad)], 2, function(x) length(x[!is.na(x) & x > 0.1]))
isn10 <- data.frame(info_snp_neigh_10pct)
isn10$sample_name <- names(info_snp_neigh_10pct)

summary_metrics <- merge(merge(merge(metrics, pen), isn1), isn10)
write.csv(summary_metrics, "summary_metrics_with_assay_noise.csv", row.names = F, quote = F)
