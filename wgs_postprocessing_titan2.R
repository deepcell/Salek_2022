#setwd(<directory with all the data>)
options(width=600)

library(scales)

cov_mb <- read.csv("genome_mb_coverage.csv", header = F, as.is = T)
names(cov_mb) <- c("cov","chrom","pos_mb","sample_name")
samples <- unique(cov_mb$sample_name)

# also read in  data from WGS14 for normalization
cov2 <- read.csv("/media/Biology/NGS/WGS/WGS14/analysis_files/genome_mb_coverage.csv", header = F, as.is = T)
names(cov2) <- c("cov","chrom","pos_mb","sample_name")
samples2 <- unique(cov2$sample_name)

normalization_samples <-
  c("GM12878-100c-Res-KAPA_S3" ,"GM12878-200c-Res-KAPA_S4",
    "GM24143-100c-Res-KAPA_S7","GM24143-200c-Res-KAPA_S8","GM24143-50c-Res-KAPA_S6",
    "GM24149-100c-Res-KAPA_S11","GM24149-200c-Res-KAPA_S12","GM24149-50c-Res-KAPA_S10",
    "GM24385-100c-Res-KAPA_S15","GM24385-200c-Res-KAPA_S16","GM24385-50c-Res-KAPA_S14",
    "HFL1-100c-Res-KAPA_S23","HFL1-200c-Res-KAPA_S24","HFL1-50c-Res-KAPA_S22")

cov_mb <- rbind(cov_mb, cov2[cov2$sample_name %in% normalization_samples,])
samples <- unique(cov_mb$sample_name)

# for each sample, determine the internal scaled cov_mb_scaled relative to mean
# then determine overall bias dist
# then scale each coverage by bias (or adjust by PCs for overkill)

cov_mb$sample_mean <- 0
for (samp in samples) {
  cov_mb$sample_mean[cov_mb$sample_name == samp] = mean(cov_mb$cov[cov_mb$sample_name == samp])
}
cov_mb$scaled_cov <- with(cov_mb, cov/sample_mean)
#print(cov_mb[!duplicated(cov_mb[, c("sample_name")]), ])

# this is the per-mb bin coverage median across all the reference samples
cov_bias_median <- with(cov_mb[cov_mb$sample_name %in% normalization_samples,], aggregate(scaled_cov, list(chrom, pos_mb), median))
names(cov_bias_median) <- c("chrom","pos_mb","bias_median")
cov_bias_mad  <- with(cov_mb[cov_mb$sample_name %in% normalization_samples,], aggregate(scaled_cov, list(chrom, pos_mb), mad))
names(cov_bias_mad) <- c("chrom","pos_mb","bias_mad")

cov_bias <- merge(cov_bias_median, cov_bias_mad)
cov_bias <- cov_bias[order(cov_bias$chrom, cov_bias$pos_mb),]
cov_bias$bias_cv <- with(cov_bias, bias_mad/bias_median)
cov_mb <- merge(cov_mb, cov_bias)
cov_mb <- cov_mb[order(cov_mb$sample_name, cov_mb$chrom, cov_mb$pos_mb),]
cov_mb$loc <- 1
for (i in 2:nrow(cov_mb)) {
  cov_mb$loc[i] <- ifelse(cov_mb$sample_name[i] == cov_mb$sample_name[i-1], cov_mb$loc[i-1] + 1, 1)
}
# cv >= 20% bins were removed
# exploratory analysis on the distribution of normalized coverage ad cv
cov_mb <- cov_mb[cov_mb$bias_cv < 0.2,]

# restrict to autosomes  for now. Don't have data for Y in particular
cov_mb <- cov_mb[!(cov_mb$chrom %in% c("chrEBV","chrM","chrY","chrX" )),]
cov_mb$chrom_num <- -1
cov_mb$chrom_num <- ifelse(cov_mb$chrom == "chrX", 23, as.numeric(sub("chr","",cov_mb$chrom)))
cov_mb <- cov_mb[order(cov_mb$sample_name, cov_mb$chrom_num, cov_mb$pos_mb),]
cov_mb$loc <- 1
for (i in 2:nrow(cov_mb)) {
  cov_mb$loc[i] <- ifelse(cov_mb$sample_name[i] == cov_mb$sample_name[i-1], cov_mb$loc[i-1] + 1, 1)
}

# now build all of the samples, one by one. Hope CNVs show up!
cov_mb$scaled_scaled_cov <- with(cov_mb, scaled_cov / bias_median)

pdf("sample_cnv_profiles_normalized.pdf")
for (samp in samples[!(samples %in% normalization_samples)]) {
  with(cov_mb[cov_mb$sample_name == samp,],
       plot(loc, scaled_scaled_cov, pch = 19, col = alpha("darkblue",0.3), cex = 0,
            main = samp, ylim = c(0,3), xlab = "genomic position (Mb)", ylab = "scaled coverage"))
  for (chrom_num in 1:23) {
    with(cov_mb[cov_mb$sample_name == samp & cov_mb$chrom_num == chrom_num,],
         points(loc, scaled_cov/bias_median, pch = 19, cex = 0.5, col = alpha(ifelse(chrom_num %% 2 == 0, "darkblue","darkred"),0.3)))
  }
  abline(h = seq(0,3,0.5), lty = "dotted", col = "grey")
}
dev.off()

# get chromosome scaling for crude aneuploidy detection
cov_chrom <- with(cov_mb, aggregate(scaled_scaled_cov, list(sample_name, chrom, chrom_num), median))
names(cov_chrom) <- c("sample_name","chrom","chrom_num","median_scaled_cov")
write.csv(cov_chrom, "chromosomal_scaled_coverage.csv", quote = F, row.names = F)

#***********************************
gc_stats <- read.csv("gc_stats.csv", header = F, as.is = T)
names(gc_stats) <- c("sample_name","gcf","gc_percentile","depth_10","depth_25","depth_50","depth_75","depth_90")
sample_means <- with(gc_stats, aggregate(depth_50, list(sample_name), mean))
names(sample_means) <- c("sample_name", "sample_mean")
sample_gcf40_median <- with(gc_stats[gc_stats$gcf == 40.0,], aggregate(depth_50, list(sample_name), mean))
names(sample_gcf40_median) <- c("sample_name", "sample_mean")
gc_stats <- merge(gc_stats, sample_gcf40_median)
gc_stats$scaled_depth50 <- with(gc_stats, depth_50 / sample_mean)
#***********************************

samples <- unique(gc_stats$sample_name)
pdf("gc_bias_plots.pdf")
for (samp in samples) {
  gci <- gc_stats[grep(samp, gc_stats$sample_name),]
  gci <- gci[order(gci$gcf),]
  with(gci, plot(gc_percentile, scaled_depth50, xlab = "GC percentile", ylab = "scaled median depth", type = "l",
                 ylim = c(0,3), main = samp))
  abline(h = seq(0,3,0.5), lty = "dotted", col = "grey")
}
dev.off()

dups <- read.csv("duplicate_stats.csv", header = F, as.is = T)
names(dups) <- c("sample_name","total_reads","unique_reads")
mismatches <- read.csv("mismatch_stats.csv", header = F, as.is = T)
names(mismatches) <- c("sample_name","mismatch_rate")
wga_stats <- merge(dups, mismatches)
wga_stats$dup_rate <- with(wga_stats, 1 - unique_reads/total_reads)

gc_stats$start <- 0
gc_stats$end <- 100
for (i in 1:nrow(gc_stats)) {
  if (i > 1 && gc_stats$sample_name[i] == gc_stats$sample_name[i-1]) {
    gc_stats$start[i] <- (gc_stats$gc_percentile[i-1] + gc_stats$gc_percentile[i])/2
  }
  if (i < nrow(gc_stats) && gc_stats$sample_name[i] == gc_stats$sample_name[i+1]) {
    gc_stats$end[i] <-  (gc_stats$gc_percentile[i+1] + gc_stats$gc_percentile[i])/2
  }
}

gc_stats$diff <- gc_stats$end - gc_stats$start

# now compute fraction less than 0.5 in each sample
gc_cov_lt_half <- with(gc_stats[gc_stats$scaled_depth50 < 0.5,], aggregate(diff, list(sample_name), function (x) sum(x)/100))
names(gc_cov_lt_half) <- c("sample_name","fraction_lt_p5")
wga_stats <- merge(wga_stats, gc_cov_lt_half)
write.csv(wga_stats, "wga_stats.csv", row.names = F, quote = F)

# get  coverage histograms
cov <- read.csv("genome_mb_coverage.csv", header=F, as.is = T)
names(cov) <- c("coverage","chrom","pos_Mb","sample_name")
sample_names <- sort(unique(cov$sample_name))
sample_names <- sample_names[!(sample_names == "Undetermined_S0")]

pdf("genome_coverage_histograms.pdf")
for (samp in sample_names[!(samp %in% normalization_samples)]) {
  d <- cov[cov$sample_name == samp,]
  d$locus <- 1:nrow(d)
  sample_median <- median(d$coverage)
  with(d, plot(locus, coverage, cex=0, col = "darkred", main = samp, xlab = "position (Mb)", ylim = c(0,3*sample_median)))
  chroms <- unique(d$chrom)
  for (i in 1:length(chroms)) {
    dc <- d[d$chrom == chroms[i],]
    color <- ifelse(i %% 2 == 0, "darkblue","darkred")
    with(dc, lines(locus, coverage, col = color, ylim = c(0,4*sample_median)))
  }
  abline(h = seq(0,5000,200), lty = "dotted", col = "grey")
  abline(h = sample_median, lwd = 2, lty = "dashed")
}
dev.off()

# get a pivot of the coverage data for easy analysis downstream
cov_mb <- read.csv("genome_mb_coverage.csv", header = F, as.is = T)
names(cov_mb) <- c("cov","chrom","pos_mb","sample_name")
samples <- unique(cov_mb$sample_name)

cov_mb_pivot <- unique(cov_mb[,2:3])
for (samp in samples) {
  cov_samp <- cov_mb[cov_mb$sample_name == samp,c("chrom","pos_mb","cov")]
  names(cov_samp)[3] <- samp
  cov_mb_pivot <- merge(cov_mb_pivot, cov_samp, all.x = T)
}
cov_mb_pivot[is.na(cov_mb_pivot)] <- 0
write.csv(cov_mb_pivot, "cov_mb_pivot.csv", quote = F, row.names = F)
