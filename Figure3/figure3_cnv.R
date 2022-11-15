library(scales)
setwd("./")

# read in DTC pre-sorted and sorted data
cov_mb1 <- read.csv("7_DTC_1mb.csv", header = F, as.is = T)
names(cov_mb1) <- c("cov","chrom","pos_mb","sample_name")
samples1 <- unique(cov_mb1$sample_name)

# also read in 12878 data
cov2_resolve <- read.csv("7_GM12878_1mb.csv", header = F, as.is = T)
names(cov2_resolve) <- c("cov","chrom","pos_mb","sample_name")
samples2_resolve <- unique(cov2_resolve$sample_name)

normalization_samples_resolve <- samples2_resolve <- samples2_resolve[grep("^12878-[0-9]+cell-Reso", samples2_resolve)]

cov_mb <- rbind(cov_mb1, cov2_resolve[cov2_resolve$sample_name %in% normalization_samples_resolve,])
samples <- unique(cov_mb$sample_name[cov_mb$sample_name != "Undetermined_S0"])

# for each sample, determine the internal scaled cov_mb_scaled relative to mean
# then determine overall bias dist
# then scale each coverage by bias (or adjust by PCs for overkill)

cov_mb$sample_mean = 0
for (samp in samples) {
  cov_mb$sample_mean[cov_mb$sample_name == samp] = mean(cov_mb$cov[cov_mb$sample_name == samp])
}
cov_mb$scaled_cov <- with(cov_mb, cov/sample_mean)

cov_bias_median_resolve <- with(cov_mb[cov_mb$sample_name %in% normalization_samples_resolve,], aggregate(scaled_cov, list(chrom, pos_mb), median))
names(cov_bias_median_resolve) <- c("chrom","pos_mb","bias_median_resolve")
cov_bias_mad_resolve  <- with(cov_mb[cov_mb$sample_name %in% normalization_samples_resolve,], aggregate(scaled_cov, list(chrom, pos_mb), mad))
names(cov_bias_mad_resolve) <- c("chrom","pos_mb","bias_mad_resolve")
cov_bias_resolve <- merge(cov_bias_median_resolve, cov_bias_mad_resolve)

cov_bias_resolve <- cov_bias_resolve[order(cov_bias_resolve$chrom, cov_bias_resolve$pos_mb),]
cov_bias_resolve$bias_cv_resolve <- with(cov_bias_resolve, bias_mad_resolve/bias_median_resolve)

cov_mb <- merge(cov_mb, cov_bias_resolve)

cov_mb <- cov_mb[order(cov_mb$sample_name, cov_mb$chrom, cov_mb$pos_mb),]
cov_mb$loc <- 1
for (i in 2:nrow(cov_mb)) {
  cov_mb$loc[i] <- ifelse(cov_mb$sample_name[i] == cov_mb$sample_name[i-1], cov_mb$loc[i-1] + 1, 1)
}
cov_mb <- cov_mb[cov_mb$bias_cv < 0.2,]

# restrict to autosomes + X for now. Don't have data for Y in particular
cov_mb <- cov_mb[!(cov_mb$chrom %in% c("chrEBV","chrM","chrY" )),]
cov_mb$chrom_num <- -1
cov_mb$chrom_num <- ifelse(cov_mb$chrom == "chrX", 23, as.numeric(sub("chr","",cov_mb$chrom)))
cov_mb <- cov_mb[order(cov_mb$sample_name, cov_mb$chrom_num, cov_mb$pos_mb),]
cov_mb$loc <- 1
for (i in 2:nrow(cov_mb)) {
  cov_mb$loc[i] <- ifelse(cov_mb$sample_name[i] == cov_mb$sample_name[i-1], cov_mb$loc[i-1] + 1, 1)
}

# first show 1 Mb-level data for 12878 samples
pdf("normalized_coverage_data_12878.pdf")
with(cov_mb[cov_mb$sample_name %in% normalization_samples_resolve,], 
     plot(loc, scaled_cov/bias_median_resolve, pch = 19, col = alpha("darkblue",0.3), cex = 0,
          main = "GM12878 normalized coverage data\nResolvex", ylim = c(0,3), xlab = "genomic position (Mb)", ylab = "scaled coverage"))
for (chrom_num in 1:23) {
  with(cov_mb[cov_mb$sample_name %in% normalization_samples_resolve & cov_mb$chrom_num == chrom_num,],
       points(loc, scaled_cov/bias_median_resolve, pch = 19, cex = 0.5, col = alpha(ifelse(chrom_num %% 2 == 0, "darkblue","darkred"),0.3)))
}
abline(h = seq(0,3,0.5), lty = "dotted", col = "grey")

dev.off()


# now build all of the samples, one by one. Hope CNVs show up!

pdf("DTC_sample_cnv_profiles.pdf")
for (samp in samples) {
  with(cov_mb[cov_mb$sample_name == samp,], 
       plot(loc, scaled_cov/bias_median_resolve, pch = 19, col = alpha("darkblue",0.3), cex = 0,
            main = samp, ylim = c(0,3), xlab = "genomic position (Mb)", ylab = "scaled coverage"))
  for (chrom_num in 1:23) {
    with(cov_mb[cov_mb$sample_name == samp & cov_mb$chrom_num == chrom_num,],
         points(loc, scaled_cov/bias_median_resolve, pch = 19, cex = 0.5, col = alpha(ifelse(chrom_num %% 2 == 0, "darkblue","darkred"),0.3)))
  }
  abline(h = seq(0,3,0.5), lty = "dotted", col = "grey")
}
dev.off()
