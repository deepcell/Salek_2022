library(ggplot2)

saf <- read.csv("snp_allele_fractions.csv", check.names = F, as.is = T)
saf$input <- sub("^.*-","",saf$sample_name)
saf$fraction <- as.numeric(sub("-.*","", saf$sample_name))/100
saf <- saf[saf$sample_qc == "PASS" & saf$assay_qc == "PASS" & saf$input == "250pg",]
saf$af[saf$total < 50] <- NA
saf$role <- "mixture"
saf$role[saf$sample_name %in% unique(saf$sample_name[saf$fraction == 0])] <- "base"
saf$role[saf$sample_name %in% unique(saf$sample_name[saf$fraction == 1])] <- "additive"
sample_means <- with(saf, aggregate(total, list(sample_name), mean))
names(sample_means) <- c("sample_name","sample_mean")
saf <- merge(saf, sample_means)
saf$scaled_total <- with(saf, total/sample_mean)


md <- c()
sample_names <- unique(saf$sample_name)
for (samp in sample_names) {
  role <- saf$role[saf$sample_name == samp][1]
  di <- saf[saf$sample_name == samp,c("CHROM","POS","scaled_total","af")]
  names(di)[3] <- paste("total", samp, role, sep = "_")
  names(di)[4] <- paste("af", samp, role, sep = "_")
  if (is.null(md)) {
    md <- di
  } else {
    md <- merge(md, di)
  }
}

# genotype two components
md$af_base <- apply(as.matrix(md[,names(md)[grep("af.*base",names(md))]]), 1, mean, na.rm = T)
md$total_base <- apply(as.matrix(md[,names(md)[grep("total.*base",names(md))]]), 1, mean, na.rm = T)
md$af_additive <- apply(as.matrix(md[,names(md)[grep("af.*additive",names(md))]]), 1, mean, na.rm = T)
md$total_additive <- apply(as.matrix(md[,names(md)[grep("total.*additive",names(md))]]), 1, mean, na.rm = T)
md$genmaj <- with(md, ifelse(af_base < 0.05, 0, ifelse(af_base > 0.95, 1, 0.5)))
#md$genadd <- with(md, ifelse(af_additive < 0.05, 0, ifelse(af_additive > 0.95, 1, 0.5)))
md$genadd <- with(md, ifelse(af_additive < 0.05, 0, ifelse(af_additive > 0.95, 1, af_additive)))
mixture_samples <- names(md)[grep("af.*mixture",names(md))]

af_estimate <- c()
for (samp in mixture_samples) {
  di <- md[,c(samp, "af_base", "total_base", "af_additive","total_additive","genmaj","genadd")]
  names(di)[1] <- "af_mixture"
  mean_af_het <- mean(c(di$af_mixture[di$genmaj == 0.5 & di$genadd == 0], 1 - di$af_mixture[di$genmaj == 0.5 & di$genadd == 1]), na.rm = T)
  max_var     <-  var(c(di$af_mixture[di$genmaj == 0.5 & di$genadd == 0], 1 - di$af_mixture[di$genmaj == 0.5 & di$genadd == 1]), na.rm = T)
  N <- round(mean_af_het * (1 - mean_af_het) / max_var,0)
  di_inf <- di[!(di$genmaj - md$genadd == 0),]

  sumlogp_vec <- c()
  fraction <- seq(0.001,0.999,0.001)
  ml_frac <- ml_ll <- NA
  for (frac in fraction) {
    af_expected <- with(di_inf, (total_base * genmaj * (1-frac) + total_additive * genadd * frac)/(total_base * (1-frac) + total_additive * frac) )
    logp <- dbinom(round(di_inf$af_mixture * N, 0), N, af_expected, log = T)
    sumlogp <- sum(logp)
    sumlogp_vec <- c(sumlogp_vec, sumlogp)

    if (is.na(ml_frac)) {
      ml_frac <- frac
      ml_ll <- sumlogp
    } else {
      if (sumlogp > ml_ll) {
        ml_frac <- frac
        ml_ll <- sumlogp
      }
    }
  }
  af_estimate <- c(af_estimate, ml_frac)
}

mixture_sample <- mixture_samples
mixture_estimates <- data.frame(mixture_sample, af_estimate, stringsAsFactors = F)
mixture_estimates$fraction <- as.numeric(sub("-.*", "", sub("af_","",mixture_estimates$mixture_sample)))/100
ggplot(mixture_estimates, aes(x=fraction, y=af_estimate)) + geom_point(size = 2, shape = 23) + geom_abline(slope = 1, intercept = 0, linetype = "dashed") + labs(x = "mixture fraction", y = "estimated fraction")
