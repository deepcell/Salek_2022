# This script is for the blood spike-ins analysis used in ExtendedFigure7
# Custom genotype A549 (hypotriploid genome). To estimate G (genotypes), N (effective GEs), e (assay error) in a MLE sense

setwd("./")
saf <- read.csv("snp_allele_fractions.csv", check.names = F, as.is = T)
safi <- saf[saf$sample_name == "B-A549" & saf$assay_qc == "PASS" & saf$total > 50,]

allowed_genotypes <- c(0, 1/3, 1/2, 2/3, 1)
j <- 0
N_prev <- 10

# initial estimates of N and e
N <- 100
e <- 0.01
max_iter <- 1000
N_estimates <- c(N)

# iteratively estimate N and the genotypes
iter <- 0
while (N_prev != N & iter < max_iter) {
  # get likelihoods for all allowed genotypes
  iter <- iter + 1
  likelihoods <- c()
  for (genotype in allowed_genotypes) {
    if (genotype == 0) {
      l <- dbinom(round(safi$na / safi$total * N, 0), N, e)
    } else if (genotype == 1) {
      l <- dbinom(round(safi$nr / safi$total * N, 0), N, e)
    } else {
      l <- dbinom(as.integer(round(safi$af * N, 0)), N, genotype)
    }
    likelihoods <- cbind(likelihoods, l)
  }
  max_likelihood <- apply(likelihoods, 1, max)
  
  for (i in 1:nrow(safi)) {
    safi$G[i] <- allowed_genotypes[which(likelihoods[i,] == max_likelihood[i])]
  }
  
  # re-estimate e and N
  e <- mean(c(safi$af[safi$G == 0], 1 - safi$af[safi$G == 1]))
  N_prev <- N
  N <- round(1/with(safi[safi$G*(1-safi$G) > 0,], mean((af - G)^2 / (G*(1-G)))),0)
  N_estimates <- c(N_estimates, N)
}

###########################################################################
# Estimate mixture proportions using an expectation maximization algorithm

mixture_sample_1 <- "S-AJ-81820-T1"
mixture_sample_2 <- "S-AJ-81820-A1"
safm1 <- saf[saf$sample_name == mixture_sample_1 & saf$assay_qc == "PASS" & saf$total > 50,c("CHROM","POS","af")]
names(safm1)[names(safm1) == "af"] <- mixture_sample_1
safm2 <- saf[saf$sample_name == mixture_sample_2 & saf$assay_qc == "PASS" & saf$total > 50,c("CHROM","POS","af")]
names(safm2)[names(safm2) == "af"] <- mixture_sample_2
safp <- merge(merge(safi[,c("CHROM","POS","af","G")], safm1), safm2)
names(safp)[names(safp) == "G"] <- "G_A549"
names(safp)[names(safp) == "af"] <- "af_A549"
max_iter <- 1000

for (mix_samp in c(mixture_sample_1, mixture_sample_2)) {
  # Set all A549 genotypes to 0.5 as an initial guess
  safp$G <- 0.5
  fprev <- iter <- 0

  # estimate f - the fraction of the *known* component, with a linear regression
  f <- lm(safp[,mix_samp] ~ safp$G + safp$G_A549 - 1)$coefficients["safp$G_A549"]
  f_iter <- c(f)

  # initial estimate of error e from homozygous clusters at either corner. This will be better estimated iteratively in EM
  e <- mean(c(safp[safp$G_A549 == 0 & safp[,mix_samp] < 0.05,mix_samp], 1 - safp[safp$G_A549 == 1 & safp[,mix_samp] > 0.95, mix_samp]))
  
  # EM loop: 
  # iterate f -> estimate N from var(allele fraction | G, f) 
  # -> Likelihoods of each allowed genotype of unknown component
  # -> estimate G from max likelihood at current f
  # -> estimate f from linear regression 
  
  while (abs(fprev - f) > 1e-4 & iter < max_iter) {
    iter <- iter + 1
    
    # find expected allele fraction at current estimate of f for all three genotypes
    e0 <- safp$G_A549 * f
    e1 <- safp$G_A549 * f + (1-f)/2
    e2 <- safp$G_A549 * f + (1-f)
    
    # compute scaled squared deviation from expected allele fraction
    d0 <- (safp[,mix_samp] - e0)^2/ (e0*(1-e0))
    d1 <- (safp[,mix_samp] - e1)^2/ (e1*(1-e1))
    d2 <- (safp[,mix_samp] - e2)^2/ (e2*(1-e2))
    
    # estimate N from variance, restricting t
    N  <- round(1/ mean(c(d1[safp$G == 0.5 & safp$G_A549 %in% c(0,1)], d2[safp$G == 1 & safp$G_A549 == 0], d0[safp$G == 0 & safp$G_A549 == 1])),0)
    
    # compute binomial likelihoods, considering error when expected proportion is 0.
    l0 <- ifelse(e0 == 0, dbinom(round(safp[,mix_samp] * N,0), N, e), dbinom(round(safp[,mix_samp] * N,0), N, e0))
    l1 <- dbinom(round(safp[,mix_samp] * N,0), N, e1)
    l2 <- ifelse(e2 == 1, dbinom(N - round(safp[,mix_samp] * N,0), N, e), dbinom(round(safp[,mix_samp] * N,0), N, e2))
    safp$G <- ifelse(l0 > l1, ifelse(l0 > l2, 0, 1), ifelse(l1 > l2, 0.5, 1))
    e <- mean(c(safp[safp$G_A549 == 0 & safp$G == 0 ,mix_samp], 1 - safp[safp$G_A549 == 1 & safp$G == 1, mix_samp]))
    
    fprev <- f
    f <- lm(safp[,mix_samp] ~ safp$G + safp$G_A549 - 1)$coefficients["safp$G_A549"]
    f_iter <- c(f_iter, f)
  }
  
  # make plot for the best set of genotypes and f
  jpeg(paste0(mix_samp,"_redone.jpg"))
  # Uncomment below to plot in svg
  # svg(paste0(mix_samp,".svg"))
  plot(safp$G_A549, safp[,mix_samp], 
       cex = 0, xlab = "genotype allele fraction in A549", ylab = "allele fraction in sorted cells", 
       main = "")
  abline(h = seq(0,1,0.2),lty = "dotted", col = "grey")
  abline(v = seq(0,1,0.2),lty = "dotted", col = "grey")
  lines(c(0,1),c(0, f), lty = "dashed")
  lines(c(0,1),c(1 - f, 1), lty = "dashed")
  lines(c(0,1),c((1-f)/2, (1+f)/2), lty = "dashed")
  points(safp[safp$G == 0, "G_A549"], safp[safp$G == 0, mix_samp], col = alpha("darkred",0.6), pch = 19)
  points(safp[safp$G == 0.5, "G_A549"], safp[safp$G == 0.5, mix_samp], col = alpha("blue",0.6), pch = 18)
  points(safp[safp$G == 1, "G_A549"], safp[safp$G == 1, mix_samp], col = alpha("darkgreen",0.6), pch = 17)
  dev.off()
  
  jpeg(paste0(mix_samp,"_raw_redone.jpg"))
  # Uncomment below to plot in svg
  # svg(paste0(mix_samp,"_raw.svg"))
  plot(safp$af_A549, safp[,mix_samp], 
       cex = 0, xlab = "allele fraction in A549", ylab = "allele fraction in sorted cells", 
       main = "")
  abline(h = seq(0,1,0.2),lty = "dotted", col = "grey")
  abline(v = seq(0,1,0.2),lty = "dotted", col = "grey")
  lines(c(0,1),c(0, f), lty = "dashed")
  lines(c(0,1),c(1 - f, 1), lty = "dashed")
  lines(c(0,1),c((1-f)/2, (1+f)/2), lty = "dashed")
  points(safp[safp$G == 0, "af_A549"], safp[safp$G == 0, mix_samp], col = alpha("darkred",0.6), pch = 19)
  points(safp[safp$G == 0.5, "af_A549"], safp[safp$G == 0.5, mix_samp], col = alpha("blue",0.6), pch = 18)
  points(safp[safp$G == 1, "af_A549"], safp[safp$G == 1, mix_samp], col = alpha("darkgreen",0.6), pch = 17)
  dev.off()
  names(safp)[names(safp) == "G"] <- paste0("G_", mix_samp)
}

write.csv(safp, "blood_spike_ins.csv", row.names = F)
