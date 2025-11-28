
# R version 3.5.0
setwd("~/path_to_repo/code/")
set.seed(100)

library(openxlsx) # v4.1.4
library(plyr) # v1.8.4
library(InPosition) # v0.12.7.1
library(TExPosition) # v2.6.10.1
library(R.matlab) # v3.6.2
library(RColorBrewer) # v1.1-2

# Prepare data
fluxes <- read.xlsx("../data/results/end_pFBA_fluxes.xlsx")
DoE <- read.xlsx("../data/raw/DoE data.xlsx")
DoE$BR <- NULL
DoE$IPTG <- NULL
DoE <- DoE[c(2:5, 7:24),]
quant.fluxes <- as.matrix(fluxes)
class(quant.fluxes) <- "numeric"
num_samples <- dim(quant.fluxes)[1]
num_fluxes <- dim(quant.fluxes)[2]
cat.DoE <- makeNominalData(DoE)
num_DoEs <- dim(cat.DoE)[2]



## PLSCA on permuted data
num_permutations <- 999
lambdas <- matrix(NA, num_permutations, 12)
for (i in c(1:num_permutations)) {
  permuted_fluxes <- quant.fluxes
  for (j in c(1:dim(quant.fluxes)[2])) { # for each data dimension
    permutation_idx <- sample.int(dim(quant.fluxes)[1], replace = FALSE)
    permuted_fluxes[, j] <- quant.fluxes[permutation_idx, j] # shuffle the data values of each axis
  }
  # Rescale fluxes
  center.normed.fluxes <- expo.scale(permuted_fluxes)
  escofier.transformed.fluxes <- cbind((1-center.normed.fluxes)/2, (1+center.normed.fluxes)/2)
  res <- tepPLSCA(escofier.transformed.fluxes, cat.DoE, graphs = FALSE)
  lambdas[i, ] <- res$TExPosition.Data$eigs
}

## PLSCA on original data
# Rescale fluxes
center.normed.fluxes <- expo.scale(quant.fluxes)
escofier.transformed.fluxes <- cbind((1-center.normed.fluxes)/2, (1+center.normed.fluxes)/2)
# Mixed-modality simple PLSCA w/ strictly quantitative SNPs.
PLSCA.res <- tepPLSCA(escofier.transformed.fluxes, cat.DoE, graphs = FALSE)

## Permutation test to evaluate PLSCA components
p_lambdas <- (colSums(lambdas >= t(replicate(num_permutations, PLSCA.res$TExPosition.Data$eigs))) + 1) / (num_permutations + 1);



# Define function to get data mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
# Select fluxes with enough variability
idx <- matrix(TRUE, num_fluxes, 1)
for (i in c(1:num_fluxes)) {
  if (sum(quant.fluxes[, i] != getmode(quant.fluxes[, i])) == 1) {
    idx[i] <- FALSE
  }
}
quant.fluxes.1 <- quant.fluxes[, idx]
# Select DoE parameters with enough variability
idx <- matrix(TRUE, num_DoEs, 1)
for (i in c(1:num_DoEs)) {
  if (sum(cat.DoE[, i]) <= 1) {
    idx[i] <- FALSE
  }
}
cat.DoE.1 <- cat.DoE[, idx]

## Jackknife sampling
bootstrap_flux_scores <- array(0, c(2*dim(quant.fluxes.1)[2], 2, num_samples))
bootstrap_DoE_scores <- array(0, c(dim(cat.DoE.1)[2], 2, num_samples))
for (i in c(1:num_samples)) {
  bootstrap_idx <- array(TRUE, c(num_samples, 1))
  bootstrap_idx[i] <- FALSE
  bootstrap_fluxes <- quant.fluxes.1[bootstrap_idx, ]
  bootstrap_DoE <- cat.DoE.1[bootstrap_idx, ]
  # Rescale fluxes
  center.normed.fluxes <- expo.scale(bootstrap_fluxes)
  escofier.transformed.fluxes <- cbind((1-center.normed.fluxes)/2, (1+center.normed.fluxes)/2)
  res <- tepPLSCA(escofier.transformed.fluxes, bootstrap_DoE, graphs = FALSE)
  bootstrap_flux_scores[, , i] <- res$TExPosition.Data$fi[, 1:2]
  bootstrap_DoE_scores[, , i] <- res$TExPosition.Data$fj[, 1:2]
}

## Bootstrap ratio test
BSR.res.fluxes <- boot.ratio.test(bootstrap_flux_scores, critical.value = 2)
BSR.res.DoE <- boot.ratio.test(bootstrap_DoE_scores, critical.value = 2)



## Save all
save.image("../data/results/PLSCA_results.RData")

## Export data frames for MATLAB
writeMat('../data/results/PLSCA_results.mat', fluxes_factor_scores=PLSCA.res$TExPosition.Data$fi, 
         fluxes_contributions=PLSCA.res$TExPosition.Data$ci,
         DoE_factor_scores=PLSCA.res$TExPosition.Data$fj, 
         DoE_contributions=PLSCA.res$TExPosition.Data$cj,
         explained=PLSCA.res$TExPosition.Data$t,
         p_lambdas=p_lambdas,
         fluxes_BSRs=BSR.res.fluxes$boot.ratios,
         DoE_BSRs=BSR.res.DoE$boot.ratios,
         rxns=colnames(quant.fluxes),
         DoE_variables=colnames(cat.DoE),
         bootstrap_rxns=colnames(quant.fluxes.1),
         bootstrap_DoE_variables=colnames(cat.DoE.1))

