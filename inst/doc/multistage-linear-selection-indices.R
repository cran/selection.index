## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup_data---------------------------------------------------------------
library(selection.index)

# Estimate phenotypic and genotypic covariance matrices for the 3 traits
# The traits are Yield, PlantHeight, DaysToMaturity
traits <- c("Yield", "PlantHeight", "DaysToMaturity")
pmat <- phen_varcov(maize_pheno[, traits], maize_pheno$Environment, maize_pheno$Genotype)
gmat <- gen_varcov(maize_pheno[, traits], maize_pheno$Environment, maize_pheno$Genotype)

# Matrix limits for Stage 1 (Traits 1 to 2)
P1 <- pmat[1:2, 1:2]
G1 <- gmat[1:2, 1:2]

# Complete Matrices for Stage 2
P <- pmat
C <- gmat

# Economic weights for the 3 traits
weights <- c(10, -5, -2)

## ----mlpsi_example------------------------------------------------------------
# We apply a selection proportion of 10% (0.10) per stage.
mlpsi_res <- mlpsi(
  P1 = P1, P = P, G1 = G1, C = C,
  wmat = weights,
  selection_proportion = 0.1
)

# Stage 1 metrics
mlpsi_res$summary_stage1

# Stage 2 metrics
mlpsi_res$summary_stage2

## ----mrlpsi_example-----------------------------------------------------------
# We constrain PlantHeight (Trait 2) at Stage 1
C1 <- matrix(0, nrow = 2, ncol = 1)
C1[2, 1] <- 1

# We constrain PlantHeight (Trait 2) at Stage 2
C2 <- matrix(0, nrow = 3, ncol = 1)
C2[2, 1] <- 1

mrlpsi_res <- mrlpsi(
  P1 = P1, P = P, G1 = G1, C = C,
  wmat = weights,
  C1 = C1, C2 = C2,
  selection_proportion = 0.1
)

# Observe that Expected Gain (E) for PlantHeight is approximately 0
mrlpsi_res$summary_stage1

## ----mppg_lpsi_example--------------------------------------------------------
# Target specific proportional gains
d1 <- c(2, 1) # Yield gains twice as much as PlantHeight at stage 1
d2 <- c(3, 1, 0.5) # Desired proportions at stage 2

mppg_res <- mppg_lpsi(
  P1 = P1, P = P, G1 = G1, C = C,
  wmat = weights,
  d1 = d1, d2 = d2,
  selection_proportion = 0.1
)

# Observe the Expected Gain (E) in the resulting summary stats aligns with d1 proportions
mppg_res$summary_stage1

## ----setup_genomic------------------------------------------------------------
set.seed(42)
reliability <- 0.7 # Simulated genomic prediction reliability

Gamma1 <- reliability * G1
Gamma <- reliability * C
A1 <- reliability * G1
A <- C[, 1:2] # n x n1 covariance mapping

## ----mlgsi_example------------------------------------------------------------
mlgsi_res <- mlgsi(
  Gamma1 = Gamma1, Gamma = Gamma, A1 = A1, A = A,
  C = C, G1 = G1, P1 = P1,
  wmat = weights,
  selection_proportion = 0.1
)

mlgsi_res$summary_stage1

## ----mrlgsi_example-----------------------------------------------------------
mrlgsi_res <- mrlgsi(
  Gamma1 = Gamma1, Gamma = Gamma, A1 = A1, A = A,
  C = C, G1 = G1, P1 = P1,
  wmat = weights,
  C1 = C1, C2 = C2,
  selection_proportion = 0.1
)

mrlgsi_res$summary_stage2

## ----mppg_lgsi_example--------------------------------------------------------
mppg_lgsi_res <- mppg_lgsi(
  Gamma1 = Gamma1, Gamma = Gamma, A1 = A1, A = A,
  C = C, G1 = G1, P1 = P1,
  wmat = weights,
  d1 = d1, d2 = d2,
  selection_proportion = 0.1
)

mppg_lgsi_res$summary_stage1

