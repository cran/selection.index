## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(selection.index)

# Load the built-in maize phenotype dataset
data("maize_pheno")

# Extract the traits of interest
traits <- c("Yield", "PlantHeight", "DaysToMaturity")

# Calculate Genotypic (gmat) and Phenotypic (pmat) covariance matrices
gmat <- gen_varcov(maize_pheno[, traits], maize_pheno$Genotype, maize_pheno$Block)
pmat <- phen_varcov(maize_pheno[, traits], maize_pheno$Genotype, maize_pheno$Block)

# Define economic weights for the three traits
wmat <- weight_mat(data.frame(Trait = traits, Weight = c(1, -1, 1)))

## ----rlpsi_example------------------------------------------------------------
# Restrict trait 2 (PHT) to have ZERO expected genetic gain
rlpsi_res <- rlpsi(
  pmat = pmat,
  gmat = gmat,
  wmat = wmat,
  restricted_traits = c(2)
)

# View the summary and coefficients
print(rlpsi_res$summary)

# View the expected genetic gains (Delta_G)
print(rlpsi_res$Delta_G)

## ----ppg_lpsi_example---------------------------------------------------------
# Specify the desired proportions
k_proportions <- c(2, 1, 1)

# Calculate the PPG-LPSI
ppg_res <- ppg_lpsi(pmat = pmat, gmat = gmat, k = k_proportions, wmat = wmat)

# View the expected genetic gains
print(ppg_res$Delta_G)

## ----dg_lpsi_example----------------------------------------------------------
# Explicit vector of desired absolute genetic gains
desired_gains <- c(5, -2, 1)

# Calculate DG-LPSI
dg_res <- dg_lpsi(pmat = pmat, gmat = gmat, d = desired_gains)

# Check the achieved proportional genetic gains
print(dg_res$Delta_G)

# The DG-LPSI also calculates implied Smith-Hazel economic weights
print(dg_res$implied_weights_normalized)

## ----manual_c_matrix----------------------------------------------------------
# Manually restrict traits #1 and #3:
# Create a 3x2 constraint matrix
C_matrix <- matrix(
  c(
    1, 0, 0, # Restrict trait 1
    0, 0, 1
  ), # Restrict trait 3
  nrow = 3, ncol = 2
)

# Pass directly to RLPSI
rlpsi_manual <- rlpsi(pmat = pmat, gmat = gmat, wmat = wmat, C = C_matrix)
print(rlpsi_manual$Delta_G)

