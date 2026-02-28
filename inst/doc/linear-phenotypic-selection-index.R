## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----matrices-----------------------------------------------------------------
library(selection.index)

# Load the synthetic phenotypic multi-environment dataset
data("maize_pheno")

# In maize_pheno: Traits are columns 4:6.
# Genotypes are in column 1, and Block/Replication is in column 3.
gmat <- gen_varcov(data = maize_pheno[, 4:6], genotypes = maize_pheno[, 1], replication = maize_pheno[, 3])
pmat <- phen_varcov(data = maize_pheno[, 4:6], genotypes = maize_pheno[, 1], replication = maize_pheno[, 3])

## ----weights------------------------------------------------------------------
# Define the economic weights for the 3 continuous traits
# (e.g., Yield, PlantHeight, DaysToMaturity)
weights <- c(10, -5, -5)

## ----lpsi---------------------------------------------------------------------
# Calculate the Optimal Combinatorial Linear Phenotypic Selection Index (LPSI)
index_results <- lpsi(
  ncomb = 3,
  pmat = pmat,
  gmat = gmat,
  wmat = as.matrix(weights),
  wcol = 1
)

## ----gains--------------------------------------------------------------------
# View the top combinatorial indices, including their selection response (R_A)
head(index_results)

# Extract the phenotypic selection scores to strategically rank the parental candidates
# using the top evaluated combinatorial index
scores <- predict_selection_score(
  index_results,
  data = maize_pheno[, 4:6],
  genotypes = maize_pheno[, 1]
)

# View the top performing candidates designated for the next breeding cycle
head(scores)

## ----marker_data, eval=FALSE--------------------------------------------------
# # Load the associated synthetic genomic dataset (500 SNPs for the 100 genotypes)
# data("maize_geno")
# 
# # Calculate the marker-assisted index combining our matrices and raw SNP profiles
# marker_index_results <- lmsi(
#   pmat = pmat,
#   gmat = gmat,
#   marker_scores = maize_geno,
#   wmat = weights
# )
# 
# summary(marker_index_results)

## ----base_index---------------------------------------------------------------
# Calculate the Base Index and automatically compare its efficiency to the LPSI
base_results <- base_index(
  pmat = pmat,
  gmat = gmat,
  wmat = weights,
  compare_to_lpsi = TRUE
)

# Observe the expected genetic gains and efficiency comparison
base_results$summary

## ----heritability-------------------------------------------------------------
# Extract the top combinatorial index results
top_index <- index_results[1, ]

# h^2_I: Heritability of the optimal index
top_index$hI2

# \rho_HI: Correlation between the LPSI and the true underlying Net Genetic Merit
top_index$rHI

