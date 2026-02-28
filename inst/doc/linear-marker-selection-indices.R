## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(selection.index)

# Load the built-in maize phenotype and genotype datasets
data("maize_pheno")
data("maize_geno")

# 1. Prepare Phenotypic Data
# We select three traits for our demonstration
traits <- c("Yield", "PlantHeight", "DaysToMaturity")
phen_mat <- as.matrix(maize_pheno[, traits])

# Calculate Genotypic (gmat) and Phenotypic (pmat) covariance matrices
gmat <- gen_varcov(maize_pheno[, traits], maize_pheno$Genotype, maize_pheno$Block)
pmat <- phen_varcov(maize_pheno[, traits], maize_pheno$Genotype, maize_pheno$Block)

# Aggregate phenotypic data by Genotype to match marker data dimensions
agg_pheno <- aggregate(maize_pheno[, traits], by = list(Genotype = maize_pheno$Genotype), FUN = mean)
phen_mat <- as.matrix(agg_pheno[, traits])

# 2. Prepare Genotypic Data
# Extract only the marker columns (dropping the Genotype ID column)
marker_mat <- as.matrix(maize_geno[, -1])

# 3. Define Economic Weights
# Suppose we want to improve Yield strongly, decrease Plant Height, and lightly increase Days to Maturity
wmat <- weight_mat(data.frame(Trait = traits, Weight = c(2, -1, 0.5)))

## ----calculate_marker_scores--------------------------------------------------
# Calculate marker scores via a simple Ridge Regression
# (We use MASS::lm.ridge or glmnet in practice, but for the vignette we approximate
# the prediction using the marker matrix directly)
library(MASS)

marker_scores <- matrix(0, nrow = nrow(phen_mat), ncol = ncol(phen_mat))
colnames(marker_scores) <- colnames(phen_mat)

for (i in seq_len(ncol(phen_mat))) {
  # Fit ridge regression to handle multicollinearity among markers
  # Note: A lambda must be chosen carefully in real scenarios
  fit <- lm.ridge(phen_mat[, i] ~ marker_mat, lambda = 10)
  # Calculate predicted marker scores
  marker_scores[, i] <- scale(marker_mat, center = fit$xm, scale = fit$scales) %*% fit$coef + fit$ym
}

## ----lmsi_example-------------------------------------------------------------
# Calculate the LMSI
lmsi_res <- lmsi(
  phen_mat = phen_mat,
  marker_scores = marker_scores,
  pmat = pmat,
  gmat = gmat,
  wmat = wmat
)

# View the LMSI coefficient summary
print(lmsi_res$summary)

# Display the expected genetic gains per trait
print(lmsi_res$Delta_G)

## ----gw_lmsi_example----------------------------------------------------------
# We use all 500 markers in the maize_geno matrix
# Since p (500) approaches n (600), we supply a ridge regularization lambda
# to ensure the covariance matrices are invertible.
gw_lmsi_res <- gw_lmsi(
  marker_mat = marker_mat,
  trait_mat = phen_mat,
  gmat = gmat,
  wmat = wmat,
  lambda = 0.05
)

# Display the phenotypic weighting coefficients
print(gw_lmsi_res$b_y)

# We can also check the condition number to observe the matrix stability
print(gw_lmsi_res$condition_number)

# And summarize the overall index statistics
print(gw_lmsi_res$summary)

