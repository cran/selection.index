## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup_data---------------------------------------------------------------
library(selection.index)
library(MASS) # For robust ridge regression implementation (lm.ridge)

# 1. Load the built-in synthetic maize datasets
data("maize_pheno")
data("maize_geno")

# Ensure markers are numeric matrices
X <- as.matrix(maize_geno)

# Ensure pheno lines up perfectly with geno
# (Note: maize_pheno contains 4 repetitions per genotype; we take the means)
Y_agg <- mean_performance(
  data = maize_pheno[, c("Yield", "PlantHeight", "DaysToMaturity")],
  genotypes = maize_pheno$Genotype,
  replications = maize_pheno$Block
)

# Sort both datasets to ensure identical ordering
Y_agg <- Y_agg[order(Y_agg$Genotypes), ]
X <- X[order(rownames(X)), ]

# Ensure overlapping genotypes
common_genotypes <- intersect(Y_agg$Genotypes, rownames(X))
Y_agg <- Y_agg[Y_agg$Genotypes %in% common_genotypes, ]
X <- X[rownames(X) %in% common_genotypes, ]

# Select only multi-traits relevant to breeding
# Yield, PlantHeight, DaysToMaturity
Y <- as.matrix(Y_agg[, c("Yield", "PlantHeight", "DaysToMaturity")])

## ----calculate_gebv-----------------------------------------------------------
# Simulate Genomic Estimated Breeding Values (GEBVs) using Ridge Regression
# In best practice, you'd use cross-validation or a separate training/testing set
# We use lambda = 100 to handle the p >> n dimensionality problem
gebv_mat <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
colnames(gebv_mat) <- colnames(Y)

for (i in seq_len(ncol(Y))) {
  # Fit a ridge regression model for trait `i` using markers `X`
  model_ridge <- lm.ridge(Y[, i] ~ X, lambda = 100)

  # Predict values using coef() which correctly un-scales the coefficients: Include intercept
  betas <- coef(model_ridge)
  intercept <- betas[1]
  beta <- betas[-1] # Exclude intercept

  # Calculate the GEBV for the trait
  gebv_mat[, i] <- intercept + (X %*% beta)
}

head(gebv_mat, 3)

## ----calculated_lgsi----------------------------------------------------------
# 2. Compute Trait Covariances
pmat <- cov(Y) # Phenotypic Covariance (Approximation)

# In optimal practice, the true Genomic Covariance Matrix (\Gamma) is
# estimated using Restricted Maximum Likelihood (REML) utilizing both general
# phenotypic and genotypic information. Here, we simulate \Gamma as proportional
# to the phenotypic variance (assuming a high heritability correlation).
gmat <- pmat * 0.4 # Approximate Genotypic Covariance (Approximating \Gamma)

# 3. Define Economic Weights
weights <- data.frame(
  Trait = c("Yield", "PlantHeight", "DaysToMaturity"),
  Weight = c(5, -0.1, -0.1)
)
wmat <- weight_mat(weights)

# 4. Calculate Linear Genomic Selection Index (LGSI)
# For the testing population where we only use genomic values
lgsi_result <- lgsi(
  gebv_mat = gebv_mat,
  gmat = gmat,
  wmat = wmat
)

# Output Summary
lgsi_result$summary

## ----calculate_clgsi----------------------------------------------------------
clgsi_result <- clgsi(
  phen_mat = Y, # Observed phenotypic data
  gebv_mat = gebv_mat, # Genomic Estimated Breeding Values
  pmat = pmat, # Expected Phenotypic traits covariance
  gmat = gmat, # Expected Genotypic traits covariance
  wmat = wmat
)

clgsi_result$summary

