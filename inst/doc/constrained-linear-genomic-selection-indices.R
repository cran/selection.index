## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(selection.index)
library(MASS) # For Ridge regression

# Load the maize datasets
data("maize_pheno", package = "selection.index")
data("maize_geno", package = "selection.index")

# Extract only the traits we care about to prevent missing value logic failures
traits <- c("Yield", "PlantHeight", "DaysToMaturity")
maize_pheno_clean <- na.omit(maize_pheno[, c("Genotype", "Block", traits)])

# Calculate mean performance for phenotypic data
# Providing only the traits to `data` to accurately generate the table
pheno_means_full <- mean_performance(
  data = maize_pheno_clean[, traits],
  genotypes = maize_pheno_clean$Genotype,
  replications = maize_pheno_clean$Block,
  method = "Mean"
)

# Remove the trailing 9 rows composed of summary text-stats automatically generated
pheno_means <- pheno_means_full[1:(nrow(pheno_means_full) - 9), ]

# Extract numerical phenotype matrix
Y <- as.matrix(pheno_means[, traits])
mode(Y) <- "numeric"
rownames(Y) <- pheno_means$Genotypes

# Format genotype matrix
X <- as.matrix(maize_geno)

# Ensure matching dimensions and align the matrices
common_genotypes <- intersect(rownames(Y), rownames(X))

# Filter and sort matrices to match genotypes sequence exactly
Y_filtered <- Y[common_genotypes, ]
X_filtered <- X[common_genotypes, ]

# Calculate GEBVs using a simple Ridge Regression model per trait
gebvs_sim <- matrix(0, nrow = nrow(Y_filtered), ncol = ncol(Y_filtered))
colnames(gebvs_sim) <- colnames(Y_filtered)
rownames(gebvs_sim) <- rownames(Y_filtered)
lambda_ridge <- 0.1

for (j in seq_len(ncol(Y_filtered))) {
  # Fit ridge regression to simulate marker effects
  model_ridge <- lm.ridge(Y_filtered[, j] ~ X_filtered, lambda = lambda_ridge)
  beta_ridge <- coef(model_ridge)[-1]

  # Predict GEBVs
  gebvs_sim[, j] <- X_filtered %*% matrix(beta_ridge, ncol = 1)
}

# Define Covariance Matrices
Gamma <- cov(gebvs_sim) # Genomic Covariance Matrix estimated from simulated GEBVs
P_matrix <- cov(Y_filtered) # Phenotypic Covariance Matrix
C_matrix <- Gamma # Genetic Covariance (Assuming Gamma captures genetic additive cov)

# Define vector of economic weights for Yield, PlantHeight, DaysToMaturity
w <- matrix(c(5, -0.1, -0.1), ncol = 1)

## ----rlgsi_example------------------------------------------------------------
# Define the restriction matrix U
# We restrict traits 2 (PlantHeight) and 3 (DaysToMaturity).
# Trait 1 (Yield) is left unrestricted.
U_rlgsi <- matrix(c(
  0, 0, # Yield unrestricted
  1, 0, # PlantHeight restricted
  0, 1 # DaysToMaturity restricted
), nrow = 3, byrow = TRUE)

# Calculate RLGSI
rlgsi_res <- rlgsi(Gamma = Gamma, wmat = w, U = U_rlgsi)

# View the resulting index coefficients and Expected Genetic Gain
cat("RLGSI Coefficients (beta_RG):\n")
print(rlgsi_res$b)

cat("\nExpected Genetic Gain per Trait:\n")
print(rlgsi_res$Summary[, "Delta_H", drop = FALSE])

## ----ppglgsi_example----------------------------------------------------------
# Define the restriction matrix U for PPG
# Restrict traits 1 (Yield) and 2 (PlantHeight) to predetermined scalars
U_ppg <- matrix(c(
  1, 0, # Yield restricted
  0, 1, # PlantHeight restricted
  0, 0 # DaysToMaturity unrestricted
), nrow = 3, byrow = TRUE)

# Define the predetermined values vector 'd'
d_vec <- matrix(c(7.0, -3.0), ncol = 1)

# Calculate PPG-LGSI
ppg_res <- ppg_lgsi(Gamma = Gamma, d = d_vec, wmat = w, U = U_ppg)

cat("PPG-LGSI Coefficients (beta_PG):\n")
print(ppg_res$b)

## ----crlgsi_example-----------------------------------------------------------
# 1. Provide Combined Matrices
# The combined Phenotypic-Genomic covariance matrix T_C
T_C <- rbind(
  cbind(P_matrix, Gamma),
  cbind(Gamma, Gamma)
)

# The combined Genetic-Genomic covariance matrix Psi_C
Psi_C <- rbind(
  cbind(C_matrix, Gamma),
  cbind(Gamma, Gamma)
)

# 2. Define Restriction Matrix U_C
# Note: For combined indices of `t` traits, U_C spans exactly 2t rows.
# Restricting Trait 1 (Yield) on both empirical and genomic parameters:
U_crlgsi <- matrix(0, nrow = 6, ncol = 2)
U_crlgsi[1, 1] <- 1 # Trait 1 phenotypic component
U_crlgsi[4, 2] <- 1 # Trait 1 genomic component

# 3. Calculate CRLGSI
crlgsi_res <- crlgsi(T_C = T_C, Psi_C = Psi_C, wmat = w, U = U_crlgsi)

cat("CRLGSI Combined Coefficients (beta_CR):\n")
print(crlgsi_res$b)

## ----cppglgsi_example---------------------------------------------------------
# Target Yield using predetermined combined proportions d
d_combined <- matrix(c(7.0, 3.5), ncol = 1)

# Calculate CPPG-LGSI dynamically
cppg_res <- cppg_lgsi(T_C = T_C, Psi_C = Psi_C, d = d_combined, wmat = w, U = U_crlgsi)

cat("CPPG-LGSI Coefficients (beta_CP):\n")
print(cppg_res$b)

## ----summary_comparison-------------------------------------------------------
comparison_df <- data.frame(
  Index = c("RLGSI", "PPG-LGSI", "CRLGSI", "CPPG-LGSI"),
  Selection_Response = c(
    rlgsi_res$R,
    ppg_res$R,
    crlgsi_res$R,
    cppg_res$R
  ),
  Overall_Genetic_Advance = c(
    rlgsi_res$GA,
    ppg_res$GA,
    crlgsi_res$GA,
    cppg_res$GA
  )
)

print(comparison_df)

