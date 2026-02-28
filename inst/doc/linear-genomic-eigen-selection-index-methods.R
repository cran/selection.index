## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup_data---------------------------------------------------------------
library(selection.index)

# Load standard phenotype dataset
data(maize_pheno)

# Define traits and design variables
traits <- c("Yield", "PlantHeight", "DaysToMaturity")
env_col <- "Environment"
genotype_col <- "Genotype"

# Phenotypic variance-covariance matrix (P)
pmat <- phen_varcov(maize_pheno[, traits], maize_pheno[[genotype_col]], maize_pheno[[env_col]])

# Genetic variance-covariance matrix (G)
gmat <- gen_varcov(maize_pheno[, traits], maize_pheno[[genotype_col]], maize_pheno[[env_col]])

# For the sake of demonstration within this vignette, we simulate the required molecular/genomic variance components:
set.seed(42)

# Simulate Gamma: Covariance between phenotypes and GEBVs
Gamma <- gmat * 0.85

# Molecular matrices for MESIM
S_M <- gmat * 0.75 # Covariance between phenotypic values and marker scores
S_Mg <- gmat * 0.70 # Covariance between genotypic values and marker scores
S_var <- gmat * 0.80 # Variance-covariance of marker scores

# Genomic matrices for GW-ESIM
G_M <- gmat * 0.82 # Covariance between true genotypic values and marker values
M <- gmat * 0.90 # Variance-covariance matrix of markers

## ----mesim_demo---------------------------------------------------------------
mes_index <- mesim(pmat, gmat, S_M, S_Mg, S_var)
summary(mes_index)

## ----gesim_demo---------------------------------------------------------------
ges_index <- gesim(pmat, gmat, Gamma)
summary(ges_index)

## ----gw_esim_demo-------------------------------------------------------------
gw_index <- gw_esim(pmat, gmat, G_M, M)
summary(gw_index)

## ----rgesim_demo--------------------------------------------------------------
# Restrict the second trait (PlantHeight)
U_mat <- matrix(c(0, 1, 0), nrow = 1)
rges_index <- rgesim(pmat, gmat, Gamma, U_mat)
summary(rges_index)

## ----ppg_gesim_demo-----------------------------------------------------------
# Desired genetic gain proportions: 1 for Yield, 0.5 for DaysToMaturity, 0 for others
d <- c(1, 0, 0.5)
ppg_ges_index <- ppg_gesim(pmat, gmat, Gamma, d)
summary(ppg_ges_index)

