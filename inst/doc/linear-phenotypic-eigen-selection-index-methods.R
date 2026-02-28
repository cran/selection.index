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

## ----esim_example-------------------------------------------------------------
# Compute the linear phenotypic eigen selection index
esim_res <- esim(
  pmat = pmat,
  gmat = gmat,
  selection_intensity = 2.063
)

# View the summary
print(esim_res$summary)

# View expected genetic gains per trait (Delta_G)
print(esim_res$Delta_G)

## ----resim_example------------------------------------------------------------
# We restrict PlantHeight (the second trait in our data frame)
resim_res <- resim(
  pmat = pmat,
  gmat = gmat,
  restricted_traits = c(2),
  selection_intensity = 2.063
)

# Expected genetic gains per trait
print(resim_res$Delta_G)

## ----ppg_esim_example---------------------------------------------------------
# Provide the vector d of desired predetermined proportional gains
# The indices represent non-zero targets. For exactly corresponding target magnitudes, we provide d.
d_vector <- c(0, 0.5, -1) # e.g. proportional mapping

# We use the ppg_esim function
ppgesim_res <- ppg_esim(
  pmat = pmat,
  gmat = gmat,
  d = d_vector,
  selection_intensity = 2.063
)

# View the proportionality of the genetic gains
print(ppgesim_res$Delta_G)

