## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(selection.index)

# To ensure reproducibility
set.seed(42)

## ----define_markers-----------------------------------------------------------
n_traits <- 3
n_loci <- 50 # Number of segregating sites / markers

# Generate random base QTL effects for the markers across the 3 traits
# Negative correlation infused between trait 1 and 2
qtl_eff <- matrix(rnorm(n_loci * n_traits), nrow = n_loci, ncol = n_traits)
qtl_eff[, 2] <- -0.5 * qtl_eff[, 1] + 0.5 * qtl_eff[, 2]

# Define heritabilities and corresponding environmental variance
heritabilities <- c(0.2, 0.5, 0.5)

# Simulate base genetic variance to deduce correct environmental variance noise
base_gv <- apply(qtl_eff, 2, var) * n_loci
env_var <- base_gv * (1 - heritabilities) / heritabilities

## ----define_weights-----------------------------------------------------------
# Equal economic trait weighting
weights <- c(1, 1, 1)

# Constraint matrix for RLPSI/RESIM: Constrain Trait 1
U_mat <- matrix(0, nrow = 3, ncol = 1)
U_mat[1, 1] <- 1

## ----run_sim------------------------------------------------------------------
# Run the stochastic selection (may take a moment)
sim_results <- simulate_selection_cycles(
  n_cycles = 5,
  n_individuals = 200,
  n_loci = n_loci,
  n_traits = n_traits,
  qtl_effects = qtl_eff,
  heritability = heritabilities,
  economic_weights = weights,
  selection_proportion = 0.25, # Select upper 25% progeny
  restricted_traits = 1
)

## ----interpret_gain-----------------------------------------------------------
# Expected: Because Trait 1 was constrained via the U_mat for the RLPSI metric,
# its expected generational gain should stabilize at 0.
print(sim_results$rlpsi_gain)

## ----interpret_lpsi_gain------------------------------------------------------
print(sim_results$lpsi_gain)

## ----interpret_var------------------------------------------------------------
# Observe the diminishing variance arrays for the LPSI evaluations
print(sim_results$lpsi_var)

