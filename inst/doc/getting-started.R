## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(selection.index)

# Load the built-in phenotypic dataset
data("seldata")

# Inspect the structure of the dataset
head(seldata)

## ----define_weights-----------------------------------------------------------
# Define economic weights for the 7 traits of interest
weights <- c(10, 8, 6, 4, 2, 1, 1)

# Calculate genotypic and phenotypic variance-covariance matrices
# Traits: columns 3:9, Genotypes: column 2, Replication: column 1
gmat <- gen_varcov(data = seldata[, 3:9], genotypes = seldata[, 2], replication = seldata[, 1])
pmat <- phen_varcov(data = seldata[, 3:9], genotypes = seldata[, 2], replication = seldata[, 1])

## ----calculate_index----------------------------------------------------------
# Calculate the combinatorial selection index for all 7 traits
index_results <- lpsi(
  ncomb = 7,
  pmat = pmat,
  gmat = gmat,
  wmat = as.matrix(weights),
  wcol = 1
)

## ----view_results-------------------------------------------------------------
# View the calculated index metrics for our 7-trait combination
head(index_results)

# Extract the final selection scores to rank the genotypes
scores <- predict_selection_score(
  index_results,
  data = seldata[, 3:9],
  genotypes = seldata[, 2]
)

# View the top ranked genotypes based on their selection scores
head(scores)

