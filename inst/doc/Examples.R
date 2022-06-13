## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(selection.index)
d<- seldata # Manually generated data for analysis which is included in package

## -----------------------------------------------------------------------------
w<- weight # Weights assigned to the traits also include in package

## -----------------------------------------------------------------------------
gmat<- gen.varcov(data = d[,3:9], genotypes = d$treat, replication = d$rep)
print(gmat)

## -----------------------------------------------------------------------------
pmat<- phen.varcov(data = d[,3:9], genotypes = d$treat, replication = d$rep)
print(pmat)

## -----------------------------------------------------------------------------
GAY<- gen.advance(phen_mat = pmat[1,1], gen_mat = gmat[1,1],
                  weight_mat = w[1,2])
print(GAY)

## -----------------------------------------------------------------------------
comb.indices(ncomb = 1, pmat = pmat, gmat = gmat, wmat = w[,-1], wcol = 1, GAY = GAY)

## -----------------------------------------------------------------------------
rcomb.indices(ncomb = 1, i = 1, pmat = pmat, gmat = gmat, wmat = w[,-1], wcol = 1, GAY = GAY)

