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
wmat<- weight.mat(w)
wmat

## -----------------------------------------------------------------------------
GA1<- sel.index(ID = 1, phen_mat = pmat[1,1], gen_mat = gmat[1,1],
                weight_mat = w[1,2])
print(GA1)

## -----------------------------------------------------------------------------
GAY<- GA1[[3]]

## -----------------------------------------------------------------------------
si_ew<- comb.indices(ncomb = 2, pmat = pmat, gmat = gmat, wmat = wmat, wcol = 1, GAY = 1.7694)
si_ew

## -----------------------------------------------------------------------------
si_h2<- comb.indices(ncomb = 2, pmat = pmat, gmat = gmat, wmat = wmat, wcol = 2, GAY = 1.7694)
si_h2

## -----------------------------------------------------------------------------
b_ew<- c(si_ew$b.1[1], si_ew$b.2[1])
sel.score.rank(data = d[,3:4], bmat = b_ew, genotype = d$treat)

## -----------------------------------------------------------------------------
b_h2<- c(si_h2$b.1[1], si_h2$b.2[1])
sel.score.rank(data = d[,3:4], bmat = b_h2, genotype = d$treat)

