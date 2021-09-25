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
si<- list()
si[[1]]<- sel.index(ID = 1, phen_mat = pmat[1,1], gen_mat = gmat[1,1],
                    weight_mat = w[1,2], GAY = GAY)
si[[2]]<- sel.index(ID = 2, phen_mat = pmat[2,2], gen_mat = gmat[2,2],
                    weight_mat = w[2,2], GAY = GAY)
si[[3]]<- sel.index(ID = 3, phen_mat = pmat[3,3], gen_mat = gmat[3,3],
                    weight_mat = w[3,2], GAY = GAY)
si[[4]]<- sel.index(ID = 4, phen_mat = pmat[4,4], gen_mat = gmat[4,4],
                    weight_mat = w[4,2], GAY = GAY)
si[[5]]<- sel.index(ID = 5, phen_mat = pmat[5,5], gen_mat = gmat[5,5],
                    weight_mat = w[5,2], GAY = GAY)
si[[6]]<- sel.index(ID = 6, phen_mat = pmat[6,6], gen_mat = gmat[6,6],
                    weight_mat = w[6,2], GAY = GAY)
si[[7]]<- sel.index(ID = 7, phen_mat = pmat[7,7], gen_mat = gmat[7,7],
                    weight_mat = w[7,2], GAY = GAY)

## -----------------------------------------------------------------------------
b<- si[[2]][[2]]
sel.score.rank(data = d[,3], bmat = b, genotype = d$treat)

## -----------------------------------------------------------------------------
comb.indices(ncomb = 1, pmat = pmat, gmat = gmat, wmat = w[,-1], wcol = 1, GAY = GAY)

## -----------------------------------------------------------------------------
rcomb.indices(ncomb = 1, i = 1, pmat = pmat, gmat = gmat, wmat = w[,-1], wcol = 1, GAY = GAY)

