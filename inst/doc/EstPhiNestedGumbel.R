## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(VecDep)

## -----------------------------------------------------------------------------
q = 4
dim = c(2,2)
n = 1000

## -----------------------------------------------------------------------------
hac = gethac(dim,c(2,3,4),type = 1)

## -----------------------------------------------------------------------------
sample = suppressWarnings(HAC::rHAC(n,hac))

## -----------------------------------------------------------------------------
estphi(sample,dim,list("hac",type = 1,M = 10000),function(t){t * log(t)})


