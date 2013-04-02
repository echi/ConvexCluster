## This script tests the proximal mappings in the prox.f90 module

# type = 1 : proj_L2                                                                                                                                                                                
# type = 2 : prox_L1                                                                                                                                                                                
# type = 3 : prox_L2                                                                                                                                                                                
# type = 4 : proj_L1                                                                                                                                                                                
# type = 5 : proj_Linfinity                                                                                                                                                                         
# type = 6 : proj_L12                                                                                                                                                                               
# type = 7 : prox_L12

## Log
#  2013.02.25: Types 1, 2, and 3 are coded.

rm(list=ls())
library(testthat)
if (is.loaded('test_prox')) {
  dyn.unload('test-prox.so')
}
dyn.load('test-prox.so')
source('prox.R')

prox = function(x,tau,type) {
  x = as.double(x)
  n = as.integer(length(x))
  type = as.integer(type)
  tau = as.double(tau)
  sol = .Fortran('test_prox',x,y=double(n),n=n,tau=tau,type=type)
  return(sol$y)
}

context("Proximal Mappings")

test_that("proj_L2 is working", {
  set.seed(12345)
  m = 1e3
  nTrials = 100
  X = matrix(rnorm(m*nTrials),m,nTrials)
  Tau = rexp(nTrials)
  type = 1
  for (i in 1:nTrials) {
    x = X[,i]
    tau = Tau[i]
    expect_that(prox(x,tau,type),equals(prox_R(x,tau,type)))
  }
})

test_that("prox_L1 is working", {
  set.seed(12345)
  m = 1e3
  nTrials = 100
  X = matrix(rnorm(m*nTrials),m,nTrials)
  Tau = rexp(nTrials)
  type = 2
  for (i in 1:nTrials) {
    x = X[,i]
    tau = Tau[i]
    expect_that(prox(x,tau,type),equals(prox_R(x,tau,type)))
  }
})

test_that("prox_L2 is working", {
  set.seed(12345)
  m = 1e3
  nTrials = 100
  X = matrix(rnorm(m*nTrials),m,nTrials)
  Tau = rexp(nTrials)
  type = 3
  for (i in 1:nTrials) {
    x = X[,i]
    tau = Tau[i]
    expect_that(prox(x,tau,type),equals(prox_R(x,tau,type)))
  }
})


test_that("proj_Linf is working", {
  set.seed(12345)
  m = 1e3
  nTrials = 100
  X = matrix(rnorm(m*nTrials),m,nTrials)
  Tau = rexp(nTrials)
  type = 4
  for (i in 1:nTrials) {
    x = X[,i]
    tau = Tau[i]
    expect_that(prox(x,tau,type),equals(prox_R(x,tau,type)))
  }
})


test_that("proj_L1 is working", {
  set.seed(12345)
  m = 1e3
  nTrials = 100
  X = matrix(rnorm(m*nTrials),m,nTrials)
  Tau = rexp(nTrials)
  type = 5
  for (i in 1:nTrials) {
    x = X[,i]
    tau = Tau[i]
    expect_that(prox(x,tau,type),equals(prox_R(x,tau,type)))
  }
})