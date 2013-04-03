## This script runs some unit tests for the convex_clustering algorithms in the admm.f90 module
##

rm(list=ls())
library(testthat)
source('R/Functions/admm.R')
source('R/Functions/cluster_path_preprocess.R')

if (is.loaded('test_prox')) {
  dyn.unload('Fortran/test-prox.so')
}
dyn.load('Fortran/test-prox.so')

test_file(path='R/Tests/test-admm.R')