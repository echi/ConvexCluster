## This script tests the updates in the ama.f90 module

rm(list=ls())
library(testthat)
source('ama_loss.R')
source('cluster_path_preprocess.R')

context("Loss functions")

test_that("loss_primal is working", {
  ## Create random problem
  set.seed(12345)
  nTrials = 100
  q = 10
  p = 20
  for (i in 1:nTrials) {
    U = matrix(rnorm(q*p),q,p)
    X = matrix(rnorm(q*p),q,p)
    w = kernel_weights(X,0)
    w = knn_weights(w,3,p)
    edge_info = compactify_edges(w,p)
    w = w[w>0]
    gamma = rexp(1)
    
    ix = edge_info$ix
    type = 2

    expect_that(loss_primalF(X,U,gamma,w,ix,type), equals(loss_primalR(X,U,gamma,w,ix)))
  }
})

test_that("loss_dual is working", {
  ## Create random problem
  set.seed(12345)
  nTrials = 100
  q = 10
  p = 20
  for (i in 1:nTrials) {
    X = matrix(rnorm(q*p),q,p)
    w = kernel_weights(X,0)
    w = knn_weights(w,3,p)
    edge_info = compactify_edges(w,p)
    
    M1 = edge_info$M1
    M2 = edge_info$M2
    s1 = edge_info$s1
    s2 = edge_info$s2
    ix = edge_info$ix
    nK = length(which(w>0))
    Lambda = matrix(rnorm(q*nK),q,nK)
    
    loss_dualF(X,Lambda,ix,M1,M2,s1,s2)
      
    expect_that(loss_dualF(X,Lambda,ix,M1,M2,s1,s2), equals(loss_dualR(X,Lambda,ix,M1,M2,s1,s2)))
  }
})