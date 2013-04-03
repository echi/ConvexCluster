## This script tests the updates in the ama.f90 module

rm(list=ls())
library(testthat)
source('ama_updates.R')
source('cluster_path_preprocess.R')

context("Parameter Updates")

test_that("update_U is working", {
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
    nK = length(which(w>0))
    Lambda = matrix(rnorm(q*nK),q,nK)
  
    expect_that(update_UF(X,Lambda,M1,M2,s1,s2), equals(update_UR(X,Lambda,M1,M2,s1,s2)))
  }
})

test_that("update_V is working", {
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
    nK = length(w)
    nu = rexp(1)
    gamma = rexp(1)
    Lambda = matrix(rnorm(q*nK),q,nK)
    ix = edge_info$ix

    expect_that(update_VF(U,Lambda,w,gamma,nu,ix), equals(update_VR(U,Lambda,w,gamma,nu,ix)))
  }
})


test_that("update_Lambda is working", {
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
    nK = length(w)
    nu = rexp(1)
    Lambda = matrix(rnorm(q*nK),q,nK)
    ix = edge_info$ix
    V = matrix(rnorm(q*nK),q,nK)
    
    expect_that(update_LambdaF(Lambda,U,V,nu,ix), equals(update_LambdaR(Lambda,U,V,nu,ix)))
  }
})
  
test_that("update_LambdaDP is working", {
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
    nK = length(w)
    nu = rexp(1)
    Lambda = matrix(rnorm(q*nK),q,nK)
    ix = edge_info$ix
    gamma = rexp(1)    
    expect_that(update_Lambda_combinedF(Lambda,U,w,gamma,nu,ix), equals(update_Lambda_combinedR(Lambda,U,w,gamma,nu,ix)))
  }
})