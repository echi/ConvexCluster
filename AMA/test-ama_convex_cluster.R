## This script tests the convex_clustering algorithms in the ama.f90 module

rm(list=ls())
setwd("/Users/ericchi/Dropbox/Work/Research/00_Active/SONCluster/Code/AMA")
library(testthat)
source('ama_loss.R')
source('ama_updates.R')
source('ama_convex_cluster.R')
source('cluster_path_preprocess.R')

context("convex_clustering functions")

test_that("Convex clustering algorithms are working", {
  ## Create random problem
  set.seed(12345)
  nTrials = 100
  q = 10
  p = 20
  nu = 10
  gamma = 0.5
  for (i in 1:nTrials) {
    X = matrix(rnorm(q*p),q,p)
    w = kernel_weights(X,0)
    w = knn_weights(w,3,p)
    edge_info = compactify_edges(w,p)
    ix = edge_info$ix
    M1 = edge_info$M1
    M2 = edge_info$M2
    s1 = edge_info$s1
    s2 = edge_info$s2
    ix = edge_info$ix
    nK = length(which(w>0))    
    w = w[w>0]
    Lambda = matrix(rnorm(q*nK),q,nK)    
#    gamma = rexp(1)
    type = 2
    max_iter = 1e6
    tol = 1e-13
    expect_that(convex_cluster_backtrack(X,Lambda,ix,M1,M2,s1,s2,w,
                                         gamma,nu,max_iter=max_iter,tol=tol)$U,
                equals(convex_cluster_fista_backtrack(X,Lambda,ix,M1,M2,s1,s2,w,
                                                      max_iter=max_iter,tol=tol,gamma,nu)$U))
  }
})