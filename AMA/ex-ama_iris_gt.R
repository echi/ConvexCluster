## This script tests the convex_clustering algorithms in the ama.f90 module
## on iris data

rm(list=ls())
setwd("/Users/ericchi/Dropbox/Work/Research/00_Active/SONCluster/Code/AMA")
library(testthat)
library(ggplot2)
source('ama_loss.R')
source('ama_updates.R')
source('ama_convex_cluster.R')
source('cluster_path_preprocess.R')
library(clusterpath)

## Run clusterpath.l2
kgamma = 1e-2
system.time({path <- clusterpath.l2(iris[,1:4],gamma=kgamma)})
gamma_seq = unique(path$lambda)

## Create longer gamma sequence to take advantage of warm starts; insert 100 extra interpoints
nGamma = length(gamma_seq)
gamma_seq_exp = gamma_seq[1]
nInsert = 100
for (i in 2:nGamma) {
  tmp = seq(gamma_seq[i-1],gamma_seq[i],length.out=nInsert)
  gamma_seq_exp = c(gamma_seq_exp,tmp[-1])
}

## Get ground truth
X = t(as.matrix(iris[,1:4]))
p = ncol(X)

w = kernel_weights(X,kgamma)
nu = (2/length(w))*0.99
tol = 1e-11
system.time({ccp = convex_cluster_path_fixed(X,w,gamma_seq_exp,nu=nu,tol=tol)})
save(list=c('ccp'),file='ex-ama_iris_ground_truth.rda')
