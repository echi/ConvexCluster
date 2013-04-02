## This script tests the convex_clustering algorithms in the ama.f90 module
## on alon data

rm(list=ls())
setwd("/Users/ericchi/Dropbox/Work/Research/00_Active/SONCluster/Code/AMA")
library(testthat)
library(ggplot2)
source('ama_loss.R')
source('ama_updates.R')
source('ama_convex_cluster.R')
source('cluster_path_preprocess.R')
library(clusterpath)

alon_path = "/Users/ericchi/Dropbox/Work/Research/00_Active/SONCluster/Code/Examples/Alon"

X = as.matrix(log10(read.table(paste(alon_path,"alon.txt",sep="/"))))

## Run clusterpath.l2
kgamma = 1e-6
system.time({path <- clusterpath.l2(iris[,1:4],gamma=kgamma)})
gamma_seq = unique(path$lambda)

## Get ground truth
X = t(X)
p = ncol(X)

w = kernel_weights(X,kgamma)
nu = (2/length(w))*0.99
tol = 1e-1
system.time({ccp = convex_cluster_path_fixed(X,w,gamma_seq,nu=nu,tol=tol)})
ix = compactify_edges(w,p)$ix
nGamma = ccp$nGamma
Uout = consolidate_U(ccp$UHx,ccp$VHx,ix,nGamma)

## Run sloppy code
X = t(as.matrix(iris[,1:4]))
p = ncol(X)

df.paths = data.frame(x=c(),y=c(), group=c(),gamma=c(),k=c())