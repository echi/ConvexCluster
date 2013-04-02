
rm(list=ls())
setwd("/Users/ericchi/Dropbox/Work/Research/00_Active/SONCluster/Code/AMA")
library(testthat)
source('ama_loss.R')
source('ama_updates.R')
source('ama_convex_cluster.R')
source('cluster_path_preprocess.R')
source('../ex_two_moons.R')
X = two_moons(n=200,p=2)$X
p = ncol(X)

df.paths = data.frame(x=c(),y=c(), group=c(),gamma=c(),k=c())

## Trial 1
kgamma = 0
knbs = 10
w = kernel_weights(X,kgamma)
w = knn_weights(w,knbs,p)
gamma_seq = seq(0.0,0.63, length.out=100)
tol = 1e-1
ix = compactify_edges(w,p)$ix
A = matrix(0,nrow(ix),p)
for (i in 1:nrow(ix)) {
  A[i,ix[,1]] = 1
  A[i,ix[,2]] = -1
}
nu = 0.99*2/(svd(A)$d[1]**2)
system.time({ccp = convex_cluster_path_fixed(X,w,gamma_seq,nu=nu,tol=tol)})