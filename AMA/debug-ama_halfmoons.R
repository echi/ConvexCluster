## Debugging L-infinity on half-moons


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


## Trial 3
kgamma = 0.5
knbs = 50
w = kernel_weights(X,kgamma)
w = knn_weights(w,knbs,p)
gamma_seq = seq(0.0,0.16, length.out=1000)
tol = 1e-1
nu = 10
gamma = gamma_seq[145]
for (igamma in 1:length(gamma_seq)) {
  gamma = gamma_seq[igamma]
  system.time({ccp = convex_cluster_path(X,w,gamma,nu=nu,tol=tol,type=3)})
}
