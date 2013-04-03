## This script tests the convex_clustering algorithms in the ama.f90 module
## on the two moons data

rm(list=ls())
setwd("/Users/ericchi/Dropbox/Work/Research/01_Submitted/SONCluster/Code/AMA")
library(testthat)
source('ama_loss.R')
source('ama_updates.R')
source('ama_convex_cluster.R')
source('cluster_path_preprocess.R')
library(clusterpath)
source('../ex_two_moons.R')

#######################################
## Part B: Timing comparisons        ##
#######################################

nRuns = 10

times_hok_l2 = double(nRuns)
times_hok_l1 = double(nRuns)
times_ama_slow_l2 = double(nRuns)
times_ama_slow_l1 = double(nRuns)
times_ama_fast_l2 = double(nRuns)
times_ama_fast_l1 = double(nRuns)

X = two_moons(n=200,p=2)$X
p = ncol(X)
q = nrow(X)
#****************************#
#* Run Hocking algorithm    *#
#****************************#
kgamma = 1e-2
for (i in 1:nRuns) {
  times_hok_l2[i] = system.time({path_l2 <- clusterpath.l2(x=t(X),gamma=kgamma,lambda=0.006,maxit=10000,verbose=FALSE)})[3] 
  print(paste0("Hocking L2: Completed ",i))
}
gamma_seq_l2 = unique(path_l2$lambda)

for (i in 1:nRuns) {
  times_hok_l1[i] = system.time({path_l1 <- clusterpath.l1.general(t(X),gamma=kgamma,lambda=0.006,maxit=10000)})[3]
  print(paste0("Hocking L1: Completed ",i))
}
gamma_seq_l1 = unique(path_l1$lambda)
#****************************#
#* Run AMA L2 - Fast        *#
#****************************#

#kgamma = 0.5
#knbs = 10
w = kernel_weights(X,kgamma)
#w = knn_weights(w,knbs,p)
#nK = length(which(w>0))
#ix = vec2tri(1:(p*(p-1)/2),p)
#D = matrix(0,p,p)
#innz = which(w>0)
#for (i in 1:length(innz)) {
#  D[ix[innz[i],1],ix[innz[i],2]] = w[innz[i]]
#}
#gamma_seq = seq(0.0,6, length.out=1000)
type = 2
tol = 1e-5
nu = 10
for (i in 1:nRuns) {
  times_ama_fast_l2[i] = system.time({halfmoons_ama_fast_l2 = convex_cluster_path_fista_backtrack(X,w,gamma_seq_l2,nu=nu,tol=tol,type=type)})[3]
  print(paste0("AMA Fast L2: Completed ",i))
}
#****************************#
#* Run AMA L1 - Fast        *#
#****************************#
type = 1
for (i in 1:nRuns) {
  times_ama_fast_l1[i] = system.time({halfmoons_ama_fast_l1 = convex_cluster_path_fista_backtrack(X,w,gamma_seq_l1,nu=nu,tol=tol,type=type)})[3]
  print(paste0("AMA Fast L1: Completed ",i))
}
#****************************#
#* Run AMA L2 - Slow        *#
#****************************#
type = 2
for (i in 1:nRuns) {
  times_ama_slow_l2[i] = system.time({halfmoons_ama_slow_l2 = convex_cluster_path_backtrack(X,w,gamma_seq_l2,nu=nu,tol=tol,type=type)})[3]
  print(paste0("AMA Slow L2: Completed ",i))
}
#****************************#
#* Run AMA L1 - Slow        *#
#****************************#
type = 1
for (i in 1:nRuns) {
  times_ama_slow_l1[i] = system.time({halfmoons_ama_slow_l1 = convex_cluster_path_backtrack(X,w,gamma_seq_l1,nu=nu,tol=tol,type=type)})[3]
  print(paste0("AMA Slow L1: Completed ",i))
}
#****************************#
#* Compute Objective        *#
#* Relative difference: L2  *#
#****************************#
nGamma = length(gamma_seq_l2)
obj_hok_l2 = double(nGamma)
obj_ama_l2 = double(nGamma)
ix = vec2tri(1:(p*(p-1)/2),p)
type = 2
for (iGamma in 1:nGamma) {
  gamma = gamma_seq_l2[iGamma]
  obj_hok_l2[iGamma] = loss_primalF(X, t(path_l2[(p*(iGamma-1)+1):(p*iGamma),1:q]),gamma=gamma,w=w,ix=ix,type=type)
  obj_ama_l2[iGamma] = loss_primalF(X, halfmoons_ama_fast_l2$UHx[[iGamma]],gamma=gamma,w=w,ix=ix,type=2)
}
rel_ama2hok_l2 =  (obj_hok_l2[-1]-obj_ama_l2[-1])/obj_hok_l2[-1]
#****************************#
#* Compute Objective        *#
#* Relative difference: L1  *#
#****************************#
nGamma = length(gamma_seq_l1)
obj_hok_l1 = double(nGamma)
obj_ama_l1 = double(nGamma)
ix = vec2tri(1:(p*(p-1)/2),p)
type = 1
for (iGamma in 1:nGamma) {
  gamma = gamma_seq_l1[iGamma]
  obj_hok_l1[iGamma] = loss_primalF(X, t(path_l1[(p*(iGamma-1)+1):(p*iGamma),1:q]),gamma=gamma,w=w,ix=ix,type=type)
  obj_ama_l1[iGamma] = loss_primalF(X, halfmoons_ama_fast_l1$UHx[[iGamma]],gamma=gamma,w=w,ix=ix,type=type)
}
rel_ama2hok_l1 =  (obj_hok_l1[-1]-obj_ama_l1[-1])/obj_hok_l1[-1]

save(list=c("times_hok_l2","times_hok_l1","gamma_seq_l2","gamma_seq_l1",
            "times_ama_fast_l2","times_ama_fast_l1","times_ama_slow_l2","times_ama_slow_l1",
            "obj_hok_l2","obj_hok_l1","obj_ama_l2","obj_ama_l1","rel_ama2hok_l2","rel_ama2hok_l1"),
     file="halfmoons_hocking_ama.rda")