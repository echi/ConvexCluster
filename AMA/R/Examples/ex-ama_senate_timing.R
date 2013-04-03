## This script tests the convex_clustering algorithms in the ama.f90 module
## on senate data

rm(list=ls())
setwd("/Users/ericchi/Dropbox/Work/Research/01_Submitted/SONCluster/Code/AMA")
library(testthat)
library(ggplot2)
source('ama_loss.R')
source('ama_updates.R')
source('ama_convex_cluster.R')
source('cluster_path_preprocess.R')
library(clusterpath)

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

senate = read.table('../Examples/senate_homals.csv',sep=",",head=TRUE)
X = as.matrix(senate[,-c(1:2)])
X = t(scale(X,center=TRUE,scale=FALSE))
dupix = which(duplicated(t(X)))
X = X[,-dupix]
p = ncol(X)
q = nrow(X)
#****************************#
#* Run Hocking algorithm    *#
#****************************#
kgamma = 1e-2
for (i in 1:nRuns) {
  times_hok_l2[i] = system.time({path_l2 <- clusterpath.l2(t(X),gamma=kgamma)})[3] 
  print(paste0("Hocking L2: Completed ",i))
}
gamma_seq_l2 = unique(path_l2$lambda)

for (i in 1:nRuns) {
  times_hok_l1[i] = system.time({path_l1 <- clusterpath.l1.general(t(X),gamma=kgamma)})[3]
  print(paste0("Hocking L1: Completed ",i))
}
gamma_seq_l1 = unique(path_l1$lambda)
#****************************#
#* Run AMA FAST: L2         *#
#****************************#
w = kernel_weights(X,kgamma)
tol = 1e-5
nu = 10
type = 2
for (i in 1:nRuns) {
  times_ama_fast_l2[i] = system.time({senate_ama_fast_l2 = convex_cluster_path_fista_backtrack(X,w,gamma_seq_l2,nu=nu,tol=tol,type=type)})[3]
  print(paste0("AMA Fast L2: Completed ",i))
}
#****************************#
#* Run AMA FAST: L1         *#
#****************************#
type=1
for (i in 1:nRuns) {
  times_ama_fast_l1[i] = system.time({senate_ama_fast_l1 = convex_cluster_path_fista_backtrack(X,w,gamma_seq_l1,nu=nu,tol=tol,type=type)})[3]
  print(paste0("AMA Fast L1: Completed ",i))
}
#****************************#
#* Run AMA SLOW: L2         *#
#****************************#
type = 2
for (i in 1:nRuns) {
  times_ama_slow_l2[i] = system.time({senate_ama_slow_l2 = convex_cluster_path_backtrack(X,w,gamma_seq_l2,nu=nu,tol=tol,type=type)})[3]
  print(paste0("AMA Slow L2: Completed ",i))
}
#****************************#
#* Run AMA SLOW: L1         *#
#****************************#
type = 1
for (i in 1:nRuns) {
  times_ama_slow_l1[i] = system.time({senate_ama_slow_l1 = convex_cluster_path_backtrack(X,w,gamma_seq_l1,nu=nu,tol=tol,type=type)})[3]
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
  obj_ama_l2[iGamma] = loss_primalF(X, senate_ama_fast_l2$UHx[[iGamma]],gamma=gamma,w=w,ix=ix,type=2)
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
  obj_ama_l1[iGamma] = loss_primalF(X, senate_ama_fast_l1$UHx[[iGamma]],gamma=gamma,w=w,ix=ix,type=type)
}
rel_ama2hok_l1 =  (obj_hok_l1[-1]-obj_ama_l1[-1])/obj_hok_l1[-1]

save(list=c("times_hok_l2","times_hok_l1","gamma_seq_l2","gamma_seq_l1",
            "times_ama_fast_l2","times_ama_fast_l1","times_ama_slow_l2","times_ama_slow_l1",
            "obj_hok_l2","obj_hok_l1","obj_ama_l2","obj_ama_l1","rel_ama2hok_l2","rel_ama2hok_l1"),
     file="senate_hocking_ama.rda")