## This script tests the convex_clustering algorithms in the admm.f90 module
## on halfmoons data

# rm(list=ls())
# setwd("/Users/ericchi/Dropbox/Work/Research/00_Active/SONCluster/Code/ADMM")
# library(ggplot2)
# source('admm.R')
# source('admm_convex_cluster.R')
# source('cluster_path_preprocess.R')
# library(clusterpath)
# source('../ex_two_moons.R')

#######################################
## Part B: Time trials               ##
#######################################

nRuns = 10

times_admm_slow_l2 = double(nRuns)
times_admm_slow_l1 = double(nRuns)
times_admm_fast_l2 = double(nRuns)
times_admm_fast_l1 = double(nRuns)

X = two_moons(n=200,p=2)$X
p = ncol(X)
q = nrow(X)
#****************************#
#* Run Hocking algorithm    *#
#* to get gamma sequences;  *#
#* the timings are done in  *#
#* ex-ama_senate_itiming.R  *#
#****************************#
kgamma = 1e-2
system.time({path_l2 <- clusterpath.l2(x=t(X),gamma=kgamma,lambda=0.006,maxit=10000,verbose=FALSE)})[3] 
gamma_seq_l2 = unique(path_l2$lambda)

system.time({path_l1 <- clusterpath.l1.general(t(X),gamma=kgamma,lambda=0.006,maxit=10000)})[3]
gamma_seq_l1 = unique(path_l1$lambda)
#****************************#
#* Run ADMM L2 - Slow       *#
#****************************#
w = kernel_weights(X,kgamma)
type = 2
tol = 5e-2
nu = 1
for (i in 1:nRuns) {
  times_admm_slow_l2[i] = system.time({halfmoons_admm_slow_l2 = convex_cluster_pathB_admm(X,w,gamma_seq_l2,nu=nu,tol=tol,type=type)})[3]
  print(paste0("ADMM slow L2: Completed ",i)) 
}
#****************************#
#* Run ADMM L2 - Fast       *#
#****************************#
for (i in 1:nRuns) {
  times_admm_fast_l2[i] = system.time({halfmoons_admm_fast_l2 = convex_cluster_path_admm_acc(X,w,gamma_seq_l2,nu=nu,tol=tol,type=type)})[3]
  print(paste0("ADMM fast L2: Completed ",i)) 
}
#****************************#
#* Run ADMM L1 - Slow       *#
#****************************#
type = 1
for (i in 1:nRuns) {
  times_admm_slow_l1[i] = system.time({halfmoons_admm_slow_l1 = convex_cluster_pathB_admm(X,w,gamma_seq_l1,nu=nu,tol=tol,type=type)})[3]
  print(paste0("ADMM slow L1: Completed ",i)) 
}
#****************************#
#* Run ADMM L1 - Fast       *#
#****************************#
for (i in 1:nRuns) {
  times_admm_fast_l1[i] = system.time({halfmoons_admm_fast_l1 = convex_cluster_path_admm_acc(X,w,gamma_seq_l1,nu=nu,tol=tol,type=type)})[3]
  print(paste0("ADMM fast L1: Completed ",i)) 
}
#######################################
## Loss function (Primal)             #
#######################################
if (is.loaded('loss_primal')) {
  dyn.unload('../AMA/Fortran/ama.so')
}
dyn.load('../AMA/Fortran/ama.so')
loss_primalF = function(X,U,gamma,w,ix,type) {
  dimU = dim(U)
  p = as.integer(dimU[2])
  q = as.integer(dimU[1])
  nK = as.integer(nrow(ix))
  gamma = as.double(gamma)  
  storage.mode(X) = "double"
  storage.mode(U) = "double"
  w = as.double(w)
  storage.mode(ix) = "integer"
  type = as.integer(type)
  sol = .Fortran('loss_primal',X=X,U=U,gamma=gamma,ix=ix,p=p,q=q,nK=nK,w=w,output=double(1),type=type)
  return(sol$output)
}
#****************************#
#* Compute Objective        *#
#* Relative difference: L2  *#
#****************************#
nGamma = length(gamma_seq_l2)
obj_admm_l2 = double(nGamma)
obj_hok_l2 = double(nGamma)
ix = vec2tri(1:(p*(p-1)/2),p)
type = 2
for (iGamma in 1:nGamma) {
  gamma = gamma_seq_l2[iGamma]
  obj_hok_l2[iGamma] = loss_primalF(X, t(path_l2[(p*(iGamma-1)+1):(p*iGamma),1:q]),gamma=gamma,w=w,ix=ix,type=type)
  obj_admm_l2[iGamma] = loss_primalF(X=X,U=halfmoons_admm_fast_l2$UHx[[iGamma]],gamma=gamma,w,ix,type)
  print(iGamma)
}
rel_admm2hok_l2 = (obj_hok_l2[-1]-obj_admm_l2[-1])/obj_hok_l2[-1]
#****************************#
#* Compute Objective        *#
#* Relative difference: L1  *#
#****************************#
nGamma = length(gamma_seq_l1)
obj_admm_l1 = double(nGamma)
obj_hok_l1 = double(nGamma)
type = 1
for (iGamma in 1:nGamma) {
  gamma = gamma_seq_l1[iGamma]
  obj_hok_l1[iGamma] = loss_primalF(X, t(path_l1[(p*(iGamma-1)+1):(p*iGamma),1:q]),gamma=gamma,w=w,ix=ix,type=type)
  obj_admm_l1[iGamma] = loss_primalF(X=X,U=halfmoons_admm_fast_l1$UHx[[iGamma]],gamma=gamma,w,ix,type)
  print(iGamma)
}

rel_admm2hok_l1 = (obj_hok_l1[-1]-obj_admm_l1[-1])/obj_hok_l1[-1]
save(list=c("gamma_seq_l2","gamma_seq_l1",
            "times_admm_fast_l2","times_admm_fast_l1","times_admm_slow_l2","times_admm_slow_l1",
            "obj_hok_l2","obj_hok_l1","obj_admm_l2","obj_admm_l1","rel_admm2hok_l2","rel_admm2hok_l1"),
     file="halfmoons_hocking_admm.rda")