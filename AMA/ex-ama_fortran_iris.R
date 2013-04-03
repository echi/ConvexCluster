## Iris Data for Fortran only implementation
rm(list=ls())
setwd("/Users/ericchi/Dropbox/Work/Research/01_Submitted/SONCluster/Code/AMA")
source('ama_loss.R')
source('ama_updates.R')
source('ama_convex_cluster.R')
source('cluster_path_preprocess.R')
library(clusterpath)

gamma = 0.13
eta = 2
nu = 10
type = 2
max_iter = 1e2
tol = 1e-3

## X
X = t(iris[,1:4])
X = matrix(rnorm(1000**2),1000,1000)
write.table(t(X),"X.dat", sep=" ", row.names=FALSE, col.names=FALSE)
q = ncol(X)
p = nrow(X)

## w
kgamma = 1e-2
w = kernel_weights(t(X),kgamma)
nK = length(which(w > 0))
write.table(w,"w.dat",sep=" ", row.names=FALSE, col.names=FALSE)

## Lambda
Lambda = matrix(0,q,nK)
write.table(t(Lambda),"Lambda.dat",sep=" ", row.names=FALSE, col.names=FALSE)

edge_info = compactify_edges(w,p)
ix = edge_info$ix  
M1 = edge_info$M1
M2 = edge_info$M2
s1 = edge_info$s1
s2 = edge_info$s2
mix1 = nrow(M1)
mix2 = nrow(M2)
write.table(t(ix),"ix.dat", sep=" ", row.names=FALSE, col.names=FALSE)
write.table(t(M1),"M1.dat", sep=" ", row.names=FALSE, col.names=FALSE)
write.table(t(M2),"M2.dat", sep=" ", row.names=FALSE, col.names=FALSE)
write.table(s1,"s1.dat", sep=" ", row.names=FALSE, col.names=FALSE)
write.table(s2,"s2.dat", sep=" ", row.names=FALSE, col.names=FALSE)
write.table(c(q,p,nK,mix1,mix2,gamma,eta,nu,type,max_iter,tol),"problem_size.dat", sep=" ", row.names=FALSE, col.names=FALSE)

