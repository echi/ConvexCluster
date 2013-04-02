## This script tests the convex_clustering algorithms in the ama.f90 module
## on the two moons data

rm(list=ls())
setwd("/Users/ericchi/Dropbox/Work/Research/00_Active/SONCluster/Code/ADMM")
library(testthat)
source('admm.R')
source('admm_convex_cluster.R')
source('cluster_path_preprocess.R')

source('../ex_two_moons.R')
X = two_moons(n=200,p=2)$X
p = ncol(X)
kgamma = 0.5
knbs = 10
w = kernel_weights(X,kgamma)
w = knn_weights(w,knbs,p)

gamma_seq = seq(0.0,6,length.out=1000)
gamma_seq = c(gamma_seq,seq(6.01,8.2, length.out=1000))
tol = 1e-3
nu = 1

system.time({ccp = convex_cluster_pathB(X,w,gamma_seq[1:5],nu=nu,tol=tol,type=2)})
system.time({ccp_acc = convex_cluster_path_acc(X,w,gamma_seq,nu=nu,tol=tol,type=2)})


## Visualize results
X_data = as.data.frame(t(X[1:2,]))
colnames(X_data) = c("x","y")

nGamma = ccp$nGamma
x = double(nGamma)
y = double(nGamma)
data_plot = ggplot(data=X_data,aes(x=x,y=y))
data_plot = data_plot + geom_point(size=2)
for (i in 1:p) {
  for (j in 1:nGamma) {
    x[j] = ccp$UHx[[j]][1,i]
    y[j] = ccp$UHx[[j]][2,i]
  }
  df = data.frame(x=x, y=y)
  data_plot = data_plot + geom_path(data=df,aes(x=x, y=y),colour='blue',alpha=0.5)
}
data_plot + theme_bw() + theme(legend.position = "none") 

df.paths = data.frame(x=c(),y=c(), group=c(),gamma=c(),k=c())
for (i in 1:p) {
  for (j in 1:nGamma) {
    x[j] = ccp$UHx[[j]][1,i]
    y[j] = ccp$UHx[[j]][2,i]
  }
  df = data.frame(x=x, y=y)
  df$group = i
  df$gamma = kgamma
  df$k = knbs
  df.paths = rbind(df.paths,df)
}
data_plot = ggplot(data=df.paths,aes(x=x,y=y))
data_plot + geom_path(aes(group=group))