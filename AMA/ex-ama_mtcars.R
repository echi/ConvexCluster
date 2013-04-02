## This script tests the convex_clustering algorithms in the ama.f90 module
## on senate data

rm(list=ls())
setwd("/Users/ericchi/Dropbox/Work/Research/00_Active/SONCluster/Code/AMA")
library(testthat)
library(ggplot2)
source('ama_loss.R')
source('ama_updates.R')
source('ama_convex_cluster.R')
source('cluster_path_preprocess.R')
library(clusterpath)

data(mtcars)

X = t(as.matrix(scale(mtcars)))
kgamma = 1e-1
p = ncol(X)

df.paths = data.frame(x=c(),y=c(), group=c())

knbs = 5
w = kernel_weights(X,kgamma)
w = knn_weights(w,knbs,p)
gamma_seq = seq(0.0,6, length.out=2000)
tol = 1e-3
nu = 10
system.time({ccp = convex_cluster_path(X,w,gamma_seq,nu=nu,tol=tol,type=1)})

ix = compactify_edges(w,p)$ix
nGamma = ccp$nGamma
x = double(nGamma)
y = double(nGamma)
Z = X
for (j in 1:nGamma) {
  Z = cbind(Z,ccp$UHx[[j]])    
}
svdX = svd(X)
svdZ = svd(Z)
pc = svdZ$u[,1:2,drop=FALSE]
pc = svdX$u[,1:2,drop=FALSE]

pc.df = as.data.frame(t(pc)%*%X)

for (j in 1:nGamma) {
  pcs = t(pc)%*%ccp$UHx[[j]]
  x = pcs[1,]
  y = pcs[2,]
  df = data.frame(x=pcs[1,], y=pcs[2,], group=1:p)
  df.paths = rbind(df.paths,df)
}

X_data = as.data.frame(t(X)%*%pc)
colnames(X_data) = c("x","y")
data_plot = ggplot(data=df.paths,aes(x=x,y=y))
data_plot = data_plot + geom_path(aes(group=group),colour='grey60',alpha=0.5)
#data_plot = data_plot + geom_text(data=X_data,aes(x=x,y=y,label=colnames(X)))
data_plot = data_plot + geom_point(data=X_data,aes(x=x,y=y))
data_plot = data_plot + xlab('Principal Component 1') + ylab('Principal Component 2')
data_plot + xlim(c(1.75,2.75)) + ylim(c(-0.75,-0.2))
data_plot + theme_bw()


