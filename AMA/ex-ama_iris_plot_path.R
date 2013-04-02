## This script plots the convex_clustering algorithms in the ama.f90 module
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

## Run sloppy code
X = t(as.matrix(iris[,1:4]))
p = ncol(X)


#knbs = 10
w = kernel_weights(X,kgamma)
#w = knn_weights(w,knbs,p)
gamma_seq = seq(0,0.0425,length.out=1000)
#gamma_seq = seq(0,1,length.out=1000)
tol = 1e-3
nu = 10
system.time({ccpB = convex_cluster_path(X,w,gamma_seq,nu=nu,tol=tol)})

nGamma = ccpB$nGamma
x = double(nGamma)
y = double(nGamma)
Z = X
for (j in 1:nGamma) {
  Z = cbind(Z,ccpB$UHx[[j]])
}
Zs = t(scale(t(Z),center=TRUE,scale=FALSE))
svdZ = svd(Zs)
pc = svdZ$u[,1:2,drop=FALSE]
pc.df = as.data.frame(t(pc)%*%X)

df.paths = data.frame(x=c(),y=c(), group=c())
for (j in 1:nGamma) {
  pcs = t(pc)%*%ccpB$UHx[[j]]
  x = pcs[1,]
  y = pcs[2,]
  df = data.frame(x=pcs[1,], y=pcs[2,], group=1:p)
  df.paths = rbind(df.paths,df)
}

X_data = as.data.frame(t(X)%*%pc)
colnames(X_data) = c("x","y")
X_data$Species = iris[,5]
data_plot = ggplot(data=df.paths,aes(x=x,y=y))
data_plot = data_plot + geom_path(aes(group=group),colour='grey60',alpha=0.5)
data_plot = data_plot + geom_point(data=X_data,aes(x=x,y=y,colour=Species,shape=Species))
data_plot = data_plot + xlab('Principal Component 1') + ylab('Principal Component 2')
data_plot + theme_bw()
filename = "iris_l2_path.pdf"
golden_ratio = 1.61803398875
height = 7
ggsave(filename=filename,width=golden_ratio*height,height=height)