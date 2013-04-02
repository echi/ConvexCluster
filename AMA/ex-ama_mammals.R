## This script tests the convex_clustering algorithms in the ama.f90 module
## on mammals data

rm(list=ls())
setwd("/Users/ericchi/Dropbox/Work/Research/00_Active/SONCluster/Code/AMA")
library(testthat)
library(ggplot2)
source('ama_loss.R')
source('ama_updates.R')
source('ama_convex_cluster.R')
source('cluster_path_preprocess.R')
library(clusterpath)

mammals = read.table('../Examples/mammals_homals.csv',sep=",",head=TRUE)

X = t(as.matrix(mammals[,-1]))
colnames(X) = mammals[,1]
X = X[,-which(duplicated(t(X)))]
p = ncol(X)

## Parameter choices for different scenarios
# L1 norm
## scenario = 1
# L2 norm
## scenario = 2
# L-infinity norm
 scenario = 3

if (scenario==1) {
  type = 1
  knbs = 5
  kgamma = 0.5
  w = kernel_weights(X,kgamma)
  w = knn_weights(w,knbs,p)
  gamma_seq = seq(0.0,40, length.out=1000)
  tol = 1e-3
  nu = 10
} else if (scenario==2) {
  type = 2
  knbs = 5
  kgamma = 0.5
  w = kernel_weights(X,kgamma)
  w = knn_weights(w,knbs,p)
  gamma_seq = seq(0.0,43, length.out=1000)
  tol = 1e-3
  nu = 10
} else if (scenario==3) {
  type = 3
  knbs = 5
  kgamma = 0.5
  w = kernel_weights(X,kgamma)
  w = knn_weights(w,knbs,p)
  gamma_seq = seq(0.0,1.24, length.out=1000)
  tol = 1e-1
  ix = compactify_edges(w,p)$ix
  A = matrix(0,nrow(ix),p)
  for (i in 1:nrow(ix)) {
    A[i,ix[,1]] = 1
    A[i,ix[,2]] = -1
  }
  nu = (2/svd(t(A)%*%A)$d[1])*0.99
}
system.time({ccp = convex_cluster_path(X,w,gamma_seq,nu=nu,tol=tol,type=type)})
#system.time({ccp = convex_cluster_path_fixed(X,w,gamma_seq,nu=nu,tol=tol,type=type)})

ix = compactify_edges(w,p)$ix
nGamma = ccp$nGamma
x = double(nGamma)
y = double(nGamma)
#Z = X
#for (j in 1:nGamma) {
#  Z = cbind(Z,ccp$UHx[[j]])    
#}
#svdZ = svd(Z)
#pc = svdZ$u[,1:2,drop=FALSE]
svdX = svd(X)
pc = svdX$u[,1:2,drop=FALSE]

pc.df = as.data.frame(t(pc)%*%X)

df.paths = data.frame(x=c(),y=c(), group=c())
for (j in 1:nGamma) {
  pcs = t(pc)%*%ccp$UHx[[j]]
  x = pcs[1,]
  y = pcs[2,]
  df = data.frame(x=pcs[1,], y=pcs[2,], group=1:p)
  df.paths = rbind(df.paths,df)
}
set.seed(12345)
X_data = as.data.frame(t(X)%*%pc)
colnames(X_data) = c("x","y")
data_plot = ggplot(data=df.paths,aes(x=x,y=y))
data_plot = data_plot + geom_path(aes(group=group),colour='grey30',alpha=0.5)
data_plot = data_plot + geom_text(data=X_data,aes(x=x,y=y,label=colnames(X)), position=position_jitter(h=0.08,w=0.08))
#data_plot = data_plot + geom_text(data=X_data,aes(x=x,y=y,label=colnames(X)))
data_plot = data_plot + geom_point(data=X_data,aes(x=x,y=y),size=1.5)
data_plot = data_plot + xlab('Principal Component 1') + ylab('Principal Component 2')
data_plot + theme_bw()
filename = "mammals.pdf"
golden_ratio = 1.61803398875
height = 12
ggsave(filename=filename,width=golden_ratio*height,height=height)
