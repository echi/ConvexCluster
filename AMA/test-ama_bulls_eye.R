## This script tests the convex_clustering algorithms in the ama.f90 module
## on the bulls eye data

rm(list=ls())
setwd("/Users/ericchi/Dropbox/Work/Research/01_Submitted/SONCluster/Code/AMA")
library(testthat)
source('ama_loss.R')
source('ama_updates.R')
source('ama_convex_cluster.R')
source('cluster_path_preprocess.R')

source('../ex_bulls_eye.R')
X = bulls_eye(nc=50,np=450)$X
p = ncol(X)
w = kernel_weights(X,1)
w = knn_weights(w,5,p)

gamma_seq = seq(0.0, 50, length.out=500)
tol = 1e-1
nu = 10

system.time({ccp = convex_cluster_path(X,w,gamma_seq,nu=nu,tol=tol)})

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