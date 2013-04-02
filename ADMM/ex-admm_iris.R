## This script tests the convex_clustering algorithms in the admm.f90 module
## on iris data

rm(list=ls())
setwd("/Users/ericchi/Dropbox/Work/Research/00_Active/SONCluster/Code/ADMM")
library(ggplot2)
source('admm.R')
source('admm_convex_cluster.R')
source('cluster_path_preprocess.R')
library(clusterpath)

#######################################
## Part A: Create plots              ##
#######################################

X = as.matrix(iris[,1:4])
X = t(scale(X,center=TRUE,scale=FALSE)) #<- Center the data.
dupix = which(duplicated(t(X))) #<- Remove duplicates
X = X[,-dupix]
p = ncol(X)
#****************************#
# Get Principal Components   #
#****************************#
svdX = svd(X)
pc = svdX$u[,1:2,drop=FALSE]
pc.df = as.data.frame(t(pc)%*%X)
#****************************#
#  ADMM parameters           #
#****************************#
tol = 1e-3
nu = 1
#****************************#
#* Run Parameter Set A      *#
#****************************#
kgamma = 0
w = kernel_weights(X,kgamma)
type = 2 #<- Use L2 norm
#****************************#
#* Create gamma sequence    *#
#* for Set A                *#
#****************************#
gamma_seq =c()
gamma_seq[1] = 0
gamma_seq[2] = 0.01
i = 3
repeat {
  gam = 1.01*gamma_seq[i-1]
  if (gam > 0.0275) {
    break
  }
  gamma_seq[i] = gam
  i = i + 1
}
gamma_seq = gamma_seq[1:(i-1)]
#****************************#
#* Run ADMM on Set A        *#
#****************************#
system.time({ccp = convex_cluster_path_acc(X,w,gamma_seq,nu=nu,tol=tol,type=type)})
nGamma = ccp$nGamma
ix = vec2tri(k=1:(p*(p-1)/2),p)
Uout = consolidate_U(ccp$UHx,ccp$VHx,ix,nGamma)
#Uout = ccp$UHx
#****************************#
#* Collect Path Data Frame  *#
#****************************#
df.paths = data.frame(x=c(),y=c(), group=c())
for (j in 1:nGamma) {
  pcs = t(pc)%*%Uout[[j]]
  x = pcs[1,]
  y = pcs[2,]
  df = data.frame(x=pcs[1,], y=pcs[2,], group=1:p, parameter_set="A")
  df.paths = rbind(df.paths,df)
}
#****************************#
#* Run Parameter Set B      *#
#****************************#
kgamma = 4
knbs = 5
w = kernel_weights(X,kgamma)
w = knn_weights(w,knbs,p)
type = 2 #<- Use L2 norm
#****************************#
#* Create gamma sequence    *#
#* for Set B                *#
#****************************#
gamma_seq =c()
gamma_seq[1] = 0
gamma_seq[2] = 1.01
i = 3
repeat {
  gam = 1.01*gamma_seq[i-1]
  if (gam > 21) {
    break
  }
  gamma_seq[i] = gam
  i = i + 1
}
gamma_seq = gamma_seq[1:(i-1)]
#****************************#
#* Run ADMM on Set B        *#
#****************************#
system.time({ccp = convex_cluster_path_acc(X,w,gamma_seq,nu=nu,tol=tol,type=type)})
nGamma = ccp$nGamma
ix = vec2tri(k=1:(p*(p-1)/2),p)
Uout = consolidate_U(ccp$UHx,ccp$VHx,ix,nGamma)
#Uout = ccp$UHx
#****************************#
#* Collect Path Data Frame  *#
#****************************#
for (j in 1:nGamma) {
  pcs = t(pc)%*%Uout[[j]]
  x = pcs[1,]
  y = pcs[2,]
  df = data.frame(x=pcs[1,], y=pcs[2,], group=1:p, parameter_set="B")
  df.paths = rbind(df.paths,df)
}
#****************************#
#* Plot Results             *#
#****************************#
X_data = as.data.frame(t(X)%*%pc)
colnames(X_data) = c("x","y")
X_data$Species = iris[-dupix,5]
data_plot = ggplot(data=df.paths,aes(x=x,y=y))
data_plot = data_plot + geom_path(aes(group=group),colour='grey60',alpha=0.5)
data_plot = data_plot + geom_point(data=X_data,aes(x=x,y=y,colour=Species,shape=Species))
data_plot = data_plot + facet_grid(.~parameter_set)
data_plot = data_plot + xlab('Principal Component 1') + ylab('Principal Component 2')
data_plot + theme_bw()
#****************************#
#* Save Results             *#
#****************************#
filename = "iris_paths_l2.pdf"
golden_ratio = 2
height = 9
ggsave(filename=filename,width=golden_ratio*height,height=height)
