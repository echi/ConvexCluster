## This script tests the convex_clustering algorithms in the admm.f90 module
## on mammals data

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
mammals = read.table('../Examples/mammals_homals.csv',sep=",",head=TRUE)
X = as.matrix(mammals[,-1])
X = t(scale(X,center=TRUE,scale=FALSE))
colnames(X) = mammals[,1]
X = X[,-which(duplicated(t(X)))]
p = ncol(X)
#****************************#
# Get Principal Components   #
#****************************#
svdX = svd(X)
pc = svdX$u[,1:2,drop=FALSE]
pc.df = as.data.frame(t(pc)%*%X)
#****************************#
#  Parameter choices for     # 
#  different scenarios       #
#                            #
#  L1 norm: scenario = 1     #
#  L2 norm: scenario = 2     #
#  L2 norm: scenario = 3     #
#****************************#
scenario = 2
if (scenario==1) {
  type = 1
  knbs = 5
  kgamma = 0.5
  w = kernel_weights(X,kgamma)
  w = knn_weights(w,knbs,p)
  gamma_seq = seq(0.0,40, length.out=1000)
  tol = 1e-3
  nu = 1
} else if (scenario==2) {
  type = 2
  knbs = 5
  kgamma = 0.5
  w = kernel_weights(X,kgamma)
  w = knn_weights(w,knbs,p)
  gamma_seq = seq(0.0,43, length.out=1000)
  tol = 1e-3
  nu = 1
} else if (scenario==3) {
  type = 3
  knbs = 5
  kgamma = 0.5
  w = kernel_weights(X,kgamma)
  w = knn_weights(w,knbs,p)
  gamma_seq = seq(0.0,1.24, length.out=1000)
  tol = 1e-3
  nu = 1
}
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
  if (gam > 0.2) {
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
#****************************#
#* Collect Path Data Frame  *#
#****************************#
df.paths = data.frame(x=c(),y=c(), group=c(), parmater_set=c())
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
kgamma = 0.5
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
  if (gam > 45) {
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
set.seed(123)
X_data = as.data.frame(t(X)%*%pc)
colnames(X_data) = c("x","y")
X_data$Name = colnames(X)
data_plot = ggplot(data=df.paths,aes(x=x,y=y))
data_plot = data_plot + geom_path(aes(group=group),colour='grey30',alpha=0.5)
data_plot = data_plot + geom_text(data=X_data,aes(x=x,y=y,label=Name), position=position_jitter(h=0.125,w=0.125))
data_plot = data_plot + geom_point(data=X_data,aes(x=x,y=y),size=1.5)
data_plot = data_plot + facet_grid(.~parameter_set)
data_plot = data_plot + xlab('Principal Component 1') + ylab('Principal Component 2')
data_plot + theme_bw()
#****************************#
#* Save Results             *#
#****************************#
filename = "mammals_paths_l2.pdf"
golden_ratio = 2
height = 9
ggsave(filename=filename,width=golden_ratio*height,height=height)