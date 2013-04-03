## This script tests the convex_clustering algorithms in the ama.f90 module
## on senate data

# rm(list=ls())
# setwd("/Users/ericchi/Dropbox/Work/Research/01_Submitted/SONCluster/Code/AMA")
# library(testthat)
# library(ggplot2)
# source('ama_loss.R')
# source('ama_updates.R')
# source('ama_convex_cluster.R')
# source('cluster_path_preprocess.R')
# library(clusterpath)

#######################################
## Part A: Create plots              ##
#######################################
senate = read.table('../Examples/senate_homals.csv',sep=",",head=TRUE)
X = as.matrix(senate[,-c(1:2)])
X = t(scale(X,center=TRUE,scale=FALSE))
dupix = which(duplicated(t(X)))
X = X[,-dupix]
p = ncol(X)
#****************************#
# Create Duplicates Lists    #
#****************************#
nSenators = nrow(senate)
unique_votes = unique(senate[,-c(1:2)])
nUnique = nrow(unique_votes)
myList = integer(nSenators)
for (i in 1:nSenators) {
  vote_i = senate[i,-c(1:2)]
  for (j in 1:nUnique) {
    vote_j = unique_votes[j,]
    if (all(vote_i==vote_j)) {
      myList[i] = as.integer(rownames(unique_votes)[j])
      break
    }
  }
}
Y=table(myList)
gnames = as.character(senate$X[-dupix])
gnames[as.integer(names(which(Y>3)))] = paste("Group",1:length(which(Y>3)))
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
  knbs = 15
  kgamma = 0.5
  w = kernel_weights(X,kgamma)
  w = knn_weights(w,knbs,p)
  gamma_seq = seq(0.0,6, length.out=2000)
  tol = 1e-3
  nu = 1
} else if (scenario==2) {
  type = 2
  knbs = 15
  kgamma = 0.5
  w = kernel_weights(X,kgamma)
  w = knn_weights(w,knbs,p)
  gamma_seq = seq(0.0,15, length.out=1000)
  tol = 1e-3
  nu = 1
} else if (scenario==3) {
  type = 2
  kgamma = 0
  w = kernel_weights(X,kgamma)
  gamma_seq = seq(0.0,15, length.out=1000)
  tol = 1e-3
  nu = 1
}


ix = compactify_edges(w,p)$ix
nGamma = length(gamma_seq)
x = double(nGamma)
y = double(nGamma)
Z = X
  for (j in 1:nGamma) {
    Z = cbind(Z,ccp$UHx[[j]])    
}

svdZ = svd(Z)
pc = svdZ$u[,1:2,drop=FALSE]
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
#X_data$Party = senate$Party
data_plot = ggplot(data=df.paths,aes(x=x,y=y))
data_plot = data_plot + geom_path(aes(group=group),colour='grey60',alpha=0.5)
#data_plot = data_plot + geom_point(data=X_data,aes(x=x,y=y,colour=Party,shape=Party))
data_plot = data_plot + xlab('Principal Component 1') + ylab('Principal Component 2')
data_plot + theme_bw()