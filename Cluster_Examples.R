## Classic Clustering

rm(list=ls())
library(ggplot2)

make_simple_clusters = function(p=100, sigma=0.1, seed=12345) {
  mu1 = matrix(0,1,2)
  mu2 = matrix(c(0.75,1),1,2)
  mu3 = matrix(c(-0.1,0.25),1,2)
  mu4 = matrix(c(0.25,-0.1),1,2)
  mu5 = matrix(c(0.8,0.4),1,2)
  
  
  n = floor(p/5)
  X = matrix(0,2,p)
  X[,1:n] = matrix(mu1,2,n) + sigma*matrix(rnorm(2*n),2,n)
  X[,(n+1):(2*n)] = matrix(mu2,2,n) +sigma* matrix(rnorm(2*n),2,n)
  X[,(2*n+1):(3*n)] = matrix(mu3,2,n) +sigma* matrix(rnorm(2*n),2,n)
  X[,(3*n+1):(4*n)] = matrix(mu4,2,n) +sigma* matrix(rnorm(2*n),2,n)
  X[,(4*n+1):(5*n)] = matrix(mu5,2,n) +sigma* matrix(rnorm(2*n),2,n)  
  return(X)
}

X = as.data.frame(t(make_simple_clusters(sigma=0.03, seed=1234)))
colnames(X) = c("x","y")
q = ggplot(data=X, aes(x=x,y=y))
q = q + geom_point(size=4)

set.seed(123)
theta = runif(n=5,min=0,max=2*pi)
r = runif(n=5,min=0,max=0.4)
centers = data.frame(x=r*cos(theta)+0.5,y=r*sin(theta)+0.5)

Lloyd_sol = kmeans(X,centers=centers, algorithm="Lloyd")
X_Lloyd = X
X_Lloyd$cluster = factor(Lloyd_sol$cluster)
q = ggplot(data=X_Lloyd,aes(x=x,y=y))
q = q + geom_point(aes(shape=cluster, colour=cluster), size=4) + scale_colour_brewer(palette="Set1")
q = q + geom_point(data=centers, aes(x=x,y=y), colour='black', size=4) + geom_point(data=centers, colour="grey100", size = 2)
q + theme_bw() + theme(legend.position = "none") 
ggsave(filename="../Manuscript/kmeans.pdf",width=8,height=8)

## Perform convex clustering
X = t(X)
q = nrow(X)
p = ncol(X)
nK = p*(p-1)/2
k = seq(1,nK,1)
w = kernel_weights(X,1)
gamma = 0.01
nu = (2/nK)*0.99

Lambda = matrix(0,q,nK)

system.time({sol = convex_cluster(X,Lambda,w,nu,gamma,tol=1e-5,max_iter=1e4)})

## Visualize the loss
library(reshape2)
loss_data = data.frame(Iteration=1:sol$iterations, primal=sol$primal, dual=sol$dual)
loss_data = melt(loss_data, id.vars=c('Iteration'))

loss_plot = ggplot(data=loss_data, aes(x=Iteration, y=value))
loss_plot = loss_plot + geom_point(aes(colour=variable, shape=variable))
loss_plot + theme_bw()

## Visualize the data + shrunken centroids
centers = data.frame(x=sol$U[1,],y=sol$U[2,])
X_data = as.data.frame(t(X))
colnames(X_data) = c("x","y")
data_plot = ggplot(data=X_data,aes(x=x,y=y))
data_plot = data_plot + geom_point(size=4)
data_plot = data_plot + geom_point(data=centers, aes(x=x,y=y), size=4, colour='red')
data_plot + theme_bw() + theme(legend.position = "none") 

## Example: Sparse centers

X1 = make_2d_sparse_clusters(shift=-15,n=10,seed=123)
X2 = make_2d_sparse_clusters(shift=15,n=10,seed=1234)

X = as.data.frame(rbind(X1,X2))
rm(list=c("X1","X2"))
colnames(X) = c("x","y")
kmeans_sol = kmeans(X,2)

X_kmeans = X
X_kmeans$cluster = factor(kmeans_sol$cluster)
q = ggplot(data=X_kmeans,aes(x=x,y=y))
q = q + geom_point(aes(shape=cluster, colour=cluster), size=4) + scale_colour_brewer(palette="Set1")
q + theme_bw() + theme(legend.position = "none") 
