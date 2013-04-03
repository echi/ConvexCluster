## Clusterpath preprocessing
library(compiler)

temp = function(i,j,p) {
  return(p*(i-1) - i*(i-1)/2 + j -i)
}
tri2vec = cmpfun(temp)

temp = function(k,p) {
  i = ceiling(0.5*(2*p-1 - sqrt((2*p-1)^2 - 8*k)))
  j = k - p*(i-1) + i*(i-1)/2 + i
  return(as.matrix(cbind(i,j)))
}
vec2tri = cmpfun(temp)

temp = function(w,p) {
  sizes1 = double(p)
  sizes2 = double(p)
  
  P = vec2tri(which(w > 0),p)
  M1 = matrix(0,nrow(P),p)
  M2 = matrix(0,nrow(P),p)
  
  for (i in 1:p) {
    group1 = which(P[,1] == i)
    sizes1[i] = length(group1)
    if (sizes1[i] > 0) {
      M1[1:sizes1[i],i] = group1
    }
    group2 = which(P[,2] == i)
    sizes2[i] = length(group2)
    if (sizes2[i] > 0) {
      M2[1:sizes2[i],i] = group2
    }
  }
  
  M1 = M1[1:max(sizes1),,drop=FALSE]
  M2 = M2[1:max(sizes2),,drop=FALSE]
  
  return(list(ix=P,M1=M1,M2=M2,s1=sizes1,s2=sizes2))
}
compactify_edges = cmpfun(temp)

temp = function(X,mu) {
  p = ncol(X)
  k = 1
  w = matrix(0,p*(p-1)/2,1)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      w[k] = exp(-mu*norm(as.matrix(X[,i,drop=FALSE]-X[,j,drop=FALSE]),'f')^2)
      k = k+1
    }
  }
  return(weights=w)
}
kernel_weights = cmpfun(temp)

temp = function(w,k,p) {
  i = 1
  neighbors = tri2vec(i,(i+1):p,p)
  keep = neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
  for (i in 2:(p-1)) {
    group_A = tri2vec(i,(i+1):p,p)
    group_B = tri2vec(1:(i-1),i,p)
    neighbors = c(group_A,group_B)
    knn = neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
    keep = union(knn,keep)
  }
  i = p
  neighbors = tri2vec(1:(i-1),i,p)
  knn = neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
  keep = union(knn,keep)
  w[-keep] = 0
  return(w)
}
knn_weights = cmpfun(temp)
rm(temp)