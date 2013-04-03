## These functions calculate the primal and dual loss for AMA
library(compiler)
#######################################
## R versions                         #
#######################################
temp = function(X,U,gamma,w,ix) {
  data_fidelity = 0.5*(norm(as.matrix(X-U),'f')^2)
  center_distances = U[,ix[,1],drop=FALSE] - U[,ix[,2],drop=FALSE]
  penalty = sum(w*apply(center_distances,2,FUN=function(x) {norm(as.matrix(x),'2')}))
  return(data_fidelity + gamma*penalty)
}
loss_primalR = cmpfun(temp)

temp = function(X,Lambda,ix,M1,M2,s1,s2) {
  q = nrow(X)
  p = ncol(X)
  nK = ncol(Lambda)
  L = matrix(0,q,p)
  for (i in 1:p) {
    if (s1[i] > 0) {
      jx = M1[1:s1[i],i]
      L[,i] = apply(Lambda[,jx,drop=FALSE],1,sum)
    }
    if (s2[i] > 0) {
      jx = M2[1:s2[i],i]
      L[,i] = L[,i] - apply(Lambda[,jx,drop=FALSE],1,sum)
    }
  }
  first_term = 0.5*norm(L,'f')^2
  second_term = sum((X[,ix[,1]] - X[,ix[,2]])*Lambda)
  return(-first_term-second_term)
}
loss_dualR = cmpfun(temp)

rm(temp)

#######################################
## R wrappers around Fortran versions #
#######################################
if (is.loaded('loss_primal')) {
  dyn.unload('Fortran/ama.so')
}
dyn.load('Fortran/ama.so')

loss_primalF = function(X,U,gamma,w,ix,type) {
  dimU = dim(U)
  p = as.integer(dimU[2])
  q = as.integer(dimU[1])
  nK = as.integer(nrow(ix))
  gamma = as.double(gamma)  
  storage.mode(X) = "double"
  storage.mode(U) = "double"
  w = as.double(w)
  storage.mode(ix) = "integer"
  type = as.integer(type)
  sol = .Fortran('loss_primal',X=X,U=U,gamma=gamma,ix=ix,p=p,q=q,nK=nK,w=w,output=double(1),type=type)
  return(sol$output)
}

loss_dualF = function(X,Lambda,ix,M1,M2,s1,s2) {
  dimX = dim(X)
  p = as.integer(dimX[2])
  q = as.integer(dimX[1])
  nK = as.integer(ncol(Lambda))
  storage.mode(X) = "double"
  storage.mode(Lambda) = "double"
  storage.mode(ix) = "integer"
  storage.mode(M1) = "integer"
  storage.mode(M2) = "integer"
  s1 = as.integer(s1)
  s2 = as.integer(s2)
  mix1 = as.integer(nrow(M1))
  mix2 = as.integer(nrow(M2))
  sol = .Fortran('loss_dual',X=X,Lambda=Lambda,ix=ix,p=p,q=q,nK=nK,s1=s1,s2=s2,M1=M1,M2=M2,mix1=mix1,mix2=mix2,output=double(1))
  return(sol$output)
}