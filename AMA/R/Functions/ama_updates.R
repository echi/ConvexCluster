## These functions perform various variable updates for AMA
source('prox.R')
library(compiler)
#######################################
## R versions                         #
#######################################
temp = function(X,Lambda,M1,M2,s1,s2) {
  p = ncol(X)
  U = X
  for (i in 1:p) {
    if (s1[i] > 0) {
      ix = M1[1:s1[i],i]
      U[,i] = U[,i] + apply(Lambda[,ix,drop=FALSE],1,sum)
    }
    if (s2[i] > 0) {
      ix = M2[1:s2[i],i]
      U[,i] = U[,i] - apply(Lambda[,ix,drop=FALSE],1,sum)
    }
  }
  return(U)
}
update_UR = cmpfun(temp)
  
temp = function(U,Lambda,w,gamma,nu,ix) {
  q = nrow(U)
  nK = ncol(Lambda)
  V = matrix(0,q,nK)
  for (kk in 1:nK) {
    i = ix[kk,1]
    j = ix[kk,2]
    z = U[,i] - U[,j] - (1/nu)*Lambda[,kk]
    V[,kk] = prox_L2(z, w[kk]*gamma/nu)
  }
  return(V)
}
update_VR = cmpfun(temp)

temp = function(Lambda,U,V,nu,ix) {
  Z = U[,ix[,1],drop=FALSE] - U[,ix[,2],drop=FALSE] - V
  return(Lambda - nu*Z)
}
update_LambdaR = cmpfun(temp)

temp = function(Lambda,U,w,gamma,nu,ix) {
  V = update_VR(U,Lambda,w,gamma,nu,ix)
  return(update_LambdaR(Lambda,U,V,nu,ix))
}
update_Lambda_combinedR = cmpfun(temp)

rm(temp)

#######################################
## R wrappers around Fortran versions #
#######################################
if (is.loaded('update_u')) {
  dyn.unload('ama.so')
}
dyn.load('ama.so')

update_UF = function(X,Lambda,M1,M2,s1,s2) {
  p = as.integer(ncol(X))
  q = as.integer(nrow(X))
  nK = as.integer(ncol(Lambda))
  mix1 = as.integer(nrow(M1))
  mix2 = as.integer(nrow(M2))
  storage.mode(M1) = "integer"
  storage.mode(M2) = "integer"
  U = matrix(0,q,p)
  storage.mode(X) = "double"
  storage.mode(U) = "double"
  storage.mode(Lambda) = "double"
  s1 = as.integer(s1)
  s2 = as.integer(s2)
  
  sol = .Fortran('update_u',X=X,Lambda=Lambda,U=U,M1=M1,M2=M2,s1=s1,s2=s2,mix1=mix1,mix2=mix2,p=p,q=q,nK=nK)
  return(sol$U)
}

update_VF = function(U,Lambda,w,gamma,nu,ix) {
  dimU = dim(U)
  p = as.integer(dimU[2])
  q = as.integer(dimU[1])
  nK = as.integer(ncol(Lambda))
  nu = as.double(nu)
  gamma = as.double(gamma)
  storage.mode(U) = "double"
  storage.mode(Lambda) = "double"
  V = matrix(0,q,nK)
  storage.mode(V) = "double"
  storage.mode(ix) = "integer"
  sol = .Fortran('update_v',U=U,Lambda=Lambda,V=V,w=w,gamma=gamma,nu=nu,ix=ix,q=q,p=p,nK=nK)
  return(sol$V)
}

update_LambdaF = function(Lambda,U,V,nu,ix) {
  dimU = dim(U)
  p = as.integer(dimU[2])
  q = as.integer(dimU[1])
  nK = as.integer(ncol(Lambda))
  nu = as.double(nu)
  storage.mode(U) = "double"
  storage.mode(Lambda) = "double"
  storage.mode(V) = "double"
  storage.mode(ix) = "integer"
  sol = .Fortran('update_lambda',Lambda=Lambda,U=U,V=V,nu=nu,ix=ix,q=q,p=p,nK=nK)
  return(sol$Lambda)
}

#subroutine update_LambdaDP(Lambda,U,nu,gamma,ix,q,p,nK,w)

update_Lambda_combinedF = function(Lambda,U,w,gamma,nu,ix) {
  dimU = dim(U)
  p = as.integer(dimU[2])
  q = as.integer(dimU[1])
  nK = as.integer(ncol(Lambda))
  nu = as.double(nu)
  gamma = as.double(gamma)
  w = as.double(w)
  storage.mode(U) = "double"
  storage.mode(Lambda) = "double"
  storage.mode(ix) = "integer"
  sol = .Fortran('update_lambdadp',Lambda=Lambda,U=U,nu=nu,gamma=gamma,ix=ix,q=q,p=p,nK=nK,w=w)
  return(sol$Lambda)
}
