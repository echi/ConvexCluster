## These functions perform various variable updates for ADMM
source('prox.R')
library(compiler)
#######################################
## R versions                         #
#######################################

temp = function(X,U,gamma,w,ix,type) {
  data_fidelity = 0.5*(norm(as.matrix(X-U),'F')^2)
  center_distances = U[,ix[,1],drop=FALSE] - U[,ix[,2],drop=FALSE]
  if (type==1) {
    penalty = sum(w*apply(center_distances,2,FUN=function(x) {norm(as.matrix(x),'1')}))
  } else if (type==2) {
    penalty = sum(w*apply(center_distances,2,FUN=function(x) {norm(as.matrix(x),'F')}))
  } else if (type==3) {
    penalty = sum(w*apply(center_distances,2,FUN=function(x) {norm(as.matrix(x),'I')}))
  }
  return(data_fidelity + gamma*penalty)
}
loss_primalR = cmpfun(temp)

temp = function(X,V,Lambda,M1,M2,s1,s2,nu) {
  p = ncol(X)
  q = nrow(X)
  U = X
  xi = p*nu
  xi = xi/(1+xi)
  xbar = matrix(apply(X,1,mean),q,1)
  for (i in 1:p) {
    if (s1[i] > 0) {
      ix = M1[1:s1[i],i]
      U[,i] = U[,i] + apply(Lambda[,ix,drop=FALSE],1,sum) + nu*apply(V[,ix,drop=FALSE],1,sum)
    }
    if (s2[i] > 0) {
      ix = M2[1:s2[i],i]
      U[,i] = U[,i] - apply(Lambda[,ix,drop=FALSE],1,sum) - nu*apply(V[,ix,drop=FALSE],1,sum)
    }
  }
  U = (1-xi)*U + xi*xbar%*%t(matrix(1,p,1))
  return(U)
}
update_UR = cmpfun(temp)

temp = function(U,Lambda,w,gamma,nu,ix,type) {
  q = nrow(U)
  nK = ncol(Lambda)
  V = matrix(0,q,nK)
  for (kk in 1:nK) {
    i = ix[kk,1]
    j = ix[kk,2]
    z = U[,i] - U[,j] - (1/nu)*Lambda[,kk]
    V[,kk] = prox_R(z, w[kk]*gamma/nu,type)
  }
  return(V)
}
update_VR = cmpfun(temp)

temp = function(Lambda,U,V,nu,ix) {
  Z = U[,ix[,1],drop=FALSE] - U[,ix[,2],drop=FALSE] - V
  return(Lambda - nu*Z)
}
update_LambdaR = cmpfun(temp)

temp = function(U,V,ix) {
  return(max(abs(U[,ix[,1]] - U[,ix[,2]] - V)))
}
residual_primalR = cmpfun(temp)

temp = function(V,V_old,M1,M2,s1,s2,nu,p) {
  q = nrow(V)
  dual_residual = matrix(0,q,p)
  residual = -1
  for (i in 1:p) {
    if (s1[i] > 0) {
      ix = M1[1:s1[i],i]
      D = V[,ix,drop=FALSE]-V_old[,ix,drop=FALSE]
      dual_residual[,i] = dual_residual[,i] + apply(D,1,sum)
    }
    if (s2[i] > 0) {
      ix = M2[1:s2[i],i]
      D = V[,ix,drop=FALSE]-V_old[,ix,drop=FALSE]
      dual_residual[,i] = dual_residual[,i] - apply(D,1,sum)
    }
    residual = max(residual, max(abs(nu*dual_residual)))
  }
  return(residual)
}
residual_dualR = cmpfun(temp)

rm(temp)

#######################################
## R wrappers around Fortran versions #
#######################################
if (is.loaded('update_u')) {
  dyn.unload('admm.so')
}
dyn.load('admm.so')

update_UF = function(X,V,Lambda,M1,M2,s1,s2,nu) {
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
  storage.mode(V) = "double"
  storage.mode(Lambda) = "double"
  s1 = as.integer(s1)
  s2 = as.integer(s2)
  xbar = as.double(apply(X,1,mean),q,1)
  nu = as.double(nu)
#  subroutine update_U(X,Lambda,U,V,xbar,M1,M2,s1,s2,mix1,mix2,p,q,nK,nu) 
  sol = .Fortran('update_u',X=X,Lambda=Lambda,U=U,V=V,xbar=xbar,
                 M1=M1,M2=M2,s1=s1,s2=s2,mix1=mix1,mix2=mix2,p=p,q=q,nK=nK,nu)
  return(sol$U)
}

update_VF = function(U,Lambda,w,gamma,nu,ix,type) {
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
  type = as.integer(type)
#  subroutine update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK,type) 
  sol = .Fortran('update_v',U=U,Lambda=Lambda,V=V,w=w,gamma=gamma,nu=nu,ix=ix,q=q,p=p,nK=nK,type=type)
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

residual_primalF = function(U,V,ix) {
  dimU = dim(U)
  p = as.integer(dimU[2])
  q = as.integer(dimU[1])
  nK = as.integer(nrow(ix))
  storage.mode(U) = "double"
  storage.mode(V) = "double"
  storage.mode(ix) = "integer"
  sol = .Fortran('residual_primal',U=U,V=V,ix=ix,p=p,q=q,nK=nK,residual=double(1))
  return(sol$residual)
}

residual_dualF = function(V,V_old,M1,M2,s1,s2,nu,p) {
  p = as.integer(p)
  q = as.integer(nrow(V))
  nK = as.integer(ncol(V))
  storage.mode(V) = "double"
  storage.mode(V_old) = "double"
  storage.mode(M1) = "integer"
  storage.mode(M2) = "integer"
  mix1 = as.integer(nrow(M1))
  mix2 = as.integer(nrow(M2))
  storage.mode(s1) = "integer"
  storage.mode(s2) = "integer"
  nu = as.double(nu)
  sol = .Fortran('residual_dual',V=V,V_old=V_old,M1=M1,M2=M2,s1=s1,s2=s2,mix1=mix1,mix2=mix2,
                 p=p,q=q,nK=nK,nu=nu,residual=double(1))
  return(sol$residual)
}