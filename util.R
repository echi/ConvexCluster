## Functions to move back and forth between tuple and single indices
#
# p - number of points
# q - dimension of points
# k - index of kth pair. There are 1:q*(q-1)/2 pairs.

temp = function(X,mu) {
  p = ncol(X)
  k = 1
  w = matrix(0,p*(p-1)/2,1)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      w[k] = exp(-mu*norm(as.matrix(X[,i,drop=FALSE]-X[,j,drop=FALSE]))^2)
      k = k+1
    }
  }
  return(weights=w)
}
kernel_weights = cmpfun(temp)

temp = function(X,w,norm_type="2") {
  p = ncol(X)
  Xbar = as.matrix(apply(X,1,mean))
  Xc = t(scale(t(X),center=TRUE,scale=FALSE))
  k = 1
  denom = 0
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      denom = denom + w[k]*norm(X[,i,drop=FALSE] - X[,j,drop=FALSE],type=norm_type)
      k = k + 1
    }
  }
  return(gamma=(norm(Xc,'F')^2)/denom)
}
gamma_max = cmpfun(temp)

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


idx = vec2tri(k,p)
M = matrix(0,p,p)

for (kk in 1:nrow(idx)) {
  M[idx[kk,1],idx[kk,2]] = kk
}


## Update U

temp = function(X, Lambda,gamma1) {
  dimU = dim(X)
  q = dimU[1]
  p = dimU[2]
  U = matrix(0,q,p)
  for (ii in 1:p) {
    if (ii > 1) {
      ix = tri2vec(seq(1,ii-1,1),ii,p)
      L1 = as.matrix(apply(Lambda[,ix,drop=FALSE],1,sum))
    } else {
      L1 = 0
    }
    if (ii < p) {
      ix = tri2vec(ii,seq(ii+1,p,1),p)
      L2 = as.matrix(apply(Lambda[,ix,drop=FALSE],1,sum))
    } else {
      L2 = 0
    }
#    U[,ii] = X[,ii] - L1 + L2
    U[,ii] = prox_L1(X[,ii] - L1 + L2,gamma1)
  }
  return(U)
}
update_U = cmpfun(temp)

temp = function(X, Lambda,gamma1) {
  dimU = dim(X)
  q = dimU[1]
  p = dimU[2]
  U = matrix(0,q,p)
  ii = integer(1)
  U = foreach(ii=1:p, .combine=cbind) %dopar% {
    if (ii > 1) {
      ix = tri2vec(seq(1,ii-1,1),ii,p)
      L1 = as.matrix(apply(Lambda[,ix,drop=FALSE],1,sum))
    } else {
      L1 = 0
    }
    if (ii < p) {
      ix = tri2vec(ii,seq(ii+1,p,1),p)
      L2 = as.matrix(apply(Lambda[,ix,drop=FALSE],1,sum))
    } else {
      L2 = 0
    }
    prox_L1(X[,ii] - L1 + L2, gamma1)
  }
  return(U)
}
update_U = cmpfun(temp)

## Update_V

temp = function(U,Lambda,w,gamma,nu) {
  p = ncol(U)
  q = nrow(U)
  nK = p*(p-1)/2
  V = matrix(0,q,nK)
  ix = vec2tri(1:nK,p)
  kk = integer(1)
  for (kk in 1:nK) {
    i = ix[kk,1]
    j = ix[kk,2]
    z = U[,i,drop=FALSE] - U[,j,drop=FALSE] - (1/nu)*Lambda[,kk,drop=FALSE]
    V[,kk] = prox_L2(z,w[kk]*gamma/nu)
#    if (any(is.na(V[,kk]))) {
#      print(paste("kk=",kk))
#    }
#    V[,kk] = prox_L1(z,w[kk]*gamma/nu)   
  }
  return(V)
}
update_V = cmpfun(temp)

temp = function(U,Lambda,w,gamma,nu) {
  p = ncol(U)
  q = nrow(U)
  nK = p*(p-1)/2
  V = matrix(0,q,nK)
  ix = vec2tri(1:nK,p)
  kk = integer(1)
  V = foreach(kk=1:nK, .combine=cbind) %dopar% {
    i = ix[kk,1]
    j = ix[kk,2]
    z = U[,i,drop=FALSE] - U[,j,drop=FALSE] - (1/nu)*Lambda[,kk,drop=FALSE]
    prox_L2(z,w[kk]*gamma/nu)
  }
  return(V)
}
update_V = cmpfun(temp)

## Update_Lambda

temp = function(Lambda,U,V,nu) {
  nK = ncol(Lambda)
  p = ncol(U)
  ix = vec2tri(1:nK,p)
  i = ix[,1]
  j = ix[,2]
  Z = U[,i,drop=FALSE] - U[,j,drop=FALSE] - V
  Lambda = Lambda - nu*Z
#  for (kk in 1:nK) {
#    i = ix[kk,1]
#    j = ix[kk,2]
#    z = U[,i,drop=FALSE] - U[,j,drop=FALSE] - V[,kk,drop=FALSE]
#    Lambda[,kk] = Lambda[,kk] - nu*z
#  }
  return(Lambda)
}
update_Lambda = cmpfun(temp)


temp = function(X,U,gamma1,gamma,w) {
  p = ncol(X)
  nK = p*(p-1)/2
  ix = vec2tri(1:nK,p)
  i = ix[,1]
  j = ix[,2]
  penalty = apply(U[,i,drop=FALSE] - U[,j,drop=FALSE],2,FUN=function(x) {return(norm(as.matrix(x),'f'))})
  penalty1 = apply(U,2,FUN=function(x) {return(norm(as.matrix(x),'1'))})
  return(0.5*norm(as.matrix(X-U),'F')^2 + gamma*sum(w*penalty) + gamma1*sum(penalty1))
}
loss_primal = cmpfun(temp)

temp = function(X,Lambda) {
  dimX = dim(X)
  q = dimX[1]
  p = dimX[2]
  U = matrix(0,q,p)
  first_term = 0
  for (ii in 1:p) {
    if (ii > 1) {
      ix = tri2vec(seq(1,ii-1,1),ii,p)
      L1 = as.matrix(apply(Lambda[,ix,drop=FALSE],1,sum))
    } else {
      L1 = 0
    }
    if (ii < p) {
      ix = tri2vec(ii,seq(ii+1,p,1),p)
      L2 = as.matrix(apply(Lambda[,ix,drop=FALSE],1,sum))
    } else {
      L2 = 0
    }
    first_term = first_term + norm(as.matrix(L2 - L1),'2')^2
  }
  k = 1
  second_term = 0
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      second_term = second_term + t(X[,i,drop=FALSE] - X[,j,drop=FALSE])%*%Lambda[,k,drop=FALSE]
      k = k + 1
    }
  }
  return(-0.5*first_term - second_term)
}
loss_dual = cmpfun(temp)

temp = function(X,Lambda,gamma,gamma1=0,max_iter=1e3,tol=1e-4) {
  obj_p = double(max_iter)
  obj_d = double(max_iter)
  duality_gap = double(max_iter)
  alpha_old = 1
  Lambda_old = Lambda
  for (iter in 1:max_iter) {
    U = update_U(X,Lambda,gamma1)
    V = update_V(U,Lambda,w,gamma,nu)
    Lambda = update_Lambda(Lambda,U,V,nu)
    alpha = 0.5*(1 + sqrt(1 + 4*alpha_old^2))
    Lambda = Lambda + ((alpha_old-1)/alpha)*(Lambda-Lambda_old)

    obj_p[iter] = loss_primal(X,U,gamma1,gamma,w)
    obj_d[iter] = loss_dual(X,Lambda)
    duality_gap[iter] = obj_p[iter] - obj_d[iter]
  if (abs(duality_gap[iter]) < tol) {
    break
  }
#    if (norm(Lambda-Lambda_old)/(norm(Lambda_old)+1) < tol) {
#      break
#    }
    Lambda_old = Lambda
    alpha_old = alpha
  }
  obj_d = obj_d[1:iter]
  obj_p = obj_p[1:iter]
  return(list(U=U,V=V,Lambda=Lambda,primal=obj_p,dual=obj_d,iterations=iter))
}
convex_cluster = cmpfun(temp)

temp = function(X, Lambda, gamma_max) {
  
}
convex_cluster_path = cmpfun(temp)