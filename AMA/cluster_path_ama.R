dyn.load('ama.so')
library(compiler)
library(igraph)

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

update_VF = function(U,Lambda,w,gamma,nu,ix) {
  dimU = dim(U)
  p = as.integer(dimU[2])
  q = as.integer(dimU[1])
  nK = as.integer(p*(p-1)/2)
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

#subroutine update_Ut(X,Lambda,U,M1,M2,s1,s2,mix1,mix2,p,q,nK)
update_UF = function(X,Lambda,w) {
  dimU = dim(X)
  p = as.integer(dimU[2])
  q = as.integer(dimU[1])
  edge_info = compactify_edges(w,p)
  L = Lambda[,which(w > 0),drop=FALSE]
  
  nK = as.integer(ncol(L))
  U = matrix(0,q,p)  
  storage.mode(U) = "double"
  storage.mode(L) = "double"
  
  M1 = edge_info$M1
  M2 = edge_info$M2
  s1 = as.integer(edge_info$s1)
  s2 = as.integer(edge_info$s2)
  mix1 = as.integer(nrow(M1))
  mix2 = as.integer(nrow(M2))
  storage.mode(M1) = "integer"
  storage.mode(M2) = "integer"
  sol = .Fortran('update_u',X=X,Lambda=L,U=U,M1=M1,M2=M2,s1=s1,s2=s2,mix1=mix1,mix2=mix2,p=p,q=q,nK=nK)
  return(sol$U)
}

loss_dualF = function(X,Lambda,w) {
  dimX = dim(X)
  p = as.integer(dimX[2])
  q = as.integer(dimX[1])
  edge_info = compactify_edges(w,p)
  L = Lambda[,which(w > 0),drop=FALSE]
  
  ix = vec2tri(which(w > 0),p)
  storage.mode(ix) = "integer"
  nK = as.integer(ncol(L))
  storage.mode(X) = "double"
  storage.mode(L) = "double"
  M1 = edge_info$M1
  M2 = edge_info$M2
  s1 = as.integer(edge_info$s1)
  s2 = as.integer(edge_info$s2)
  mix1 = as.integer(nrow(M1))
  mix2 = as.integer(nrow(M2))
  storage.mode(M1) = "integer"
  storage.mode(M2) = "integer"
  sol = .Fortran('loss_dualt',X=X,Lambda=L,ix=ix,p=p,q=q,nK=nK,s1=s1,s2=s2,M1=M1,M2=M2,mix1=mix1,mix2=mix2,output=double(1))
  return(sol$output)
}

loss_primalF = function(X,U,gamma1,gamma,w,ix) {
  dimU = dim(U)
  p = as.integer(dimU[2])
  q = as.integer(dimU[1])
  nK = as.integer(nrow(ix))
  gamma = as.double(gamma)  
  storage.mode(X) = "double"
  storage.mode(U) = "double"
  w = as.double(w)
  storage.mode(ix) = "integer"
  sol = .Fortran('loss_primal',X=X,U=U,gamma=gamma,ix=ix,p=p,q=q,nK=nK,w=w,output=double(1))
  return(sol$output)
}

convex_clusterF = function(X,Lambda,ix,w,gamma,s1,s2,M1,M2,nu=1,eta=2,max_iter=1e5,tol=1e-4) {
  dimX = dim(X)
  p = as.integer(dimX[2])
  q = as.integer(dimX[1])
  nK = as.integer(nrow(ix))
  gamma = as.double(gamma)
  max_iter = as.integer(max_iter)
  tol = as.double(tol)
  storage.mode(X) = "double"
  w = as.double(w)
  storage.mode(Lambda) = "double"
  storage.mode(ix) = "integer"
  U = matrix(0,q,p)
  storage.mode(U) = "double"
  V = matrix(0,q,nK)
  storage.mode(V) = "double"
  nu = as.double(nu)
  eta = as.double(eta)
  primal = matrix(0,max_iter,1)
  primal = as.double(primal)
  dual = matrix(0,max_iter,1)
  dual = as.double(dual)
  
  s1 = as.integer(s1)
  s2 = as.integer(s2)
  storage.mode(M1) = "integer"
  storage.mode(M2) = "integer"
  mix1 = as.integer(nrow(M1))
  mix2 = as.integer(nrow(M2))
  type = as.integer(2)
  
  sol = .Fortran('convex_cluster',X=X,Lambda=Lambda,U=U,V=V,q=q,p=p,nK=nK,
                 ix=ix,w=w,gamma=gamma,nu=nu,eta=eta,
                 s1=s1,s2=s2,M1=M1,M2=M2,mix1=mix1,mix2=mix2,
                 primal=primal,dual=dual,
                 max_iter=max_iter,iter=integer(1),tol=tol,type=type)
  return(list(U=sol$U,V=sol$V,Lambda=sol$Lambda,primal=sol$primal[1:(sol$iter)],
              dual=sol$dual[1:(sol$iter)],iterations=sol$iter,nu=sol$nu))  
}

temp = function(X, w, gamma_seq,nu0=1,tol=1e-4) {
  nGamma = length(gamma_seq)
  p = ncol(X)
  q = nrow(X)
  edge_info = compactify_edges(w,p)
#  nK = length(w)
  nK = length(which(w > 0))
  Lambda = matrix(0,q,nK)
  list_U = vector(mode="list",length=nGamma)
  list_V = vector(mode="list",length=nGamma)
  list_Lambda = vector(mode="list",length=nGamma)
  nu = nu0
#  ix = vec2tri(1:nK,p)
  ix = edge_info$ix

  M1 = edge_info$M1
  M2 = edge_info$M2
  s1 = edge_info$s1
  s2 = edge_info$s2
  
  for (ig in 1:nGamma) {
    gamma = gamma_seq[ig]      
    cc = convex_clusterF(X,Lambda,ix,w[which(w>0)],gamma,s1=s1,s2=s2,M1=M1,M2=M2,nu=nu,tol=tol)
      
    nu = cc$nu
    Lambda = cc$Lambda
    list_U[[ig]] = cc$U
    list_V[[ig]] = cc$V
    list_Lambda[[ig]] = Lambda 
    print(paste0("iters: ", cc$iter,"| primal: ", cc$primal[cc$iter],
                 "| dual: ",cc$dual[cc$iter],
                 "| gap: ", cc$primal[cc$iter]-cc$dual[cc$iter]))
    print(paste("Completed gamma",ig))
    if (norm(cc$V,'1')==0) {
      print('Single cluster')
      break
    }
  }
  return(list(UHx=list_U,VHx=list_V,LambdaHx=list_Lambda,nGamma=ig))
}
convex_cluster_path = cmpfun(temp)
  
temp = function(V,p,ix,tol=1e-10) {
  cList = which(apply(V,2,FUN=function(x) {norm(as.matrix(x),'1')}) <= tol)
  if (length(cList) > 0) {
    edgeList = ix[cList,,drop=FALSE]
    g = graph(t(edgeList), p, directed=FALSE)
    return(clusters(g))
  } else {
    return(list(membership=1:p,csize=double(p)+1,no=p))
  }
}
get_clusters = cmpfun(temp)

temp = function(VList,p,ix) {
  ngamma = length(VList)
  nClusters = integer(ngamma)
  for (i in 1:ngamma) {
    V = VList[[i]]
    if (length(V) > 0) {
      nClusters[i] = get_clusters(VList[[i]],p,ix)$no
#      print(paste('i = ',i))
    } else {
      break
    }
  }
  return(nClusters[1:i-1])
}
number_clusters = cmpfun(temp)

rm(list="temp")
