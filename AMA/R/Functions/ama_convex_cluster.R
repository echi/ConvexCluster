## These functions perform various versions of convex clustering via AMA
library(compiler)
library(igraph)

#######################################
## R wrappers around Fortran versions #
#######################################
source('cluster_path_preprocess.R')

if (is.loaded('convex_cluster')) {
  dyn.unload('ama.so')
}
dyn.load('ama.so')

convex_cluster = function(X,Lambda,ix,M1,M2,s1,s2,w,gamma,nu,eta=2,type=2,max_iter=1e2,tol=1e-4) {
  q = as.integer(nrow(X))
  p = as.integer(ncol(X))
  nK = as.integer(ncol(Lambda))
  mix1 = as.integer(nrow(M1))
  mix2 = as.integer(nrow(M2))
  storage.mode(X) = "double"
  storage.mode(Lambda) = "double"
  U = matrix(0,q,p)
  storage.mode(U) = "double"
  V = matrix(0,q,nK)
  storage.mode(V) = "double"
  storage.mode(ix) = "integer"
  storage.mode(M1) = "integer"
  storage.mode(M2) = "integer"
  s1 = as.integer(s1)
  s2 = as.integer(s2)
  w = as.double(w)
  gamma = as.double(gamma)
  nu = as.double(nu)
  eta = as.double(eta)
  type = as.integer(type)
  max_iter = as.integer(max_iter)
  tol = as.double(tol)
  primal = double(max_iter)
  dual = double(max_iter)
  sol = .Fortran('convex_cluster',X=X,Lambda=Lambda,U=U,V=V,q=q,p=p,nK=nK,ix=ix,w=w,gamma=gamma,nu=nu,
                 eta=eta,s1=s1,s2=s2,M1=M1,M2=M2,mix1=mix1,mix2=mix2,primal=primal,dual=dual,
                 max_iter=max_iter,iter=integer(1),tol=tol,type=type)
  return(list(U=sol$U,V=sol$V,Lambda=sol$Lambda,nu=sol$nu,
              primal=sol$primal[1:sol$iter],dual=sol$dual[1:sol$iter],iter=sol$iter))
}

convex_cluster_fista_fixed = function(X,Lambda,ix,M1,M2,s1,s2,w,gamma,nu,eta=2,type=2,max_iter=1e2,tol=1e-4) {
  q = as.integer(nrow(X))
  p = as.integer(ncol(X))
  nK = as.integer(ncol(Lambda))
  mix1 = as.integer(nrow(M1))
  mix2 = as.integer(nrow(M2))
  storage.mode(X) = "double"
  storage.mode(Lambda) = "double"
  U = matrix(0,q,p)
  storage.mode(U) = "double"
  V = matrix(0,q,nK)
  storage.mode(V) = "double"
  storage.mode(ix) = "integer"
  storage.mode(M1) = "integer"
  storage.mode(M2) = "integer"
  s1 = as.integer(s1)
  s2 = as.integer(s2)
  w = as.double(w)
  gamma = as.double(gamma)
  nu = as.double(nu)
  eta = as.double(eta)
  type = as.integer(type)
  max_iter = as.integer(max_iter)
  tol = as.double(tol)
  primal = double(max_iter)
  dual = double(max_iter)
  sol = .Fortran('convex_cluster_fista_fixed',X=X,Lambda=Lambda,U=U,V=V,q=q,p=p,nK=nK,ix=ix,w=w,gamma=gamma,nu=nu,
                 eta=eta,s1=s1,s2=s2,M1=M1,M2=M2,mix1=mix1,mix2=mix2,primal=primal,dual=dual,
                 max_iter=max_iter,iter=integer(1),tol=tol,type=type)
  return(list(U=sol$U,V=sol$V,Lambda=sol$Lambda,nu=sol$nu,
              primal=sol$primal[1:sol$iter],dual=sol$dual[1:sol$iter],iter=sol$iter))
}

convex_cluster_backtrack = function(X,Lambda,ix,M1,M2,s1,s2,w,gamma,nu,eta=2,type=2,max_iter=1e2,tol=1e-4) {
  q = as.integer(nrow(X))
  p = as.integer(ncol(X))
  nK = as.integer(ncol(Lambda))
  mix1 = as.integer(nrow(M1))
  mix2 = as.integer(nrow(M2))
  storage.mode(X) = "double"
  storage.mode(Lambda) = "double"
  U = matrix(0,q,p)
  storage.mode(U) = "double"
  V = matrix(0,q,nK)
  storage.mode(V) = "double"
  storage.mode(ix) = "integer"
  storage.mode(M1) = "integer"
  storage.mode(M2) = "integer"
  s1 = as.integer(s1)
  s2 = as.integer(s2)
  w = as.double(w)
  gamma = as.double(gamma)
  nu = as.double(nu)
  eta = as.double(eta)
  type = as.integer(type)
  max_iter = as.integer(max_iter)
  tol = as.double(tol)
  primal = double(max_iter)
  dual = double(max_iter)
  sol = .Fortran('convex_cluster_backtrack',X=X,Lambda=Lambda,U=U,V=V,q=q,p=p,nK=nK,ix=ix,w=w,gamma=gamma,nu=nu,
                 eta=eta,s1=s1,s2=s2,M1=M1,M2=M2,mix1=mix1,mix2=mix2,primal=primal,dual=dual,
                 max_iter=max_iter,iter=integer(1),tol=tol,type=type)
  return(list(U=sol$U,V=sol$V,Lambda=sol$Lambda,nu=sol$nu,
              primal=sol$primal[1:sol$iter],dual=sol$dual[1:sol$iter],iter=sol$iter))
}

convex_cluster_fista_backtrack = function(X,Lambda,ix,M1,M2,s1,s2,w,gamma,nu,eta=2,type=2,max_iter=1e2,tol=1e-4) {
  q = as.integer(nrow(X))
  p = as.integer(ncol(X))
  nK = as.integer(ncol(Lambda))
  mix1 = as.integer(nrow(M1))
  mix2 = as.integer(nrow(M2))
  storage.mode(X) = "double"
  storage.mode(Lambda) = "double"
  U = matrix(0,q,p)
  storage.mode(U) = "double"
  V = matrix(0,q,nK)
  storage.mode(V) = "double"
  storage.mode(ix) = "integer"
  storage.mode(M1) = "integer"
  storage.mode(M2) = "integer"
  s1 = as.integer(s1)
  s2 = as.integer(s2)
  w = as.double(w)
  gamma = as.double(gamma)
  nu = as.double(nu)
  eta = as.double(eta)
  type = as.integer(type)
  max_iter = as.integer(max_iter)
  tol = as.double(tol)
  primal = double(max_iter)
  dual = double(max_iter)
  sol = .Fortran('convex_cluster_fista_backtrack',X=X,Lambda=Lambda,U=U,V=V,q=q,p=p,nK=nK,ix=ix,w=w,gamma=gamma,nu=nu,
                 eta=eta,s1=s1,s2=s2,M1=M1,M2=M2,mix1=mix1,mix2=mix2,primal=primal,dual=dual,
                 max_iter=max_iter,iter=integer(1),tol=tol,type=type)
  return(list(U=sol$U,V=sol$V,Lambda=sol$Lambda,nu=sol$nu,
              primal=sol$primal[1:sol$iter],dual=sol$dual[1:sol$iter],iter=sol$iter))
}

temp = function(X, w, gamma_seq,nu0=10,tol=1e-4,max_iter=1e6,type=2) {
  nGamma = length(gamma_seq)
  p = ncol(X)
  q = nrow(X)
  edge_info = compactify_edges(w,p)
  nK = length(which(w > 0))
  Lambda = matrix(0,q,nK)
  list_U = vector(mode="list",length=nGamma)
  list_V = vector(mode="list",length=nGamma)
  list_Lambda = vector(mode="list",length=nGamma)
  nu = nu0
  ix = edge_info$ix  
  M1 = edge_info$M1
  M2 = edge_info$M2
  s1 = edge_info$s1
  s2 = edge_info$s2
  iter_vec = integer(nGamma)
  for (ig in 1:nGamma) {
    gamma = gamma_seq[ig]      
    cc = convex_cluster_backtrack(X,Lambda,ix,M1,M2,s1,s2,w[w>0],gamma,nu,max_iter=1e4,tol=tol,type=type)
#    cc = convex_cluster_fista_backtrack(X,Lambda,ix,M1,M2,s1,s2,w[w>0],gamma,nu,max_iter=1e4,tol=tol,type=type)    
    iter_vec[ig] = cc$iter
    nu = cc$nu
    Lambda = cc$Lambda
    list_U[[ig]] = cc$U
    list_V[[ig]] = cc$V
    list_Lambda[[ig]] = Lambda 
#        print(paste0("iters: ", cc$iter,"| primal: ", cc$primal[cc$iter],
#                     "| dual: ",cc$dual[cc$iter],
#                     "| gap: ", cc$primal[cc$iter]-cc$dual[cc$iter]))
#        print(paste("Completed gamma",ig))
#        if (norm(cc$V,'1')==0) {
#          print('Single cluster')
#          break
#        }
  }
  return(list(UHx=list_U,VHx=list_V,LambdaHx=list_Lambda,nGamma=ig,iters=iter_vec))
}
convex_cluster_path_backtrack = cmpfun(temp)

temp = function(X, w, gamma_seq,nu0=10,tol=1e-4,max_iter=1e6,type=2) {
  nGamma = length(gamma_seq)
  p = ncol(X)
  q = nrow(X)
  edge_info = compactify_edges(w,p)
  nK = length(which(w > 0))
  Lambda = matrix(0,q,nK)
  list_U = vector(mode="list",length=nGamma)
  list_V = vector(mode="list",length=nGamma)
  list_Lambda = vector(mode="list",length=nGamma)
  nu = nu0
  ix = edge_info$ix  
  M1 = edge_info$M1
  M2 = edge_info$M2
  s1 = edge_info$s1
  s2 = edge_info$s2
  iter_vec = integer(nGamma)
  for (ig in 1:nGamma) {
    gamma = gamma_seq[ig]      
    cc = convex_cluster_fista_backtrack(X,Lambda,ix,M1,M2,s1,s2,w[w>0],gamma,nu,max_iter=1e4,tol=tol,type=type)    
    iter_vec[ig] = cc$iter
    nu = cc$nu
    Lambda = cc$Lambda
    list_U[[ig]] = cc$U
    list_V[[ig]] = cc$V
    list_Lambda[[ig]] = Lambda 
    print(paste0("iters: ", cc$iter,"| primal: ", cc$primal[cc$iter],
                 "| dual: ",cc$dual[cc$iter],
                 "| gap: ", cc$primal[cc$iter]-cc$dual[cc$iter]))
    print(paste("Completed gamma",ig))
    #        if (norm(cc$V,'1')==0) {
    #          print('Single cluster')
    #          break
    #        }
  }
  return(list(UHx=list_U,VHx=list_V,LambdaHx=list_Lambda,nGamma=ig,iters=iter_vec))
}
convex_cluster_path_fista_backtrack = cmpfun(temp)

temp = function(X, w, gamma_seq,nu0=10,tol=1e-4,max_iter=1e6,type=2) {
  nGamma = length(gamma_seq)
  p = ncol(X)
  q = nrow(X)
  edge_info = compactify_edges(w,p)
  nK = length(which(w > 0))
  Lambda = matrix(0,q,nK)
  list_U = vector(mode="list",length=nGamma)
  list_V = vector(mode="list",length=nGamma)
  list_Lambda = vector(mode="list",length=nGamma)
  nu = nu0
  ix = edge_info$ix  
  M1 = edge_info$M1
  M2 = edge_info$M2
  s1 = edge_info$s1
  s2 = edge_info$s2
  iter_vec = integer(nGamma)  
  for (ig in 1:nGamma) {
    gamma = gamma_seq[ig]      
    cc = convex_cluster_fista_backtrack(X,Lambda,ix,M1,M2,s1,s2,w[w>0],gamma,nu,max_iter=max_iter,tol=tol,type=type)
#    cc = convex_cluster_backtrack(X,Lambda,ix,M1,M2,s1,s2,w[w>0],gamma,nu,max_iter=1e4,tol=tol,type=type)
    iter_vec[ig] = cc$iter    
    nu = cc$nu
    Lambda = cc$Lambda
    list_U[[ig]] = cc$U
    list_V[[ig]] = cc$V
    list_Lambda[[ig]] = Lambda 
    print(paste0("iters: ", cc$iter,"| primal: ", cc$primal[cc$iter],
                 "| dual: ",cc$dual[cc$iter],
                 "| gap: ", cc$primal[cc$iter]-cc$dual[cc$iter]))
    print(paste("Completed gamma",ig))
#    if (norm(cc$V,'1')==0) {
#      print('Single cluster')
#      break
#    }
  }
  return(list(UHx=list_U,VHx=list_V,LambdaHx=list_Lambda,nGamma=ig,iters=iter_vec))
}
convex_cluster_path = cmpfun(temp)

temp = function(X, w, gamma_seq,nu0=10,tol=1e-4,max_iter=1e6,type=2) {
  nGamma = length(gamma_seq)
  p = ncol(X)
  q = nrow(X)
  edge_info = compactify_edges(w,p)
  nK = length(which(w > 0))
  Lambda = matrix(0,q,nK)
  list_U = vector(mode="list",length=nGamma)
  list_V = vector(mode="list",length=nGamma)
  list_Lambda = vector(mode="list",length=nGamma)
  nu = nu0
  ix = edge_info$ix  
  M1 = edge_info$M1
  M2 = edge_info$M2
  s1 = edge_info$s1
  s2 = edge_info$s2
  
  for (ig in 1:nGamma) {
    gamma = gamma_seq[ig]      
#    cc = convex_cluster_fista_backtrack(X,Lambda,ix,M1,M2,s1,s2,w[w>0],gamma,nu,max_iter=1e5,tol=tol)
#    cc = convex_cluster(X,Lambda,ix,M1,M2,s1,s2,w,gamma,nu,eta=2,type=2,max_iter=1e5,tol=tol)
    cc = convex_cluster_fista_fixed(X,Lambda,ix,M1,M2,s1,s2,w,gamma,nu,eta=2,type=type,max_iter=max_iter,tol=tol)
    
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
convex_cluster_path_fixed = cmpfun(temp)

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

temp = function(Ulist,Vlist,ix,nU) {
  dimU = dim(as.matrix(Ulist[[1]]))
  q = dimU[1]
  p = dimU[2]
  U_out = vector(mode="list",length=nU)
  for (j in 1:nU) {
    U = as.matrix(Ulist[[j]])
    V = as.matrix(Vlist[[j]])
    cluster_info = get_clusters(V,p,ix)
    nClusters = cluster_info$no
    U_out[[j]] = matrix(0,q,p)
    for (i in 1:nClusters) {
      ix_cluster = which(cluster_info$membership == i)
      U_out[[j]][,ix_cluster] = apply(U[,ix_cluster,drop=FALSE],1,mean)
    }
    print(paste(j))
  }
  return(U_out)
}
consolidate_U = cmpfun(temp)

rm(list="temp")
