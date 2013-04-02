## These functions perform various versions of convex clustering via ADMM
library(compiler)
library(igraph)

#######################################
## R wrappers around Fortran versions #
#######################################
source('cluster_path_preprocess.R')

if (is.loaded('convex_cluster_admm')) {
  dyn.unload('admm.so')
}
dyn.load('admm.so')

convex_cluster_admm = function(X,Lambda,V,ix,M1,M2,s1,s2,w,gamma,nu=1,type=2,max_iter=1e2,tol=1e-4) {
  q = as.integer(nrow(X))
  p = as.integer(ncol(X))
  nK = as.integer(ncol(Lambda))
  mix1 = as.integer(nrow(M1))
  mix2 = as.integer(nrow(M2))
  storage.mode(X) = "double"
  storage.mode(Lambda) = "double"
  U = matrix(0,q,p)
  storage.mode(U) = "double"
  storage.mode(V) = "double"
  storage.mode(ix) = "integer"
  storage.mode(M1) = "integer"
  storage.mode(M2) = "integer"
  s1 = as.integer(s1)
  s2 = as.integer(s2)
  w = as.double(w)
  gamma = as.double(gamma)
  nu = as.double(nu)
  type = as.integer(type)
  max_iter = as.integer(max_iter)
  tol = as.double(tol)
  primal = double(max_iter)
  dual = double(max_iter)
  
#  subroutine convex_cluster(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,s1,s2,M1,M2,mix1,mix2,primal,dual,max_iter,iter,tol,type)
  
  sol = .Fortran('convex_cluster_admm',X=X,Lambda=Lambda,U=U,V=V,q=q,p=p,nK=nK,ix=ix,w=w,gamma=gamma,nu=nu,
                 s1=s1,s2=s2,M1=M1,M2=M2,mix1=mix1,mix2=mix2,primal=primal,dual=dual,
                 max_iter=max_iter,iter=integer(1),tol=tol,type=type)
  return(list(U=sol$U,V=sol$V,Lambda=sol$Lambda,nu=sol$nu,
              primal=sol$primal[1:sol$iter],dual=sol$dual[1:sol$iter],iter=sol$iter))
}
  
convex_clusterB_admm = function(X,Lambda,V,ix,w,gamma,nu=1,type=2,max_iter=1e2,tol=1e-4) {

  q = as.integer(nrow(X))
  p = as.integer(ncol(X))
  nK = as.integer(ncol(Lambda))
  gamma = as.double(gamma)
  max_iter = as.integer(max_iter)
  tol = as.double(tol)
  storage.mode(X) = "double"
  w = as.double(w)
  storage.mode(Lambda) = "double"
  storage.mode(V) = "double"
  storage.mode(ix) = "integer"
  U = matrix(0,q,p)
  storage.mode(U) = "double"
  nu = as.double(nu)
  type = as.integer(type)
  primal = double(max_iter)
  dual = double(max_iter)
  
#  subroutine convex_clusterB(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,primal,dual,max_iter,iter,tol,type)
  
  sol = .Fortran('convex_clusterB_admm',X=X,Lambda=Lambda,U=U,V=V,q=q,p=p,nK=nK,ix=ix,w=w,gamma=gamma,nu=nu,
                 primal=primal,dual=dual,max_iter=max_iter,iter=integer(1),tol=tol,type=type)
  return(list(U=sol$U,V=sol$V,Lambda=sol$Lambda,nu=sol$nu,
              primal=sol$primal[1:sol$iter],dual=sol$dual[1:sol$iter],iter=sol$iter))
}

convex_cluster_admm_acc = function(X,Lambda,V,ix,w,gamma,nu=1,type=2,max_iter=1e2,tol=1e-4) {
  
  q = as.integer(nrow(X))
  p = as.integer(ncol(X))
  nK = as.integer(ncol(Lambda))
  gamma = as.double(gamma)
  max_iter = as.integer(max_iter)
  tol = as.double(tol)
  storage.mode(X) = "double"
  w = as.double(w)
  storage.mode(Lambda) = "double"
  storage.mode(V) = "double"
  storage.mode(ix) = "integer"
  U = matrix(0,q,p)
  storage.mode(U) = "double"
  nu = as.double(nu)
  type = as.integer(type)
  primal = double(max_iter)
  dual = double(max_iter)
  
#  subroutine convex_cluster_acc(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,primal,dual,max_iter,iter,tol,type)
  
  sol = .Fortran('convex_cluster_admm_acc',X=X,Lambda=Lambda,U=U,V=V,q=q,p=p,nK=nK,ix=ix,w=w,gamma=gamma,nu=nu,
                 primal=primal,dual=dual,max_iter=max_iter,iter=integer(1),tol=tol,type=type)
  return(list(U=sol$U,V=sol$V,Lambda=sol$Lambda,nu=sol$nu,
              primal=sol$primal[1:sol$iter],dual=sol$dual[1:sol$iter],iter=sol$iter))
}

temp = function(X, w, gamma_seq,nu=1,tol=1e-4,max_iter=1e6,type=2) {
  nGamma = length(gamma_seq)
  p = ncol(X)
  q = nrow(X)
  nK = p*(p-1)/2
  Lambda = matrix(0,q,nK)
  V = matrix(0,q,nK)
  list_U = vector(mode="list",length=nGamma)
  list_V = vector(mode="list",length=nGamma)
  list_Lambda = vector(mode="list",length=nGamma)
  ix = vec2tri(1:nK,p)  
      
  for (ig in 1:nGamma) {
    gamma = gamma_seq[ig]      
    cc = convex_clusterB_admm(X,Lambda,V,ix,w,gamma,nu=nu,type=type,max_iter=max_iter,tol=tol)    
    Lambda = cc$Lambda
    V = cc$V
    list_U[[ig]] = cc$U
    list_V[[ig]] = V
    list_Lambda[[ig]] = Lambda 
#    print(paste0("iters: ", cc$iter,"| primal: ", cc$primal[cc$iter],
#                 "| dual: ",cc$dual[cc$iter],
#                 "| max: ", max(cc$primal[cc$iter],cc$dual[cc$iter])))
#    print(paste("Completed gamma",ig))
#    if (norm(cc$V,'1')==0) {
#      print('Single cluster')
#      break
#    }
  }
  return(list(UHx=list_U,VHx=list_V,LambdaHx=list_Lambda,nGamma=ig))
}
convex_cluster_pathB_admm = cmpfun(temp)

temp = function(X, w, gamma_seq,nu=1,tol=1e-4,max_iter=1e6,type=2) {
  nGamma = length(gamma_seq)
  p = ncol(X)
  q = nrow(X)
  nK = p*(p-1)/2
  Lambda = matrix(0,q,nK)
  V = matrix(0,q,nK)
  list_U = vector(mode="list",length=nGamma)
  list_V = vector(mode="list",length=nGamma)
  list_Lambda = vector(mode="list",length=nGamma)
  ix = vec2tri(1:nK,p)  
  
  for (ig in 1:nGamma) {
    gamma = gamma_seq[ig]      
    cc = convex_cluster_admm_acc(X,Lambda,V,ix,w,gamma,nu=nu,type=type,max_iter=max_iter,tol=tol)    
    Lambda = cc$Lambda
    V = cc$V
    list_U[[ig]] = cc$U
    list_V[[ig]] = V
    list_Lambda[[ig]] = Lambda 
#    print(paste0("iters: ", cc$iter,"| primal: ", cc$primal[cc$iter],
#                 "| dual: ",cc$dual[cc$iter],
#                 "| max: ", max(cc$primal[cc$iter],cc$dual[cc$iter])))
#    print(paste("Completed gamma",ig))
#    if (norm(cc$V,'1')==0) {
#      print('Single cluster')
#      break
#    }
  }
  return(list(UHx=list_U,VHx=list_V,LambdaHx=list_Lambda,nGamma=ig))
}
convex_cluster_path_admm_acc = cmpfun(temp)

temp = function(X, w, gamma_seq,nu=1,tol=1e-4,max_iter=1e6,type=2) {
  nGamma = length(gamma_seq)
  p = ncol(X)
  q = nrow(X)
  edge_info = compactify_edges(w,p)
  nK = length(which(w > 0))
  Lambda = matrix(0,q,nK)
  V = matrix(0,q,nK)
  list_U = vector(mode="list",length=nGamma)
  list_V = vector(mode="list",length=nGamma)
  list_Lambda = vector(mode="list",length=nGamma)
  ix = edge_info$ix  
  M1 = edge_info$M1
  M2 = edge_info$M2
  s1 = edge_info$s1
  s2 = edge_info$s2
  
  for (ig in 1:nGamma) {
    gamma = gamma_seq[ig]      
    cc = convex_cluster_admm(X,Lambda,V,ix,M1,M2,s1,s2,w[w>0],gamma,nu=nu,type=type,max_iter=max_iter,tol=tol)    
    Lambda = cc$Lambda
    V = cc$V
    list_U[[ig]] = cc$U
    list_V[[ig]] = V
    list_Lambda[[ig]] = Lambda 
    print(paste0("iters: ", cc$iter,"| primal: ", cc$primal[cc$iter],
                 "| dual: ",cc$dual[cc$iter],
                 "| max: ", max(cc$primal[cc$iter],cc$dual[cc$iter])))
    print(paste("Completed gamma",ig))
    if (norm(cc$V,'1')==0) {
      print('Single cluster')
      break
    }
  }
  return(list(UHx=list_U,VHx=list_V,LambdaHx=list_Lambda,nGamma=ig))
}
convex_cluster_path_admm = cmpfun(temp)

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
