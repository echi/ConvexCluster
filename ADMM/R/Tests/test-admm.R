## This script tests the updates and residual calculations in the admm.f90 module

context("Parameter Updates")

test_that("Proximal Mappings are correct", {
  n = 10
  set.seed(12345)
  nTrials = 100
  for (i in 1:nTrials) {
    x = rnorm(n)
    tau = rgamma(1,shape=1)
    type = 1
    expect_that(prox_F(x,tau,type),equals(prox_R(x,tau,type)))
    type = 2
    expect_that(prox_F(x,tau,type),equals(prox_R(x,tau,type)))
    type = 3
    expect_that(prox_F(x,tau,type),equals(prox_R(x,tau,type)))
  }
})

test_that("residual_primal is working", {
  ## Create random problem
  set.seed(12345)
  nTrials = 100
  q = 10
  p = 20
  for (i in 1:nTrials) {
    X = matrix(rnorm(q*p),q,p)
    w = kernel_weights(X,0)
    w = knn_weights(w,3,p)
    ix = compactify_edges(w,p)$ix
    nK = nrow(ix)
    U = matrix(rnorm(q*p),q,p)
    V = matrix(rnorm(q*nK),q,nK)
    expect_that(residual_primalF(U,V,ix), equals(residual_primalR(U,V,ix)))
  }    
})

test_that("residual_dual is working", {
  ## Create random problem
  set.seed(12345)
  nTrials = 100
  q = 10
  p = 20
  for (i in 1:nTrials) {
    X = matrix(rnorm(q*p),q,p)
    w = kernel_weights(X,0)
    w = knn_weights(w,3,p)
    edge_info = compactify_edges(w,p)
    
    M1 = edge_info$M1
    M2 = edge_info$M2
    s1 = edge_info$s1
    s2 = edge_info$s2
    nK = length(which(w>0))
    V = matrix(rnorm(q*nK),q,nK) 
    V_old = matrix(rnorm(q*nK),q,nK)
    nu = 1
    expect_that(residual_dualF(V,V_old,M1,M2,s1,s2,nu,p), equals(residual_dualR(V,V_old,M1,M2,s1,s2,nu,p)))
  }    
})
  
test_that("update_U is working", {
  ## Create random problem
  set.seed(12345)
  nTrials = 100
  q = 10
  p = 20
  for (i in 1:nTrials) {
    X = matrix(rnorm(q*p),q,p)
    w = kernel_weights(X,0)
    w = knn_weights(w,3,p)
    edge_info = compactify_edges(w,p)
    
    M1 = edge_info$M1
    M2 = edge_info$M2
    s1 = edge_info$s1
    s2 = edge_info$s2
    nK = length(which(w>0))
    Lambda = matrix(rnorm(q*nK),q,nK)
    V = matrix(rnorm(q*nK),q,nK)      
    nu = 1
    expect_that(update_UF(X,V,Lambda,M1,M2,s1,s2,nu), equals(update_UR(X,V,Lambda,M1,M2,s1,s2,nu)))
  }
})

test_that("update_V is working", {
  ## Create random problem
  set.seed(12345)
  nTrials = 100
  q = 10
  p = 20
  for (i in 1:nTrials) {
    U = matrix(rnorm(q*p),q,p)
    X = matrix(rnorm(q*p),q,p)
    w = kernel_weights(X,0)
    w = knn_weights(w,3,p)
    edge_info = compactify_edges(w,p)
    w = w[w>0]
    nK = length(w)
    nu = rexp(1)
    nu = 1
    gamma = rexp(1)
    Lambda = matrix(rnorm(q*nK),q,nK)
    ix = edge_info$ix
    type = 1
    expect_that(update_VF(U,Lambda,w,gamma,nu,ix,type), equals(update_VR(U,Lambda,w,gamma,nu,ix,type)))
    type = 2
    expect_that(update_VF(U,Lambda,w,gamma,nu,ix,type), equals(update_VR(U,Lambda,w,gamma,nu,ix,type)))
    type = 3
    expect_that(update_VF(U,Lambda,w,gamma,nu,ix,type), equals(update_VR(U,Lambda,w,gamma,nu,ix,type)))    
  }
})

test_that("update_Lambda is working", {
  ## Create random problem
  set.seed(12345)
  nTrials = 100
  q = 10
  p = 20
  for (i in 1:nTrials) {
    U = matrix(rnorm(q*p),q,p)
    X = matrix(rnorm(q*p),q,p)
    w = kernel_weights(X,0)
    w = knn_weights(w,3,p)
    edge_info = compactify_edges(w,p)
    w = w[w>0]
    nK = length(w)
    nu = rexp(1)
    Lambda = matrix(rnorm(q*nK),q,nK)
    ix = edge_info$ix
    V = matrix(rnorm(q*nK),q,nK)
    
    expect_that(update_LambdaF(Lambda,U,V,nu,ix), equals(update_LambdaR(Lambda,U,V,nu,ix)))
  }
})