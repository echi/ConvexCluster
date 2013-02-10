## Test case 0
make_simple_clusters = function(q=2, p=100, sigma=0.1, seed=12345) {
  mu1 = matrix(0,q,1)
  mu2 = matrix(1,q,1)

  X = matrix(0,q,p)
  X[,1:(p/2)] = matrix(mu1,q,p/2) + sigma*matrix(rnorm(q*p/2),q,p/2)
  X[,(p/2 + 1):p] = matrix(mu2,q,p/2) +sigma* matrix(rnorm(q*p/2),q,p/2)
  return(X)
}

make_2d_sparse_clusters = function(shift=0, n=500, sigma=2, seed=12345) {
  x = rnorm(n=n)
  y = 2*x + sigma*rnorm(n=n)
  return(X=cbind(x+shift,y))
}