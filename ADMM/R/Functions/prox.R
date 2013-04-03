## These functions perform various common proximal mappings.

# type = 1 : prox_L1                                                                                                                                                                                
# type = 2 : prox_L2                                                                                                                                                                                
# type = 3 : prox_L-infinity                                                                                                                                                                                
# type = 4 : proj_L-infinity                                                                                                                                                                                
# type = 5 : proj_L2                                                                                                                                                                        
# type = 6 : proj_L1                                                                                                                                                                               

library(compiler)

## Proximal mapping: h(x) = tau * ||x||_2
temp = function(x, tau) {
  lv = norm(as.matrix(x),type="f")
  if (lv == 0) {
    px = x
  } else {
    px = max(0, 1-tau/lv)*x
  }
  return(px)
}
prox_L2 = cmpfun(temp)

## Proximal mapping: h(x) = tau * ||x||_1
temp = function(x,tau) {
  return(sign(x)*sapply(abs(x)-tau,FUN=function(x) {max(x,0)}))
}
prox_L1 = cmpfun(temp)

## Projection onto L2 ball of radius tau
temp = function(x,tau) {
  lv = norm(as.matrix(x),'f')
  u = x
  if (lv > tau) {
    u = (tau/lv)*x
  }
  return(u)
}
proj_L2 = cmpfun(temp)

## Projection onto L-infinty ball of radius tau
temp = function(x,tau) {
  u = as.matrix(abs(x))
  cbind(u,tau)
  u = apply(cbind(u,tau),1,min)
  u = u*sign(x)
  return(u)
}
proj_Linf = cmpfun(temp)

## Projection onto a simplex of size z
temp = function(v,z) {
  mu = sort(v, decreasing=T)
  
  partial_sum <- mu[1]
  rho = 1
  for (j in 2:length(v)) {
    partial_sum <- partial_sum + mu[j]
    if (j*mu[j] - partial_sum + z <= 0) {
      rho = j - 1
      break
    }
    rho = j
  }
  theta <- mean(mu[1:rho]) - z/rho
  w <- v - theta
  w[w <= 0] = 0
  return(w)
}
project_to_simplex = cmpfun(temp)

## Projection onto L1 ball of radius tau
temp = function(x,tau) {
  w = project_to_simplex(abs(x),tau)
  return(w*sign(x))
}
proj_L1 = cmpfun(temp)

## Proximal mapping of L-infinity norm
temp = function(x,tau) {
  return(x - proj_L1(x,tau))
}
prox_Linf = cmpfun(temp)

temp = function(x,tau,type) {
  if (type == 1) {
    return(prox_L1(x,tau))
  } else if (type==2) {
    return(prox_L2(x,tau))
  } else if (type==3) {
    return(prox_Linf(x,tau))
  } else if (type==4) {
    return(proj_Linf(x,tau))  
  } else if (type==5) {
    return(proj_L2(x,tau))    
  } else if (type==6) {
    return(proj_L1(x,tau))
  }
}
prox_R = cmpfun(temp)
rm(temp)

prox_F = function(x,tau,type) {
  n = as.integer(length(x))
  x = as.double(x)
  tau = as.double(tau)
  type = as.integer(type)
  sol = .Fortran('test_prox',vector_in=x,n=n,vector_out=double(n),tau=tau,type=type)
  return(sol$vector_out)
#  subroutine prox(vector_in,n,vector_out,tau,type)
  
}