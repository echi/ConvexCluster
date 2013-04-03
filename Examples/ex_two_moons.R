## Two moons data
library(ggplot2)
library(compiler)

temp = function(n=100,p=500,var=0.01,seed=12345) {
  set.seed(seed)
  class = double(n)
  class[(n/2+1):n] = 1
  x = seq(0,pi,length.out=n/2)
  M = cbind(cos(x), sin(x))
  M = rbind(M,cbind(cos(x)+1,0.5-sin(x)))
  M = cbind(M, matrix(0,n,p-2))
  M = M + sqrt(var)*matrix(rnorm(n*p),n,p)
  return(list(X=t(M),class=class))
}
two_moons = cmpfun(temp)

plot_two_moons = function(M,class) {
## Visualize the data
  tm = data.frame(x=M[1,],y=M[2,],class=factor(class))
  q = ggplot(data=tm, aes(x=x,y=y))
  q + geom_point(aes(colour=class)) + theme_bw()
}