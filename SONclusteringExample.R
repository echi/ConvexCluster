#
# SONclustering
#
# Date: October 19, 2012
# Author: Eric Chi
#

#
# Demo example: Really easy two cluster problem
#
n = 50
p = 2
set.seed(12345)
X = matrix(0,n,p)
X[1:25,] = 0
X[26:50,] = c(1,1)
X = X + 0.125*matrix(rnorm(n*p),n,p)
plot(X[,1],X[,2])