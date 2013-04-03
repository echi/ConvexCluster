## This script tests the convex_clustering algorithms in the ama.f90 module
## on the two moons data

X = two_moons(n=200,p=2)$X
p = ncol(X)

df.paths = data.frame(x=c(),y=c(), group=c(),gamma=c(),k=c())

## Trial 1
kgamma = 0
knbs = 10
w = kernel_weights(X,kgamma)
w = knn_weights(w,knbs,p)
gamma_seq = seq(0.0,0.63, length.out=1000)
tol = 1e-3
nu = 10
system.time({ccp = convex_cluster_path(X,w,gamma_seq,nu=nu,tol=tol,type=2)})
ix = compactify_edges(w,p)$ix
nGamma = ccp$nGamma
Uout = consolidate_U(ccp$UHx,ccp$VHx,ix,nGamma)

x = double(nGamma)
y = double(nGamma)
for (i in 1:p) {
  for (j in 1:nGamma) {
    U = Uout[[j]]
    x[j] = U[1,i]
    y[j] = U[2,i]
  }
  df = data.frame(x=x, y=y)
  df$group = i
  df$gamma = kgamma
  df$k = knbs
  df.paths = rbind(df.paths,df)
  print(paste(i))
}

## Trial 2
kgamma = 0
knbs = 50
w = kernel_weights(X,kgamma)
w = knn_weights(w,knbs,p)
gamma_seq = seq(0.0,0.96, length.out=1000)
gamma_seq = seq(0.0,0.027,length.out=1000)
tol = 1e-3
nu = 10
system.time({ccp = convex_cluster_path(X,w,gamma_seq,nu=nu,tol=tol,type=2)})
ix = compactify_edges(w,p)$ix
nGamma = ccp$nGamma
Uout = consolidate_U(ccp$UHx,ccp$VHx,ix,nGamma)
x = double(nGamma)
y = double(nGamma)
for (i in 1:p) {
  for (j in 1:nGamma) {
    U = Uout[[j]]
    x[j] = U[1,i]
    y[j] = U[2,i]    
#    x[j] = ccp$UHx[[j]][1,i]
#    y[j] = ccp$UHx[[j]][2,i]
  }
  df = data.frame(x=x, y=y)
  df$group = i
  df$gamma = kgamma
  df$k = knbs
  df.paths = rbind(df.paths,df)
}

## Trial 3
kgamma = 0.5
knbs = 50
w = kernel_weights(X,kgamma)
w = knn_weights(w,knbs,p)
gamma_seq = seq(0.0,0.16, length.out=1000)
tol = 1e-3
nu = 10
system.time({ccp = convex_cluster_path(X,w,gamma_seq,nu=nu,tol=tol,type=2)})
nGamma = ccp$nGamma
ix = compactify_edges(w,p)$ix
Uout = consolidate_U(ccp$UHx,ccp$VHx,ix,nGamma)
x = double(nGamma)
y = double(nGamma)
for (i in 1:p) {
  for (j in 1:nGamma) {
#    x[j] = ccp$UHx[[j]][1,i]
#    y[j] = ccp$UHx[[j]][2,i]
    U = Uout[[j]]
    x[j] = U[1,i]
    y[j] = U[2,i]        
  }
  df = data.frame(x=x, y=y)
  df$group = i
  df$gamma = kgamma
  df$k = knbs
  df.paths = rbind(df.paths,df)
}

## Trial 4
kgamma = 0.5
knbs = 10
w = kernel_weights(X,kgamma)
w = knn_weights(w,knbs,p)
gamma_seq = seq(0.0,6, length.out=1000)
tol = 1e-3
nu = 10
system.time({ccp1 = convex_cluster_path(X,w,gamma_seq,nu=nu,tol=tol,type=2)})
nGamma = ccp1$nGamma
ix = compactify_edges(w,p)$ix
Uout = consolidate_U(ccp1$UHx,ccp1$VHx,ix,nGamma)
x = double(nGamma)
y = double(nGamma)
for (i in 1:p) {
  for (j in 1:nGamma) {
    U = Uout[[j]]
    x[j] = U[1,i]
    y[j] = U[2,i]          
#    x[j] = ccp1$UHx[[j]][1,i]
#    y[j] = ccp1$UHx[[j]][2,i]
  }
  df = data.frame(x=x, y=y)
  df$group = i
  df$gamma = kgamma
  df$k = knbs
  df.paths = rbind(df.paths,df)
}

tol = 2.5e-3
nu = 10

gamma_seq = seq(6.01,8.2, length.out=1000)

system.time({ccp2 = convex_cluster_path(X,w,gamma_seq,nu=nu,tol=tol)})
nGamma = ccp2$nGamma
Uout = consolidate_U(ccp2$UHx,ccp2$VHx,ix,nGamma)
x = double(nGamma)
y = double(nGamma)
for (i in 1:p) {
  for (j in 1:nGamma) {
    U = Uout[[j]]
    x[j] = U[1,i]
    y[j] = U[2,i]            
#    x[j] = ccp2$UHx[[j]][1,i]
#    y[j] = ccp2$UHx[[j]][2,i]
  }
  df = data.frame(x=x, y=y)
  df$group = i
  df$gamma = kgamma
  df$k = knbs
  df.paths = rbind(df.paths,df)
}
X_data = as.data.frame(t(X[1:2,]))
colnames(X_data) = c("x","y")
data_plot = ggplot(data=df.paths,aes(x=x,y=y))
data_plot = data_plot + geom_path(aes(group=group),colour='grey60',alpha=0.5)
data_plot = data_plot + facet_grid(k~gamma)
data_plot = data_plot + geom_point(data=X_data,aes(x=x,y=y))
data_plot + theme_bw() + theme(legend.position = "none")
filename = "halfmoons.pdf"
golden_ratio = 1.61803398875
height = 7
ggsave(filename=filename,width=golden_ratio*height,height=height)

## Visualize results
# X_data = as.data.frame(t(X[1:2,]))
# colnames(X_data) = c("x","y")
# 
# nGamma = ccp$nGamma
# x = double(nGamma)
# y = double(nGamma)
# data_plot = ggplot(data=X_data,aes(x=x,y=y))
# data_plot = data_plot + geom_point(size=2)
# for (i in 1:p) {
#   for (j in 1:nGamma) {
#     x[j] = ccp$UHx[[j]][1,i]
#     y[j] = ccp$UHx[[j]][2,i]
#   }
#   df = data.frame(x=x, y=y)
#   data_plot = data_plot + geom_path(data=df,aes(x=x, y=y),colour='grey60',alpha=0.5)
# }
# nGamma = ccp2$nGamma
# x = double(nGamma)
# y = double(nGamma)
# for (i in 1:p) {
#   for (j in 1:nGamma) {
#     x[j] = ccp2$UHx[[j]][1,i]
#     y[j] = ccp2$UHx[[j]][2,i]
#   }
#   df = data.frame(x=x, y=y)
#   data_plot = data_plot + geom_path(data=df,aes(x=x, y=y),colour='grey60',alpha=0.5)
# }
# data_plot + theme_bw() + theme(legend.position = "none") 
# 
# df.paths = data.frame(x=c(),y=c(), group=c(),gamma=c(),k=c())
# for (i in 1:p) {
#   for (j in 1:nGamma) {
#     x[j] = ccp1$UHx[[j]][1,i]
#     y[j] = ccp1$UHx[[j]][2,i]
#   }
#   df = data.frame(x=x, y=y)
#   df$group = i
#   df$gamma = kgamma
#   df$k = knbs
#   df.paths = rbind(df.paths,df)
# }
# data_plot = ggplot(data=df.paths,aes(x=x,y=y))
# data_plot + geom_path(aes(group=group),colour='grey60',alpha=0.5) + theme_bw() + theme(legend.position = "none")