## This script tests the convex_clustering algorithms in the ama.f90 module
## on iris data

#######################################
## Part B: Timing comparisons        ##
#######################################

nRuns = 10

times_hok_l2 = double(nRuns)
times_hok_l1 = double(nRuns)
times_ama_fast_l2 = double(nRuns)
times_ama_fast_l1 = double(nRuns)
times_ama_slow_l2 = double(nRuns)
times_ama_slow_l1 = double(nRuns)

X = as.matrix(iris[,1:4])
X = t(scale(X,center=TRUE,scale=FALSE))
p = ncol(X)
q = nrow(X)
#****************************#
#* Run Hocking algorithm    *#
#****************************#
kgamma = 1e-2
for (i in 1:nRuns) {
  times_hok_l2[i] = system.time({path_l2 <- clusterpath.l2(t(X),gamma=kgamma)})[3] 
  print(paste0("Hocking L2: Completed ",i))
}
gamma_seq_l2 = unique(path_l2$lambda)

for (i in 1:nRuns) {
  times_hok_l1[i] = system.time({path_l1 <- clusterpath.l1.general(t(X),verbose=0,gamma=kgamma)})[3]
  print(paste0("Hocking L1: Completed ",i))
}
gamma_seq_l1 = unique(path_l1$lambda)

#****************************#
#* Run AMA FAST: L2         *#
#****************************#
w = kernel_weights(X,kgamma)
tol = 1e-5
nu = 10
type = 2
for (i in 1:nRuns) {
  times_ama_fast_l2[i] = system.time({iris_ama_fast_l2 = convex_cluster_path_fista_backtrack(X,w,gamma_seq_l2,nu=nu,tol=tol,type=type)})[3]
  print(paste0("AMA Fast L2: Completed ",i))
}
#****************************#
#* Run AMA FAST: L1         *#
#****************************#
type=1
for (i in 1:nRuns) {
  times_ama_fast_l1[i] = system.time({iris_ama_fast_l1 = convex_cluster_path_fista_backtrack(X,w,gamma_seq_l1,nu=nu,tol=tol,type=type)})[3]
  print(paste0("AMA Fast L1: Completed ",i))
}
#****************************#
#* Run AMA SLOW: L2         *#
#****************************#
type = 2
for (i in 1:nRuns) {
  times_ama_slow_l2[i] = system.time({iris_ama_slow_l2 = convex_cluster_path_backtrack(X,w,gamma_seq_l2,nu=nu,tol=tol,type=type)})[3]
  print(paste0("AMA Slow L2: Completed ",i))
}
#****************************#
#* Run AMA SLOW: L1         *#
#****************************#
type = 1
for (i in 1:nRuns) {
  times_ama_slow_l1[i] = system.time({iris_ama_slow_l1 = convex_cluster_path_backtrack(X,w,gamma_seq_l1,nu=nu,tol=tol,type=type)})[3]
  print(paste0("AMA Slow L1: Completed ",i))
}
#****************************#
#* Plot results for         *#
#* verification purposes    *#
#* Use ADMM code for plots  *#
#* included in paper.       *#
#****************************#
#ix = compactify_edges(w,p)$ix
#nGamma = iris_ama_fast_l1$nGamma
##Uout = consolidate_U(iris_ama_fast_l1$UHx,iris_ama_fast_l1$VHx,ix,nGamma)
#Uout = iris_ama_slow_l2$UHx
#df.paths = data.frame(x=c(),y=c(), group=c())
#nGamma = iris_ama_fast_l1$nGamma
#svdX = svd(X)
#pc = svdX$u[,1:2,drop=FALSE]
#pc.df = as.data.frame(t(pc)%*%X)
#
#x = double(nGamma)
#y = double(nGamma)
#for (j in 1:nGamma) {
#  pcs = t(pc)%*%Uout[[j]]
#  x = pcs[1,]
#  y = pcs[2,]
#  df = data.frame(x=pcs[1,], y=pcs[2,], group=1:p)
#  df.paths = rbind(df.paths,df)
#}
#
#X_data = as.data.frame(t(X)%*%pc)
#colnames(X_data) = c("x","y")
#X_data$Species = iris[,5]
#data_plot = ggplot(data=df.paths,aes(x=x,y=y))
#data_plot = data_plot + geom_path(aes(group=group),colour='grey60',alpha=0.5)
#data_plot = data_plot + geom_point(data=X_data,aes(x=x,y=y,colour=Species,shape=Species))
#data_plot = data_plot + xlab('Principal Component 1') + ylab('Principal Component 2')
#data_plot + theme_bw()
#****************************#
#* Compute Objective        *#
#* Relative difference: L2  *#
#****************************#
nGamma = length(gamma_seq_l2)
obj_hok_l2 = double(nGamma)
obj_ama_l2 = double(nGamma)
ix = vec2tri(1:(p*(p-1)/2),p)
type = 2
for (iGamma in 1:nGamma) {
  gamma = gamma_seq_l2[iGamma]
  obj_hok_l2[iGamma] = loss_primalF(X, t(path_l2[(p*(iGamma-1)+1):(p*iGamma),1:q]),gamma=gamma,w=w,ix=ix,type=type)
  obj_ama_l2[iGamma] = loss_primalF(X, iris_ama_fast_l2$UHx[[iGamma]],gamma=gamma,w=w,ix=ix,type=2)
}
rel_ama2hok_l2 =  (obj_hok_l2[-1]-obj_ama_l2[-1])/obj_hok_l2[-1]
#****************************#
#* Compute Objective        *#
#* Relative difference: L1  *#
#****************************#
nGamma = length(gamma_seq_l1)
obj_hok_l1 = double(nGamma)
obj_ama_l1 = double(nGamma)
ix = vec2tri(1:(p*(p-1)/2),p)
type = 1
for (iGamma in 1:nGamma) {
  gamma = gamma_seq_l1[iGamma]
  obj_hok_l1[iGamma] = loss_primalF(X, t(path_l1[(p*(iGamma-1)+1):(p*iGamma),1:q]),gamma=gamma,w=w,ix=ix,type=type)
  obj_ama_l1[iGamma] = loss_primalF(X, iris_ama_fast_l1$UHx[[iGamma]],gamma=gamma,w=w,ix=ix,type=type)
}
rel_ama2hok_l1 =  (obj_hok_l1[-1]-obj_ama_l1[-1])/obj_hok_l1[-1]

save(list=c("times_hok_l2","times_hok_l1","gamma_seq_l2","gamma_seq_l1",
            "times_ama_fast_l2","times_ama_fast_l1","times_ama_slow_l2","times_ama_slow_l1",
            "obj_hok_l2","obj_hok_l1","obj_ama_l2","obj_ama_l1","rel_ama2hok_l2","rel_ama2hok_l1"),
     file="iris_hocking_ama.rda")

## Compute relative error in objective and centroid values for the various algorithms
#library(ggplot2)
#library(reshape2)

#load('iris_ground_truth.rda')

#gamma_seq_exp = gamma_seq[1]
#nInsert = 100
#for (i in 2:nGamma) {
#  tmp = seq(gamma_seq[i-1],gamma_seq[i],length.out=nInsert)
#  gamma_seq_exp = c(gamma_seq_exp,tmp[-1])
#}

#obj_rel_hok = (obj_hok - obj_true) / obj_true
#obj_rel_ama = (obj_ama - obj_true) / obj_true
#compare_iris1 = data.frame(ActiveSet=obj_rel_hok[-1], AMA=obj_rel_ama[-1], gamma=gamma_seq[-1])
#compare_iris1 = melt(compare_iris1,id.vars="gamma")
#compare_iris1$type = "Objective"

#error_ama = double(nGamma)
#error_amaB = double(nGamma)
#error_hok = double(nGamma)
#gix = match(gamma_seq,gamma_seq_exp)
#for (i in 1:nGamma) {
#  error_ama[i] = norm(as.matrix(ccp$UHx[[gix[i]]]) - as.matrix(ccpB$UHx[[i]]),'f')/norm(as.matrix(ccp$UHx[[gix[i]]]),'f')
#  error_hok[i] = norm(as.matrix(ccp$UHx[[gix[i]]]) - t(as.matrix(path[(p*(i-1)+1):(p*i),1:4])),'f')/norm(as.matrix(ccp$UHx[[gix[i]]]),'f')
#}
#compare_iris2 = data.frame(ActiveSet=error_hok[-1], AMA=error_ama[-1], gamma=gamma_seq[-1])
#compare_iris2 = melt(compare_iris2,id.vars="gamma")
#compare_iris2$type = "Centroids"
#compare_iris = rbind(compare_iris1, compare_iris2)
#q = ggplot(data=compare_iris, aes(x=gamma,y=value))
#q = q + geom_point(aes(colour=variable, shape=variable)) + scale_shape(name="Algorithms",labels=c("Hocking et. al.","AMA"))
#q = q + scale_colour_manual(name="Algorithms",labels=c("Hocking et. al.","AMA"), values = c(rgb(241/255, 163/255, 64/255), rgb(153/255, 142/255, 195/255)))
#q = q + facet_grid(.~ type, scales="free_y")
#q + theme_bw() + scale_x_log10(breaks=c(0.01,0.01710339,0.03225100),labels=signif(c(0.01,0.01710339,0.03225100),digits=3)) + scale_y_log10() + ylab("Relative Error") + xlab(expression(gamma))
#filename = "iris_l2.pdf"
#golden_ratio = 1.61803398875
#height = 7
#ggsave(filename=filename,width=golden_ratio*height,height=height)

## Plot Timing Comparisons
#times = data.frame(times = times_hok, Method="Hok")
#times = rbind(times, data.frame(times=times_ama, Method="AMA"))
#times$method = factor(times$Method)
#q = ggplot(data=times,aes(x=Method,y=times))
#q + geom_boxplot()