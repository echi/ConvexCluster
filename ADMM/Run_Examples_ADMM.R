## This script tests the convex_clustering algorithms in the admm.f90 module
## on four data sets.
##
## Note these scripts only reproduce timing results and figures, except for the halfmoons figure.

rm(list=ls())
source('R/Functions/admm.R')
source('R/Functions/admm_convex_cluster.R')
source('R/Functions/cluster_path_preprocess.R')
library(ggplot2)
library(clusterpath)

#######################################
## Example 1: Two Moons              ##
#######################################

source('../Examples/ex_two_moons.R')
source('R/Examples/ex-admm_halfmoons_timing.R')

#######################################
## Example 2: Iris                   ##
#######################################

source('R/Examples/ex-admm_iris.R')
source('R/Examples/ex-admm_iris_timing.R')

#######################################
## Example 3: Senate                 ##
#######################################

source('R/Examples/ex-admm_senate.R')
source('R/Examples/ex-admm_senate_timing.R')

#######################################
## Example 4: Mammals                ##
#######################################

source('R/Examples/ex-admm_mammals.R')
source('R/Examples/ex-admm_mammals_timing.R')