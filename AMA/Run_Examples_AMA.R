## This script tests the convex_clustering algorithms in the ama.f90 module
## on four data sets.
##
## Note these scripts only reproduce timing results (except the halfmoons example).
## See ADMM/Run_Examples.R for R code that generates figures.

rm(list=ls())
source('R/Functions/ama_loss.R')
source('R/Functions/ama_updates.R')
source('R/Functions/ama_convex_cluster.R')
source('R/Functions/cluster_path_preprocess.R')
library(ggplot2)
library(clusterpath)

#######################################
## Example 1: Two Moons              ##
#######################################

source('../Examples/ex_two_moons.R')
source('R/Examples//ex-ama_halfmoons.R')
source('R/Examples//ex-ama_halfmoons_timing.R')

#######################################
## Example 2: Iris                   ##
#######################################

source('R/Examples/ex-ama_iris_timing.R')

#######################################
## Example 3: Senate                 ##
#######################################

source('R/Examples/ex-ama_senate_timing.R')

#######################################
## Example 4: Mammals                ##
#######################################

source('R/Examples/ex-ama_mammals_timing.R')