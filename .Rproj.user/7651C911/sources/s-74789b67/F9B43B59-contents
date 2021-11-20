library(Rcpp)
library(fabricatr)
library(mltools)
library(onehot)
library(mgcv)
library(tidyverse)
library(foreach)
library(doParallel)
library(kableExtra)
library(boot)
library(tensorr)
library(reshape2)
library(dvmisc)
library(hash)
source("Utils.R")
source("pfilter.R")
sourceCpp("Cov_struct_func.cpp")
setwd("~/FABBrainConnectome/FAB_corr")
options(warn=-1)


func_sim_boot <- function(data_info, R=2000){

  data_raw <- data_info$data_raw
  data_ref <- data_info$data_ref
  ttl <- data_info$ttl
  Sigma <- data_info$Sigma
  zero_idx <- data_info$zero_idx
  null_idx_cbind <- data_info$null_idx_cbind
  alter_idx_cbind <- data_info$alter_idx_cbind
  count <- data_info$count
  num_groups <- data_info$num_groups
  num_null <- data_info$num_null
  n <- data_info$n
  ref_n <- data_info$ref_n
  p <- data_info$p
  num_in_group <- data_info$num_in_group
  
  get_cor <- function(data, i){
    sub_data <- data[i, ]
    return(atanh(cor(sub_data)))
  }
  
  F_r <- atanh(cor(data_raw, use = "pairwise.complete.obs"))
  F_r_ref <- atanh(cor(data_ref, use = "pairwise.complete.obs"))
  
  cores=detectCores()
  bootobject <- boot(data = data_ref, statistic = get_cor, R=R, parallel = "multicore", ncpus = cores)
  
  bootsample <- array(0, c(R, p, p))
  for(i in 1:R){
    bootsample[i, ,] = matrix(bootobject$t[i, ], nrow=p)
  }
  
  rm(bootobject)
  
  group_idx <- NULL
  for(i in 2:p){
    for(j in 1:(i-1)){
      group_idx <- rbind(group_idx, c(i, j, F_r_ref[i, j]))
    }
  }
  group_idx <- data.frame(group_idx)
  
  colnames(group_idx) <- c("row_idx", "col_idx", "F_r_value")
  group_idx <- group_idx[order(group_idx$F_r_value), ]
  group_idx$group_idx <- as.numeric(split_quantile(x = 1:ttl, type = num_groups))
  
  
  FAB_ps <- matrix(0, p, p)
  UMPU_ps <- matrix(0, p, p)
  
  
  for(idx in 1:num_groups){
    group_subset <- subset(group_idx, group_idx==idx)
    group_mtx <- tensor_to_group_mtx(bootsample, group_subset)
    group_Cor_F_r <- cor(t(group_mtx))
    group_F_r <- F_r[cbind(group_subset$row_idx, group_subset$col_idx)]
    group_p <- dim(group_subset)[1]
    group_result <- getFABnUMPUpvals(group_Cor_F_r, group_F_r, n, group_p)
    
    for(i in 1:group_p){
      row_idx = group_subset$row_idx[i]
      col_idx = group_subset$col_idx[i]
      UMPU_ps[row_idx, col_idx] = group_result[1, i]
      UMPU_ps[col_idx, row_idx] = group_result[1, i]
      FAB_ps[row_idx, col_idx] = group_result[2, i]
      FAB_ps[col_idx, row_idx] = group_result[2, i]
    }
    if(idx %% 500 == 0){
      print(idx)
    }
  }
  
  results <- hash()
  results[["UMPU_ps"]]<- UMPU_ps
  results[["FAB_ps"]]<- FAB_ps
  results[["Sigma"]]<- Sigma
  results[["null_idx_cbind"]]<- null_idx_cbind
  results[["alter_idx_cbind"]]<- alter_idx_cbind
  
  
  return(results)  
}

