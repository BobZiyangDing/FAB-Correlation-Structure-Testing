pkg_test("Rcpp")
pkg_test("fabricatr")
pkg_test("mltools")
pkg_test("onehot")
pkg_test("mgcv")
pkg_test("tidyverse")
pkg_test("foreach")
pkg_test("doParallel")
pkg_test("kableExtra")
pkg_test("boot")
pkg_test("tensorr")
pkg_test("reshape2")
pkg_test("dvmisc")
pkg_test("tidyverse")
pkg_test("hash")
sourceCpp("Cov_struct_func.cpp")
source("Utils.R")
setwd("~/FABBrainConnectome/FAB_corr")
options(warn=-1)

FAB_idpt <- function(data_raw, data_ref, num_in_group = 100, need_UMPU = FALSE){
  n <- dim(data_raw)[1]
  ref_n <- dim(data_ref)[1]
  p <- dim(data_raw)[2]
  ttl <- p*(p-1)/2
  num_groups <- ttl %/% num_in_group
  
  # UMPU F_r
  F_r <- atanh(cor(data_raw, use = "pairwise.complete.obs")) #, use = "pairwise.complete.obs"
  # external reference group F_r for Z score
  F_r_ref <- atanh(cor(data_ref, use = "pairwise.complete.obs")) #, use = "pairwise.complete.obs"
  
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
  
  get_group_F_r <- function(group_subset, F_r){
    group_index_list <- group_subset %>% dplyr::select(row_idx, col_idx) %>% as.matrix()
    group_F_r <- F_r[group_index_list]
    return(group_F_r)
  }
  
  FAB_ps <- matrix(0, p, p)
  
  if(need_UMPU) 
    UMPU_ps <- matrix(0, p, p)
  
  for(idx in 1:num_groups){
    group_subset <- subset(group_idx, group_idx==idx)
    group_direct_F_r <- get_group_F_r(group_subset, F_r)
    group_indirect_F_r <- get_group_F_r(group_subset, F_r_ref)
    
    group_p <- dim(group_subset)[1]
    
    for(j in 1:group_p){
      group_indirect_info <- group_indirect_F_r[-j]
      mu_MLE <- mean(group_indirect_info)
      psi2_MLE <- group_p/(group_p-1)*var(group_indirect_info) - 1/(ref_n-3) + 0.5
      
      tilde_m <- mu_MLE
      tilde_v <- psi2_MLE
      
      fab_p <- 1 - abs(pnorm(group_direct_F_r[j]*sqrt(n-3) + 2*tilde_m*sqrt(ref_n-3)/tilde_v , 0, 1) - 
                         pnorm(-group_direct_F_r[j]*sqrt(n-3), 0, 1) )
      
      row_idx = group_subset$row_idx[j]
      col_idx = group_subset$col_idx[j]
      FAB_ps[row_idx, col_idx] = fab_p
      FAB_ps[col_idx, row_idx] = fab_p
      
      if(need_UMPU){
        umpu_p <- 1 - abs(pnorm(group_direct_F_r[j]*sqrt(n-3) , 0, 1) - 
                            pnorm(-group_direct_F_r[j]*sqrt(n-3), 0, 1) )
        UMPU_ps[row_idx, col_idx] = umpu_p
        UMPU_ps[col_idx, row_idx] = umpu_p
      }
    }
    if(idx %% 100 == 0){
      print(idx)
    }
  }
  
  results <- hash()
  if(need_UMPU){
    results[["UMPU_ps"]]<- UMPU_ps
  }
  results[["FAB_ps"]]<- FAB_ps
  
  return(results)
}



FAB_boot <- function(data_raw, data_ref, num_in_group = 50, R = 2000, need_UMPU = FALSE){
  n <- dim(data_raw)[1]
  ref_n <- dim(data_ref)[1]
  p <- dim(data_raw)[2]
  ttl <- p*(p-1)/2
  num_groups <- ttl %/% num_in_group
  
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
  
  if(need_UMPU) 
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
      if(need_UMPU)
      {
        UMPU_ps[row_idx, col_idx] = group_result[1, i]
        UMPU_ps[col_idx, row_idx] = group_result[1, i]
      }
      FAB_ps[row_idx, col_idx] = group_result[2, i]
      FAB_ps[col_idx, row_idx] = group_result[2, i]
    }
    if(idx %% 500 == 0){
      print(idx)
    }
  }
  
  results <- hash()
  if(need_UMPU){
    results[["UMPU_ps"]]<- UMPU_ps
  }
  results[["FAB_ps"]]<- FAB_ps
  
  return(results)  
}


FAB_general_boot <- function(data_raw, data_ref, num_in_group = 50, shrinkage=0.05, R = 2000, need_UMPU = FALSE, 
                             powers = c(0, 1) ){
  n <- dim(data_raw)[1]
  ref_n <- dim(data_ref)[1]
  p <- dim(data_raw)[2]
  ttl <- p*(p-1)/2
  num_groups <- ttl %/% num_in_group
  
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
  
  if(need_UMPU) 
    UMPU_ps <- matrix(0, p, p)
  
  for(idx in 1:num_groups){
    group_subset <- subset(group_idx, group_idx==idx)
    group_mtx <- tensor_to_group_mtx(bootsample, group_subset)
    group_Cor_F_r <- cor(t(group_mtx))
    group_F_r <- F_r[cbind(group_subset$row_idx, group_subset$col_idx)]
    
    if( length( powers ) )
    {
      base <- as.matrix(F_r_ref[cbind(group_subset$row_idx, group_subset$col_idx)])
      group_W <- NULL
      for( power in powers )
      {
        group_W <- cbind(group_W, base ^ power )
      }
    }
    else
    {
      print( "Empty Linking Model" )
      return()
    }

    group_p <- dim(group_subset)[1]
    group_result <- getGeneralFABnUMPUpvals(Cor_F_r = group_Cor_F_r, 
                                            F_r = group_F_r, 
                                            one = group_W, 
                                            lambda = shrinkage, 
                                            n = n, 
                                            p = group_p)
    
    for(i in 1:group_p){
      row_idx = group_subset$row_idx[i]
      col_idx = group_subset$col_idx[i]
      if(need_UMPU)
      {
        UMPU_ps[row_idx, col_idx] = group_result[1, i]
        UMPU_ps[col_idx, row_idx] = group_result[1, i]
      }
      FAB_ps[row_idx, col_idx] = group_result[2, i]
      FAB_ps[col_idx, row_idx] = group_result[2, i]
    }
    if(idx %% 500 == 0){
      print(idx)
    }
  }
  
  results <- hash()
  if(need_UMPU){
    results[["UMPU_ps"]]<- UMPU_ps
  }
  results[["FAB_ps"]]<- FAB_ps
  
  return(results)  
}
