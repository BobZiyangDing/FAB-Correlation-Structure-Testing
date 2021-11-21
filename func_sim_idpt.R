# pkg_test("Rcpp")
# pkg_test("fabricatr")
# pkg_test("mltools")
# pkg_test("onehot")
# pkg_test("mgcv")
# pkg_test("tidyverse")
# pkg_test("foreach")
# pkg_test("doParallel")
# pkg_test("kableExtra")
# pkg_test("boot")
# pkg_test("tensorr")
# pkg_test("reshape2")
# pkg_test("dvmisc")
# pkg_test("tidyverse")
# pkg_test("hash")
# source("Utils.R")
# setwd("~/FABBrainConnectome/FAB_corr")
options(warn=-1)

func_sim_idpt <- function(data_info){

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
      umpu_p <- 1 - abs(pnorm(group_direct_F_r[j]*sqrt(n-3) , 0, 1) - 
                          pnorm(-group_direct_F_r[j]*sqrt(n-3), 0, 1) )
      
      row_idx = group_subset$row_idx[j]
      col_idx = group_subset$col_idx[j]
      FAB_ps[row_idx, col_idx] = fab_p
      FAB_ps[col_idx, row_idx] = fab_p
      UMPU_ps[row_idx, col_idx] = umpu_p
      UMPU_ps[col_idx, row_idx] = umpu_p
    }
    if(idx %% 100 == 0){
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





