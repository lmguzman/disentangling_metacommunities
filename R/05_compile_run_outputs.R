library(dplyr)
library(stringr)


library(dplyr)
library(stringr)


#### all descriptive stats including time vp

rds_files <- list.files("intermediary_outputs/SS_no_VP/")

all_res <- data.frame()

for(i in 1:length(rds_files)){
  
  file1 <- readRDS(paste0("intermediary_outputs/SS_no_VP/",rds_files[i]))
  
  a_ratios <-  colMeans(file1[[3]]$abundance_ratios[,2:3])
  
  names(a_ratios) <- paste0("abundance_ratios_", c("ratio_1_0", "ratio_2_0"))
  
  o_ratios <-  colMeans(file1[[3]]$occupancy_ratios[,2:3])
  
  names(o_ratios) <- paste0("occupancy_ratios_", c("ratio_1_0", "ratio_2_0"))
  
  g_ratios <-  colMeans(file1[[3]]$gamma_ratios[,2:3])
  
  names(g_ratios) <- paste0("gamma_ratios_", c("ratio_1_0", "ratio_2_0"))
  
  time_var_part <- file1[[4]]$part$indfrac$Adj.R.square
  
  time_var_part[time_var_part < 0] <- 0
  
  names(time_var_part) <- c("time_var_part_Space", "time_var_part_Environment", "time_var_part_Time", "time_var_part_space_environ", "time_var_part_environ_time", "time_var_part_time_space", "time_var_part_all", "time_var_part_residual")
  
  params <- unlist(str_split(rds_files[i], "_"))
  
  res <- c(rep=params[2], dispersal = params[4], sig_niche = params[7], alpha = params[9], unlist(file1[[1]]), file1[[2]], a_ratios, o_ratios, g_ratios, time_var_part)
  
  all_res <- bind_rows(all_res, res)
}

write.csv(all_res, "results/desc_stat_time_vp.csv")




## space vp ##


library(dplyr)
library(stringr)

rds_files <- list.files("intermediary_outputs/vp/")

all_res <- data.frame()

for(i in 1:length(rds_files)){
  
  file1 <- readRDS(paste0("intermediary_outputs/vp/",rds_files[i]))
  
  space_var_part <- file1$part$indfrac$Adj.R.squared
  
  space_var_part[space_var_part < 0] <- 0
  
  names(space_var_part) <- c("space_var_part_Space", "space_var_part_both", "space_var_part_Environment", "space_var_part_Residual")
  
  params <- unlist(str_split(rds_files[i], "_"))
  
  res <- c(rep=params[2], dispersal = params[4], sig_niche = params[7], alpha = params[9], space_var_part)
  
  all_res <- bind_rows(all_res, res)
}

write.csv(all_res, "results/vp.csv")




### HMSC ###

library(dplyr)
library(stringr)

rds_files <- list.files("intermediary_outputs/HMSC_res/")

all_res <- data.frame()

for(i in 1:length(rds_files)){
  
  file1 <- readRDS(paste0("intermediary_outputs/HMSC_res/",rds_files[i]))
  

  if(!is.list(file1)) next
  
  hmsc_r2 <- mean(file1[[1]]$SR2, na.rm = TRUE)
  
  names(hmsc_r2) <- c("hmsc_r2")
  
  R2 <- file1[[1]]$SR2
  
  R2_1 <- R2
  
  R2_1[R2_1<1] <- 1
  
  hmsc_var_part_raw <- rowMeans(file1[[2]]$vals)
  
  names(hmsc_var_part_raw) <- paste0(names(hmsc_var_part_raw), "_raw")
  
  hmsc_var_part_std <- rowMeans(sapply(1:length(file1[[1]]$SR2), FUN = function(x) file1[[2]]$vals[,x]*R2[x]), na.rm = TRUE)
  
  names(hmsc_var_part_std) <- paste0(names(hmsc_var_part_std), "_std")
  
  file1_3_1 <- file1[[3]][[1]]$mean
  
  p1a <- file1_3_1[upper.tri(file1_3_1)]
  
  supportLevel = 0.95
  
  POS_space_time = (sum(file1[[3]][[1]]$support>(supportLevel))-file1[[6]])/(file1[[6]]*file1[[6]]-file1[[6]])
  
  NEG_space_time = (sum(file1[[3]][[1]]$support<(1-supportLevel)))/(file1[[6]]*file1[[6]]-file1[[6]])
  
  ## mean of positives and mean of negatives, mean absolute on all the things 
  space_time_associations <- c(summary(p1a), var = var(p1a), prop_high = sum(p1a>0)/length(p1a), prop_lower = sum(p1a<0)/length(p1a),
                               prop_zero = sum(p1a==0)/length(p1a), mean_abs = mean(abs(p1a)), 
                               mean_pos = mean(p1a[p1a>0]), mean_neg = mean(p1a[p1a<0]))
  
  names(space_time_associations) <- paste("hmsc_space_time_ass", names(space_time_associations), sep = "_")
  
  file1_3_2 <- file1[[3]][[2]]$mean
  
  s1a <- file1_3_2[upper.tri(file1_3_2)]
  
  POS_space = (sum(file1[[3]][[2]]$support>(supportLevel))-file1[[6]])/(file1[[6]]*file1[[6]]-file1[[6]])
  
  NEG_space = (sum(file1[[3]][[2]]$support<(1-supportLevel)))/(file1[[6]]*file1[[6]]-file1[[6]])
  
  space_associations <- c(summary(s1a), var = var(s1a), prop_high = sum(s1a>0)/length(s1a), prop_lower = sum(s1a<0)/length(s1a),
                          prop_zero = sum(s1a==0)/length(s1a), mean_abs = mean(abs(s1a)), 
                          mean_pos = mean(s1a[s1a>0]), mean_neg = mean(s1a[s1a<0]))
  
  names(space_associations) <- paste("hmsc_space_ass", names(space_associations), sep = "_")
  
  file1_3_3 <- file1[[3]][[3]]$mean
  
  t1a <- file1_3_3[upper.tri(file1_3_3)]
  
  POS_time = (sum(file1[[3]][[3]]$support>(supportLevel))-file1[[6]])/(file1[[6]]*file1[[6]]-file1[[6]])
  
  NEG_time = (sum(file1[[3]][[3]]$support<(1-supportLevel)))/(file1[[6]]*file1[[6]]-file1[[6]])
  
  
  time_associations <- c(summary(t1a), var = var(t1a), prop_high = sum(t1a>0)/length(t1a), prop_lower = sum(t1a<0)/length(t1a),
                         prop_zero = sum(t1a==0)/length(t1a), mean_abs = mean(abs(t1a)), 
                         mean_pos = mean(t1a[t1a>0]), mean_neg = mean(t1a[t1a<0]))
  
  
  names(time_associations) <- paste("hmsc_time_ass", names(time_associations), sep = "_")
  
  params <- unlist(str_split(rds_files[i], "_"))
  
  MEA_space_time <- mean(file1[[4]][[1]][[1]])
  MEA_space <- mean(file1[[4]][[2]][[1]])
  MEA_time <- mean(file1[[4]][[3]][[1]])
  
  SS_space_time <- mean(file1[[4]][[1]][[1]]>0)
  SS_space <- mean(file1[[4]][[2]][[1]]>0)
  SS_time <- mean(file1[[4]][[3]][[1]]>0)
  
  res <- c(rep=params[2], dispersal = params[4], sig_niche = params[7], alpha = params[9],
           hmsc_r2, hmsc_var_part_raw, hmsc_var_part_std, space_time_associations, space_associations, time_associations,
           POS_space_time = POS_space_time, NEG_space_time = NEG_space_time, POS_space = POS_space, 
           NEG_space = NEG_space, POS_time = POS_time, NEG_time = NEG_time, MEA_space_time = MEA_space_time,
           MEA_space = MEA_space, MEA_time= MEA_time, SS_space_time = SS_space_time, 
           SS_space = SS_space, SS_time = SS_time)
  

  all_res <- bind_rows(all_res, res)
  print(i)
}

write.csv(all_res, "results/hmsc.csv")




### compile gamma diversity ###

n_gamma <- data.frame()

for(r in 1:15){
  
  load(paste0("simulation_outputs/outputfile_",r,".RData"))
  
  r1_gamma <- model.df %>% 
    group_by(dispersal, sig_niche, alpha, Species) %>% 
    summarise(n()) %>% 
    ungroup() %>% 
    group_by(dispersal, sig_niche, alpha) %>% 
    summarise(n_sp = n()) %>% 
    mutate(rep = r) %>% 
    ungroup()
  
  n_gamma <- bind_rows(n_gamma, r1_gamma)
  
}

write.csv(n_gamma, "results/n_gamma.csv")



