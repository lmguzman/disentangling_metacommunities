library(dplyr)
library(stringr)

rds_files <- list.files("intermediary_outputs/sensitivity")

all_res <- data.frame()

for(i in 1:length(rds_files)){
  
  file1 <- readRDS(paste0("intermediary_outputs/sensitivity/",rds_files[i]))
  
  for(r in 1:1000){
    
    a_ratios <-  colMeans(file1[[r]][[2]][[3]]$abundance_ratios[,2:3])
    
    names(a_ratios) <- paste0("abundance_ratios_", c("ratio_1_0", "ratio_2_0"))
    
    o_ratios <-  colMeans(file1[[r]][[2]][[3]]$occupancy_ratios[,2:3])
    
    names(o_ratios) <- paste0("occupancy_ratios_", c("ratio_1_0", "ratio_2_0"))
    
    g_ratios <-  colMeans(file1[[r]][[2]][[3]]$gamma_ratios[,2:3])
    
    names(g_ratios) <- paste0("gamma_ratios_", c("ratio_1_0", "ratio_2_0"))
    
    time_var_part <- file1[[r]][[2]][[4]]$part$indfrac$Adj.R.square
    
    time_var_part[time_var_part < 0] <- 0
    
    names(time_var_part) <- c("time_var_part_Space", "time_var_part_Environment", "time_var_part_Time", "time_var_part_space_environ", "time_var_part_environ_time", "time_var_part_time_space", "time_var_part_all", "time_var_part_residual")
    
    params <- unlist(str_split(rds_files[i], "_"))
    
    res <- c(r = r, time=params[2], patch = params[4], unlist(file1[[r]][[2]][[1]]), file1[[r]][[2]][[2]], a_ratios, o_ratios, g_ratios, time_var_part)
    
    all_res <- bind_rows(all_res, res)
  }
  
}

write.csv(all_res, "results/sensitivity.csv")

