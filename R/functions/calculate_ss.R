library(vegan)
library(dplyr)
library(tidyr)
library(adespatial)
library(reshape2)
library(purrr)
library(ape)
library(gtools)
library(circular)
source("sTURN_simulation_runs/functions/Movement_functions.R")
source("sTURN_simulation_runs/functions/pcoa.mod.R")

summary_stats <- function(Model.output, sub.df){
  
  #Model.output <- clean_data$community_long 
  
  #sub.df<-clean_data$community
  
  #sub.df$time <- as.numeric(sub.df$time)
  
  sub.com<-sub.df %>%
    ungroup() %>% 
    dplyr::select(-time, -patch) 
  
  sub.com$`a` <- 0
  
  if(nrow(sub.com) == 1){
    N <- renyi(sub.com, scales = 0:2, hill = TRUE)
  } else{
    N <- colMeans(renyi(sub.com, scales = 0:2, hill = TRUE), na.rm = TRUE)
  }
  
  J <- log(exp(vegan::diversity(sub.com, index = 'shannon')/vegan::diversity(sub.com, index = 'invsimpson'))) %>% mean()
  
  
  alpha<- sub.com %>% 
    renyi(scales = 0:2,hill = TRUE)
  
  alpha[is.na(alpha)]<-0
  
  alpha$time<-sub.df$time
  alpha$patch<-sub.df$patch
  
  alpha_mean<-data.frame(alpha) %>% 
    ungroup() %>%
    dplyr::select(-patch) %>% 
    filter(time == max(time)) %>% 
    summarise_all(mean,na.rm=TRUE) %>% 
    ungroup() %>% 
    dplyr::select(-time)
  
  gamma<-Model.output %>% 
    filter(time==max(time)) %>% 
    group_by(time,species) %>%
    dplyr::summarise(abundance = sum(abundance)) %>% 
    spread(key = species, value = abundance) %>% 
    ungroup() %>% 
    dplyr::select(-time) %>% 
    renyi(scales = 0:2,hill = TRUE)
  
  beta<-gamma/alpha_mean
  
  alpha_time_mean<-data.frame(alpha) %>% 
    group_by(patch) %>% 
    summarise_all(mean,na.rm=TRUE) %>% 
    ungroup() %>% 
    dplyr::select(-patch,-time)
  
  gamma_time<-Model.output %>% 
    group_by(patch,species) %>%
    dplyr::summarise(abundance = sum(abundance)) %>% 
    spread(key = species, value = abundance, fill = 0) %>% 
    ungroup() %>% 
    dplyr::select(-patch) %>% 
    renyi(scales = 0:2,hill = TRUE)
  
  gamma_time[is.na(gamma_time)]<-0
  
  beta_time<-colMeans(gamma_time/alpha_time_mean, na.rm=TRUE)
  beta_space<-colMeans(beta,na.rm=TRUE)
  
  a_CV<- Model.output %>% 
    group_by(patch, time) %>% 
    dplyr::summarise(biomass = sum(abundance)) %>% 
    ungroup() %>% 
    group_by(patch) %>% 
    dplyr::summarise(a_CV = sd(biomass,na.rm=TRUE)/mean(biomass,na.rm=TRUE)) %>% 
    ungroup() %>% 
    dplyr::summarise(a_CV = mean(a_CV,na.rm=TRUE))
  
  g_CV<-Model.output %>% 
    group_by(time) %>% 
    dplyr::summarise(biomass = sum(abundance)) %>% 
    dplyr::summarise(g_CV = sd(biomass,na.rm=TRUE)/mean(biomass, na.rm = TRUE))
  
  min_max_abun <- Model.output %>%
    group_by(species) %>% 
    dplyr::summarise(mean_abundance = mean(abundance)) %>% 
    filter(mean_abundance > 0) %>% 
    summarise(max_abun = max(mean_abundance), min_abun =  min(mean_abundance)) %>% 
    mutate(ratio_max_min = max_abun/min_abun)
  
  prop_patches <- Model.output %>% 
    group_by(species, patch) %>%
    dplyr::summarise(mean_abundance = mean(abundance)) %>% 
    filter(mean_abundance  > 0) %>% 
    dplyr::summarise(n_patches = n()) %>% 
    mutate(prop_patches = n_patches/max(Model.output$patch)) %>% 
    summarise(mean_prop_patches = mean(prop_patches), min_prop_patches = min(prop_patches), max_prop_patches = max(prop_patches))
  
  beta_div_comp_patch <- Model.output %>% 
    group_by(time,species) %>%
    dplyr::summarise(abundance = sum(abundance)) %>% 
    spread(key = species, value = abundance, fill = 0) %>% 
    ungroup() %>% 
    dplyr::select(-time) %>% 
    beta.div.comp(coef = "J", quant = TRUE, save.abc = FALSE)
  
  
  beta_div_comp_time <- Model.output %>% 
    group_by(patch,species) %>%
    dplyr::summarise(abundance = sum(abundance)) %>% 
    spread(key = species, value = abundance, fill = 0) %>% 
    ungroup() %>% 
    dplyr::select(-patch) %>% 
    beta.div.comp(coef = "J", quant = TRUE, save.abc = FALSE)

  
  output <- c('J' = J, beta_space, beta_time, a_CV, g_CV, min_max_abun, prop_patches, 
              beta_div_comp_patch$part, beta_div_comp_time$part)
  names(output)<-c("J","B_s_0","B_s_1","B_s_2","B_t_0","B_t_1","B_t_2", "alpha_bio_cv", "gamma_bio_cv", 
                   "max_abun", "min_abun", "ratio_max_min", "mean_prop_patches", "min_prop_patches", "max_prop_patches", 
                   "BDtotal_patch", "Repl_patch", "RichDif_patch", "Repl_BDtotal_patch", "RichDi_BDtotal_patch",
                   "BDtotal_time", "Repl_time", "RichDif_time", "Repl_BDtotal_time", "RichDi_BDtotal_time")
  return(output)
}





PCoA_trajectories <- function(Model.output, sub.df, tot_env){
  
  #Model.output <- clean_data$community 
  
  #sub.df<-clean_data$community
 
  #sub.df$time <- as.numeric(sub.df$time)
  
  sub.com<-sub.df %>%
    ungroup() %>% 
    dplyr::select(-time, -patch) 
  
  pcoa1<-pcoa.mod(vegdist(sub.com,method = "bray"))
  
  pcoaVects<-c(1:2)
  if(dim(pcoa1$vectors)[2]<2){
    pcoaVects<-1:dim(pcoa1$vectors)[2]
  }
  
  time_all <- sub.df%>% 
    dplyr::select(time, patch)
  
  env_all <-tot_env %>% 
    rename(env1 = env) %>% 
    mutate(env2 = 0)
  
  env_all$time <- as.numeric(env_all$time)
  
  p_all <-data.frame(pcoa1$vectors[,pcoaVects], time_all)
  
  mvt_properties <- list()
  
  env_list <- list()
  
  patch_vec <-p_all$patch %>%unique()
  
  for(i in patch_vec){
    
    p1 <- p_all %>% 
      filter(patch == i) %>% 
      dplyr::select(-time, -patch)
    
     if(nrow(p1) < 3){next}    
  
     if(nrow(unique(p1)) < 2){next}
    


    time <- time_all %>% 
      filter(patch == i) %>% 
      dplyr::select(time)
    
    env.mat <- env_all %>% 
      filter(patch == i)
    
    env_dist <- melt(as.matrix(vegdist(env.mat, "euclidean")) , varnames = c("Time1", "Time2")) %>%
      dplyr::rename(env.dist = value) %>%
      mutate(one_step = ifelse(Time2 == Time1+1, TRUE, FALSE)) %>%
      filter(one_step == TRUE)
    
    ##pcoa_environment<-pcoa(vegdist(env.mat, "euclidean"))
    ## IF there are more than 2 environmental variables, use pcoa of environment
    
    env_list[[i]] <- cbind(env.mat, time)
    
    change_in_xx <- env.mat$env1[-length(env.mat$env1)]-env.mat$env1[-1]
  
    if(any(change_in_xx == 0)){
      change_in_xx[which(change_in_xx == 0)] <- 0.00001
    }
    
    environmental_angle <- anglefun(change_in_xx, env.mat$env2[-length(env.mat$env2)]-env.mat$env2[-1],as.deg = TRUE)
    
     change_in_xx_p1 <- p1$Axis.1[-length(p1$Axis.1)]-p1$Axis.1[-1]
    
    if(any(change_in_xx_p1 == 0)){
      change_in_xx_p1[which(change_in_xx_p1 == 0)] <- 0.00001
    }
  
    abs_angle <- anglefun(change_in_xx_p1, p1$Axis.2[-length(p1$Axis.2)]-p1$Axis.2[-1], as.deg = TRUE)
    
    autocor_angle <- as.vector(acf(abs_angle, na.action = na.pass, plot = FALSE)$acf)
    
    rel_angle <- rel.angle(abs_angle)
    
    step_len <- step_length(p1$Axis.1, p1$Axis.2)
    
    environmental_step_len <- step_length(env.mat$env1, env.mat$env2)
    
    net_disp <- net_displacement(p1$Axis.1, p1$Axis.2)
    
    gross_disp <- round(cumsum(step_len),2)
    
    step_env_disp <- step_len[-length(step_len)]/env_dist$env.dist
    step_env_disp[is.infinite(step_env_disp)]<-NA
    
    net_env_std_disp <- sqrt(net_disp[-1])/env_dist$env.dist
    net_env_std_disp[is.infinite(net_env_std_disp)]<-NA
    
    mvt_properties[[i]] <- list(max_net_disp = max(sqrt(net_disp), na.rm=T),
                                net_disp_final = sqrt(net_disp[length(net_disp)]),
                                mean_net_disp = mean(sqrt(net_disp), na.rm=T),
                                max_step = max(step_len, na.rm=T),
                                #min_step = round(min(step_len, na.rm=T),2),
                                sd_step = sd(step_len, na.rm=T),
                                sd_env_disp = sd(step_env_disp, na.rm=T),
                                max_env_disp = max(step_env_disp, na.rm=T),
                                min_env_disp = min(step_env_disp, na.rm=T), 
                                sd_env_disp = sd(net_env_std_disp, na.rm=T),
                                max_env_disp = max(net_env_std_disp, na.rm=T),
                                min_env_disp = min(net_env_std_disp, na.rm=T),
                                abs_angle = mean(as.circular(abs_angle), na.rm = T),
                                abs_angle_var = angular.variance(abs_angle, na.rm = TRUE),
                                #rel_angle = round(mean(rel_angle, na.rm = T),2),
                                auto_cor_angle1 = autocor_angle[2],
                                angle_correlation = cor.circular(environmental_angle, abs_angle),
                                step_correlation = tryCatch(cor(environmental_step_len, step_len, use = "complete.obs"), error=function(e) NA),
                                step_correlation = tryCatch(cor(environmental_step_len[-length(environmental_step_len)], step_len[-1], use = "complete.obs"), error=function(e) NA)
    )
  }
  
  
  output <- mvt_properties %>% 
    map_df(~ as.data.frame(.x)) %>% 
    summarise_all(funs(mean = mean, sd = sd, min = min, max = max),  na.rm = TRUE)
  
  return(output)
}


cors_ss <- function(sub.df, tot_env, xy){

  x <- c(sub.df[, "time", drop = FALSE])$time
  
  x <- as.numeric(x - min(x))
  
  composition_distance <- sqrt(vegdist(sub.df[, -c(1,2)], method = 'bray'))
  env_distance <- vegdist(scale(tot_env[, -c(1,2)]), method = 'euclidean')
  spatial_distance <- vegdist(xy, method = 'euclidean')
  temporal_distance <- vegdist(x, method = 'euclidean')
  
  cor_comp_env <- cor(composition_distance, env_distance) 
  cor_comp_space <- cor(composition_distance, spatial_distance)
  cor_comp_time <- cor(c(composition_distance), temporal_distance)
  cor_env_space <- cor(env_distance, spatial_distance)
  cor_env_time <- cor(env_distance, temporal_distance)
  
  cor_output <- c("cor_comp_env"= cor_comp_env, "cor_comp_space" = cor_comp_space, 
                  "cor_comp_time" = cor_comp_time, "cor_env_space" = cor_env_space, 
                  "cor_env_time" = cor_env_time)
  
  return(cor_output)
}


abun_occupancy <- function(sub.df){
  
  #sub.df<-clean_data$community
  
  #sub.df$time <- as.numeric(sub.df$time)

  
  sub.df$`a` <- 0
  
  tp <- unite(sub.df ,"time_patch", c("time", "patch"), sep = "_")
  
  tp_l <- dplyr::select(tp, -time_patch)
  
  abun_hill <- apply(tp_l, 1, FUN = function(x) renyi(x, scales = 0:2, hill = TRUE))
  abun_r1 <- apply(abun_hill, 2, FUN = function(x) x[2]/x[1])
  abun_r2 <- apply(abun_hill, 2, FUN = function(x) x[3]/x[1])
  
  
  abun_ratios <- data.frame(time_patch = tp$time_patch, ratio_1.1 = abun_r1, ratio_2.2 = abun_r2) %>% 
    separate(time_patch, c("time", "patch"), sep = "_") %>% 
    group_by(time) %>% 
    summarise(ratio1_0 = mean(ratio_1.1), ratio2_0 = mean(ratio_2.2))
  
  sub_p_t <-dplyr::select(sub.df, -time, -patch)
  
  sub_p_t[sub_p_t > 0] <- 1
  
  sub_p_t$time <- sub.df$time
  
  ocup_mat <- sub_p_t %>% 
    group_by(time) %>% 
    summarise_all(sum) %>% 
    ungroup() %>% 
    dplyr::select(-time)
  
  ocup_hill <- apply(ocup_mat, 1, FUN = function(x) renyi(x, scales = 0:2, hill = TRUE))
  ocup_r1 <- apply(ocup_hill, 2, FUN = function(x) x[2]/x[1])
  ocup_r2 <- apply(ocup_hill, 2, FUN = function(x) x[3]/x[1])
  
  occup_ratios <- data.frame(time = unique(sub.df$time), ratio1_0 = ocup_r1, ratio2_0 = ocup_r2)
  
  
  sumed_abun <- sub.df %>% 
    dplyr::select(-patch) %>% 
    group_by(time) %>% 
    summarise_all(sum) %>% 
    ungroup() %>% 
    dplyr::select(-time)
  
  gamma_hill <- apply(sumed_abun, 1, FUN = function(x) renyi(x, scales = 0:2, hill = TRUE))
  gamma_r1 <- apply(gamma_hill, 2, FUN = function(x) x[2]/x[1])
  gamma_r2 <- apply(gamma_hill, 2, FUN = function(x) x[3]/x[1])
  
  gamma_ratios <- data.frame(time = unique(sub.df$time), ratio1_0 = gamma_r1, ratio2_0 = gamma_r2)
  
  all_ratios <- list(abundance_ratios = abun_ratios, occupancy_ratios = occup_ratios, gamma_ratios = gamma_ratios)
  
  return(all_ratios)  
}


variation_partitioning_space <- function(sub.df,  mem_to_use, tot_env){
  
  community_only <- sub.df %>% 
    filter(time == max(time)) %>% 
    dplyr::select(-time, -patch)
  
  env_only <- tot_env %>% 
    filter(time == max(time)) %>% 
    dplyr::select(-time, -patch)
  
  mem_final <- bind_cols(dplyr::select(tot_env, time, patch), mem_to_use) %>% 
    filter(time == max(time)) %>% 
    dplyr::select(-time, -patch)
  
  Comhel<-decostand(community_only,"hellinger")
  
  v_p <- varpart(Comhel,  mem_final, env_only)
  
  return(v_p)
}



variation_partitioning_space_time <- function(sub.df,  mem_to_use, tot_env){
  
  community_only <- sub.df %>% 
    dplyr::select(-time, -patch)
  
  env_only <- tot_env %>% 
    dplyr::select(-time, -patch)
  
  mem_final <-mem_to_use
  
  no.timepoints <- max(sub.df$time)
  
  aem1<-aem.time(no.timepoints, w = NULL, moran = FALSE)
  
  aems<-aem1$aem
  
  order_aems <- sub.df %>% 
    left_join(data.frame(time= 1:no.timepoints, aems)) %>% 
    dplyr::select(starts_with("X"))
  
  Comhel<-decostand(community_only,"hellinger")
  
  v_p <- varpart(Comhel,  mem_final, env_only, order_aems)
  
  return(v_p)
}

