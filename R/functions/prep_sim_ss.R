all_ss <- function(Model.sel, land){

if(nrow(Model.sel)  < 2) {res <- c(space = NA, time = NA, environment = NA)} 

  mem_init <- dbmem(as.data.frame(land[,c("x","y")]),MEM.autocor = "positive",silent=TRUE)
  mem_all <- cbind(X = land$X, as.data.frame(mem_init))

  
  # species data
  sub.df<-Model.sel %>% 
    dplyr::select(env, time = Time, patch = Patch, N, Species) %>% 
    spread(key = Species, value = N, fill = 0) %>% as.data.frame() 
  
  sub.df <- sub.df[,colSums(sub.df) !=0]
  
  sub.df2<-Model.sel %>% 
    group_by(Patch) %>% 
    dplyr::summarise(presence=sum(N)) %>% 
    filter(presence>0) %>% 
    rename(patch = Patch)

 if(nrow(sub.df2) == 0){res <- c(space = NA, time = NA, environment = NA)}
  
  sub.df<-left_join(sub.df2,sub.df)
  
  # predictors
  sub.com<-sub.df %>%
    ungroup() %>% 
    dplyr::select(-time, -patch)
  
  f1<-sub.df %>% 
    dplyr::select(patch,time,env) %>% 
    left_join(land, by = c('patch' = "X")) %>% 
    left_join(mem_all, by = c('patch' = "X"))
  
  Model.output <- Model.sel %>% 
    rename(abundance = N, time = Time, patch = Patch, species = Species)

 sub.df_ss<-left_join(sub.df2,sub.df) %>% 
    dplyr::select(-env, -presence)
  
  env_tit <- dplyr::select(f1, patch, time, env)
  
  xy <- as.data.frame(f1[,c("x","y")])
  
  mem_to_use <- dplyr::select(f1, starts_with('MEM'))
  
  ss_div <- summary_stats(Model.output, sub.df_ss)
  
 # pcoa <- PCoA_trajectories(Model.output, sub.df_ss, env_tit)
  
  cors_general <- cors_ss(sub.df_ss, env_tit, xy)
  
  sad <- abun_occupancy(sub.df_ss)
 
 # varpart_space <- variation_partitioning_space(sub.df_ss,  mem_to_use, env_tit)
  
  varpart_space_time <- variation_partitioning_space_time(sub.df_ss,  mem_to_use, env_tit)
  
  res <- list(ss_div,  cors_general, sad, varpart_space_time) 

return(res) 


}





var_p <- function(Model.sel, land){

if(nrow(Model.sel)  < 2) {res <- c(space = NA, time = NA, environment = NA)}

  mem_init <- dbmem(as.data.frame(land[,c("x","y")]),MEM.autocor = "positive",silent=TRUE)
  mem_all <- cbind(X = land$X, as.data.frame(mem_init))


  # species data
  sub.df<-Model.sel %>%
    dplyr::select(env, time = Time, patch = Patch, N, Species) %>%
    spread(key = Species, value = N, fill = 0) %>% as.data.frame()

  sub.df <- sub.df[,colSums(sub.df) !=0]

  sub.df2<-Model.sel %>%
    group_by(Patch) %>%
    dplyr::summarise(presence=sum(N)) %>%
    filter(presence>0) %>%
    rename(patch = Patch)

 if(nrow(sub.df2) == 0){res <- c(space = NA, time = NA, environment = NA)}

  sub.df<-left_join(sub.df2,sub.df)

  # predictors
  sub.com<-sub.df %>%
    ungroup() %>%
    dplyr::select(-time, -patch)

  f1<-sub.df %>%
    dplyr::select(patch,time,env) %>%
    left_join(land, by = c('patch' = "X")) %>%
    left_join(mem_all, by = c('patch' = "X"))

Model.output <- Model.sel %>%
    rename(abundance = N, time = Time, patch = Patch, species = Species)

 sub.df_ss<-left_join(sub.df2,sub.df) %>%
    dplyr::select(-env, -presence)

  env_tit <- dplyr::select(f1, patch, time, env)

  xy <- as.data.frame(f1[,c("x","y")])

  mem_to_use <- dplyr::select(f1, starts_with('MEM'))

varpart_space <- variation_partitioning_space(sub.df_ss,  mem_to_use, env_tit)

return(varpart_space)
}
