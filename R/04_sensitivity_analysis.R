library(dplyr)
library(tidyr)
library(Hmsc)
library(plotrix)
library(adespatial)
library(data.table)
library(som.nn)

# R script

# Passing arguments to R from command lines
args = commandArgs(trailingOnly=TRUE)
output.file <- args[1]

# try to get SGE_TASK_ID from submit script, otherwise fall back to 1
task.id = as.integer(Sys.getenv("SGE_TASK_ID", "1"))

r <- 1
landscape <- read.csv(paste0("simulation_outputs/landscape_",r,".csv"))
load(paste0("simulation_outputs/outputfile_",r,".RData"))

# load parameter table
load("R/combs.RData")

# create parameter object for each task
cur_combr <- combs[task.id,]

source("functions/calculate_ss.R")

source("functions/prep_sim_ss.R")

mid_d <- unique(round(model.df$dispersal,5))[8]
mid_s <- unique(round(model.df$sig_niche,5))[6]

run_one <- model.df[which(round(model.df$dispersal, 5) == mid_d & round(model.df$sig_niche, 5) == mid_s &  model.df$alpha == "0.5"), ]

multi_res <-list()

for(i in 1:1000){
  
  patches <- sample(max(run_one$Patch),cur_combr$Patch)
  
  landscape <- filter(all_land, X %in% patches)
  
  Model.sel <- filter(run_one, Time <= cur_combr$Time, Patch %in% patches)
  
  res <- all_ss(Model.sel, landscape)
  
  all_res <- list(patches, res)
  
  multi_res[[i]] <- all_res
}


saveRDS(multi_res,file=paste0('intermediary_outputs/sensitivity/time_', cur_combr$Time, '_patch_', cur_combr$Patch,'_.rds'))
