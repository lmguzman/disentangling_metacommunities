print("started script")

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
print(paste("task.id =", task.id))

# load parameter table
load("R/pars.Rdata")

head(pars)

# create parameter object for each task
p <- pars[task.id,]

print(paste("p",p))

# Model code

source("R/functions/HMSC_fun.R")

r <- p$r

print(paste("r",r))
landscape <- read.csv(paste0("simulation_outputs/landscape_",r,".csv"))
load(paste0("simulation_outputs/outputfile_",r,".RData"))

run_one <- model.df[which(round(model.df$dispersal, 5) == round(p$dispersal, 5) & round(model.df$sig_niche, 5) == round(p$sig_niche, 5) &  model.df$alpha == p$alpha), ]

run_one <- run_one[which(run_one$Time %in% seq(1,max(run_one$Time),3)),]

sp_to_remove <- run_one %>% 
  group_by(Species, Time) %>% 
  summarise(n()) %>% 
  ungroup() %>% 
  group_by(Species) %>% 
  summarise(n_sp = n()) %>% 
  filter(n_sp == 1)

run_one <- filter(run_one, !(Species %in% sp_to_remove$Species))

landscape2 <- landscape[,-1]

hms_res <- hmsc_fun(run_one, landscape2)
  
# Save output to a .Rdata file

saveRDS(hms_res,file=paste0('intermediary_outputs/HMSC_res/rep_', r, '_dispersal_', p$dispersal, '_sig_niche_', p$sig_niche, '_alpha_', p$alpha, '_.rds'))



