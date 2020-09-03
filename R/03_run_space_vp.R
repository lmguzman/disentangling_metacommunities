library(dplyr)
library(tidyr)
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

# load parameter table
load("R/pars.Rdata")

# create parameter object for each task
p <- pars[task.id,]

# Model code

source("functions/calculate_ss.R")

source("functions/prep_sim_ss.R")

r <- p$r

landscape <- read.csv(paste0("simulation_outputs/landscape_",r,".csv"))
load(paste0("simulation_outputs/outputfile_",r,".RData"))

Model.sel <- model.df[which(round(model.df$dispersal, 5) == round(p$dispersal, 5) & round(model.df$sig_niche, 5) == round(p$sig_niche, 5) &  model.df$alpha == p$alpha), ]

res <- var_p(Model.sel, landscape)
  
# Save output to a .Rdata file

saveRDS(res,file=paste0('intermediary_outputs/vp/rep_', r, '_dispersal_', p$dispersal, '_sig_niche_', p$sig_niche, '_alpha_', p$alpha, '_.rds'))



