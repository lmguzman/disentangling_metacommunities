library(dplyr)
library(randomForest)
library(ggplot2)
library(tibble)
library(tidyr)
library(viridis)
library(cowplot)
library(caret)
library(mlbench)
library(e1071)

n_gamma <- read.csv('results/n_gamma.csv', stringsAsFactors = FALSE) %>% 
  select(-X) %>% 
  mutate(alpha = ifelse(alpha == 'patch_dynamics', 'patch', alpha))

n_space_var <-  read.csv("results/vp.csv") %>% 
  select(rep, dispersal, sig_niche, alpha, contains("space_var"))  %>% 
  group_by(alpha, sig_niche, dispersal) %>% 
  summarise(n_reps = n())

space_var <- read.csv("results/vp.csv") %>% 
  select(rep, dispersal, sig_niche, alpha, contains("space_var"))  %>% 
  left_join(n_space_var) %>% 
  filter(n_reps > 11) %>% 
  left_join(n_gamma) %>% 
  filter(n_sp >1) %>% 
  select(-n_reps, -n_sp) 
  


n_all_s <- read.csv("results/desc_stat_time_vp.csv") %>% 
  select(rep, dispersal, sig_niche, alpha) %>% 
  group_by(alpha, sig_niche, dispersal) %>% 
  summarise(n_reps = n())

all_ss2 <- read.csv("results/desc_stat_time_vp.csv") %>% 
  select(-X) %>% 
  left_join(n_all_s) %>% 
  filter(n_reps > 11) %>% 
  left_join(n_gamma) %>% 
  filter(n_sp >1) %>% 
  select(-n_reps, -n_sp)


n_hms_ss <- read.csv("results/hmsc.csv") %>% 
  select(rep, dispersal, sig_niche, alpha) %>% 
  group_by(alpha, sig_niche, dispersal) %>% 
  summarise(n_reps = n()) 

hms_ss <- read.csv("results/hmsc.csv") %>% 
  select(rep, dispersal, sig_niche, alpha,hmsc_r2, env_raw:`Random..Time_std`, contains('Median'), 
         contains('Mean'), contains('var'), contains('prop'), contains('POS'), contains('NEG'),
         contains('MEA'), SS_space_time, SS_space, SS_time) %>% 
  left_join(n_hms_ss) %>% 
  filter(n_reps > 11) %>% 
  left_join(n_gamma) %>% 
  filter(n_sp >1) %>% 
  select(-n_reps, -n_sp) 


all_outputs <- all_ss2  %>% 
  left_join(hms_ss, by = c("rep", "dispersal", "sig_niche", "alpha")) %>% 
  left_join(space_var, by = c("rep", "dispersal", "sig_niche", "alpha")) %>% 
  select(-max_abun , -min_abun, -ratio_max_min, -time_var_part_residual, -space_var_part_Residual)


name_df <- data.frame(original_names = names(all_outputs)[c(6:95)],
pretty_names = c("Beta diversity in space (Hill 0)", "Beta diversity in space (Hill 1)", "Beta diversity in space (Hill 2)",
  "Beta diversity in time (Hill 0)", "Beta diversity in time (Hill 1)","Beta diversity in time (Hill 2)",
  "Coefficient of variation \n in abundance at local scale", "Coefficient of variation \n in abundance at regional scale",
  "Mean proportion of \npatches occupied", "Minimum proportion of \npatches occupied", "Maximum proportion of patches occupied",
  "Total beta diversity in space", "Replacement component of beta diversity in space", "Total richness difference \ndiversity in space",
  "Total replacement diversity/ \n Total beta diversity in space", "Total richness difference diversity (or nestedness)/\n Total beta diversity in space","Total beta diversity in time", "Replacement component of beta diversity in time", "Total richness difference \ndiversity in time",
  "Total replacement diversity/\n Total beta diversity in time", "Total richness difference diversity (or nestedness)/\n Total beta diversity in time",
  "Correlation of compositional \n distance and environmental distance", "Correlation of compositional \n distance and spatial distance",
  "Correlation of compositional \n distance and temporal distance", "Correlation of environmental \ndistance and spatial distance", 
  "Correlation of environmental \ndistance and temporal distance", "Ratio of abundances at local \nscale (Hill 1/Hill 0)",
  "Ratio of abundances at local \nscale (Hill 2/Hill 0)", "Ratio of occupancies at local \nscale (Hill 1/ Hill 0)", 
  "Ratio of occupancies at local \nscale (Hill 2/ Hill 0)","Ratio of abundances at regional \nscale (Hill 1/Hill 0)",
  "Ratio of abundances at regional\n scale (Hill 2/Hill 0)", "Variation partitioning in time and \nspace (Space component)",
  "Variation partitioning in time and \nspace (Environment component)", "Variation partitioning in time and \nspace (Time component)",
  "Variation partitioning in time and \nspace (Space and environment component)", "Variation partitioning in time and \nspace (Environment and time component)",
  "Variation partitioning in time and \nspace (Space and time component)", "Variation partitioning in time and \nspace (Space, time and environment component)",
  "HMSC mean R2", "HMSC environment component raw", "HMSC space-time component raw", "HMSC space component raw", "HMSC time component raw", 
  "HMSC environment component standardized", "HMSC space-time component standardized", "HMSC space component standardized", "HMSC time component standardized",
  "HMSC median of species associations in space-time","HMSC median of species associations in space", "HMSC median of species associations in time",
  "HMSC mean of species associations in space-time", "HMSC mean of absolute species associations in space-time", "HMSC mean of positive species associations in space-time", "HMSC mean of negative species associations in space-time",
  "HMSC mean of species associations in space", "HMSC mean of absolute species associations in space", "HMSC mean of positive species associations in space", "HMSC mean of negative species associations in space",
  "HMSC mean of species associations in time", "HMSC mean of absolute species associations in time", "HMSC mean of positive species associations in time", "HMSC mean of negative species associations in time",
  "HMSC variance of species associations in space-time","HMSC variance of species associations in space", "HMSC variance of species associations in time",
  "HMSC proportion of species associations \nhigher than zero in space-time", "HMSC proportion of species associations \nlower than zero in space-time", "HMSC proportion of species associations \nequal to zero in space-time",
  "HMSC proportion of species associations \nhigher than zero in space", "HMSC proportion of species associations \nlower than zero in space", "HMSC proportion of species associations \nequal to zero in space",
  "HMSC proportion of species associations \nhigher than zero in time", "HMSC proportion of species associations \nlower than zero in time", "HMSC proportion of species associations \nequal to zero in space",
  "HMSC space-time POS", "HMSC space POS", "HMSC time POS", "HMSC space-time NEG", "HMSC space NEG", "HMSC time NEG",
  "HMSC space-time MEA", "HMSC space MEA", "HMSC time MEA", "HMSC space-time SS", "HMSC space SS", "HMSC time SS",
  "Variation partitioning in \nspace (Space component)", "Variation partitioning in \nspace (Space and environment component)",
  "Variation partitioning in \nspace (Environment component)")) %>% 
  mutate(original_names = as.character(original_names), pretty_names = as.character(pretty_names))


####### alpha ######

### space summary  ##

all_summary_data <- all_outputs %>%
  select(alpha, B_s_0:B_s_2, BDtotal_patch:RichDi_BDtotal_patch) 

nrow(all_summary_data)

all_summary_rf <- randomForest(as.factor(alpha)~ ., data = all_summary_data, importance = TRUE)



### all summary ##

all_summary_data <- all_outputs %>% 
  select(alpha, B_s_0:gamma_ratios_ratio_2_0) %>% 
  select(-contains('cor'))

nrow(all_summary_data)

all_summary_rf <- randomForest(as.factor(alpha)~ ., data = all_summary_data, importance = TRUE)



## var part space ##

var_part_data <- space_var %>% 
  select(alpha, space_var_part_Space, space_var_part_both, space_var_part_Environment) %>% 
  filter(!is.na(space_var_part_Space))
nrow(var_part_data)
var_part_orig <- randomForest(as.factor(alpha)~ ., data = var_part_data, importance = TRUE)

var_part_orig


## var part time ## 

var_part_time_data <- all_outputs %>% 
  select(alpha, starts_with("time_var_part")) %>% 
  filter(!is.na(time_var_part_Space))

nrow(var_part_time_data)
var_part_time <- randomForest(as.factor(alpha)~ ., data = var_part_time_data, importance = TRUE)


## HMSC in time ##

## without NAs

hmsc_time_data <- hms_ss %>%
  select(alpha, hmsc_r2:SS_time) %>% 
  filter(!is.na(hmsc_time_ass_mean_pos)) %>% 
  filter(!is.na(hmsc_space_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_time_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_space_time_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_space_time_ass_mean_pos)) %>% 
  filter(!is.na(hmsc_space_ass_mean_pos))

nrow(hmsc_time_data)

hmsc_time_rf <- randomForest(as.factor(alpha)~ ., data = hmsc_time_data, importance = TRUE)


## everything model
#without Nas
# without NAs
everything_data <- all_outputs %>% 
  select(alpha, B_s_0:space_var_part_Environment) %>% 
  select(-contains('cor')) %>% 
  filter(!is.na(hmsc_time_ass_mean_pos)) %>% 
  filter(!is.na(hmsc_space_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_time_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_space_time_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_space_time_ass_mean_pos)) %>% 
  filter(!is.na(hmsc_space_ass_mean_pos)) %>% 
  filter(!is.na(space_var_part_Environment))

summary(everything_data)

nrow(everything_data)

everything_rf <- randomForest(as.factor(alpha)~ ., data = everything_data, importance = TRUE)

## minimal model
## without NA

everything_data_cv <- all_outputs %>% 
  select(alpha, B_s_0:space_var_part_Environment) %>%
  select(-contains('cor')) %>% 
  filter(!is.na(hmsc_time_ass_mean_pos)) %>% 
  filter(!is.na(hmsc_space_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_time_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_space_time_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_space_time_ass_mean_pos)) %>% 
  filter(!is.na(hmsc_space_ass_mean_pos)) %>% 
  filter(!is.na(space_var_part_Environment))

everything_data_cv %>% 
  group_by(alpha) %>% summarise(min(alpha_bio_cv), mean(alpha_bio_cv), median(alpha_bio_cv))

everything_data_cv %>% 
  group_by(alpha) %>% summarise(min(RichDi_BDtotal_patch), mean(RichDi_BDtotal_patch), median(RichDi_BDtotal_patch))

# recursive feature elimitation
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results_alpha <- rfe(everything_data_cv[,-1], as.factor(everything_data_cv$alpha), sizes=c(10:50), rfeControl=control)

results_alpha$bestSubset

results_alpha$fit$importance

saveRDS(results_alpha, "results/minimal_alpha.rds")



####### dispersal ######

### space summary  ##

all_summary_data_dis <- all_outputs %>%
  select(dispersal, B_s_0:B_s_2, BDtotal_patch:RichDi_BDtotal_patch) 

nrow(all_summary_data_dis)

space_summary_rf_dis <- randomForest(dispersal~ ., data = all_summary_data_dis, importance = TRUE)


### summary that we have ##

all_summary_data_dis <- all_outputs %>% 
  select(dispersal, B_s_0:gamma_ratios_ratio_2_0) %>% 
  select(-contains('cor'))

nrow(all_summary_data_dis)

all_summary_rf_dis <- randomForest(dispersal~ ., data = all_summary_data_dis, importance = TRUE)


## var part space ##

var_part_data_dis <- space_var %>% 
  select(dispersal, space_var_part_Space, space_var_part_both, space_var_part_Environment) %>% 
  filter(!is.na(space_var_part_Space))

nrow(var_part_data_dis)

var_part_orig_dis <- randomForest(dispersal~ ., data = var_part_data_dis, importance = TRUE)

## var part time ## 

var_part_time_data_dis <- all_outputs %>% 
  select(dispersal, starts_with("time_var_part")) %>% 
  filter(!is.na(time_var_part_environ_time))

nrow(var_part_time_data_dis)

var_part_time_dis <- randomForest(dispersal~ ., data = var_part_time_data_dis, importance = TRUE)


## HMSC in time ##

hmsc_time_data_dis <- hms_ss %>%
  select(dispersal, hmsc_r2:SS_time) %>% 
  filter(!is.na(hmsc_time_ass_var)) %>% 
  filter(!is.na(hmsc_space_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_time_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_time_ass_mean_pos)) %>% 
  filter(!is.na(hmsc_space_time_ass_mean_pos)) %>% 
  filter(!is.na(hmsc_space_time_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_space_ass_mean_pos))

summary(hmsc_time_data_dis)


hmsc_time_rf_dis <- randomForest(dispersal~ ., data = hmsc_time_data_dis, importance = TRUE)



## everything model

everything_data_dis <- all_outputs %>% 
  select(dispersal, B_s_0:space_var_part_Environment) %>% 
  select(-contains('cor')) %>% 
  filter(!is.na(hmsc_time_ass_mean_pos)) %>% 
  filter(!is.na(hmsc_space_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_time_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_space_time_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_space_time_ass_mean_pos)) %>% 
  filter(!is.na(hmsc_space_ass_mean_pos)) %>% 
  filter(!is.na(space_var_part_Environment))

nrow(everything_data_dis)

everything_rf_dis <- randomForest(dispersal~ ., data = everything_data_dis, importance = TRUE)


## minimal model

# recursive feature elimitation
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results_dis <- rfe(everything_data_dis[,-1], everything_data_dis$dispersal, sizes=c(10:50), rfeControl=control)

results_dis$bestSubset

results_dis$fit$importance

saveRDS(results_dis, "results/minimal_dispersal.rds")


####### env_niche ######


### summary space only ##

all_summary_data_env <- all_outputs %>% 
  select(sig_niche, B_s_0:B_s_2, BDtotal_patch:RichDi_BDtotal_patch)

nrow(all_summary_data_env)

all_summary_rf_env <- randomForest(sig_niche~ ., data = all_summary_data_env, importance = TRUE)



### summary that we have ##

all_summary_data_env <- all_outputs %>% 
  select(sig_niche, B_s_0:gamma_ratios_ratio_2_0) %>% 
  select(-contains('cor'))

nrow(all_summary_data_env)

all_summary_rf_env <- randomForest(sig_niche~ ., data = all_summary_data_env, importance = TRUE)



## var part space ##

var_part_data_env <- space_var %>% 
  select(sig_niche, space_var_part_Space, space_var_part_both, space_var_part_Environment) %>% 
  filter(!is.na(space_var_part_Space))

nrow(var_part_data_env)
var_part_orig_env <- randomForest(sig_niche~ ., data = var_part_data_env, importance = TRUE)


## var part time ## 

var_part_time_data_env <- all_outputs %>% 
  select(sig_niche, starts_with("time_var_part")) %>% 
  filter(!is.na(time_var_part_environ_time))

nrow(var_part_time_data_env)

var_part_time_env <- randomForest(sig_niche~ ., data = var_part_time_data_env, importance = TRUE)


## HMSC in time ##

hmsc_time_data_env <- hms_ss %>% 
  select(sig_niche, hmsc_r2:SS_time) %>% 
  filter(!is.na(hmsc_time_ass_mean_pos)) %>% 
  filter(!is.na(hmsc_space_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_time_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_space_time_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_space_time_ass_mean_pos)) %>% 
  filter(!is.na(hmsc_space_ass_mean_pos)) 

hmsc_time_rf_env <- randomForest(sig_niche~ ., data = hmsc_time_data_env, importance = TRUE)


## everything model

everything_data_env <- all_outputs %>% 
  select(sig_niche, B_s_0:space_var_part_Environment) %>% 
  select(-contains('cor')) %>% 
  filter(!is.na(hmsc_time_ass_mean_pos)) %>% 
  filter(!is.na(hmsc_space_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_time_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_space_time_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_space_time_ass_mean_pos)) %>% 
  filter(!is.na(hmsc_space_ass_mean_pos)) %>% 
  filter(!is.na(space_var_part_Space))

everything_rf_env <- randomForest(sig_niche~ ., data = everything_data_env, importance = TRUE)

everything_rf_env$importance

# Env by alpha
everything_data_env <- all_outputs %>% 
  filter(alpha == 'patch') %>% 
  select(sig_niche, B_s_0:space_var_part_Environment) %>% 
  select(-contains('cor')) %>% 
  filter(!is.na(hmsc_time_ass_mean_pos)) %>% 
  filter(!is.na(hmsc_space_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_time_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_space_time_ass_mean_neg)) %>% 
  filter(!is.na(hmsc_space_time_ass_mean_pos)) %>% 
  filter(!is.na(hmsc_space_ass_mean_pos)) %>% 
  filter(!is.na(space_var_part_Space))

everything_rf_env <- randomForest(sig_niche~ ., data = everything_data_env, importance = TRUE)

#68, 51, 68, 50

ordered_importance <- everything_rf_env$importance[, 1, drop = FALSE] %>% rowSums() %>% sort(decreasing = TRUE) %>% names()

everything_rf_env$importance[, 1, drop = FALSE] %>% 
  as.data.frame() %>% 
  rownames_to_column("explanatory") %>%
  gather(key = 'interaction', value = 'importance', -explanatory) %>% 
  mutate(explanatory = factor(explanatory, levels = ordered_importance)) %>% 
  ggplot(aes(x = explanatory, y = interaction )) + geom_tile(aes(fill = importance), colour = "white") + 
  scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text = element_text(size = 10))

## minimal model

control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results_sig_niche <- rfe(everything_data_env[,-1], everything_data_env$sig_niche, sizes=c(10:50), rfeControl=control)

results_sig_niche$bestSubset

results_sig_niche$fit$importance

saveRDS(results_sig_niche, "results/minimal_sig_niche.rds")



################# Figures #############

### carpet all Figure 2 #######

results_alpha <- readRDS("results/minimal_alpha.rds")
results_dis <- readRDS("results/minimal_dispersal.rds")
results_sig_niche <- readRDS("results/minimal_sig_niche.rds")


ordered_importance_alpha <- results_alpha$fit$importance[, 1:4] %>% rowSums() %>% sort(decreasing = TRUE) %>% names()

ordered_pretty <- data.frame(ordered_importance_explanatory = ordered_importance_alpha) %>% 
  left_join(name_df, by = c("ordered_importance_explanatory" = "original_names")) %>% 
  mutate(pretty_names = as.character(pretty_names))

alpha_minimal_importance <- results_alpha$fit$importance[, 1:4] %>% 
  as.data.frame() %>% 
  rownames_to_column("explanatory") %>%
  gather(key = 'interaction', value = 'importance', -explanatory) %>%
  left_join(name_df, by = c("explanatory" = "original_names")) %>% 
  mutate(pretty_names = factor(pretty_names, levels = rev(ordered_pretty$pretty_names)))  %>% 
  ggplot(aes(x = pretty_names, y = interaction )) + geom_tile(aes(fill = importance), colour = "white") + 
  scale_fill_viridis() +  theme_cowplot() +
  theme(axis.text = element_text(size = 10),
        legend.position = 'none') +
  ylab("") + xlab("") + scale_y_discrete(labels = c("Stabilizing\ncompetition", "Mixed\ncompetition", "Equal\ncompetition", 
                                                    "CC-\ntrade-off")) +
  coord_flip()



ordered_importance_dispersal <- results_dis$fit$importance[, 1, drop = FALSE] %>% rowSums() %>% sort(decreasing = TRUE) %>% names()

ordered_pretty <- data.frame(ordered_importance_explanatory = ordered_importance_dispersal) %>% 
  left_join(name_df, by = c("ordered_importance_explanatory" = "original_names")) %>% 
  mutate(pretty_names = as.character(pretty_names))

dispersal_minimal_importance <- results_dis$fit$importance[, 1, drop = FALSE] %>% 
  as.data.frame() %>% 
  rownames_to_column("explanatory") %>%
  gather(key = 'interaction', value = 'importance', -explanatory) %>%
  left_join(name_df, by = c("explanatory" = "original_names")) %>% 
  mutate(pretty_names = factor(pretty_names, levels = rev(ordered_pretty$pretty_names)))  %>% 
  ggplot(aes(x = pretty_names, y = interaction )) + geom_tile(aes(fill = importance), colour = "white") + 
  scale_fill_viridis() + theme_cowplot() +
  theme(axis.text = element_text(size = 10),
        legend.position = 'none') +
  ylab("") + xlab("") + scale_y_discrete(labels = c("Dispersal")) +
  coord_flip()


ordered_importance_sig_niche <- results_sig_niche$fit$importance[, 1, drop = FALSE] %>% rowSums() %>% sort(decreasing = TRUE) %>% names()

ordered_pretty <- data.frame(ordered_importance_explanatory = ordered_importance_sig_niche) %>% 
  left_join(name_df, by = c("ordered_importance_explanatory" = "original_names")) %>% 
  mutate(pretty_names = as.character(pretty_names))


niche_minimal_importance <- results_sig_niche$fit$importance[, 1, drop = FALSE] %>% 
  as.data.frame() %>% 
  rownames_to_column("explanatory") %>%
  gather(key = 'interaction', value = 'importance', -explanatory)%>%
  left_join(name_df, by = c("explanatory" = "original_names")) %>% 
  mutate(pretty_names = factor(pretty_names, levels = rev(ordered_pretty$pretty_names)))  %>% 
  ggplot(aes(x = pretty_names, y = interaction )) + geom_tile(aes(fill = importance), colour = "white") + 
  scale_fill_viridis() +theme_cowplot() +
  theme(axis.text = element_text(size = 10),
        legend.position = 'none') + 
  ylab("") + xlab("") + scale_y_discrete(labels = c("Density-independent responses to abiotic conditions")) + coord_flip()


ggsave(alpha_minimal_importance, filename = "figures/alpha_carpet.jpeg", width = 8, height = 14)
ggsave(dispersal_minimal_importance, filename = "figures/dispersal_carpet.jpeg", height = 10)
ggsave(niche_minimal_importance, filename = "figures/niche_carpet.jpeg", height = 10)




niche_importance <- results_sig_niche$fit$importance[1:20, 1, drop = FALSE] %>% 
  as.data.frame() %>% 
  rownames_to_column("explanatory") %>%
  rename(Niche = `%IncMSE`) %>% 
  gather(key = 'interaction', value = 'importance', -explanatory) %>% 
  mutate(importance = importance/max(importance))

dispersal_importance <- results_dis$fit$importance[1:20, 1, drop = FALSE] %>% 
  as.data.frame() %>% 
  rownames_to_column("explanatory")%>%
  rename(Dispersal = `%IncMSE`) %>% 
  gather(key = 'interaction', value = 'importance', -explanatory)%>% 
  mutate(importance = importance/max(importance))

alpha_importance <- results_alpha$fit$importance[1:20, 1:4] %>% 
  as.data.frame() %>% 
  rownames_to_column("explanatory") %>%
  gather(key = 'interaction', value = 'importance', -explanatory) %>% 
  mutate(importance = importance/max(importance))

all_importance <- rbind(alpha_importance, dispersal_importance, niche_importance)

ordered_importance_all <- all_importance %>% 
  arrange(desc(importance)) %>% 
  select(explanatory) %>% unique() %>% 
  left_join(name_df, by = c("explanatory" = "original_names"))


all_importance_df <- all_importance %>% 
  left_join(name_df, by = c("explanatory" = "original_names")) %>% 
  mutate(pretty_names = factor(pretty_names, levels = rev(ordered_importance_all$pretty_names))) %>% 
  mutate(interaction = factor(interaction, levels = c("0.5", "1.5", "equal", "patch", "Niche","Dispersal"))) %>% 
  mutate(big_label = case_when(interaction %in% c("0.5", "1.5", "equal", "patch") ~ "Density-dependent",
                               interaction == 'Niche' ~ "Density-independent", 
                               interaction == "Dispersal" ~ "Dispersal"))
  
  
carpet_all <-  all_importance_df %>% 
  ggplot(aes(x = pretty_names, y = interaction )) + geom_tile(aes(fill = importance), colour = "white") + 
  scale_fill_viridis() + theme_cowplot() +
  theme(axis.text = element_text(size = 16),
        legend.position = 'none') +
  ylab("") + xlab("")  + coord_flip() +
  scale_y_discrete(position = 'right', 
    labels =  c("Stabilizing\ncompetition", "                 Density \n \n Mixed \ncompetition", "-Dependent \n \n Equal \ncompetition", 
                                                     "CC-\ntrade-off", "Density-Independent \n  \n  \n ", "  Dispersal \n \n \n ")) 


ggsave(carpet_all, filename = "figures/carpet_all.jpeg", width = 15, height = 18)


######### number of predictors performance ########

trellis.par.set(caretTheme())

alpha_accurancey_cv <- plot(results_alpha, type = c("g", "o"))

alpha_accurancey_cv$xlab <- "Number of predictors"

disper_R2_cv <- plot(results_dis, type = c("g", "o"))

disper_R2_cv$xlab <- "Number of predictors"

sig_niche_R2_cv <- plot(results_sig_niche, type = c("g", "o"))

sig_niche_R2_cv$xlab <- "Number of predictors"

cross_validation_n_pred <- plot_grid(alpha_accurancey_cv, disper_R2_cv, sig_niche_R2_cv, nrow = 1, labels = c("A.", "B.", "C."))

ggsave(cross_validation_n_pred, filename = "figures/cross_validation_n_pred.jpeg", width = 10)



################ best preditors with median and iqr #############

results_alpha <- readRDS("results/minimal_alpha.rds")
results_dis <- readRDS("results/minimal_dispersal.rds")
results_sig_niche <- readRDS("results/minimal_sig_niche.rds")

alpha_minimal_importance <- results_alpha$fit$importance[, 1:4, drop = FALSE] %>% rowSums() %>% sort(decreasing = TRUE) %>% names()
alpha_minimal_importance[1]

dispersal_minimal_importance <- results_dis$fit$importance[, 1, drop = FALSE] %>% rowSums() %>% sort(decreasing = TRUE) %>% names()
dispersal_minimal_importance[1]

sig_niche_minimal_importance <- results_sig_niche$fit$importance[, 1, drop = FALSE] %>% rowSums() %>% sort(decreasing = TRUE) %>% names()
sig_niche_minimal_importance[1]


important_names <- name_df %>% 
  filter(original_names %in% c(alpha_minimal_importance[1], dispersal_minimal_importance[1], sig_niche_minimal_importance[1]))


star_alpha <- data.frame(alpha = "patch", mn = 2.1, pretty_names = "Coefficient of variation \n in abundance at local scale")
star_alpha_2 <- data.frame(alpha = "patch", mn = 1.1, pretty_names = "Ratio of occupancies at local \nscale (Hill 1/ Hill 0)")
star_alpha_3 <- data.frame(alpha = "patch", mn = 1.8, pretty_names = 'Variation partitioning in time and \nspace (Environment component)')


alpha_important <- all_outputs %>% 
  select(dispersal, sig_niche, alpha, alpha_bio_cv, occupancy_ratios_ratio_1_0, time_var_part_Environment) %>% 
  gather(key = "original_names", value = value, alpha_bio_cv:time_var_part_Environment) %>% 
  left_join(important_names) %>% 
  ggplot(aes(x = alpha, y = value)) + facet_wrap(~pretty_names, nrow = 1, scales = 'free',strip.position = "left")+geom_violin(fill = 'grey') +
  geom_boxplot(width=.1, outlier.colour=NA) +
  scale_x_discrete(labels = c("Stabilizing\ncompetition", "Mixed\ncompetition", "Equal\ncompetition", 
                              "CC-\ntrade-off")) + xlab("") + theme_cowplot() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 24), 
        strip.text = element_text(size = 19),
        strip.placement = "outside") + scale_colour_viridis() + ylab(NULL) + scale_y_log10() +
  geom_point(data = star_alpha, aes(x = alpha, y = mn), size = 7, shape = 23, fill = "#FDE725FF") + 
  geom_point(data = star_alpha_2, aes(x = alpha, y = mn), size = 7, shape = 23, fill = "#481567FF")+ 
  geom_point(data = star_alpha_3, aes(x = alpha, y = mn), size = 7, shape = 23, fill = "white")


star_dispersal <- data.frame(dispersal = 0.15, mn = 1.1, pretty_names = "Ratio of occupancies at local \nscale (Hill 1/ Hill 0)")
star_dispersal_2 <- data.frame(dispersal = 0.15, mn = 0.29, pretty_names = "Variation partitioning in time and \nspace (Environment component)")
star_dispersal_3 <- data.frame(dispersal = 0.15, mn = 0.45, pretty_names = "Coefficient of variation \n in abundance at local scale")


dispersal_important <- all_outputs %>% 
  select(dispersal, sig_niche, alpha, alpha_bio_cv, occupancy_ratios_ratio_1_0, time_var_part_Environment) %>% 
  gather(key = "original_names", value = value, alpha_bio_cv:time_var_part_Environment) %>% 
  left_join(important_names) %>% 
  group_by(dispersal, pretty_names) %>% 
  summarise(mn = median(value, na.rm = TRUE), se = sd(value, na.rm = TRUE)/sqrt(n()), q1 = quantile(value,1/4, na.rm = TRUE),
            q3 = quantile(value,3/4, na.rm = TRUE)) %>% 
  ggplot() + facet_wrap(~pretty_names, nrow = 1, scales = 'free',strip.position = "left")+geom_line(aes(x = dispersal, y = mn))  +
  scale_x_log10() + geom_ribbon(aes(x = dispersal, ymin = q1, ymax = q3), alpha = 0.2) + ylab("") + xlab("Dispersal") +
  theme_cowplot() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 24), 
        strip.text = element_text(size = 19),
        strip.placement = "outside") + ylab(NULL) + 
  geom_point(data = star_dispersal, aes(x = dispersal, y = mn), size = 7, shape = 23, fill = "#FDE725FF") + 
  geom_point(data = star_dispersal_2, aes(x = dispersal, y = mn), size = 7, shape = 23, fill = "#481567FF")+ 
  geom_point(data = star_dispersal_3, aes(x = dispersal, y = mn), size = 7, shape = 23, fill = "#31688EFF")


star_niche <- data.frame(sig_niche = 4.5, mn = 1.38, pretty_names = "Coefficient of variation \n in abundance at local scale")
star_niche_2 <- data.frame(sig_niche = 4.5, mn = 0.96, pretty_names = "Ratio of occupancies at local \nscale (Hill 1/ Hill 0)")
star_niche_3 <- data.frame(sig_niche = 4.5, mn = 0.30, pretty_names = "Variation partitioning in time and \nspace (Environment component)")

niche_important <- all_outputs %>% 
  select(dispersal, sig_niche, alpha, alpha_bio_cv, occupancy_ratios_ratio_1_0, time_var_part_Environment) %>% 
  gather(key = "original_names", value = value, alpha_bio_cv:time_var_part_Environment) %>% 
  left_join(important_names) %>% 
  group_by(sig_niche, pretty_names) %>% 
  summarise(mn = median(value, na.rm = TRUE), se = sd(value, na.rm = TRUE)/sqrt(n()), q1 = quantile(value,1/4, na.rm = TRUE),
            q3 = quantile(value,3/4, na.rm = TRUE)) %>% 
  ggplot() + facet_wrap(~pretty_names, nrow = 1, scales = 'free',strip.position = "left")+geom_line(aes(x = sig_niche, y = mn))  +
  scale_x_log10() + geom_ribbon(aes(x = sig_niche, ymin = q1, ymax = q3), alpha = 0.2)+ xlab("Abiotic responses") +
  theme_cowplot() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 24), 
        strip.text = element_text(size = 19),
        strip.placement = "outside") + ylab(NULL) + 
  geom_point(data = star_niche, aes(x = sig_niche, y = mn), size = 7, shape = 23, fill = "#31688EFF") +
  geom_point(data = star_niche_2, aes(x = sig_niche, y = mn), size = 7, shape = 23, fill = "white") +
  geom_point(data = star_niche_3, aes(x = sig_niche, y = mn), size = 7, shape = 23, fill = "#FDE725FF")

important_pred <- plot_grid(alpha_important,dispersal_important, niche_important, nrow = 3)

important_pred2 <- ggdraw(important_pred) + draw_plot_label(label=c("A","B", "C", "D", "E", "F", "G", "H", "I"), 
                                         x=c(0,0.33, 0.645, 0,0.34, 0.655, 0,0.34, 0.655), y=c(1,1,1, 0.67, 0.67, 0.67, 0.35, 0.35, 0.35))


ggsave(important_pred2, filename = "figures/important_predictors.jpeg", width = 18, height = 15)




