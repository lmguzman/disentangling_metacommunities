library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(cowplot)

sensitivy <- read.csv("results/sensitivy_small.csv")

mid_d <- 0.00215
mid_s <- 4.64159

true_value <- read.csv("results/desc_stat_time_vp.csv") %>% 
  filter(rep ==1, round(dispersal, 5) == mid_d, round(sig_niche,5) == mid_s, alpha == "0.5") 



name_df <- data.frame(original_names2 = original_names[c(6:13, 17:29, 35:47)],
                      pretty_names = c("Beta diversity in space (Hill 0)", "Beta diversity in space (Hill 1)", "Beta diversity in space (Hill 2)",
                                       "Beta diversity in time (Hill 0)", "Beta diversity in time (Hill 1)","Beta diversity in time (Hill 2)",
                                       "Coefficient of variation \n in abundance at local scale", "Coefficient of variation \n in abundance at regional scale",
                                       "Mean proportion of \npatches occupied", "Minimum proportion of \npatches occupied", "Maximum proportion of patches occupied",
                                       "Total beta diversity in space", "Replacement component of\n beta diversity in space", "Total richness difference \ndiversity in space",
                                       "Total replacement diversity/ \n Total beta diversity in space", "Total richness difference diversity/\n Total beta diversity in space","Total beta diversity in time", "Replacement component of \nbeta diversity in time", "Total richness difference \ndiversity in time",
                                       "Total replacement diversity/\n Total beta diversity in time", "Total richness difference diversity/\n Total beta diversity in time",
                                       "Ratio of abundances at local \nscale (Hill 1/Hill 0)",
                                       "Ratio of abundances at local \nscale (Hill 2/Hill 0)", "Ratio of occupancies at local \nscale (Hill 1/ Hill 0)", 
                                       "Ratio of occupancies at local \nscale (Hill 2/ Hill 0)","Ratio of abundances at regional \nscale (Hill 1/Hill 0)",
                                       "Ratio of abundances at regional\n scale (Hill 2/Hill 0)", "Variation partitioning in time and \nspace (Space component)",
                                       "Variation partitioning in time and \nspace (Environment component)", "Variation partitioning in time and \nspace (Time component)",
                                       "Variation partitioning in time and \nspace (Space and environment component)", "Variation partitioning in time and \nspace (Environment and time component)",
                                       "Variation partitioning in time and \nspace (Space and time component)", "Variation partitioning in time and \nspace (Space, time and environment component)")) %>% 
  mutate(original_names2 = as.character(original_names2), pretty_names = as.character(pretty_names))




##### final sensitivy plots


sens_1 <- sensitivy %>% 
  filter(patch %in% c(4, 36, 68, 100)) %>% 
  mutate(diff_alpha_bio_cv = true_value$alpha_bio_cv- alpha_bio_cv) %>% 
  group_by(time, patch) %>% 
  summarise(mn_diff = mean(diff_alpha_bio_cv), sd_diff = sd(diff_alpha_bio_cv)) %>% 
  ggplot(aes(x = time, y = mn_diff, colour = factor(patch))) + geom_line() + 
  geom_ribbon(aes(x = time, ymin = mn_diff-sd_diff, ymax = mn_diff+sd_diff, fill = factor(patch)), alpha = 0.4) +
  ylab("Difference of true value and sampled value") +
  geom_hline(aes(yintercept = 0), linetype = 'dashed') + xlab("")+ theme_cowplot()+
  theme(legend.position = 'none',axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))


sens_2 <- sensitivy %>% 
  filter(patch %in% c(4, 36, 68, 100)) %>% 
  mutate(diff_RichDif_time = true_value$RichDif_time-RichDif_time) %>% 
  group_by(time, patch) %>% 
  summarise(mn_diff = mean(diff_RichDif_time), sd_diff = sd(diff_RichDif_time)) %>% 
  ggplot(aes(x = time, y = mn_diff, colour = factor(patch))) + geom_line() + 
  geom_ribbon(aes(x = time, ymin = mn_diff-sd_diff, ymax = mn_diff+sd_diff, fill = factor(patch)), alpha = 0.4) +
  ylab("") + xlab("Number of time points") +
  geom_hline(aes(yintercept = 0), linetype = 'dashed') + theme_cowplot()+
  theme(legend.position = 'none',axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))


sens_3 <- sensitivy %>%
  filter(patch %in% c(4, 36, 68, 100)) %>% 
  mutate(diff_mean_prop_patches = true_value$mean_prop_patches-mean_prop_patches) %>% 
  group_by(time, patch) %>% 
  summarise(mn_diff = mean(diff_mean_prop_patches), sd_diff = sd(mean_prop_patches)) %>% 
  ggplot(aes(x = time, y = mn_diff, colour = factor(patch))) + geom_line() + 
  geom_ribbon(aes(x = time, ymin = mn_diff-sd_diff, ymax = mn_diff+sd_diff, fill = factor(patch)), alpha = 0.4) +
  ylab("") +
  geom_hline(aes(yintercept = 0), linetype = 'dashed') + xlab("") + labs(fill = "Number of\n patches", 
                                                                         colour = "Number of\n patches") +
  theme_cowplot() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))

sensitivity_figure <- plot_grid(sens_1, sens_2, sens_3, nrow = 1, labels = c("A", "B","C"), rel_widths = c(0.78,0.76,1))

ggsave(sensitivity_figure, file = 'figures/sensitivity.jpeg', width = 15, height = 8)



###### half life###

### sensitivity is >300MB could not be included on github, please email if needed.

sensitivy_t <- read.csv("results/sensitivity.csv") 

original_names <- c("X", "r","time", "patch", "J", "B_s_0", "B_s_1", "B_s_2", "B_t_0", "B_t_1", "B_t_2",
                    "alpha_bio_cv", "gamma_bio_cv", "max_abun", "min_abun","ratio_max_min","mean_prop_patches",
                    "min_prop_patches","max_prop_patches", "BDtotal_patch","Repl_patch","RichDif_patch","Repl_BDtotal_patch",
                    "RichDi_BDtotal_patch","BDtotal_time", "Repl_time","RichDif_time","Repl_BDtotal_time", "RichDi_BDtotal_time",
                    "cor_comp_env", "cor_comp_space","cor_comp_time","cor_env_space","cor_env_time", "abundance_ratios_ratio_1_0",
                    "abundance_ratios_ratio_2_0","occupancy_ratios_ratio_1_0", "occupancy_ratios_ratio_2_0","gamma_ratios_ratio_1_0", 
                    "gamma_ratios_ratio_2_0","time_var_part_Space","time_var_part_Environment", "time_var_part_Time",
                    "time_var_part_space_environ", "time_var_part_environ_time", "time_var_part_time_space","time_var_part_all",
                    "time_var_part_residual")

colnames(sensitivy_t) <- original_names



#### Patch #####

mean_full_time <- sensitivy_t %>% 
  select(-X) %>% 
  filter(time == 60) %>% 
  group_by(time, patch) %>% 
  summarise_all(mean) %>% 
  ungroup()

name_df$original_names2

patch_df <- select(mean_full_time, name_df$original_names2) 

true_df <- select(true_value, name_df$original_names2) %>% 
  slice(rep(1:n(), each = 25))

patch_diff <- patch_df-true_df

patch_diff2 <- cbind(patch = mean_full_time$patch, patch_diff)

patch_decay_plot <- patch_diff2 %>%
  gather(key = "original_names2", value = 'value', -patch) %>% 
  left_join(name_df) %>% 
  ggplot(aes(x = patch, y = value)) + geom_line() + facet_wrap(~pretty_names, scales = 'free') + 
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = 'red') +
  theme_cowplot() + xlab("Number of patches") +
  theme(strip.background = element_blank()) + ylab("Difference of true value and sampled value")

patch_diff2

ggsave(patch_decay_plot, file = "figures/patch_sensitivity_decay.jpeg", height = 15, width = 25, dpi = 200)



patch_half_life <- c()
for(i in 2:ncol(patch_diff2)){
  index <- which(abs((patch_diff2[,i] - patch_diff2[1,i]/2)) == min(abs(patch_diff2[,i] - patch_diff2[1,i]/2)))
  patch_half_life <-c(patch_half_life, patch_diff2$patch[index])
}


names(patch_half_life) <- name_df$pretty_names

sort(patch_half_life) %>% View()





#### Time #####

mean_full_patch <- sensitivy_t %>% 
  select(-X) %>% 
  filter(patch == 100) %>% 
  group_by(time, patch) %>% 
  summarise_all(mean) %>% 
  ungroup()

name_df$original_names2

time_df <- select(mean_full_patch, name_df$original_names2) 

true_df <- select(true_value, name_df$original_names2) %>% 
  slice(rep(1:n(), each = 15))

time_diff <- time_df-true_df

time_diff2 <- cbind(time = mean_full_patch$time, time_diff)

time_decay_plot <- time_diff2 %>%
  gather(key = "original_names2", value = 'value', -time) %>% 
  left_join(name_df) %>% 
  ggplot(aes(x = time, y = value)) + geom_line() + facet_wrap(~pretty_names, scales = 'free') + 
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = 'red') +
  theme_cowplot() + xlab("Number of time points") +
  theme(strip.background = element_blank()) + ylab("Difference of true value and sampled value")


ggsave(time_decay_plot , file = "figures/time_sensitivity_decay.jpeg", height = 15, width = 25, dpi = 200)



time_half_life <- c()
for(i in 2:ncol(time_diff2)){
  index <- which(abs(time_diff2[,i] - time_diff2[1,i]/2) == min(abs(time_diff2[,i] - time_diff2[1,i]/2)))
  time_half_life <-c(time_half_life, time_diff2$time[min(index)])
}


names(time_half_life) <- name_df$pretty_names

sort(time_half_life) %>% View()

