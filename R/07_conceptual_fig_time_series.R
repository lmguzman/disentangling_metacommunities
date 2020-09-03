library(dplyr)
library(ggplot2)
library(cowplot)

r <- 1
load("R/pars.Rdata")
load(paste0("simulation_outputs/outputfile_",r,".RData"))

i <- 250
 
p <- pars[i,]

run_one <- model.df[which(round(model.df$dispersal, 5) == round(p$dispersal, 5) & round(model.df$sig_niche, 5) == round(p$sig_niche, 5) &  model.df$alpha == p$alpha), ]

run_one %>%
  filter(Patch < 5) %>% 
  group_by(Species, Patch) %>% 
  summarise(n()) %>% as.data.frame()

run_one %>%
  filter(Patch < 5) %>% 
  select(Patch, Time, Species, N) %>% 
  spread(key = Species, value = N, fill = 0) %>% 
  gather(key = Species, value = N, -Patch, -Time) %>% 
  mutate(p = paste("Patch", Patch)) %>% 
  ggplot(aes(x = Time, y = N, colour = as.character(Species))) + geom_line() +
  scale_colour_viridis_d() + facet_wrap(~p) +
  theme(legend.position = 'none',
        strip.background = element_rect(colour="white", fill="white")) +
  ylab("Abundance")
