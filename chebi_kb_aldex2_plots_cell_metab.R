
# Packages ----------------------------------------------------------------


library(rstatix)
library(ggsci)
library(ggthemes)
library(hrbrthemes)
library(gcookbook)
library(ggpubr)
library(ggprism)
library(forcats)
library(RColorBrewer)
library(tidyverse)
library(scales)
show_col(pal_futurama("planetexpress")(12))
futurama <- pal_futurama("planetexpress")(12)
source('./Scripts/theme_clean_rgd.R')
windowsFonts("Times" = windowsFont("Times"))


# Species_A_C -------------------------------------------------------------
group_label <- 1.5
limits <- c(-1.5,1.5)
breaks <- c(-1.5,-1,-0.5,0,0.5,1,1.5)
y.position <- -1
## aldex2_MODSUG_A_C ####

up <- read_csv(file = './Results/chebi_aldex2_MODSUG_A_C.csv') %>% 
  mutate(group_trials = as.factor('MODSUG_A_C')) %>% 
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/chebi_aldex2_MODSUG_A_C.csv') %>% 
  mutate(group_trials = as.factor('MODSUG_A_C')) %>% 
  slice_min(effect, n = 10)

middle <- down %>% mutate(species = str_replace_all(species, '--.*', '')) %>% slice_max(effect) %>% pull(species) 

MODSUG_A_C <- bind_rows(up, down) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) 

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
MODSUG_A_C_plot <- MODSUG_A_C %>% ggplot()+
  geom_col(aes(y = effect, x = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = group_label, x = middle, label = 'MODSUG', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10, face = 'italic'),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = NULL, x = NULL)+
  scale_y_continuous(breaks = breaks, limits = limits)+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(MODSUG_A_C, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'species', xmax = NULL, y.position = y.position,
             inherit.aes = F, family = 'Times')


## aldex2_LOWCHO_A_C ####
up <- read_csv(file = './Results/chebi_aldex2_LOWCHO_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWCHO_A_C')) %>% 
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/chebi_aldex2_LOWCHO_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWCHO_A_C')) %>% 
  slice_min(effect, n = 10)

middle <- down %>% mutate(species = str_replace_all(species, '--.*', '')) %>% slice_max(effect) %>% pull(species) 

LOWCHO_A_C <- bind_rows(up, down) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.105,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) 

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
LOWCHO_A_C_plot <- LOWCHO_A_C %>% ggplot()+
  geom_col(aes(y = effect, x = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = group_label, x = middle, label = 'LOWCHO', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10, face = 'italic'),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = 'Estimated ALDEx2 Effect Size', x = NULL)+
  scale_y_continuous(breaks = breaks, limits = limits)+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(LOWCHO_A_C, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'species', xmax = NULL, y.position = y.position,
             inherit.aes = F, family = 'Times')



## aldex2_LOWSUG_A_C ####
up <- read_csv(file = './Results/chebi_aldex2_LOWSUG_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWSUG_A_C')) %>% 
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/chebi_aldex2_LOWSUG_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWSUG_A_C')) %>% 
  slice_min(effect, n = 10)

middle <- down %>% mutate(species = str_replace_all(species, '--.*', '')) %>% slice_max(effect) %>% pull(species) 

LOWSUG_A_C <- bind_rows(up, down) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) 

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
LOWSUG_A_C_plot <- LOWSUG_A_C %>% ggplot()+
  geom_col(aes(y = effect, x = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = group_label, x = middle, label = 'LOWSUG', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10, face = 'italic'),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = NULL, x = 'Species')+
  scale_y_continuous(breaks = breaks, limits = limits)+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(LOWSUG_A_C, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'species', xmax = NULL, y.position = y.position,
             inherit.aes = F, family = 'Times')



## Species_A_C Arrange -----------------------------------------------------------

ggarrange(MODSUG_A_C_plot, LOWSUG_A_C_plot, LOWCHO_A_C_plot, ncol = 1, nrow = 3, align = 'v')

ggsave(filename = './Figures/chebi_aldex2_effect_Species_all_groups_A_C_2.tiff',
       device = 'tiff',
       width = 15,
       height = 24,
       units = 'cm',
       dpi = 600)

# Genus A B -------------------------------------------------------------------

## aldex2_MODSUG_A_B ####

up <- read_csv(file = './Results/chebi_aldex2_G_MODSUG_A_B.csv') %>% 
  mutate(group_trials = as.factor('MODSUG_A_B')) %>% 
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/chebi_aldex2_G_MODSUG_A_B.csv') %>% 
  mutate(group_trials = as.factor('MODSUG_A_B')) %>% 
  slice_min(effect, n = 10)

middle <- down %>% mutate(species = str_replace_all(species, '--.*', '')) %>% slice_max(effect) %>% pull(species) 

MODSUG_A_B <- bind_rows(up, down) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) 

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
MODSUG_A_B_plot <- MODSUG_A_B %>% ggplot()+
  geom_col(aes(y = effect, x = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = 1.5, x = middle, label = 'MODSUG', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10, face = 'italic'),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = 'Estimated ALDEx2 Effect Size', x = 'Genus')+
  scale_y_continuous(breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), limits = c(-1.5,1.5))+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(MODSUG_A_B, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'species', xmax = NULL, y.position = -1.5,
             inherit.aes = F, family = 'Times')


## aldex2_LOWCHO_A_B ####
up <- read_csv(file = './Results/chebi_aldex2_G_LOWCHO_A_B.csv') %>% 
  mutate(group_trials = as.factor('LOWCHO_A_B')) %>% 
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/chebi_aldex2_G_LOWCHO_A_B.csv') %>% 
  mutate(group_trials = as.factor('LOWCHO_A_B')) %>% 
  slice_min(effect, n = 10)

middle <- down %>% mutate(species = str_replace_all(species, '--.*', '')) %>% slice_max(effect) %>% pull(species) 

LOWCHO_A_B <- bind_rows(up, down) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) 

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
LOWCHO_A_B_plot <- LOWCHO_A_B %>% ggplot()+
  geom_col(aes(y = effect, x = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = 1.5, x = middle, label = 'LOWCHO', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10, face = 'italic'),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = 'Estimated ALDEx2 Effect Size', x = 'Genus')+
  scale_y_continuous(breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), limits = c(-1.5,1.5))+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(LOWCHO_A_B, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'species', xmax = NULL, y.position = -1.5,
             inherit.aes = F, family = 'Times')
  

    
## aldex2_LOWSUG_A_B ####
up <- read_csv(file = './Results/chebi_aldex2_G_LOWSUG_A_B.csv') %>% 
  mutate(group_trials = as.factor('LOWSUG_A_B')) %>% 
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/chebi_aldex2_G_LOWSUG_A_B.csv') %>% 
  mutate(group_trials = as.factor('LOWSUG_A_B')) %>% 
  slice_min(effect, n = 10)

middle <- down %>% mutate(species = str_replace_all(species, '--.*', '')) %>% slice_max(effect) %>% pull(species) 

LOWSUG_A_B <- bind_rows(up, down) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) 

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
LOWSUG_A_B_plot <- LOWSUG_A_B %>% ggplot()+
  geom_col(aes(y = effect, x = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = 1.5, x = middle, label = 'LOWSUG', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10, face = 'italic'),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = 'Estimated ALDEx2 Effect Size', x = 'Genus')+
  scale_y_continuous(breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), limits = c(-1.5,1.5))+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(LOWSUG_A_B, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'species', xmax = NULL, y.position = -1.5,
             inherit.aes = F, family = 'Times')



## Genus A B Arrange -----------------------------------------------------------

ggarrange(MODSUG_A_B_plot, LOWSUG_A_B_plot, LOWCHO_A_B_plot, ncol = 1, nrow = 3, align = 'v')

ggsave(filename = './Figures/chebi_aldex2_effect_Genus_all_groups_A_B.tiff',
       device = 'tiff',
       width = 15,
       height = 24,
       units = 'cm',
       dpi = 600)


# Genus A_C -------------------------------------------------------------------

limits <- c(-2,2)
breaks <- c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2)
y.position <- -2
group_label <- 2

## aldex2_MODSUG_A_C ####

up <- read_csv(file = './Results/chebi_aldex2_G_MODSUG_A_C.csv') %>% 
  mutate(group_trials = as.factor('MODSUG_A_C')) %>% 
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/chebi_aldex2_G_MODSUG_A_C.csv') %>% 
  mutate(group_trials = as.factor('MODSUG_A_C')) %>% 
  slice_min(effect, n = 10)

middle <- down %>% mutate(species = str_replace_all(species, '--.*', '')) %>% slice_max(effect) %>% pull(species) 

MODSUG_A_C <- bind_rows(up, down) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) 

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
MODSUG_A_C_plot <- MODSUG_A_C %>% ggplot()+
  geom_col(aes(y = effect, x = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = group_label, x = middle, label = 'MODSUG', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10, face = 'italic'),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = 'Estimated ALDEx2 Effect Size', x = 'Genus')+
  scale_y_continuous(breaks = breaks, limits = limits)+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(MODSUG_A_C, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'species', xmax = NULL, y.position = y.position,
             inherit.aes = F, family = 'Times')


## aldex2_LOWCHO_A_C ####
up <- read_csv(file = './Results/chebi_aldex2_G_LOWCHO_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWCHO_A_C')) %>% 
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/chebi_aldex2_G_LOWCHO_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWCHO_A_C')) %>% 
  slice_min(effect, n = 10)

middle <- down %>% mutate(species = str_replace_all(species, '--.*', '')) %>% slice_max(effect) %>% pull(species) 

LOWCHO_A_C <- bind_rows(up, down) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) 

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
LOWCHO_A_C_plot <- LOWCHO_A_C %>% ggplot()+
  geom_col(aes(y = effect, x = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = group_label, x = middle, label = 'LOWCHO', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10, face = 'italic'),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = 'Estimated ALDEx2 Effect Size', x = 'Genus')+
  scale_y_continuous(breaks = breaks, limits = limits)+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(LOWCHO_A_C, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'species', xmax = NULL, y.position = y.position,
             inherit.aes = F, family = 'Times')



## aldex2_LOWSUG_A_C ####
up <- read_csv(file = './Results/chebi_aldex2_G_LOWSUG_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWSUG_A_C')) %>% 
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/chebi_aldex2_G_LOWSUG_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWSUG_A_C')) %>% 
  slice_min(effect, n = 10)

middle <- down %>% mutate(species = str_replace_all(species, '--.*', '')) %>% slice_max(effect) %>% pull(species) 

LOWSUG_A_C <- bind_rows(up, down) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) 

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
LOWSUG_A_C_plot <- LOWSUG_A_C %>% ggplot()+
  geom_col(aes(y = effect, x = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = group_label, x = middle, label = 'LOWSUG', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10, face = 'italic'),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = 'Estimated ALDEx2 Effect Size', x = 'Genus')+
  scale_y_continuous(breaks = breaks, limits = limits)+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(LOWSUG_A_C, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'species', xmax = NULL, y.position = y.position,
             inherit.aes = F, family = 'Times')



## Genus A_C Arrange -----------------------------------------------------------

ggarrange(MODSUG_A_C_plot, LOWSUG_A_C_plot, LOWCHO_A_C_plot, ncol = 1, nrow = 3, align = 'v')

ggsave(filename = './Figures/chebi_aldex2_effect_Genus_all_groups_A_C.tiff',
       device = 'tiff',
       width = 15,
       height = 24,
       units = 'cm',
       dpi = 600)

# LOWCHO_A_C Only

## Genus_LOWCHO_A_C_only ####

LOWCHO_A_C <- read_csv(file = './Results/chebi_aldex2_G_LOWCHO_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWCHO_A_C')) %>% 
  filter(wi.eBH < 0.1) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) 

colours <- c("#101E9D","#46A6CE","#EEA327")

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
# LOWCHO_A_C_plot <- 
  LOWCHO_A_C %>% ggplot()+
  geom_col(aes(y = effect, x = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#EEA327", high = '#101E9D', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  theme(text = element_text(family = 'Helvetica',
                            colour = 'black',
                            size = 10, 
                            face = 'bold'),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10, face = 'italic'),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = 'Estimated ALDEx2 Effect Size', x = 'Genus')+
  scale_y_continuous(breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), limits = c(-1.75,1.7))+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(LOWCHO_A_C, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'species', xmax = NULL, y.position = -1.75,
             inherit.aes = F, family = 'Helvetica')

ggsave(filename = './Figures/chebi_aldex2_effect_Genus_LOWCHO_A_C_allsig_cell_metab.tiff',
       device = 'tiff',
       width = 15,
       height = 24,
       units = 'cm',
       dpi = 600)

# Phylum A B -------------------------------------------------------------------
## aldex2_MODSUG_A_B ####

MODSUG_A_B <- read_csv(file = './Results/chebi_aldex2_P_MODSUG_A_B.csv') %>% 
  mutate(group_trials = as.factor('MODSUG_A_B')) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  filter(species != 'Chordata') %>% 
  
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) 

middle <- MODSUG_A_B %>% mutate(species = str_replace_all(species, '--.*', '')) %>% arrange(effect) %>% slice(7) %>% pull(species)

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
MODSUG_A_B_plot <- MODSUG_A_B %>% ggplot()+
  geom_col(aes(y = effect, x = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = 0.75, x = middle, label = 'MODSUG', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = 'Estimated ALDEx2 Effect Size', x = 'Phylum')+
  scale_y_continuous(breaks = c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75), limits = c(-0.75,0.75))+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(MODSUG_A_B, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'species', xmax = NULL, y.position = -0.75,
             inherit.aes = F, family = 'Times')
MODSUG_A_B_plot

## aldex2_LOWCHO_A_B ####
LOWCHO_A_B <- read_csv(file = './Results/chebi_aldex2_P_LOWCHO_A_B.csv') %>% 
  mutate(group_trials = as.factor('LOWCHO_A_B')) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  filter(species != 'Chordata') %>% 
  
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) 

middle <- LOWCHO_A_B %>% mutate(species = str_replace_all(species, '--.*', '')) %>% arrange(effect) %>% slice(7) %>% pull(species)

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
LOWCHO_A_B_plot <- LOWCHO_A_B %>% ggplot()+
  geom_col(aes(y = effect, x = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = 0.75, x = middle, label = 'LOWCHO', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = 'Estimated ALDEx2 Effect Size', x = 'Phylum')+
  scale_y_continuous(breaks = c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75), limits = c(-0.75,0.75))+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(LOWCHO_A_B, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'species', xmax = NULL, y.position = -0.75,
             inherit.aes = F, family = 'Times')
LOWCHO_A_B_plot




## aldex2_LOWSUG_A_B ####
LOWSUG_A_B <- read_csv(file = './Results/chebi_aldex2_P_LOWSUG_A_B.csv') %>% 
  mutate(group_trials = as.factor('LOWSUG_A_B')) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  filter(species != 'Chordata') %>% 
  
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) 

middle <- LOWSUG_A_B %>% mutate(species = str_replace_all(species, '--.*', '')) %>% arrange(effect) %>% slice(7) %>% pull(species)

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
LOWSUG_A_B_plot <- LOWSUG_A_B %>% ggplot()+
  geom_col(aes(y = effect, x = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = 0.75, x = middle, label = 'LOWSUG', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = 'Estimated ALDEx2 Effect Size', x = 'Phylum')+
  scale_y_continuous(breaks = c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75), limits = c(-0.75,0.75))+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(LOWSUG_A_B, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'species', xmax = NULL, y.position = -0.75,
             inherit.aes = F, family = 'Times')
LOWSUG_A_B_plot



## Phylum A B Arrange -----------------------------------------------------------

ggarrange(MODSUG_A_B_plot, LOWSUG_A_B_plot, LOWCHO_A_B_plot, ncol = 1, nrow = 3, align = 'v')

ggsave(filename = './Figures/chebi_aldex2_effect_Phylum_all_groups_A_B_corrections.tiff',
       device = 'tiff',
       width = 15,
       height = 24,
       units = 'cm',
       dpi = 600)





# Phylum A C -------------------------------------------------------------------
## aldex2_MODSUG_A_C ####

MODSUG_A_C <- read_csv(file = './Results/chebi_aldex2_P_MODSUG_A_C.csv') %>% 
  mutate(group_trials = as.factor('MODSUG_A_C')) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  filter(species != 'Chordata') %>% 
  
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) 

middle <- MODSUG_A_C %>% mutate(species = str_replace_all(species, '--.*', '')) %>% arrange(effect) %>% slice(7) %>% pull(species)

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
MODSUG_A_C_plot <- MODSUG_A_C %>% ggplot()+
  geom_col(aes(y = effect, x = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = 0.75, x = middle, label = 'MODSUG', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = 'Estimated ALDEx2 Effect Size', x = 'Phylum')+
  scale_y_continuous(breaks = c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75), limits = c(-0.75,0.75))+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(MODSUG_A_C, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'species', xmax = NULL, y.position = -0.75,
             inherit.aes = F, family = 'Times')
MODSUG_A_C_plot

## aldex2_LOWCHO_A_C ####
LOWCHO_A_C <- read_csv(file = './Results/chebi_aldex2_P_LOWCHO_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWCHO_A_C')) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  filter(species != 'Chordata') %>% 
  
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) 

middle <- LOWCHO_A_C %>% mutate(species = str_replace_all(species, '--.*', '')) %>% arrange(effect) %>% slice(7) %>% pull(species)

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
LOWCHO_A_C_plot <- LOWCHO_A_C %>% ggplot()+
  geom_col(aes(y = effect, x = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = 0.75, x = middle, label = 'LOWCHO', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = 'Estimated ALDEx2 Effect Size', x = 'Phylum')+
  scale_y_continuous(breaks = c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75), limits = c(-0.75,0.75))+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(LOWCHO_A_C, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'species', xmax = NULL, y.position = -0.75,
             inherit.aes = F, family = 'Times')
LOWCHO_A_C_plot




## aldex2_LOWSUG_A_C ####
LOWSUG_A_C <- read_csv(file = './Results/chebi_aldex2_P_LOWSUG_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWSUG_A_C')) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  filter(species != 'Chordata') %>% 
  
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) 

middle <- LOWSUG_A_C %>% mutate(species = str_replace_all(species, '--.*', '')) %>% arrange(effect) %>% slice(7) %>% pull(species)

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
LOWSUG_A_C_plot <- LOWSUG_A_C %>% ggplot()+
  geom_col(aes(y = effect, x = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = 0.75, x = middle, label = 'LOWSUG', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = 'Estimated ALDEx2 Effect Size', x = 'Phylum')+
  scale_y_continuous(breaks = c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75), limits = c(-0.75,0.75))+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(LOWSUG_A_C, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'species', xmax = NULL, y.position = -0.75,
             inherit.aes = F, family = 'Times')
LOWSUG_A_C_plot



## Phylum A C Arrange -----------------------------------------------------------

ggarrange(MODSUG_A_C_plot, LOWSUG_A_C_plot, LOWCHO_A_C_plot, ncol = 1, nrow = 3, align = 'v')

ggsave(filename = './Figures/chebi_aldex2_effect_Phylum_all_groups_A_C_corrections.tiff',
       device = 'tiff',
       width = 15,
       height = 24,
       units = 'cm',
       dpi = 600)
