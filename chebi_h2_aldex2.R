# Packages ####
library(zCompositions)
library(vegan)
library(ggsci)
library(beepr)
library(ALDEx2)
library(rstatix)
library(ggsci)
library(ggthemes)
library(hrbrthemes)
library(gcookbook)
library(ggpubr)
library(forcats)
library(RColorBrewer)
library(tidyverse)
library(scales)
windowsFonts("Times" = windowsFont("Times"))
source('./Scripts/theme_clean_rgd.R')

# Data ####
source('./Scripts/chebi_h2_load_data.R')
samples$id_2
data_1$id_2
rownames(vegan)

# Between trials ----------------------------------------------------------

# ALDEx2 requires `abund_table` of raw reads, taxa as rows, samples as columns.
# Also requires `groups`, factor variable of groups in correct order.

# Generic setup, only need to edit group and trials ####
samples_filtered <- samples %>% 
  filter(group == 'LOWCHO') %>% # Edit group
  filter(trial == 'A'| trial == 'C') %>% # Edit trials
  filter(participant != '12', participant != '14', participant != '21', participant != '56', participant != '51', participant != '52', participant != '8') %>% # Filter only one trial
  arrange(id_2)

samples_filtered$group
samples_filtered$id_2
data_filtered <- data_1 %>%
  semi_join(samples_filtered, by = 'id_2') %>% 
  arrange(id_2) %>% 
  column_to_rownames(var = 'id_2') %>% 
  t() %>% 
  as.data.frame() 

df <- as.data.frame(sapply(data_filtered, as.integer))
rownames(df) <- rownames(data_filtered)

groups <- samples_filtered %>% 
  pull(trial) %>% 
  as.factor()

colnames(df)
colnames(data_filtered)
samples_filtered$id_2
samples_filtered$group
groups

## Step 1: Generate Instances of the Centred Log-Ratio Transformed Values Using the Function aldex.clr().#### 
vdr <- aldex.clr(df,
                 groups,
                 mc.samples = 1000,
                 verbose = TRUE)
## Step 2: Perform the Welch’s t and Wilcoxon Rank Sum Test Using aldex.ttest().#### 
vdr_t <- aldex.ttest(vdr,
                     groups,
                     paired.test = T)
## Step 3: Estimate Effect Size Using the Function aldex.effect(). ####
vdr_effect <- aldex.effect(vdr,
                           groups,
                           include.sample.summary = FALSE,
                           CI = T,
                           verbose = FALSE)
## Step 4: Merge all Data into One Object and Make a Data Frame for Result Viewing and Downstream Analysis.#### 
vdr_all <- data.frame(vdr_t, vdr_effect) %>% 
  rownames_to_column(var = 'pathway')
## Check for significant results#### 
sig_by_both <- which(vdr_all$we.ep < 0.05 & vdr_all$wi.ep < 0.05)
## Significant after B-H correction#### 
sig_by_both_fdr <- which(vdr_all$we.eBH < 0.25 & vdr_all$wi.eBH < 0.25)
## Plot results#### 
plot(vdr_all$diff.win, vdr_all$diff.btw, pch=19, cex=0.3, col=rgb(0,0,0,0.3), xlab="Dispersion", ylab="Difference")
points(vdr_all$diff.win[sig_by_both], vdr_all$diff.btw[sig_by_both], pch=19, cex=0.7, col=rgb(0,0,1,0.5))
points(vdr_all$diff.win[sig_by_both_fdr], vdr_all$diff.btw[sig_by_both_fdr], pch=19, cex=0.7, col=rgb(1,0,0,1))
abline(0,1,lty=2)
abline(0,-1,lty=2) ; beep(8)
## Save results ####
# write_csv(vdr_all, file = './Results/chebi_G_aldex2_MODSUG_B_C.csv')
vdr_all_edit <- vdr_all %>% 
  mutate(effect = diff.btw/diff.win)

write_csv(vdr_all_edit, file = './Results/chebi_G_unstr_aldex2_LOWCHO_A_C.csv')


# Effect Size Plots -------------------------------------------------------------------

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

# Updated - ggarrange -----------------------------------------------------

limits <- c(-1,1)
breaks <- c(-1,-0.5,0,0.5,1)
y.position <- -2
group_label <- 1

## aldex2_MODSUG_A_C ####

up <-
  read_csv(file = './Results/chebi_G_unstr_aldex2_MODSUG_A_C.csv') %>% 
  mutate(group_trials = as.factor('MODSUG_A_C')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>%
  
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/chebi_G_unstr_aldex2_MODSUG_A_C.csv') %>% 
  mutate(group_trials = as.factor('MODSUG_A_C')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>%
  
  slice_min(effect, n = 10)

middle <- down %>% slice_max(effect) %>% pull(pathway) 

MODSUG_A_C <- bind_rows(up, down) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>%
  mutate(pathway = as.factor(pathway)) %>%
  mutate(pathway = fct_reorder(.f = pathway, .x = effect, .desc = F)) 

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
MODSUG_A_C_plot <- 
  MODSUG_A_C %>% 
  ggplot()+
  geom_col(aes(y = effect, x = pathway, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = pathway, label = wi.eBH, group = pathway), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = group_label, x = middle, label = 'MODSUG', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = NULL, x = NULL)+
  scale_y_continuous(breaks = breaks, limits = limits)+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(MODSUG_A_C, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'pathway', xmax = NULL, y.position = y.position,
             inherit.aes = F, family = 'Times')


## aldex2_LOWCHO_A_C ####
up <- read_csv(file = './Results/chebi_G_unstr_aldex2_LOWCHO_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWCHO_A_C')) %>% 
  # mutate(pathway = str_replace_all(pathway, ':.*', '')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>%
  
  
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/chebi_G_unstr_aldex2_LOWCHO_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWCHO_A_C')) %>% 
  # mutate(pathway = str_replace_all(pathway, ':.*', '')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>%
  
  
  slice_min(effect, n = 10)

middle <- down %>% slice_max(effect) %>% pull(pathway) 

LOWCHO_A_C <- bind_rows(up, down) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  # mutate(pathway = str_replace_all(pathway, ':.*', '')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>%
  
  mutate(pathway = as.factor(pathway)) %>%
  mutate(pathway = fct_reorder(.f = pathway, .x = effect, .desc = F)) 

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
LOWCHO_A_C_plot <- LOWCHO_A_C %>% ggplot()+
  geom_col(aes(y = effect, x = pathway, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = pathway, label = wi.eBH, group = pathway), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = group_label, x = middle, label = 'LOWCHO', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = 'Estimated ALDEx2 Effect Size', x = NULL)+
  scale_y_continuous(breaks = breaks, limits = limits)+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(LOWCHO_A_C, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'pathway', xmax = NULL, y.position = y.position,
             inherit.aes = F, family = 'Times')



## aldex2_LOWSUG_A_C ####
up <- read_csv(file = './Results/chebi_G_unstr_aldex2_LOWSUG_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWSUG_A_C'),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>% 
  # mutate(pathway = str_replace_all(pathway, ':.*', '')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', '')) %>%
  
  
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/chebi_G_unstr_aldex2_LOWSUG_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWSUG_A_C'),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>% 
  # mutate(pathway = str_replace_all(pathway, ':.*', '')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', '')) %>%
  
  
  slice_min(effect, n = 10)

middle <- down %>% slice_max(effect) %>% pull(pathway) 

LOWSUG_A_C <- bind_rows(up, down) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  # mutate(pathway = str_replace_all(pathway, ':.*', '')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>%
  
  mutate(pathway = as.factor(pathway)) %>%
  mutate(pathway = fct_reorder(.f = pathway, .x = effect, .desc = F)) 

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
LOWSUG_A_C_plot <- LOWSUG_A_C %>% ggplot()+
  geom_col(aes(y = effect, x = pathway, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = pathway, label = wi.eBH, group = pathway), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = group_label, x = middle, label = 'LOWSUG', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = NULL, x = 'Pathway')+
  scale_y_continuous(breaks = breaks, limits = limits)+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(LOWSUG_A_C, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'pathway', xmax = NULL, y.position = y.position,
             inherit.aes = F, family = 'Times')



## A_C Arrange -----------------------------------------------------------

ggarrange(MODSUG_A_C_plot, LOWSUG_A_C_plot, LOWCHO_A_C_plot, ncol = 1, nrow = 3, align = 'v')

ggsave(filename = './Figures/chebi_aldex2_effect_G_all_groups_A_C_arranged_2.tiff',
       device = 'tiff',
       width = 15,
       height = 24,
       units = 'cm',
       dpi = 600)


# NEW_AB_Arranged ---------------------------------------------------------

## aldex2_MODSUG_A_B ####

up <- read_csv(file = './Results/chebi_G_unstr_aldex2_MODSUG_A_B.csv') %>% 
  mutate(group_trials = as.factor('MODSUG_A_B')) %>% 
  # mutate(pathway = str_replace_all(pathway, ':.*', '')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>%
  
  
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/chebi_G_unstr_aldex2_MODSUG_A_B.csv') %>% 
  mutate(group_trials = as.factor('MODSUG_A_B')) %>% 
  # mutate(pathway = str_replace_all(pathway, ':.*', '')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>%
  
  
  slice_min(effect, n = 10)

middle <- down %>% slice_max(effect) %>% pull(pathway) 

MODSUG_A_B <- bind_rows(up, down) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  # mutate(pathway = str_replace_all(pathway, ':.*', '')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>%
  
  mutate(pathway = as.factor(pathway)) %>%
  mutate(pathway = fct_reorder(.f = pathway, .x = effect, .desc = F)) 

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
MODSUG_A_B_plot <- MODSUG_A_B %>% ggplot()+
  geom_col(aes(y = effect, x = pathway, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = pathway, label = wi.eBH, group = pathway), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = group_label, x = middle, label = 'MODSUG', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = NULL, x = NULL)+
  scale_y_continuous(breaks = breaks, limits = limits)+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(MODSUG_A_B, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'pathway', xmax = NULL, y.position = y.position,
             inherit.aes = F, family = 'Times')


## aldex2_LOWCHO_A_B ####
up <- read_csv(file = './Results/chebi_G_unstr_aldex2_LOWCHO_A_B.csv') %>% 
  mutate(group_trials = as.factor('LOWCHO_A_B')) %>% 
  # mutate(pathway = str_replace_all(pathway, ':.*', '')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>%
  
  
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/chebi_G_unstr_aldex2_LOWCHO_A_B.csv') %>% 
  mutate(group_trials = as.factor('LOWCHO_A_B')) %>% 
  # mutate(pathway = str_replace_all(pathway, ':.*', '')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>%
  
  
  slice_min(effect, n = 10)

middle <- down %>% slice_max(effect) %>% pull(pathway) 

LOWCHO_A_B <- bind_rows(up, down) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  # mutate(pathway = str_replace_all(pathway, ':.*', '')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>%
  
  mutate(pathway = as.factor(pathway)) %>%
  mutate(pathway = fct_reorder(.f = pathway, .x = effect, .desc = F)) 

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
LOWCHO_A_B_plot <- LOWCHO_A_B %>% ggplot()+
  geom_col(aes(y = effect, x = pathway, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = pathway, label = wi.eBH, group = pathway), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = group_label, x = middle, label = 'LOWCHO', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = 'Estimated ALDEx2 Effect Size', x = NULL)+
  scale_y_continuous(breaks = breaks, limits = limits)+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(LOWCHO_A_B, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'pathway', xmax = NULL, y.position = y.position,
             inherit.aes = F, family = 'Times')



## aldex2_LOWSUG_A_B ####
up <- read_csv(file = './Results/chebi_G_unstr_aldex2_LOWSUG_A_B.csv') %>% 
  mutate(group_trials = as.factor('LOWSUG_A_B')) %>% 
  # mutate(pathway = str_replace_all(pathway, ':.*', '')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>%
  slice_max(effect, n = 10)

down <- read_csv(file = './Results/chebi_G_unstr_aldex2_LOWSUG_A_B.csv') %>% 
  mutate(group_trials = as.factor('LOWSUG_A_B')) %>% 
  # mutate(pathway = str_replace_all(pathway, ':.*', '')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>%
  
  
  slice_min(effect, n = 10)

middle <- down %>% slice_max(effect) %>% pull(pathway) 

LOWSUG_A_B <- bind_rows(up, down) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1), symbols = c("****", "***", "**", "*", "")) %>% 
  # mutate(pathway = str_replace_all(pathway, ':.*', '')) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP')) %>%
  
  mutate(pathway = as.factor(pathway)) %>%
  mutate(pathway = fct_reorder(.f = pathway, .x = effect, .desc = F)) 

# With Wilcoxon BH corrected sig * <0.1, ** <0.05
LOWSUG_A_B_plot <- LOWSUG_A_B %>% ggplot()+
  geom_col(aes(y = effect, x = pathway, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = pathway, label = wi.eBH, group = pathway), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', y = group_label, x = middle, label = 'LOWSUG', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10)
  )+
  labs(y = NULL, x = 'Pathway')+
  scale_y_continuous(breaks = breaks, limits = limits)+
  guides(fill = 'none')+
  coord_flip()+
  add_pvalue(LOWSUG_A_B, label = "wi.eBH.signif", remove.bracket = TRUE, xmin = 'pathway', xmax = NULL, y.position = y.position,
             inherit.aes = F, family = 'Times')



## A_B Arrange -----------------------------------------------------------

ggarrange(MODSUG_A_B_plot, LOWSUG_A_B_plot, LOWCHO_A_B_plot, ncol = 1, nrow = 3, align = 'v')

ggsave(filename = './Figures/chebi_aldex2_effect_G_all_groups_A_B_arranged_2.tiff',
       device = 'tiff',
       width = 15,
       height = 24,
       units = 'cm',
       dpi = 600)




# Unsure? -----------------------------------------------------------------



library(rstatix)
library(ggsci)
library(ggthemes)
library(hrbrthemes)
library(gcookbook)
library(ggpubr)
library(forcats)
library(RColorBrewer)
library(tidyverse)
library(scales)
show_col(pal_futurama("planetexpress")(12))
futurama <- pal_futurama("planetexpress")(12)
source('./Scripts/theme_clean_rgd.R')
windowsFonts("Times" = windowsFont("Times"))

## aldex2_LOWCHO_A_B ####
read_csv(file = './Results/chebi_G_unstr_aldex2_LOWCHO_A_B.csv') %>% 
  mutate(group_trials = as.factor('LOWCHO_A_B')) %>% 
  filter(effect <= -0.55 | effect >= 0.45) %>% 
  arrange(effect) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) %>% 
  
  ggplot()+
  geom_col(aes(x = effect, y = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  ggtitle('LOWCHO: Baseline - 4 Weeks')+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.x = element_text(size = 8))+
  labs(x = 'Estimated ALDEx2 Effect Size', y = 'Species')+
  scale_x_continuous(breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), limits = c(-1.5,1.5))+
  # xlim(c(-1.5,1.5))+
  guides(fill = 'none')

ggsave(filename = './Figures/chebi_aldex2_effect_LOWCHO_A_B.tiff',
       device = 'tiff',
       width = 16,
       height = 8,
       units = 'cm',
       dpi = 600)


## aldex2_MODSUG_A_B ####
read_csv(file = './Results/chebi_aldex2_MODSUG_A_B.csv') %>% 
  mutate(group_trials = as.factor('MODSUG_A_B')) %>% 
  filter(effect <= -0.39 | effect >= 0.4) %>% 
  arrange(effect) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) %>% 
  
  ggplot()+
  geom_col(aes(x = effect, y = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  ggtitle('MODSUG: Baseline - 4 Weeks')+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.x = element_text(size = 8))+
  labs(x = 'Estimated ALDEx2 Effect Size', y = 'Species')+
  scale_x_continuous(breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), limits = c(-1.5,1.5))+
  guides(fill = 'none')

ggsave(filename = './Figures/chebi_aldex2_effect_MODSUG_A_B.tiff',
       device = 'tiff',
       width = 16,
       height = 8,
       units = 'cm',
       dpi = 600)


## aldex2_LOWSUG_A_B ####
read_csv(file = './Results/chebi_aldex2_LOWSUG_A_B.csv') %>% 
  mutate(group_trials = as.factor('LOWSUG_A_B')) %>% 
  filter(effect <= -0.4 | effect >= 0.45) %>% 
  arrange(effect) %>%
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) %>% 
  
  ggplot()+
  geom_col(aes(x = effect, y = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  ggtitle('LOWSUG: Baseline - 4 Weeks')+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.x = element_text(size = 8))+
  labs(x = 'Estimated ALDEx2 Effect Size', y = 'Species')+
  scale_x_continuous(breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), limits = c(-1.5,1.5))+
  guides(fill = 'none')

ggsave(filename = './Figures/chebi_aldex2_effect_LOWSUG_A_B.tiff',
       device = 'tiff',
       width = 16,
       height = 8,
       units = 'cm',
       dpi = 600)


## aldex2_MODSUG_A_C ####
read_csv(file = './Results/chebi_aldex2_MODSUG_A_C.csv') %>% 
  mutate(group_trials = as.factor('MODSUG_A_C')) %>% 
  filter(effect <= -0.45 | effect >= 0.4) %>% 
  arrange(effect) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = T)) %>% 
  
  ggplot()+
  geom_col(aes(x = effect, y = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean()+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.x = element_text(size = 8))+
  labs(x = 'Estimated ALDEx2 Effect Size', y = 'Species')+
  guides(fill = 'none')

## aldex2_LOWSUG_A_C ####
read_csv(file = './Results/chebi_aldex2_LOWSUG_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWSUG_A_C')) %>% 
  filter(effect <= -0.4 | effect >= 0.4) %>% 
  arrange(effect) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = T)) %>% 
  
  ggplot()+
  geom_col(aes(x = effect, y = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean()+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.x = element_text(size = 8))+
  labs(x = 'Estimated ALDEx2 Effect Size', y = 'Species')+
  guides(fill = 'none')

## aldex2_LOWCHO_A_C ####
read_csv(file = './Results/chebi_aldex2_LOWCHO_A_C.csv') %>% 
  mutate(group_trials = as.factor('LOWCHO_A_C')) %>% 
  filter(effect <= -0.55 | effect >= 0.45) %>% 
  arrange(effect) %>% 
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) %>% 
  
  ggplot()+
  geom_col(aes(x = effect, y = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean()+
  theme(text = element_text(family = 'Times', colour = 'black', size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.x = element_text(size = 8))+
  labs(x = 'Estimated ALDEx2 Effect Size', y = 'Species')+
  xlim(c(-1.5,1.5))+
  guides(fill = 'none')

ggsave(filename = './Figures/chebi_aldex2_effect_LOWCHO_A_C.tiff',
       device = 'tiff',
       width = 16,
       height = 12,
       units = 'cm',
       dpi = 600)


## aldex2_MODSUG_B_C #### 
# Error - all effect sizes the same so plot not possible
# read_csv(file = './Results/chebi_aldex2_MODSUG_B_C.csv') %>% view() 
#   mutate(group_trials = as.factor('MODSUG_B_C')) %>% 
#   filter(effect <= -0.15 | effect >= 0.15) %>% 
#   arrange(effect) %>% 
#   mutate(species = as.factor(species)) %>%
#   mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) %>% 
#   
#   ggplot()+
#   geom_col(aes(x = effect, y = species, fill = effect), colour = 'black')+
#   scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
#   # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
#   # scale_fill_viridis_b(option = 'plasma')+
#   theme_clean()+
#   theme(text = element_text(family = 'Times', colour = 'black', size = 10),
#         plot.title = element_text(size = 10),
#         plot.background = element_rect(color = "white"),
#         axis.text.x = element_text(size = 8))+
#   labs(x = 'Estimated ALDEx2 Effect Size', y = 'Species')+
#   # xlim(c(-1.2,1.2))+
#   guides(fill = 'none')

## aldex2_LOWSUG_B_C ####
read_csv(file = './Results/chebi_aldex2_LOWSUG_B_C.csv') %>% 
  mutate(group_trials = as.factor('LOWSUG_B_C')) %>% 
  filter(effect <= -0.42 | effect >= 0.42) %>% 
  arrange(effect) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) %>% 
  
  ggplot()+
  geom_col(aes(x = effect, y = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean()+
  theme(text = element_text(family = 'Times', colour = 'black', size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.x = element_text(size = 8))+
  labs(x = 'Estimated ALDEx2 Effect Size', y = 'Species')+
  xlim(c(-1.2,1.2))+
  guides(fill = 'none')

## aldex2_LOWCHO_B_C ####
read_csv(file = './Results/chebi_aldex2_LOWCHO_B_C.csv') %>% 
  mutate(group_trials = as.factor('LOWCHO_B_C')) %>% 
  filter(effect <= -0.41 | effect >= 0.35) %>% 
  arrange(effect) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) %>% 
  
  ggplot()+
  geom_col(aes(x = effect, y = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean()+
  theme(text = element_text(family = 'Times', colour = 'black', size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.x = element_text(size = 8))+
  labs(x = 'Estimated ALDEx2 Effect Size', y = 'Species')+
  xlim(c(-1.2,1.2))+
  guides(fill = 'none')

## Faceted_All_A_B --------------------------------------------

# Edit
aldex2_MODSUG_A_B <- read_csv(file = './Results/chebi_aldex2_MODSUG_A_B.csv')
  
aldex2_MODSUG_A_B %>% 
  arrange(effect) %>% 
  slice(1:10,2989:2998) %>% view()
  


aldex2_MODSUG_A_B <- read_csv(file = './Results/chebi_aldex2_MODSUG_A_B.csv') %>% 
  mutate(group_trials = 'MODSUG') %>% 
  arrange(effect) %>% 
  slice(1:10,2989:2998)
aldex2_LOWSUG_A_B <- read_csv(file = './Results/chebi_aldex2_LOWSUG_A_B.csv') %>% 
  mutate(group_trials = 'LOWSUG') %>% 
  arrange(effect) %>% 
  slice(1:10,2989:2998)
aldex2_LOWCHO_A_B <- read_csv(file = './Results/chebi_aldex2_LOWCHO_A_B.csv') %>% 
  mutate(group_trials = 'LOWCHO') %>% 
  arrange(effect) %>% 
  slice(1:10,2989:2998)

aldex2_all_A_B <- bind_rows(aldex2_MODSUG_A_B, aldex2_LOWSUG_A_B, aldex2_LOWCHO_A_B) %>% 
  mutate(group_trials = factor(group_trials, levels = c('MODSUG','LOWSUG','LOWCHO')))

aldex2_all_A_B %>%
  # filter(effect <= -0.4 | effect >= 0.4) %>%
  arrange(effect) %>%
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) %>%

  ggplot()+
  geom_col(aes(x = effect, y = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(face = 'italic'))+
  labs(x = 'Estimated ALDEx2 Effect Size', y = 'Species')+
  scale_x_continuous(breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), limits = c(-1.5,1.5))+
  guides(fill = 'none')+
  facet_wrap(facets = vars(group_trials),
            nrow = 3, ncol = 1,
            drop = T,
            scales = 'free',
            strip.position = 'right')
  # facet_grid(rows = vars(group_trials),
  #            drop = T, 
  #            scales = 'free')

ggsave(filename = './Figures/chebi_aldex2_effect_all_groups_A_B.tiff',
       device = 'tiff',
       width = 16,
       height = 24,
       units = 'cm',
       dpi = 600)

## Faceted_All_A_C --------------------------------------------


aldex2_MODSUG_A_C <- read_csv(file = './Results/chebi_aldex2_MODSUG_A_C.csv') %>% 
  mutate(group_trials = 'MODSUG') %>% 
  arrange(effect) %>% 
  slice(1:10,2989:2998)
aldex2_LOWSUG_A_C <- read_csv(file = './Results/chebi_aldex2_LOWSUG_A_C.csv') %>% 
  mutate(group_trials = 'LOWSUG') %>% 
  arrange(effect) %>% 
  slice(1:10,2989:2998)
aldex2_LOWCHO_A_C <- read_csv(file = './Results/chebi_aldex2_LOWCHO_A_C.csv') %>% 
  mutate(group_trials = 'LOWCHO') %>% 
  arrange(effect) %>% 
  slice(1:10,2989:2998)

aldex2_all_A_C <- bind_rows(aldex2_MODSUG_A_C, aldex2_LOWSUG_A_C, aldex2_LOWCHO_A_C) %>% 
  mutate(group_trials = factor(group_trials, levels = c('MODSUG','LOWSUG','LOWCHO')))

aldex2_all_A_C %>%
  # filter(effect <= -0.4 | effect >= 0.4) %>%
  arrange(effect) %>%
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) %>% view()
  
  ggplot()+
  geom_col(aes(x = effect, y = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(face = 'italic'))+
  labs(x = 'Estimated ALDEx2 Effect Size', y = 'Species')+
  scale_x_continuous(breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), limits = c(-1.5,1.5))+
  guides(fill = 'none')+
  facet_wrap(facets = vars(group_trials),
             nrow = 3, ncol = 1,
             drop = T,
             scales = 'free',
             strip.position = 'right')
# facet_grid(rows = vars(group_trials),
#            drop = T, 
#            scales = 'free')

ggsave(filename = './Figures/chebi_aldex2_effect_all_groups_A_C.tiff',
       device = 'tiff',
       width = 16,
       height = 24,
       units = 'cm',
       dpi = 600)

## Faceted_All_B_C --------------------------------------------


aldex2_MODSUG_B_C <- read_csv(file = './Results/chebi_aldex2_MODSUG_B_C.csv') %>% 
  mutate(group_trials = 'MODSUG') %>% 
  arrange(effect) %>% 
  slice(1:10,2989:2998)
aldex2_LOWSUG_B_C <- read_csv(file = './Results/chebi_aldex2_LOWSUG_B_C.csv') %>% 
  mutate(group_trials = 'LOWSUG') %>% 
  arrange(effect) %>% 
  slice(1:10,2989:2998)
aldex2_LOWCHO_B_C <- read_csv(file = './Results/chebi_aldex2_LOWCHO_B_C.csv') %>% 
  mutate(group_trials = 'LOWCHO') %>% 
  arrange(effect) %>% 
  slice(1:10,2989:2998)

aldex2_all_B_C <- bind_rows(aldex2_MODSUG_B_C, aldex2_LOWSUG_B_C, aldex2_LOWCHO_B_C) %>% 
  mutate(group_trials = factor(group_trials, levels = c('MODSUG','LOWSUG','LOWCHO')))

aldex2_all_B_C %>%
  # filter(effect <= -0.4 | effect >= 0.4) %>%
  arrange(effect) %>%
  mutate(species = str_replace_all(species, '--.*', '')) %>% 
  mutate(species = as.factor(species)) %>%
  mutate(species = fct_reorder(.f = species, .x = effect, .desc = F)) %>%
  
  ggplot()+
  geom_col(aes(x = effect, y = species, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  # geom_text(aes(x = effect, y = species, label = wi.eBH, group = species), size = 2, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(face = 'italic'))+
  labs(x = 'Estimated ALDEx2 Effect Size', y = 'Species')+
  scale_x_continuous(breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), limits = c(-1.5,1.5))+
  guides(fill = 'none')+
  facet_wrap(facets = vars(group_trials),
             nrow = 3, ncol = 1,
             drop = T,
             scales = 'free',
             strip.position = 'right')
# facet_grid(rows = vars(group_trials),
#            drop = T, 
#            scales = 'free')

ggsave(filename = './Figures/chebi_aldex2_effect_all_groups_B_C.tiff',
       device = 'tiff',
       width = 16,
       height = 24,
       units = 'cm',
       dpi = 600)


# Between groups ----------------------------------------------------------
## Data ####
source('./Scripts/chebi_G_load_data.R')
samples <- samples %>% 
  arrange(participant, trial)

# ALDEx2 requires `abund_table` of raw reads, taxa as rows, samples as columns.
# Also requires `groups`, factor variable of groups in correct order.

## Generic setup, only need to edit group and trials ####
samples_filtered <- samples %>% 
  
  # filter(group == 'MODSUG' | group == 'LOWSUG') %>% # Edit groups
  # filter(group == 'MODSUG' | group == 'LOWCHO') %>% # Edit groups
  filter(group == 'LOWSUG' | group == 'LOWCHO') %>% # Edit groups
  
  # filter(trial == 'A') %>% # Edit trial
  # filter(trial == 'B') %>% # Edit trial
  filter(trial == 'C') %>% # Edit trial
  arrange(id_2)

samples_filtered$group
samples_filtered$id_2
data_filtered <- data_1 %>% 
  semi_join(samples_filtered, by = 'id_2') %>% 
  arrange(id_2) %>% 
  column_to_rownames(var = 'id_2') %>% 
  t() %>% 
  as.data.frame() 

df <- as.data.frame(sapply(data_filtered, as.integer))
rownames(df) <- rownames(data_filtered)

colnames(data_filtered)
groups <- samples_filtered %>% 
  # mutate(group = factor(group, levels = c('MODSUG','LOWSUG'))) %>% # Edit group levels
  # mutate(group = factor(group, levels = c('MODSUG','LOWCHO'))) %>% # Edit group levels
  mutate(group = factor(group, levels = c('LOWSUG','LOWCHO'))) %>% # Edit group levels
  pull(group) %>% 
  as.factor()
groups

## Step 1: Generate Instances of the Centred Log-Ratio Transformed Values Using the Function aldex.clr().#### 
vdr <- aldex.clr(df,
                 groups,
                 mc.samples = 1000,
                 verbose = TRUE)
## Step 2: Perform the Welch’s t and Wilcoxon Rank Sum Test Using aldex.ttest().#### 
vdr_t <- aldex.ttest(vdr,
                     groups,
                     paired.test = F)
## Step 3: Estimate Effect Size Using the Function aldex.effect(). ####
vdr_effect <- aldex.effect(vdr,
                           groups,
                           verbose = T, 
                           CI = T)
## Step 4: Merge all Data into One Object and Make a Data Frame for Result Viewing and Downstream Analysis.#### 
vdr_all <- data.frame(vdr_t, vdr_effect) %>% 
  rownames_to_column(var = 'pathway')
# Check for significant results#### 
sig_by_both <- which(vdr_all$we.ep < 0.05 & vdr_all$wi.ep < 0.05)
# Significant after B-H correction#### 
sig_by_both_fdr <- which(vdr_all$we.eBH < 0.2 & vdr_all$wi.eBH < 0.2)
# Plot results#### 
plot(vdr_all$diff.win, vdr_all$diff.btw, pch=19, cex=0.3, col=rgb(0,0,0,0.3), xlab="Dispersion", ylab="Difference")
points(vdr_all$diff.win[sig_by_both], vdr_all$diff.btw[sig_by_both], pch=19, cex=0.7, col=rgb(0,0,1,0.5))
points(vdr_all$diff.win[sig_by_both_fdr], vdr_all$diff.btw[sig_by_both_fdr], pch=19, cex=0.7, col=rgb(1,0,0,1))
abline(0,1,lty=2)
abline(0,-1,lty=2) 
# Save results ####
vdr_all_edit <- vdr_all %>% 
  mutate(effect = diff.btw/diff.win) ; beep(1)

write_csv(vdr_all_edit, file = './Results/ALDEx2_between/chebi_aldex2_G_LOWSUG_LOWCHO_C.csv')

# Between groups plots ----------------------------------------------------

library(ggprism)
library(rstatix)
library(ggsci)
library(ggthemes)
library(hrbrthemes)
library(gcookbook)
library(ggpubr)
library(forcats)
library(RColorBrewer)
library(tidyverse)
library(scales)
show_col(pal_futurama("planetexpress")(12))
futurama <- pal_futurama("planetexpress")(12)
source('./Scripts/theme_clean_rgd.R')
windowsFonts("Times" = windowsFont("Times"))


# MODSUG_LOWSUG -----------------------------------------------------------
## h2_MODSUG_LOWSUG_A_plot ####
up <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_MODSUG_LOWSUG_A.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
         ) %>%
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_MODSUG_LOWSUG_A.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%
  slice_min(effect, n = 10)
middle <- 
  down %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  slice_max(effect) %>% 
  pull(pathway) 
top <- up %>% 
  arrange(effect, desc = T) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  pull(pathway) 
bottom <- down %>% 
  arrange(effect) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  pull(pathway)
h2_MODSUG_LOWSUG_A <-
  bind_rows(up, down) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  )%>%
  mutate(pathway = as.factor(pathway)) %>%
  mutate(pathway = fct_reorder(.f = pathway, .x = effect, .desc = F)) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1),
                   symbols = c("****", "***", "**", "*", "")) 
h2_MODSUG_LOWSUG_A_plot <-
  h2_MODSUG_LOWSUG_A %>% 
  ggplot()+
  geom_col(aes(x = effect, y = pathway, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  geom_text(aes(x = effect, y = pathway, label = wi.eBH.signif, group = pathway),
            size = 6, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', x = 2, y = middle, label = 'Baseline', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
  )+
  labs(x = NULL, y = 'Pathway')+
  scale_x_continuous(breaks = c(-2,-1,0,1,2), limits = c(-2,2))+
  # xlim(c(-1.5,1.5))+
  guides(fill = 'none')+
  
  annotate(geom = 'text', x = 1.3, y = top, label = '\u2192 Higher\nin LOWSUG',
           family = 'Times', size = 3)+
  annotate(geom = 'text', x = -1.3, y = bottom, label = 'Higher in\nMODSUG \u2190',
           family = 'Times', size = 3)


## h2_MODSUG_LOWSUG_B_plot ####
up <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_MODSUG_LOWSUG_B.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_MODSUG_LOWSUG_B.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%
  slice_min(effect, n = 10)
middle <- 
  down %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  slice_max(effect) %>% 
  pull(pathway) 
top <- up %>% 
  arrange(effect, desc = T) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  pull(pathway) 
bottom <- down %>% 
  arrange(effect) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  pull(pathway)
h2_MODSUG_LOWSUG_B <-
  bind_rows(up, down) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%
  mutate(pathway = as.factor(pathway)) %>%
  mutate(pathway = fct_reorder(.f = pathway, .x = effect, .desc = F)) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1),
                   symbols = c("****", "***", "**", "*", "")) 
h2_MODSUG_LOWSUG_B_plot <-
  h2_MODSUG_LOWSUG_B %>% 
  ggplot()+
  geom_col(aes(x = effect, y = pathway, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  geom_text(aes(x = effect, y = pathway, label = wi.eBH.signif, group = pathway),
            size = 6, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', x = 2, y = middle, label = 'Week 4', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
  )+
  labs(x = NULL, y = 'Pathway')+
  scale_x_continuous(breaks = c(-2,-1,0,1,2), limits = c(-2,2))+
  # xlim(c(-1.5,1.5))+
  guides(fill = 'none')+
  
  annotate(geom = 'text', x = 1.3, y = top, label = '\u2192 Higher\nin LOWSUG',
           family = 'Times', size = 3)+
  annotate(geom = 'text', x = -1.3, y = bottom, label = 'Higher in\nMODSUG \u2190',
           family = 'Times', size = 3)


## h2_MODSUG_LOWSUG_C_plot ####
up <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_MODSUG_LOWSUG_C.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  ) %>%
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_MODSUG_LOWSUG_C.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  ) %>%
  slice_min(effect, n = 10)
middle <- 
  down %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  ) %>%  slice_max(effect) %>% 
  pull(pathway) 
top <- up %>% 
  arrange(effect, desc = T) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  ) %>%  pull(pathway) 
bottom <- down %>% 
  arrange(effect) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  )%>% pull(pathway)
h2_MODSUG_LOWSUG_C <-
  bind_rows(up, down) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  ) %>%
  mutate(pathway = as.factor(pathway)) %>%
  mutate(pathway = fct_reorder(.f = pathway, .x = effect, .desc = F)) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1),
                   symbols = c("****", "***", "**", "*", "")) 
h2_MODSUG_LOWSUG_C_plot <-
  h2_MODSUG_LOWSUG_C %>% 
  ggplot()+
  geom_col(aes(x = effect, y = pathway, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  geom_text(aes(x = effect, y = pathway, label = wi.eBH.signif, group = pathway),
            size = 6, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', x = 2, y = middle, label = 'Week 12', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
  )+
  labs(x = 'Estimated ALDEx2 Effect Size', y = 'Pathway')+
  scale_x_continuous(breaks = c(-2,-1,0,1,2), limits = c(-2,2))+
  # xlim(c(-1.5,1.5))+
  guides(fill = 'none')+
  
  annotate(geom = 'text', x = 1.3, y = top, label = '\u2192 Higher\nin LOWSUG',
           family = 'Times', size = 3)+
  annotate(geom = 'text', x = -1.3, y = bottom, label = 'Higher in\nMODSUG \u2190',
           family = 'Times', size = 3)


## Arrange h2_MODSUG_LOWSUG -------------------------------------------------
ggarrange(h2_MODSUG_LOWSUG_A_plot, h2_MODSUG_LOWSUG_B_plot, h2_MODSUG_LOWSUG_C_plot,
          ncol = 1, nrow = 3, align = 'v',
          labels = 'AUTO', font.label = list(family = 'Times'), vjust = 1)

ggsave(filename = './Figures/chebi_h2_MODSUG_LOWSUG_ABC_arranged.tiff',
       device = 'tiff',
       width = 15,
       height = 24,
       units = 'cm',
       dpi = 600)

# MODSUG_LOWCHO -----------------------------------------------------------
## h2_MODSUG_LOWCHO_A_plot ####
up <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_MODSUG_LOWCHO_A.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_MODSUG_LOWCHO_A.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%
  slice_min(effect, n = 10)
middle <- 
  down %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  slice_max(effect) %>% 
  pull(pathway) 
top <- up %>% 
  arrange(effect, desc = T) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  pull(pathway) 
bottom <- down %>% 
  arrange(effect) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  pull(pathway)
h2_MODSUG_LOWCHO_A <-
  bind_rows(up, down) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  )%>%
  mutate(pathway = as.factor(pathway)) %>%
  mutate(pathway = fct_reorder(.f = pathway, .x = effect, .desc = F)) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1),
                   symbols = c("****", "***", "**", "*", "")) 
h2_MODSUG_LOWCHO_A_plot <-
  h2_MODSUG_LOWCHO_A %>% 
  ggplot()+
  geom_col(aes(x = effect, y = pathway, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  geom_text(aes(x = effect, y = pathway, label = wi.eBH.signif, group = pathway),
            size = 6, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', x = 2, y = middle, label = 'Baseline', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
  )+
  labs(x = NULL, y = 'Pathway')+
  scale_x_continuous(breaks = c(-2,-1,0,1,2), limits = c(-2,2))+
  # xlim(c(-1.5,1.5))+
  guides(fill = 'none')+
  
  annotate(geom = 'text', x = 1.4, y = top, label = '\u2192\nHigher in\nLOWCHO',
           family = 'Times', size = 3)+
  annotate(geom = 'text', x = -1.4, y = bottom, label = 'Higher in\nMODSUG\n\u2190',
           family = 'Times', size = 3)


## h2_MODSUG_LOWCHO_B_plot ####
up <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_MODSUG_LOWCHO_B.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_MODSUG_LOWCHO_B.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%
  slice_min(effect, n = 10)
middle <- 
  down %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  slice_max(effect) %>% 
  pull(pathway) 
top <- up %>% 
  arrange(effect, desc = T) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  pull(pathway) 
bottom <- down %>% 
  arrange(effect) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  pull(pathway)
h2_MODSUG_LOWCHO_B <-
  bind_rows(up, down) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%
  mutate(pathway = as.factor(pathway)) %>%
  mutate(pathway = fct_reorder(.f = pathway, .x = effect, .desc = F)) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.105,1),
                   symbols = c("****", "***", "**", "*", "")) 
h2_MODSUG_LOWCHO_B_plot <-
  h2_MODSUG_LOWCHO_B %>% 
  ggplot()+
  geom_col(aes(x = effect, y = pathway, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  geom_text(aes(x = 1.3, y = pathway, label = wi.eBH.signif, group = pathway),
            size = 6, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', x = 2, y = middle, label = 'Week 4', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
  )+
  labs(x = NULL, y = 'Pathway')+
  scale_x_continuous(breaks = c(-2,-1,0,1,2), limits = c(-2,2))+
  # xlim(c(-1.5,1.5))+
  guides(fill = 'none')+
  
  annotate(geom = 'text', x = 1.4, y = top, label = '\u2192\nHigher in\nLOWCHO',
           family = 'Times', size = 3)+
  annotate(geom = 'text', x = -1.4, y = bottom, label = 'Higher in\nMODSUG\n\u2190',
           family = 'Times', size = 3)#+
  # geom_text(aes(x = 1.2, y = pathway, label = wi.eBH.signif, group = pathway),
  #           size = 6, position = position_dodge2(width = 1, reverse = T))
  

## h2_MODSUG_LOWCHO_C_plot ####
up <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_MODSUG_LOWCHO_C.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  ) %>%
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_MODSUG_LOWCHO_C.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  ) %>%
  slice_min(effect, n = 10)
middle <- 
  down %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  ) %>%  slice_max(effect) %>% 
  pull(pathway) 
top <- up %>% 
  arrange(effect, desc = T) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  ) %>%  pull(pathway) 
bottom <- down %>% 
  arrange(effect) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  )%>% pull(pathway)
h2_MODSUG_LOWCHO_C <-
  bind_rows(up, down) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  ) %>%
  mutate(pathway = as.factor(pathway)) %>%
  mutate(pathway = fct_reorder(.f = pathway, .x = effect, .desc = F)) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.105,1),
                   symbols = c("****", "***", "**", "*", "")) 
h2_MODSUG_LOWCHO_C_plot <-
  h2_MODSUG_LOWCHO_C %>% 
  ggplot()+
  geom_col(aes(x = effect, y = pathway, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  geom_text(aes(x = effect, y = pathway, label = wi.eBH.signif, group = pathway),
            size = 6, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', x = 2, y = middle, label = 'Week 12', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
  )+
  labs(x = 'Estimated ALDEx2 Effect Size', y = 'Pathway')+
  scale_x_continuous(breaks = c(-2,-1,0,1,2), limits = c(-2,2))+
  # xlim(c(-1.5,1.5))+
  guides(fill = 'none')+
  
  annotate(geom = 'text', x = 1.4, y = top, label = '\u2192\nHigher in\nLOWCHO',
           family = 'Times', size = 3)+
  annotate(geom = 'text', x = -1.4, y = bottom, label = 'Higher in\nMODSUG\n\u2190',
           family = 'Times', size = 3)

## Arrange h2_MODSUG_LOWCHO -------------------------------------------------
ggarrange(h2_MODSUG_LOWCHO_A_plot, h2_MODSUG_LOWCHO_B_plot, h2_MODSUG_LOWCHO_C_plot,
          ncol = 1, nrow = 3, align = 'v',
          labels = 'AUTO', font.label = list(family = 'Times'), vjust = 1)

ggsave(filename = './Figures/chebi_h2_MODSUG_LOWCHO_ABC_arranged_pval.tiff',
       device = 'tiff',
       width = 15,
       height = 24,
       units = 'cm',
       dpi = 600)

# LOWSUG_LOWCHO -----------------------------------------------------------

## h2_LOWSUG_LOWCHO_A_plot ####
up <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_LOWSUG_LOWCHO_A.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_LOWSUG_LOWCHO_A.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%
  slice_min(effect, n = 10)
middle <- 
  down %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  slice_max(effect) %>% 
  pull(pathway) 
top <- up %>% 
  arrange(effect, desc = T) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  pull(pathway) 
bottom <- down %>% 
  arrange(effect) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  pull(pathway)
h2_LOWSUG_LOWCHO_A <-
  bind_rows(up, down) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  )%>%
  mutate(pathway = as.factor(pathway)) %>%
  mutate(pathway = fct_reorder(.f = pathway, .x = effect, .desc = F)) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1),
                   symbols = c("****", "***", "**", "*", "")) 
h2_LOWSUG_LOWCHO_A_plot <-
  h2_LOWSUG_LOWCHO_A %>% 
  ggplot()+
  geom_col(aes(x = effect, y = pathway, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  geom_text(aes(x = effect, y = pathway, label = wi.eBH.signif, group = pathway),
            size = 6, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', x = 2, y = middle, label = 'Baseline', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
  )+
  labs(x = NULL, y = 'Pathway')+
  scale_x_continuous(breaks = c(-2,-1,0,1,2), limits = c(-2,2))+
  # xlim(c(-1.5,1.5))+
  guides(fill = 'none')+
  
  annotate(geom = 'text', x = 1.4, y = top, label = '\u2192\nHigher in\nLOWCHO',
           family = 'Times', size = 3)+
  annotate(geom = 'text', x = -1.4, y = bottom, label = 'Higher in\nLOWSUG\n\u2190',
           family = 'Times', size = 3)


## h2_LOWSUG_LOWCHO_B_plot ####
up <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_LOWSUG_LOWCHO_B.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_LOWSUG_LOWCHO_B.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%
  slice_min(effect, n = 10)
middle <- 
  down %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  slice_max(effect) %>% 
  pull(pathway) 
top <- up %>% 
  arrange(effect, desc = T) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  pull(pathway) 
bottom <- down %>% 
  arrange(effect) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%  pull(pathway)
h2_LOWSUG_LOWCHO_B <-
  bind_rows(up, down) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.')
  ) %>%
  mutate(pathway = as.factor(pathway)) %>%
  mutate(pathway = fct_reorder(.f = pathway, .x = effect, .desc = F)) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1),
                   symbols = c("****", "***", "**", "*", "")) 
h2_LOWSUG_LOWCHO_B_plot <-
  h2_LOWSUG_LOWCHO_B %>% 
  ggplot()+
  geom_col(aes(x = effect, y = pathway, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  geom_text(aes(x = effect, y = pathway, label = wi.eBH.signif, group = pathway),
            size = 6, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', x = 2, y = middle, label = 'Week 4', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
  )+
  labs(x = NULL, y = 'Pathway')+
  scale_x_continuous(breaks = c(-2,-1,0,1,2), limits = c(-2,2))+
  # xlim(c(-1.5,1.5))+
  guides(fill = 'none')+
  
  annotate(geom = 'text', x = 1.4, y = top, label = '\u2192\nHigher in\nLOWCHO',
           family = 'Times', size = 3)+
  annotate(geom = 'text', x = -1.4, y = bottom, label = 'Higher in\nLOWSUG\n\u2190',
           family = 'Times', size = 3)

## h2_LOWSUG_LOWCHO_C_plot ####
up <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_LOWSUG_LOWCHO_C.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  ) %>%
  slice_max(effect, n = 10)
down <- read_csv(file = './Results/ALDEx2_between/chebi_aldex2_h2_LOWSUG_LOWCHO_C.csv') %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  ) %>%
  slice_min(effect, n = 10)
middle <- 
  down %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  ) %>%  slice_max(effect) %>% 
  pull(pathway) 
top <- up %>% 
  arrange(effect, desc = T) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  ) %>%  pull(pathway) 
bottom <- down %>% 
  arrange(effect) %>% 
  slice(5) %>%
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  )%>% pull(pathway)
h2_LOWSUG_LOWCHO_C <-
  bind_rows(up, down) %>% 
  mutate(pathway = str_replace_all(pathway, '.*:', ''),
         pathway = str_replace_all(pathway, 'superpathway', 'SP'),
         pathway = str_replace_all(pathway, 'formation', 'f.'),
         
         pathway = str_replace_all(pathway, 'biosynthesis', 'biosyn.'),
         pathway = str_replace_all(pathway, 'containing', 'con.'),
         pathway = str_replace_all(pathway, 'acetylmuramoyl-pentapeptide', 'A-P'),
         pathway = str_replace_all(pathway, 'acetylglucosamine', 'AG'),
         pathway = str_replace_all(pathway, 'pyruvate','pyr')
  ) %>%
  mutate(pathway = as.factor(pathway)) %>%
  mutate(pathway = fct_reorder(.f = pathway, .x = effect, .desc = F)) %>% 
  add_significance("wi.eBH", cutpoints = c(0,0.001,0.01,0.05,0.1,1),
                   symbols = c("****", "***", "**", "*", "")) 
h2_LOWSUG_LOWCHO_C_plot <-
  h2_LOWSUG_LOWCHO_C %>% 
  ggplot()+
  geom_col(aes(x = effect, y = pathway, fill = effect), colour = 'black')+
  scale_fill_gradient2(low = "#C71000FF", high = '#8A4198FF', mid = "white", midpoint = 0)+
  geom_text(aes(x = effect, y = pathway, label = wi.eBH.signif, group = pathway),
            size = 6, position = position_dodge2(width = 1, reverse = T))+
  # scale_fill_viridis_b(option = 'plasma')+
  theme_clean_rgd_grid()+
  annotate(geom = 'text', x = 2, y = middle, label = 'Week 12', size = 4, family = 'Times', angle = 270)+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
  )+
  labs(x = 'Estimated ALDEx2 Effect Size', y = 'Pathway')+
  scale_x_continuous(breaks = c(-2,-1,0,1,2), limits = c(-2,2))+
  # xlim(c(-1.5,1.5))+
  guides(fill = 'none')+
  
  annotate(geom = 'text', x = 1.4, y = top, label = '\u2192\nHigher in\nLOWCHO',
           family = 'Times', size = 3)+
  annotate(geom = 'text', x = -1.4, y = bottom, label = 'Higher in\nLOWSUG\n\u2190',
           family = 'Times', size = 3)

## Arrange h2_LOWSUG_LOWCHO -------------------------------------------------
ggarrange(h2_LOWSUG_LOWCHO_A_plot, h2_LOWSUG_LOWCHO_B_plot, h2_LOWSUG_LOWCHO_C_plot,
          ncol = 1, nrow = 3, align = 'v',
          labels = 'AUTO', font.label = list(family = 'Times'), vjust = 1)

ggsave(filename = './Figures/chebi_h2_LOWSUG_LOWCHO_ABC_arranged.tiff',
       device = 'tiff',
       width = 15,
       height = 24,
       units = 'cm',
       dpi = 600)
