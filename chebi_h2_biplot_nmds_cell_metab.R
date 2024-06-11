# # Packages ####
# library(tidyverse)
# library(vegan)
# library(zCompositions)
# library(ggrepel)
# library(beepr)
# set.seed(123)
# 
# # OLD ####
# data_1 <- read_delim('./Data/kb_vegan_5_perm_2998_id')
# samples <- read.csv2('./Data/meta_chebi.csv')
# vegan <- data_1 %>% column_to_rownames(var = "id_2")
# 
# # CLR Transform
# vegan_no0 <- cmultRepl(vegan, method = "CZM", output = "p-counts")
# 
# vegan_clr <- t(apply(vegan_no0, 1, function(x){log(x) - mean(log(x))})) %>% as.data.frame()
# 
# # Generate distance matrix
# vegan_clr_dist <- vegdist(vegan_clr, method = "euclidean") %>% as.dist() # Aitchison distance matrix
# 
# starting_config <- cmdscale(vegan_clr_dist)
# 
# set.seed(2)
# 
# # Basic biplot
# scores(vegan_clr_dist) %>% 
#   as_tibble(rownames = "id") %>% 
#   inner_join(metadata, by = "id") %>%
#   filter(group == "Intervention") %>% 
#   
#   ggplot(aes(x = NMDS1,
#              y = NMDS2,
#              colour = group_trial)) +
#   # geom_point() +
#   geom_text(aes(label = id), size = 3, check_overlap = F) +
#   stat_ellipse() +
#   theme_classic()
# 


# UPDATED Robust Aitchison ------------------------------------------------

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
library(ggrepel)
library(tidyverse)
library(scales)
windowsFonts("Times" = windowsFont("Times"))
source('./Scripts/theme_clean_rgd_cell_metab.R')

# Between trials ----------------------------------------------------------


# Data ####
source('./Scripts/chebi_h2_load_data.R')
samples$combined <- paste0(samples$group,sep = "-", samples$trial)
samples$id_2
data_1$id_2
rownames(vegan)
colours <- c("#101E9D","#46A6CE","#EEA327")

# LOWCHO --------------------------------------------------------------------

samples_filtered <- samples %>% 
  filter(group == 'LOWCHO') %>% # Edit group
  # filter(trial == 'B'| trial == 'C') %>% # Edit trials
  # filter(participant != '14', participant != '8', participant != '21', participant != '56') %>% # Filter only one trial
  arrange(id_2)
samples_filtered$group
samples_filtered$id_2
data_filtered <- data_1 %>% 
  semi_join(samples_filtered, by = 'id_2') %>% 
  arrange(id_2) %>% 
  column_to_rownames(var = 'id_2') 
rownames(data_filtered)


## Run NMDS ----------------------------------------------------------------

set.seed(1)
vegan_nmds <- metaMDS(data_filtered, distance = 'robust.aitchison', trymax = 100)

scores <- scores(vegan_nmds)
stress <- vegan_nmds$stress

nmds_samples <- scores$sites %>% 
  as_tibble(rownames = "id_2") %>% 
  inner_join(samples, by = "id_2")

centroid <- nmds_samples %>% 
  group_by(combined) %>% 
  summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))


## Plot --------------------------------------------------------------------
lowcho_plot <- 
nmds_samples %>%
  # filter(group == "Control") %>%
  
  ggplot(aes(x = NMDS1,
             y = NMDS2,
             colour = trial)) +
  geom_point()+
  # geom_text(aes(label = id_2),family = 'Helvetica', size = 2.5, check_overlap = F) +
  stat_ellipse(level = 0.95) +
  geom_point(data = centroid, size = 5, shape = 21, colour = 'black', aes(fill = combined)) +
  geom_text(aes(x = 6, y = 8, label = 'Stress = 0.193'),
            family = 'Times', size = 3,inherit.aes = F)+
  
  scale_fill_futurama(labels = c('Baseline', 'Week 4', 'Week 12'))+
  scale_colour_futurama()+
  scale_x_continuous(limits = c(-8,8), breaks = c(-8,-6,-4,-2,0,2,4,6,8))+
  scale_y_continuous(limits = c(-8,8), breaks = c(-8,-6,-4,-2,0,2,4,6,8))+
  guides(colour = F)+
  labs(colour = NULL, fill = 'Trial', family = 'Times')+
  # ggtitle('LOWCHO')+
  theme_clean_rgd_no_grid()+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        legend.text = element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'),
        legend.position = 'top',
        legend.background = element_blank())

# ggsave(filename = './Figures/chebi_h2_nmds_LOWCHO.tiff',
#        device = 'tiff',
#        width = 15,
#        height = 8,
#        units = 'cm',
#        dpi = 600)

# MODSUG --------------------------------------------------------------------

samples$id_2
rownames(vegan)

samples_filtered <- samples %>% 
  filter(group == 'MODSUG') %>% # Edit group
  # filter(trial == 'B'| trial == 'C') %>% # Edit trials
  # filter(participant != '14', participant != '8', participant != '21', participant != '56') %>% # Filter only one trial
  arrange(id_2)
samples_filtered$group
samples_filtered$id_2
data_filtered <- data_1 %>% 
  semi_join(samples_filtered, by = 'id_2') %>% 
  arrange(id_2) %>% 
  column_to_rownames(var = 'id_2') 
rownames(data_filtered)


## Run NMDS ----------------------------------------------------------------

set.seed(1)
vegan_nmds <- metaMDS(data_filtered, distance = 'robust.aitchison', trymax = 100)

scores <- scores(vegan_nmds)
stress <- vegan_nmds$stress

nmds_samples <- scores$sites %>% 
  as_tibble(rownames = "id_2") %>% 
  inner_join(samples, by = "id_2")

centroid <- nmds_samples %>% 
  group_by(combined) %>% 
  summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))


## Plot --------------------------------------------------------------------

modsug_plot <- 
nmds_samples %>%
  # filter(group == "Control") %>%
  
  ggplot(aes(x = NMDS1,
             y = NMDS2,
             colour = trial)) +
  theme_clean_rgd_no_grid()+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        legend.text = element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'), 
        legend.position = 'none',
        legend.background = element_blank())+
  # geom_text(aes(label = id_2),family = 'Times', size = 2.5, check_overlap = F) +
  geom_point()+
  
  stat_ellipse(level = 0.95) +
  geom_point(data = centroid, size = 5, shape = 21, colour = 'black', aes(fill = combined)) +
  geom_text(aes(x = 6, y = 8, label = 'Stress = 0.140'),
            family = 'Times', size = 3,inherit.aes = F)+
  
  scale_fill_futurama(labels = c('Baseline', 'Week 4', 'Week 12'))+
  scale_colour_futurama()+
  scale_x_continuous(limits = c(-8,8), breaks = c(-8,-6,-4,-2,0,2,4,6,8))+
  scale_y_continuous(limits = c(-8,8), breaks = c(-8,-6,-4,-2,0,2,4,6,8))+
  guides(colour = F)+
  labs(colour = NULL, fill = 'Trial', family = 'Times')

# ggsave(filename = './Figures/chebi_h2_nmds_MODSUG.tiff',
#        device = 'tiff',
#        width = 15,
#        height = 8,
#        units = 'cm',
#        dpi = 600)

# LOWSUG --------------------------------------------------------------------

samples$id_2
rownames(vegan)

samples_filtered <- samples %>% 
  filter(group == 'LOWSUG') %>% # Edit group
  # filter(trial == 'B'| trial == 'C') %>% # Edit trials
  # filter(participant != '14', participant != '8', participant != '21', participant != '56') %>% # Filter only one trial
  arrange(id_2)
samples_filtered$group
samples_filtered$id_2
data_filtered <- data_1 %>% 
  semi_join(samples_filtered, by = 'id_2') %>% 
  arrange(id_2) %>% 
  column_to_rownames(var = 'id_2') 
rownames(data_filtered)


## Run NMDS ----------------------------------------------------------------

set.seed(1)
vegan_nmds <- metaMDS(data_filtered, distance = 'robust.aitchison', trymax = 100)

scores <- scores(vegan_nmds)
stress <- vegan_nmds$stress

nmds_samples <- scores$sites %>% 
  as_tibble(rownames = "id_2") %>% 
  inner_join(samples, by = "id_2")

centroid <- nmds_samples %>% 
  group_by(combined) %>% 
  summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))


## Plot --------------------------------------------------------------------

lowsug_plot <-
nmds_samples %>%
  # filter(group == "Control") %>%
  
  ggplot(aes(x = NMDS1,
             y = NMDS2,
             colour = trial)) +
  theme_clean_rgd_no_grid()+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        legend.text = element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'), 
        legend.position = 'none',
        legend.background = element_blank())+
  # geom_text(aes(label = id_2),family = 'Times', size = 2.5, check_overlap = F) +
  geom_point()+
  
  # geom_point(aes(fill = combined), shape = 21, colour = 'black')+
  stat_ellipse(level = 0.95) +
  geom_point(data = centroid, size = 5, shape = 21, colour = 'black', aes(fill = combined)) +
  geom_text(aes(x = 6, y = 8, label = 'Stress = 0.198'),
            family = 'Times', size = 3,inherit.aes = F)+
  
  scale_fill_futurama(labels = c('Baseline', 'Week 4', 'Week 12'))+
  scale_colour_futurama()+
  scale_x_continuous(limits = c(-8,8), breaks = c(-8,-6,-4,-2,0,2,4,6,8))+
  scale_y_continuous(limits = c(-8,8), breaks = c(-8,-6,-4,-2,0,2,4,6,8))+
  guides(colour = F)+
  labs(colour = NULL, fill = 'Trial')

# ggsave(filename = './Figures/chebi_h2_nmds_LOWSUG.tiff',
#        device = 'tiff',
#        width = 15,
#        height = 8,
#        units = 'cm',
#        dpi = 600)


# Arrange H2 NMDS Groups --------------------------------------------------

ggarrange(modsug_plot, lowsug_plot, lowcho_plot, nrow = 1, ncol = 3, common.legend = T,
          labels = list('MODSUG','LOWSUG','LOWCHO'),
          font.label = list(family = 'Times', size = 10), 
          hjust = -0.8)

ggsave(filename = './Figures/chebi_h2_nmds_groups_arranged_landscape_stress_point.tiff',
       device = 'tiff',
       width = 25.7,
       height = 13.5,
       units = 'cm',
       dpi = 600)

# Baseline --------------------------------------------------------------------

samples$id_2
rownames(vegan)

samples_filtered <- samples %>% 
  # filter(group == 'LOWCHO') %>% # Edit group
  filter(trial == 'A') %>% # Edit trials
  # filter(participant != '14', participant != '8', participant != '21', participant != '56') %>% # Filter only one trial
  arrange(id_2)
samples_filtered$group
samples_filtered$id_2
data_filtered <- data_1 %>% 
  semi_join(samples_filtered, by = 'id_2') %>% 
  arrange(id_2) %>% 
  column_to_rownames(var = 'id_2') 
rownames(data_filtered)


## Run NMDS ----------------------------------------------------------------

set.seed(1)
vegan_nmds <- metaMDS(data_filtered, distance = 'robust.aitchison', trymax = 100)

scores <- scores(vegan_nmds)
stress <- vegan_nmds$stress

nmds_samples <- scores$sites %>% 
  as_tibble(rownames = "id_2") %>% 
  inner_join(samples, by = "id_2")

centroid <- nmds_samples %>% 
  group_by(combined) %>% 
  summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2)) %>% 
  mutate(combined = factor(combined, levels = c('MODSUG-A', 'LOWSUG-A', 'LOWCHO-A')))


## Plot --------------------------------------------------------------------

baseline_plot <- 
nmds_samples %>%
  mutate(group = factor(group, levels = c('MODSUG', 'LOWSUG', 'LOWCHO'))) %>% 
  mutate(combined = factor(combined, levels = c('MODSUG-A', 'LOWSUG-A', 'LOWCHO-A'))) %>% 
  
  # filter(group == "Control") %>%
  
  ggplot(aes(x = NMDS1,
             y = NMDS2,
             colour = group)) +
  theme_clean_rgd_no_grid()+
  theme(text = element_text(family = 'Helvetica',
                            face = 'bold',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        legend.text = element_text(family = 'Helvetica'),
        legend.title = element_text(family = 'Helvetica'),
        legend.position = 'top',
        legend.background = element_blank())+
  # geom_text(aes(label = id_2),family = 'Helvetica', size = 2.5, check_overlap = F) +
  geom_point(size = 0.75)+
  
  stat_ellipse(level = 0.95) +
  geom_point(data = centroid, size = 3, shape = 21, colour = 'black', aes(fill = combined)) +
  geom_text(aes(x = 6, y = 8, label = 'Stress = 0.185'),
            family = 'Helvetica', size = 2,inherit.aes = F)+
  
  scale_fill_manual(values = colours, labels = c('MODSUG', 'LOWSUG','LOWCHO'))+
  scale_colour_manual(values = colours)+
  # scale_fill_futurama(labels = c('MODSUG', 'LOWSUG','LOWCHO'))+
  # scale_colour_futurama()+
  scale_x_continuous(limits = c(-8,8), breaks = c(-8,-6,-4,-2,0,2,4,6,8))+
  scale_y_continuous(limits = c(-8,8), breaks = c(-8,-6,-4,-2,0,2,4,6,8))+
  guides(colour = F)+
  labs(colour = NULL, fill = 'Group')

# ggsave(filename = './Figures/chebi_h2_nmds_A.tiff',
#        device = 'tiff',
#        width = 15,
#        height = 8,
#        units = 'cm',
#        dpi = 600)

# Week 4 --------------------------------------------------------------------

samples$id_2
rownames(vegan)

samples_filtered <- samples %>% 
  # filter(group == 'LOWCHO') %>% # Edit group
  filter(trial == 'B') %>% # Edit trials
  # filter(participant != '14', participant != '8', participant != '21', participant != '56') %>% # Filter only one trial
  arrange(id_2)
samples_filtered$group
samples_filtered$id_2
data_filtered <- data_1 %>% 
  semi_join(samples_filtered, by = 'id_2') %>% 
  arrange(id_2) %>% 
  column_to_rownames(var = 'id_2') 
rownames(data_filtered)


## Run NMDS ----------------------------------------------------------------

set.seed(1)
vegan_nmds <- metaMDS(data_filtered, distance = 'robust.aitchison', trymax = 100)

scores <- scores(vegan_nmds)
stress <- vegan_nmds$stress

nmds_samples <- scores$sites %>% 
  as_tibble(rownames = "id_2") %>% 
  inner_join(samples, by = "id_2")

centroid <- nmds_samples %>% 
  group_by(combined) %>% 
  summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2)) %>% 
  mutate(combined = factor(combined, levels = c('MODSUG-B', 'LOWSUG-B', 'LOWCHO-B')))%>% 
  mutate(sig = c('+','',''))


## Plot --------------------------------------------------------------------

week_4_plot <- 
nmds_samples %>%
  mutate(group = factor(group, levels = c('MODSUG', 'LOWSUG', 'LOWCHO')))%>%
  mutate(combined = factor(combined, levels = c('MODSUG-B', 'LOWSUG-B', 'LOWCHO-B'))) %>% 
  
  # filter(group == "Control") %>%
  
  ggplot(aes(x = NMDS1,
             y = NMDS2,
             colour = group)) +
  theme_clean_rgd_no_grid()+
  theme(text = element_text(family = 'Helvetica',
                            face = 'bold',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        legend.text = element_text(family = 'Helvetica'),
        legend.title = element_text(family = 'Helvetica'),
        legend.background = element_blank())+
  # geom_text(aes(label = id_2),family = 'Helvetica', size = 2.5, check_overlap = F) +
  geom_point(size = 0.75)+
  
  stat_ellipse(level = 0.95) +
  geom_point(data = centroid, size = 3, shape = 21, colour = 'black', aes(fill = combined)) +
  geom_text(aes(x = 6, y = 8, label = 'Stress = 0.156'),
            family = 'Helvetica', size = 2,inherit.aes = F)+
  geom_text(aes(x = NMDS1, y = NMDS2, label = sig),
            data = centroid,family = 'Helvetica', size = 2, inherit.aes = F)+
  
  
  scale_fill_manual(values = colours)+
  scale_colour_manual(values = colours)+
  scale_x_continuous(limits = c(-8,8), breaks = c(-8,-6,-4,-2,0,2,4,6,8))+
  scale_y_continuous(limits = c(-8,8), breaks = c(-8,-6,-4,-2,0,2,4,6,8))+
  guides(colour = F)+
  labs(colour = NULL, fill = 'Group')

# ggsave(filename = './Figures/chebi_h2_nmds_B.tiff',
#        device = 'tiff',
#        width = 15,
#        height = 8,
#        units = 'cm',
#        dpi = 600)

# Week 12 --------------------------------------------------------------------

samples$id_2
rownames(vegan)

samples_filtered <- samples %>% 
  # filter(group == 'LOWCHO') %>% # Edit group
  filter(trial == 'C') %>% # Edit trials
  # filter(participant != '14', participant != '8', participant != '21', participant != '56') %>% # Filter only one trial
  arrange(id_2)
samples_filtered$group
samples_filtered$id_2
data_filtered <- data_1 %>% 
  semi_join(samples_filtered, by = 'id_2') %>% 
  arrange(id_2) %>% 
  column_to_rownames(var = 'id_2') 
rownames(data_filtered)


## Run NMDS ----------------------------------------------------------------

set.seed(1)
vegan_nmds <- metaMDS(data_filtered, distance = 'robust.aitchison', trymax = 100)

scores <- scores(vegan_nmds)
stress <- vegan_nmds$stress

nmds_samples <- scores$sites %>% 
  as_tibble(rownames = "id_2") %>% 
  inner_join(samples, by = "id_2")

centroid <- nmds_samples %>% 
  group_by(combined) %>% 
  summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2)) %>% 
  mutate(combined = factor(combined, levels = c('MODSUG-C', 'LOWSUG-C', 'LOWCHO-C')))%>% 
  mutate(sig = c('#','',''))


## Plot --------------------------------------------------------------------

week_12_plot <-
nmds_samples %>%
  mutate(group = factor(group, levels = c('MODSUG', 'LOWSUG', 'LOWCHO')))%>%
  mutate(combined = factor(combined, levels = c('MODSUG-C', 'LOWSUG-C', 'LOWCHO-C'))) %>% 
  
  ggplot(aes(x = NMDS1,
             y = NMDS2,
             colour = group)) +
  theme_clean_rgd_no_grid()+
  theme(text = element_text(family = 'Helvetica',
                            face = 'bold',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        legend.text = element_text(family = 'Helvetica'),
        legend.title = element_text(family = 'Helvetica'),
        legend.background = element_blank())+
  # geom_text(aes(label = id_2),family = 'Helvetica', size = 2.5, check_overlap = F) +
  geom_point(size = 0.75)+
  stat_ellipse(level = 0.95) +
  geom_point(data = centroid, size = 3, shape = 21, colour = 'black', aes(fill = combined)) +
  geom_text(aes(x = 6, y = 8, label = 'Stress = 0.178'),
            family = 'Helvetica', size = 2,inherit.aes = F)+
  geom_text(aes(x = NMDS1, y = NMDS2, label = sig),
            data = centroid,family = 'Helvetica', size = 2, inherit.aes = F)+
  
  
  scale_fill_manual(values = colours)+
  scale_colour_manual(values = colours)+
  scale_x_continuous(limits = c(-8,8), breaks = c(-8,-6,-4,-2,0,2,4,6,8))+
  scale_y_continuous(limits = c(-8,8), breaks = c(-8,-6,-4,-2,0,2,4,6,8))+
  guides(colour = F)+
  labs(colour = NULL, fill = 'Group')

# ggsave(filename = './Figures/chebi_h2_nmds_C.tiff',
#        device = 'tiff',
#        width = 16,
#        height = 8,
#        units = 'cm',
#        dpi = 600)



# Arrange H2 NMDS Trials --------------------------------------------------

ggarrange(baseline_plot, week_4_plot, week_12_plot, nrow = 1, ncol = 3, common.legend = T, align = 'hv',
          labels = list('Baseline','Week 4','Week 12'),
          font.label = list(family = 'Helvetica', size = 10), 
          hjust = -1.05)

ggsave(filename = './Figures/chebi_h2_nmds_trials_arranged_landscape_stress_sig_cell_metab.tiff',
       device = 'tiff',
       width = 25.7,
       height = 13.5,
       units = 'cm',
       dpi = 600)

h2_baseline_plot <- baseline_plot
h2_week_4_plot <- week_4_plot
h2_week_12_plot <- week_12_plot

# With Pathways NOT USING ------------------------------------------------------------
# LOWCHO --------------------------------------------------------------------

samples$id_2
rownames(vegan)

samples_filtered <- samples %>% 
  filter(group == 'LOWCHO') %>% # Edit group
  # filter(trial == 'B'| trial == 'C') %>% # Edit trials
  # filter(participant != '14', participant != '8', participant != '21', participant != '56') %>% # Filter only one trial
  arrange(id_2)
samples_filtered$group
samples_filtered$id_2
data_filtered <- data_1 %>% 
  semi_join(samples_filtered, by = 'id_2') %>% 
  arrange(id_2) %>% 
  column_to_rownames(var = 'id_2') 
rownames(data_filtered)


## Run NMDS ----------------------------------------------------------------

set.seed(1)
vegan_nmds <- metaMDS(data_filtered, distance = 'robust.aitchison', trymax = 100)

scores <- scores(vegan_nmds)
scores$sites

nmds_samples <- scores$sites %>% 
  as_tibble(rownames = "id_2") %>% 
  inner_join(samples, by = "id_2")

centroid <- nmds_samples %>% 
  group_by(combined) %>% 
  summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))




## Add Species to Make Biplots ---------------------------------------------

vegan_no0 <- cmultRepl(data_filtered, method = "CZM", output = "p-counts")
vegan_clr <- t(apply(vegan_no0, 1, function(x)  {
  
  log(x) - mean(log(x))
  
}
)) %>% as.data.frame()


shared_tbl <- vegan_clr %>% 
  as_tibble(rownames = 'id_2') %>% 
  pivot_longer(-id_2)

scores <- scores(vegan_nmds)

shared_tbl %>% count(id_2)
rownames(vegan_clr)
scores$sites

nmds_positions <- scores$sites %>% 
  as_tibble(rownames = 'id_2')

nmds_shared <- inner_join(shared_tbl, nmds_positions, by = 'id_2')




cor_x <- nmds_shared %>% 
  nest(data = -name) %>% 
  mutate(cor_x = map(data,
                     ~cor.test(.x$value,
                               .x$NMDS1, 
                               method = 'spearman',
                               exact = F) %>% tidy())) %>% 
  unnest(cor_x) %>% 
  dplyr::select(name, estimate, p.value)


cor_y <- nmds_shared %>% 
  nest(data = -name) %>% 
  mutate(cor_y = map(data,
                     ~ cor.test(.x$value, .x$NMDS2, 
                                method = 'spearman',
                                exact = F) %>% tidy())) %>% 
  unnest(cor_y) %>% 
  dplyr::select(name, estimate, p.value)


correlations <- inner_join(cor_x, cor_y, by = 'name')

high_cor <- correlations %>% 
  filter(estimate.x > 0.75 | 
           estimate.y > 0.65 |
           estimate.x < -0.5 |
           estimate.y < -0.7) %>% 
  mutate(name = str_replace_all(name, '--.*', ''))
  


## Plot --------------------------------------------------------------------

nmds_samples %>%
  ggplot(aes(x = NMDS1,
             y = NMDS2,
             colour = trial)) +
  
  # geom_text(aes(label = id_2),family = 'Times', size = 3, check_overlap = F) +
  geom_point(aes(fill = combined), size = 1.5) +
  stat_ellipse(level = 0.95) +
  geom_point(data = centroid, size = 5, shape = 21, colour = 'black', aes(fill = combined)) +
  geom_segment(data = high_cor,
               aes(x = 0, xend = estimate.x*20, y = 0, yend = estimate.y*20),
               inherit.aes = F)+
  geom_text_repel(data = high_cor, family = 'Times', max.overlaps = 30, direction = 'both',
                  segment.colour = NA, point.padding = 20,
                  aes(x = estimate.x*20, y = estimate.y*20, label = name),
                  size = 2.5, inherit.aes = F)+
  scale_fill_futurama(labels = c('Baseline', 'Week 4', 'Week 12'))+
  scale_colour_futurama()+
  xlim(-20,20) +
  ylim(-20,20)+
  guides(colour = F)+
  labs(colour = NULL, fill = 'Trial', family = 'Times')+
  ggtitle('LOWCHO')+
  theme_clean_rgd_no_grid()+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        legend.text = element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'))

ggsave(filename = './Figures/chebi_nmds_biplot_LOWCHO.tiff',
       device = 'tiff',
       width = 16,
       height = 8,
       units = 'cm',
       dpi = 600)







# MODSUG --------------------------------------------------------------------

samples$id_2
rownames(vegan)

samples_filtered <- samples %>% 
  filter(group == 'MODSUG') %>% # Edit group
  # filter(trial == 'B'| trial == 'C') %>% # Edit trials
  # filter(participant != '14', participant != '8', participant != '21', participant != '56') %>% # Filter only one trial
  arrange(id_2)
samples_filtered$group
samples_filtered$id_2
data_filtered <- data_1 %>% 
  semi_join(samples_filtered, by = 'id_2') %>% 
  arrange(id_2) %>% 
  column_to_rownames(var = 'id_2') 
rownames(data_filtered)


## Run NMDS ----------------------------------------------------------------

set.seed(1)
vegan_nmds <- metaMDS(data_filtered, distance = 'robust.aitchison', trymax = 100)

scores <- scores(vegan_nmds)
scores$sites

nmds_samples <- scores$sites %>% 
  as_tibble(rownames = "id_2") %>% 
  inner_join(samples, by = "id_2")

centroid <- nmds_samples %>% 
  group_by(combined) %>% 
  summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))




## Add Species to Make Biplots ---------------------------------------------

vegan_no0 <- cmultRepl(data_filtered, method = "CZM", output = "p-counts")
vegan_clr <- t(apply(vegan_no0, 1, function(x)  {
  
  log(x) - mean(log(x))
  
}
)) %>% as.data.frame()


shared_tbl <- vegan_clr %>% 
  as_tibble(rownames = 'id_2') %>% 
  pivot_longer(-id_2)

scores <- scores(vegan_nmds)

shared_tbl %>% count(id_2)
rownames(vegan_clr)
scores$sites

nmds_positions <- scores$sites %>% 
  as_tibble(rownames = 'id_2')

nmds_shared <- inner_join(shared_tbl, nmds_positions, by = 'id_2')




cor_x <- nmds_shared %>% 
  nest(data = -name) %>% 
  mutate(cor_x = map(data,
                     ~cor.test(.x$value,
                               .x$NMDS1, 
                               method = 'spearman',
                               exact = F) %>% tidy())) %>% 
  unnest(cor_x) %>% 
  dplyr::select(name, estimate, p.value)


cor_y <- nmds_shared %>% 
  nest(data = -name) %>% 
  mutate(cor_y = map(data,
                     ~ cor.test(.x$value, .x$NMDS2, 
                                method = 'spearman',
                                exact = F) %>% tidy())) %>% 
  unnest(cor_y) %>% 
  dplyr::select(name, estimate, p.value)


correlations <- inner_join(cor_x, cor_y, by = 'name')

high_cor <- correlations %>% 
  filter(estimate.x > 0.85 | 
           estimate.y > 0.85 |
           estimate.x < -0.65 |
           estimate.y < -0.65) %>% 
  mutate(name = str_replace_all(name, '--.*', ''))



## Plot --------------------------------------------------------------------

nmds_samples %>%
  ggplot(aes(x = NMDS1,
             y = NMDS2,
             colour = trial)) +
  
  # geom_text(aes(label = id_2),family = 'Times', size = 3, check_overlap = F) +
  geom_point(aes(fill = combined), size = 1.5) +
  stat_ellipse(level = 0.95) +
  geom_point(data = centroid, size = 5, shape = 21, colour = 'black', aes(fill = combined)) +
  geom_segment(data = high_cor,
               aes(x = 0, xend = estimate.x*20, y = 0, yend = estimate.y*20),
               inherit.aes = F)+
  geom_text_repel(data = high_cor, family = 'Times', max.overlaps = 30, direction = 'both',
                  segment.colour = NA, point.padding = 2,
                  aes(x = estimate.x*20, y = estimate.y*20, label = name),
                  size = 2.5, inherit.aes = F)+
  scale_fill_futurama(labels = c('Baseline', 'Week 4', 'Week 12'))+
  scale_colour_futurama()+
  xlim(-20,20) +
  ylim(-20,21)+
  guides(colour = F)+
  labs(colour = NULL, fill = 'Trial', family = 'Times')+
  ggtitle('MODSUG')+
  theme_clean_rgd_no_grid()+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        legend.text = element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'))

ggsave(filename = './Figures/chebi_nmds_biplot_MODSUG.tiff',
       device = 'tiff',
       width = 16,
       height = 8,
       units = 'cm',
       dpi = 600)







# LOWSUG --------------------------------------------------------------------

samples$id_2
rownames(vegan)

samples_filtered <- samples %>% 
  filter(group == 'LOWSUG') %>% # Edit group
  # filter(trial == 'B'| trial == 'C') %>% # Edit trials
  # filter(participant != '14', participant != '8', participant != '21', participant != '56') %>% # Filter only one trial
  arrange(id_2)
samples_filtered$group
samples_filtered$id_2
data_filtered <- data_1 %>% 
  semi_join(samples_filtered, by = 'id_2') %>% 
  arrange(id_2) %>% 
  column_to_rownames(var = 'id_2') 
rownames(data_filtered)


## Run NMDS ----------------------------------------------------------------

set.seed(1)
vegan_nmds <- metaMDS(data_filtered, distance = 'robust.aitchison', trymax = 100)

scores <- scores(vegan_nmds)
scores$sites

nmds_samples <- scores$sites %>% 
  as_tibble(rownames = "id_2") %>% 
  inner_join(samples, by = "id_2")

centroid <- nmds_samples %>% 
  group_by(combined) %>% 
  summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))




## Add Species to Make Biplots ---------------------------------------------

vegan_no0 <- cmultRepl(data_filtered, method = "CZM", output = "p-counts")
vegan_clr <- t(apply(vegan_no0, 1, function(x)  {
  
  log(x) - mean(log(x))
  
}
)) %>% as.data.frame()


shared_tbl <- vegan_clr %>% 
  as_tibble(rownames = 'id_2') %>% 
  pivot_longer(-id_2)

scores <- scores(vegan_nmds)

shared_tbl %>% count(id_2)
rownames(vegan_clr)
scores$sites

nmds_positions <- scores$sites %>% 
  as_tibble(rownames = 'id_2')

nmds_shared <- inner_join(shared_tbl, nmds_positions, by = 'id_2')




cor_x <- nmds_shared %>% 
  nest(data = -name) %>% 
  mutate(cor_x = map(data,
                     ~cor.test(.x$value,
                               .x$NMDS1, 
                               method = 'spearman',
                               exact = F) %>% tidy())) %>% 
  unnest(cor_x) %>% 
  dplyr::select(name, estimate, p.value)


cor_y <- nmds_shared %>% 
  nest(data = -name) %>% 
  mutate(cor_y = map(data,
                     ~ cor.test(.x$value, .x$NMDS2, 
                                method = 'spearman',
                                exact = F) %>% tidy())) %>% 
  unnest(cor_y) %>% 
  dplyr::select(name, estimate, p.value)


correlations <- inner_join(cor_x, cor_y, by = 'name')

high_cor <- correlations %>% 
  filter(estimate.x > 0.65 | 
           estimate.y > 0.6 |
           estimate.x < -0.85 |
           estimate.y < -0.8) %>% 
  mutate(name = str_replace_all(name, '--.*', ''))



## Plot --------------------------------------------------------------------

nmds_samples %>%
  ggplot(aes(x = NMDS1,
             y = NMDS2,
             colour = trial)) +
  
  # geom_text(aes(label = id_2),family = 'Times', size = 3, check_overlap = F) +
  geom_point(aes(fill = combined), size = 1.5) +
  stat_ellipse(level = 0.95) +
  geom_point(data = centroid, size = 5, shape = 21, colour = 'black', aes(fill = combined)) +
  geom_segment(data = high_cor,
               aes(x = 0, xend = estimate.x*20, y = 0, yend = estimate.y*20),
               inherit.aes = F)+
  geom_text_repel(data = high_cor, family = 'Times', max.overlaps = 30, direction = 'both',
                  segment.colour = NA, point.padding = 2,
                  aes(x = estimate.x*20, y = estimate.y*20, label = name),
                  size = 2.5, inherit.aes = F)+
  scale_fill_futurama(labels = c('Baseline', 'Week 4', 'Week 12'))+
  scale_colour_futurama()+
  xlim(-20,20) +
  ylim(-20,20)+
  guides(colour = F)+
  labs(colour = NULL, fill = 'Trial', family = 'Times')+
  ggtitle('LOWSUG')+
  theme_clean_rgd_no_grid()+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        legend.text = element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'))

ggsave(filename = './Figures/chebi_nmds_biplot_LOWSUG.tiff',
       device = 'tiff',
       width = 16,
       height = 8,
       units = 'cm',
       dpi = 600)







# Baseline --------------------------------------------------------------------

samples$id_2
rownames(vegan)

samples_filtered <- samples %>% 
  # filter(group == 'LOWSUG') %>% # Edit group
  filter(trial == 'A') %>% # Edit trials
  # filter(participant != '14', participant != '8', participant != '21', participant != '56') %>% # Filter only one trial
  arrange(id_2)
samples_filtered$group
samples_filtered$id_2
data_filtered <- data_1 %>% 
  semi_join(samples_filtered, by = 'id_2') %>% 
  arrange(id_2) %>% 
  column_to_rownames(var = 'id_2') 
rownames(data_filtered)


## Run NMDS ----------------------------------------------------------------

set.seed(1)
vegan_nmds <- metaMDS(data_filtered, distance = 'robust.aitchison', trymax = 100)

scores <- scores(vegan_nmds)
scores$sites

nmds_samples <- scores$sites %>% 
  as_tibble(rownames = "id_2") %>% 
  inner_join(samples, by = "id_2")

centroid <- nmds_samples %>% 
  group_by(combined) %>% 
  summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))




## Add Species to Make Biplots ---------------------------------------------

vegan_no0 <- cmultRepl(data_filtered, method = "CZM", output = "p-counts")
vegan_clr <- t(apply(vegan_no0, 1, function(x)  {
  
  log(x) - mean(log(x))
  
}
)) %>% as.data.frame()


shared_tbl <- vegan_clr %>% 
  as_tibble(rownames = 'id_2') %>% 
  pivot_longer(-id_2)

scores <- scores(vegan_nmds)

shared_tbl %>% count(id_2)
rownames(vegan_clr)
scores$sites

nmds_positions <- scores$sites %>% 
  as_tibble(rownames = 'id_2')

nmds_shared <- inner_join(shared_tbl, nmds_positions, by = 'id_2')

cor_x <- nmds_shared %>% 
  nest(data = -name) %>% 
  mutate(cor_x = map(data,
                     ~cor.test(.x$value,
                               .x$NMDS1, 
                               method = 'spearman',
                               exact = F) %>% tidy())) %>% 
  unnest(cor_x) %>% 
  dplyr::select(name, estimate, p.value)


cor_y <- nmds_shared %>% 
  nest(data = -name) %>% 
  mutate(cor_y = map(data,
                     ~ cor.test(.x$value, .x$NMDS2, 
                                method = 'spearman',
                                exact = F) %>% tidy())) %>% 
  unnest(cor_y) %>% 
  dplyr::select(name, estimate, p.value)


correlations <- inner_join(cor_x, cor_y, by = 'name')

high_cor <- correlations %>% 
  filter(estimate.x > 0.75 | 
           estimate.y > 0.5 |
           estimate.x < -0.6 |
           estimate.y < -0.55) %>% 
  mutate(name = str_replace_all(name, '--.*', ''))



## Plot --------------------------------------------------------------------

nmds_samples %>%
  ggplot(aes(x = NMDS1,
             y = NMDS2,
             colour = group)) +
  
  # geom_text(aes(label = id_2),family = 'Times', size = 3, check_overlap = F) +
  geom_point(aes(fill = combined), size = 1.5) +
  stat_ellipse(level = 0.95) +
  geom_point(data = centroid, size = 5, shape = 21, colour = 'black', aes(fill = combined)) +
  geom_segment(data = high_cor,
               aes(x = 0, xend = estimate.x*20, y = 0, yend = estimate.y*20),
               inherit.aes = F)+
  geom_text_repel(data = high_cor, family = 'Times', max.overlaps = 30, direction = 'both',
                  segment.colour = NA, point.padding = 2,
                  aes(x = estimate.x*20, y = estimate.y*20, label = name),
                  size = 2.5, inherit.aes = F)+
  scale_fill_futurama(labels = c('LOWCHO', 'LOWSUG', 'MODSUG'))+
  scale_colour_futurama()+
  xlim(-20,20) +
  ylim(-20,20)+
  guides(colour = F)+
  labs(colour = NULL, fill = 'Group', family = 'Times')+
  ggtitle('Baseline')+
  theme_clean_rgd_no_grid()+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        legend.text = element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'))

ggsave(filename = './Figures/chebi_nmds_biplot_A.tiff',
       device = 'tiff',
       width = 16,
       height = 8,
       units = 'cm',
       dpi = 600)








# Week 4 --------------------------------------------------------------------

samples$id_2
rownames(vegan)

samples_filtered <- samples %>% 
  # filter(group == 'LOWSUG') %>% # Edit group
  filter(trial == 'B') %>% # Edit trials
  # filter(participant != '14', participant != '8', participant != '21', participant != '56') %>% # Filter only one trial
  arrange(id_2)
samples_filtered$group
samples_filtered$id_2
data_filtered <- data_1 %>% 
  semi_join(samples_filtered, by = 'id_2') %>% 
  arrange(id_2) %>% 
  column_to_rownames(var = 'id_2') 
rownames(data_filtered)


## Run NMDS ----------------------------------------------------------------

set.seed(1)
vegan_nmds <- metaMDS(data_filtered, distance = 'robust.aitchison', trymax = 100)

scores <- scores(vegan_nmds)
scores$sites

nmds_samples <- scores$sites %>% 
  as_tibble(rownames = "id_2") %>% 
  inner_join(samples, by = "id_2")

centroid <- nmds_samples %>% 
  group_by(combined) %>% 
  summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))




## Add Species to Make Biplots ---------------------------------------------

vegan_no0 <- cmultRepl(data_filtered, method = "CZM", output = "p-counts")
vegan_clr <- t(apply(vegan_no0, 1, function(x)  {
  
  log(x) - mean(log(x))
  
}
)) %>% as.data.frame()


shared_tbl <- vegan_clr %>% 
  as_tibble(rownames = 'id_2') %>% 
  pivot_longer(-id_2)

scores <- scores(vegan_nmds)

shared_tbl %>% count(id_2)
rownames(vegan_clr)
scores$sites

nmds_positions <- scores$sites %>% 
  as_tibble(rownames = 'id_2')

nmds_shared <- inner_join(shared_tbl, nmds_positions, by = 'id_2')

cor_x <- nmds_shared %>% 
  nest(data = -name) %>% 
  mutate(cor_x = map(data,
                     ~cor.test(.x$value,
                               .x$NMDS1, 
                               method = 'spearman',
                               exact = F) %>% tidy())) %>% 
  unnest(cor_x) %>% 
  dplyr::select(name, estimate, p.value)


cor_y <- nmds_shared %>% 
  nest(data = -name) %>% 
  mutate(cor_y = map(data,
                     ~ cor.test(.x$value, .x$NMDS2, 
                                method = 'spearman',
                                exact = F) %>% tidy())) %>% 
  unnest(cor_y) %>% 
  dplyr::select(name, estimate, p.value)


correlations <- inner_join(cor_x, cor_y, by = 'name')

high_cor <- correlations %>% 
  filter(estimate.x > 0.82 | 
           estimate.y > 0.6|
           estimate.x < -0.58 |
           estimate.y < -0.5) %>% 
  mutate(name = str_replace_all(name, '--.*', ''))



## Plot --------------------------------------------------------------------

nmds_samples %>%
  ggplot(aes(x = NMDS1,
             y = NMDS2,
             colour = group)) +
  
  # geom_text(aes(label = id_2),family = 'Times', size = 3, check_overlap = F) +
  geom_point(aes(fill = combined), size = 1.5) +
  stat_ellipse(level = 0.95) +
  geom_point(data = centroid, size = 5, shape = 21, colour = 'black', aes(fill = combined)) +
  geom_segment(data = high_cor,
               aes(x = 0, xend = estimate.x*20, y = 0, yend = estimate.y*20),
               inherit.aes = F)+
  geom_text_repel(data = high_cor, family = 'Times', max.overlaps = 30, direction = 'both',
                  segment.colour = NA, point.padding = 2,
                  aes(x = estimate.x*20, y = estimate.y*20, label = name),
                  size = 2.5, inherit.aes = F)+
  scale_fill_futurama(labels = c('LOWCHO', 'LOWSUG', 'MODSUG'))+
  scale_colour_futurama()+
  xlim(-20,20) +
  ylim(-20,20)+
  guides(colour = F)+
  labs(colour = NULL, fill = 'Group', family = 'Times')+
  ggtitle('Week 4')+
  theme_clean_rgd_no_grid()+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        legend.text = element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'))

ggsave(filename = './Figures/chebi_nmds_biplot_B.tiff',
       device = 'tiff',
       width = 16,
       height = 8,
       units = 'cm',
       dpi = 600)









# Week 12 --------------------------------------------------------------------

samples$id_2
rownames(vegan)

samples_filtered <- samples %>% 
  # filter(group == 'LOWSUG') %>% # Edit group
  filter(trial == 'C') %>% # Edit trials
  # filter(participant != '14', participant != '8', participant != '21', participant != '56') %>% # Filter only one trial
  arrange(id_2)
samples_filtered$group
samples_filtered$id_2
data_filtered <- data_1 %>% 
  semi_join(samples_filtered, by = 'id_2') %>% 
  arrange(id_2) %>% 
  column_to_rownames(var = 'id_2') 
rownames(data_filtered)


## Run NMDS ----------------------------------------------------------------

set.seed(1)
vegan_nmds <- metaMDS(data_filtered, distance = 'robust.aitchison', trymax = 100)

scores <- scores(vegan_nmds)
scores$sites

nmds_samples <- scores$sites %>% 
  as_tibble(rownames = "id_2") %>% 
  inner_join(samples, by = "id_2")

centroid <- nmds_samples %>% 
  group_by(combined) %>% 
  summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))




## Add Species to Make Biplots ---------------------------------------------

vegan_no0 <- cmultRepl(data_filtered, method = "CZM", output = "p-counts")
vegan_clr <- t(apply(vegan_no0, 1, function(x)  {
  
  log(x) - mean(log(x))
  
}
)) %>% as.data.frame()


shared_tbl <- vegan_clr %>% 
  as_tibble(rownames = 'id_2') %>% 
  pivot_longer(-id_2)

scores <- scores(vegan_nmds)

shared_tbl %>% count(id_2)
rownames(vegan_clr)
scores$sites

nmds_positions <- scores$sites %>% 
  as_tibble(rownames = 'id_2')

nmds_shared <- inner_join(shared_tbl, nmds_positions, by = 'id_2')

cor_x <- nmds_shared %>% 
  nest(data = -name) %>% 
  mutate(cor_x = map(data,
                     ~cor.test(.x$value,
                               .x$NMDS1, 
                               method = 'spearman',
                               exact = F) %>% tidy())) %>% 
  unnest(cor_x) %>% 
  dplyr::select(name, estimate, p.value)


cor_y <- nmds_shared %>% 
  nest(data = -name) %>% 
  mutate(cor_y = map(data,
                     ~ cor.test(.x$value, .x$NMDS2, 
                                method = 'spearman',
                                exact = F) %>% tidy())) %>% 
  unnest(cor_y) %>% 
  dplyr::select(name, estimate, p.value)


correlations <- inner_join(cor_x, cor_y, by = 'name')

high_cor <- correlations %>% 
  filter(estimate.x > 0.8 | 
           estimate.y > 0.55|
           estimate.x < -0.68 |
           estimate.y < -0.6) %>% 
  mutate(name = str_replace_all(name, '--.*', ''))



## Plot --------------------------------------------------------------------

nmds_samples %>%
  ggplot(aes(x = NMDS1,
             y = NMDS2,
             colour = group)) +
  
  # geom_text(aes(label = id_2),family = 'Times', size = 3, check_overlap = F) +
  geom_point(aes(fill = combined), size = 1.5) +
  stat_ellipse(level = 0.95) +
  geom_point(data = centroid, size = 5, shape = 21, colour = 'black', aes(fill = combined)) +
  geom_segment(data = high_cor,
               aes(x = 0, xend = estimate.x*20, y = 0, yend = estimate.y*20),
               inherit.aes = F)+
  geom_text_repel(data = high_cor, family = 'Times', max.overlaps = 30, direction = 'both',
                  segment.colour = NA, point.padding = 2,
                  aes(x = estimate.x*20, y = estimate.y*20, label = name),
                  size = 2.5, inherit.aes = F, xlim = c(-20,22))+
  scale_fill_futurama(labels = c('LOWCHO', 'LOWSUG', 'MODSUG'))+
  scale_colour_futurama()+
  xlim(-20,20) +
  ylim(-20,20)+
  guides(colour = F)+
  labs(colour = NULL, fill = 'Group', family = 'Times')+
  ggtitle('Week 12')+
  theme_clean_rgd_no_grid()+
  theme(text = element_text(family = 'Times',
                            colour = 'black',
                            size = 10),
        plot.title = element_text(size = 10),
        plot.background = element_rect(color = "white"),
        legend.text = element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'))

ggsave(filename = './Figures/chebi_nmds_biplot_C.tiff',
       device = 'tiff',
       width = 16,
       height = 8,
       units = 'cm',
       dpi = 600)







