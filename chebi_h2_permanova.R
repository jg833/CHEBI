# Packages ####
library(zCompositions)
library(tidyverse)
library(vegan)
library(ggsci)
library(RVAideMemoire)
library(beepr)

# Data ####
source('./Scripts/chebi_h2_load_data.R')
samples$combined <- paste0(samples$group,sep = "-", samples$trial)


# CLR Transform ####
vegan_no0 <- cmultRepl(vegan, method = "CZM", output = "p-counts")

vegan_clr <- t(apply(vegan_no0, 1, function(x){log(x) - mean(log(x))})) %>% as.data.frame()

# PERMANOVA ####
# All samples, group * time
set.seed(123)
adonis2(vegan_clr ~ group * trial,
        data = samples,
        method = "euclidean",
        permutations = 9999, 
        strata = samples$participant) %>% 
  print()

# samples$combined <- paste0(samples$group,sep = "_", samples$trial)

dist <- vegdist(vegan_clr, method = "euclidean")
set.seed(123)
pairwise.perm.manova(dist, samples$combined,
                     p.method = "BH",
                     nperm = 9999) 




# Test - adonis2 method = 'bray' ------------------------------------------

# set.seed(123)
# adonis2(vegan_no0 ~ group * trial,
#         data = samples,
#         method = "bray",
#         permutations = 9999, 
#         strata = samples$participant) %>% 
#   print() # Significant interaction effect - hard to interpret, try separating data by group/ trial and rerunning?
# 

# adonis2 method = 'robust.aitchison' ------------------------------------------

set.seed(123)
adonis2(vegan ~ group * trial,
        data = samples,
        method = "robust.aitchison",
        permutations = 9999, 
        strata = samples$participant) %>% 
  print()

dist <- vegdist(vegan, method = "robust.aitchison")
set.seed(123)




# Proper way - Final! -----------------------------------------------------

# Near significant effect of group. 
# Pairwise comparisons suggest LOWCHO-C is different to LOWSUG-B & C, MODSUG-A & B

samples$id_2
rownames(vegan)

vegan_dist <- vegdist(vegan, method = 'robust.aitchison')

set.seed(123)
adonis2(vegan_dist ~ group * trial,
        data = samples,
        permutations = 9999,
        strata = samples$participant, 
        pairwise = T) %>% 
  print()

set.seed(123)
pairwise.perm.manova(vegan_dist, samples$combined,
                     p.method = "BH",
                     nperm = 9999) 

## Meandist and  mrpp -----------------------------------------------------------



 # Meandist gives same results as in beta boxplot script, reassuring
meandist(dat = vegan_clr,
         grouping = samples$combined,
         distance = 'robust.aitchison',
         strata = samples$participant,
         dist = vegan_dist) %>% 
  plot()
# %>% 
#   summary()

mrpp(dat = vegan_dist, grouping = samples$combined, permutations = 9999, strata = samples$participant)



# Subsetting --------------------------------------------------------------

# Main effect of trial found for LOWCHO group only.

samples_f <- samples %>% 
  filter(group == 'LOWCHO') # Change group here

vegan_f <- data_1 %>%
  semi_join(samples_f, by = 'id_2') %>% 
  column_to_rownames(var = 'id_2')

samples_f$id_2
rownames(vegan_f)

v_f_dist <- vegdist(vegan_f, method = 'robust.aitchison')

set.seed(123)
adonis2(v_f_dist ~ trial,
        data = samples_f,
        permutations = 9999, 
        strata = samples_f$participant) %>% 
  print()


# No significant differences using pairwise comparisons, due to lack of stratification by participant?
set.seed(123)
pairwise.perm.manova(resp = v_f_dist,
                     fact = samples_f$combined,
                     p.method = "BH",
                     nperm = 99999)


# Pairwise - trials -------------------------------------------------------


# Filter by group and also trial - homemade pairwise PERMANOVA 


samples_f <- samples %>% 
  filter(group == 'LOWCHO', trial == 'A' | trial == 'C') # Change group and trial here

vegan_f <- data_1 %>%
  semi_join(samples_f, by = 'id_2') %>% 
  column_to_rownames(var = 'id_2')

samples_f$id_2
rownames(vegan_f)

v_f_dist <- vegdist(vegan_f, method = 'robust.aitchison')

set.seed(123)
adonis2(v_f_dist ~ trial,
        data = samples_f,
        permutations = 9999, 
        strata = samples_f$participant) %>% 
  print()



pvals <- c(0.7334,
           0.26,
           0.4744,
           0.4028,
           0.8708,
           0.7907,
           0.4256,
           0.2653,
           0.0691)

p.adjust(pvals, method = 'BH') 




# Pairwise - groups -------------------------------------------------------

samples_f <- samples %>% 
  filter(trial == 'C', group == 'LOWSUG' | group == 'LOWCHO') # Change group and trial here

vegan_f <- data_1 %>%
  semi_join(samples_f, by = 'id_2') %>% 
  column_to_rownames(var = 'id_2')

samples_f$id_2
rownames(vegan_f)
samples_f$combined

v_f_dist <- vegdist(vegan_f, method = 'robust.aitchison')

set.seed(123)
adonis2(v_f_dist ~ group,
        data = samples_f,
        permutations = 9999) %>% 
  print()


pvals <- c(0.3824,
           0.1981,
           0.8181,
           
           0.6071,
           0.0037,
           0.1157,
           
           0.6412,
           0.0032,
           0.0139)
p.adjust(pvals, method = 'BH') 




# Beta disper -------------------------------------------------------------

# https://rfunctions.blogspot.com/2016/08/measuring-and-comparing-beta-diversity.html

vegan_dist
samples

bd <- betadisper(vegan_dist, group = samples$combined)

bd_boxplot <- betadisper(vegan_dist, group = samples$combined) %>% 
  boxplot()

bd_pairwise <- betadisper(vegan_dist, group = samples$combined) %>% 
  permutest(permutations = 9999, strata = samples$participant, pairwise = T)

bd_pairwise$pairwise$observed %>% 
  p.adjust(method = 'BH') %>% 
  as.data.frame() %>% 
  min()

bd_pairwise$pairwise$permuted %>% 
  # p.adjust(method = 'BH') %>% 
  as.data.frame() 


bray.part



m <- betadiver(vegan_clr)
plot(m)
