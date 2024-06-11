# CHEBI Kraken/Braken Data Cleaning

# Separates files by study,
# then joins to add participant, trial, participant_trial and group columns,
# then saves as .txt files for direct import for analysis.

# packages ####
library(tidyverse)
library(PERFect)

# clean up sample lists ####
# chebi
samples <- read_csv2('./Data/chebi_vitadex_sample_list.csv')
samples_split <- strsplit(samples$id, "(?=[ABC])", perl = TRUE) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>% 
  rename(participant = V1, trial = V2, study = V3)
samples_joined <- bind_cols(samples, samples_split)
# write.csv2(samples_joined, file = './Data/chebi_vdex_sample_list_cleaned.csv')
samples_2 <- read_csv2('./Data/chebi_vdex_sample_list_cleaned.csv') %>% 
  subset(select = id:study)
samples_chebi <- samples_2 %>% filter(study == 'chebi')
samples_groups <- read_delim('./Data/chebi_sample_groups.csv', delim = ',', col_types = 'ccc')
samples_chebi <- samples_chebi %>%
  inner_join(samples_groups, by = 'participant') %>% 
  mutate(id_2 = paste0(participant,trial, sep = "")) %>% 
  dplyr::select(id:trial, id_2, group)

write_csv2(samples_chebi, file = './Data/meta_chebi.csv')


# vdex
# metadata_vdex <- read_delim('./Data/vdex_samples_list_3.csv', delim = ',', col_types = 'ccccccccccc') %>%
#   rename(code = id_2)
# samples_vdex <- samples_2 %>% filter(study == 'vitadex')
# samples_vdex_joined <- inner_join(samples_vdex, metadata_vdex, by = 'participant') %>% 
#   mutate(id_2 = paste0(participant,trial, sep = "")) %>% 
#   select(id:trial, id_2, group)


# tidy braken data ####
kb_s <- read.delim2('./Data/RD_bath_bracken_S_counts.txt',
                    header = F,
                    sep = "\t" ) %>% t()

colnames(kb_s) <- kb_s[1,]
kb_s <- kb_s[-1,] %>%
  as_tibble()

kb_s <- kb_s %>% 
  rename('id' = 'name--taxonomy_id') %>% 
  mutate(id = str_replace_all(id, "_S\\d.*", ""))

# separate braken data by study ####
# chebi
kb_s_chebi <- inner_join(samples_chebi, kb_s, by = 'id')
write_delim(kb_s_chebi, file = './Data/kb_s_chebi.txt', delim = '\t')

# vdex
kb_s_vdex <- inner_join(samples_vdex_joined, kb_s, by = 'id')
write_delim(kb_s_vdex, file = './Data/kb_s_vdex.txt', delim = '\t')
kb_s_vdex_2 <- read_delim('./Data/kb_s_vdex.txt', col_types = 'ccccc')

# present in >=5% of samples  ####
data <- read_delim('./Data/kb_s_chebi.txt', col_types = 'cdccc')

data_vegan <- data %>% 
  arrange(participant, trial) %>% 
  column_to_rownames('id_2') %>% 
  select(-id:-group)

# Keep taxa present in >X% of samples:
data_vegan_5 <- data_vegan[ ,colMeans(data_vegan > 0) >= .05] # Down from 7702 to 6247
# samples_present <- colSums(data_vegan_5 != 0)
# min(samples_present) # Shows smallest number of samples any taxon is present in, minimum is 7 for 5% threshold with 137 samples

# data_vegan_10 <- data_vegan[ ,colMeans(data_vegan > 0) >= .1] # Down from 7702 to 5513
# samples_present <- colSums(data_vegan_10 != 0) 
# min(samples_present)
# 
# data_vegan_50 <- data_vegan[ ,colMeans(data_vegan > 0) >= .5] # Down from 7702 to 1632
# samples_present <- colSums(data_vegan_50 != 0) 
# min(samples_present)

# PERFect ####
# Packages
library(tidyverse)
library(PERFect)
library(beepr)

# On 5% filtered data, drops from 6247 to 4208!
# data_vegan_5_perf_sim <- PERFect_sim(data_vegan_5)
# data_vegan_5_sim <- data_vegan_5_perf_sim$filtX %>% as_tibble()

# On 5% filtered data, drops from 6247 to 2998
set.seed(123)
data_vegan_5_perf_perm <- PERFect_perm(data_vegan_5, algorithm = "fast")
data_vegan_5_perm <- data_vegan_5_perf_perm$filtX %>% as.data.frame() %>% rownames_to_column(var = "id_2")

write_delim(data_vegan_5_perm, file = './Data/chebi_kb_vegan_5_perm_seed123_id_2', delim = '\t')

# # On unfiltered data, removes only 5 species.
# data_vegan_perfect_sim <- PERFect_sim(data_vegan)
# data_sim <- data_vegan_perfect_sim$filtX %>% as_tibble()
# 
# data_vegan_f_perm <- PERFect_perm(data_vegan_f, algorithm = "fast")
# data_vegan_f_perm_out <- data_vegan_f_perm$filtX %>% as.data.frame()
# write_delim(data_vegan_f_perm_out, file = './Data/kb_vegan_perfperm.txt', delim = '\t')
# 
# data_vegan_perm <- PERFect_perm(data_vegan, algorithm = "fast")
# data_vegan_f_perm_out <- data_vegan_perm$filtX %>% as.data.frame()


# Genus and other levels -------------------------------------------------------------------

samples <- read.csv2('./Data/meta_chebi.csv')

kb_s <- read.delim2('./Data/RD_bath_bracken_P_counts.txt',
                    header = F,
                    sep = "\t" ) %>% t()

colnames(kb_s) <- kb_s[1,]
kb_s <- kb_s[-1,] %>%
  as_tibble()

kb_s <- kb_s %>% 
  rename('id' = 'name--taxonomy_id') %>% 
  mutate(id = str_replace_all(id, "_S\\d.*", ""))

# separate braken data by study ####
# chebi
samples_chebi <- samples %>% 
  select(id, id_2)
kb_s_chebi <- inner_join(samples_chebi, kb_s, by = 'id') %>% 
  relocate(id_2, .after = id) %>% 
  select(-id) %>% 
  as_tibble() %>% 
  arrange(id_2)
  
write_delim(kb_s_chebi, file = './Data/kb_p_chebi.txt', delim = '\t')


# PERFect -----------------------------------------------------------------

library(tidyverse)
library(PERFect)

# Phylum
data_vegan <- read_delim('./Data/kb_p_chebi.txt') %>% 
  column_to_rownames(var = 'id_2')
  
data_vegan_5 <- data_vegan[ ,colMeans(data_vegan > 0) >= .05] # Down from 60 to 49

set.seed(123)
data_vegan_5_perf_perm <- PERFect_perm(data_vegan_5, algorithm = "fast")
data_vegan_5_perm <- data_vegan_5_perf_perm$filtX %>% as.data.frame() %>% rownames_to_column(var = "id_2")

write_delim(data_vegan_5_perm, file = './Data/chebi_kb_Phylum_vegan_5_perm_seed123_id_2', delim = '\t')

# Class
data_vegan <- read_delim('./Data/kb_c_chebi.txt') %>% 
  column_to_rownames(var = 'id_2')

data_vegan_5 <- data_vegan[ ,colMeans(data_vegan > 0) >= .05] # Down from 60 to 49

set.seed(123)
data_vegan_5_perf_perm <- PERFect_perm(data_vegan_5, algorithm = "fast")
data_vegan_5_perm <- data_vegan_5_perf_perm$filtX %>% as.data.frame() %>% rownames_to_column(var = "id_2")

write_delim(data_vegan_5_perm, file = './Data/chebi_kb_Class_vegan_5_perm_seed123_id_2', delim = '\t')

# Order
data_vegan <- read_delim('./Data/kb_o_chebi.txt') %>% 
  column_to_rownames(var = 'id_2')

data_vegan_5 <- data_vegan[ ,colMeans(data_vegan > 0) >= .05] # Down from 60 to 49

set.seed(123)
data_vegan_5_perf_perm <- PERFect_perm(data_vegan_5, algorithm = "fast")
data_vegan_5_perm <- data_vegan_5_perf_perm$filtX %>% as.data.frame() %>% rownames_to_column(var = "id_2")

write_delim(data_vegan_5_perm, file = './Data/chebi_kb_Order_vegan_5_perm_seed123_id_2', delim = '\t')

# Family
data_vegan <- read_delim('./Data/kb_f_chebi.txt') %>% 
  column_to_rownames(var = 'id_2')

data_vegan_5 <- data_vegan[ ,colMeans(data_vegan > 0) >= .05] # Down from 60 to 49

set.seed(123)
data_vegan_5_perf_perm <- PERFect_perm(data_vegan_5, algorithm = "fast")
data_vegan_5_perm <- data_vegan_5_perf_perm$filtX %>% as.data.frame() %>% rownames_to_column(var = "id_2")

write_delim(data_vegan_5_perm, file = './Data/chebi_kb_Family_vegan_5_perm_seed123_id_2', delim = '\t')


# Genus
data_vegan <- read_delim('./Data/kb_g_chebi.txt') %>% 
  column_to_rownames(var = 'id_2')

data_vegan_5 <- data_vegan[ ,colMeans(data_vegan > 0) >= .05] # Down from 60 to 49

set.seed(123)
data_vegan_5_perf_perm <- PERFect_perm(data_vegan_5, algorithm = "fast")
data_vegan_5_perm <- data_vegan_5_perf_perm$filtX %>% as.data.frame() %>% rownames_to_column(var = "id_2")

write_delim(data_vegan_5_perm, file = './Data/chebi_kb_Genus_vegan_5_perm_seed123_id_2', delim = '\t')