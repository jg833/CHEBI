library(tidyverse)

source('./Scripts/chebi_kb_load_data.R')
samples <- samples %>% 
  arrange(participant, trial)
samples$combined <- paste0(samples$group,sep = "-", samples$trial)

samples_id_2 <- samples %>% as_tibble()


# mrc_h2_p_abund_norm_unstr_id_2 ------------------------------------------

h2 <- read_delim('./Data/RD_bath_H3_PA_cpm_renorm_unstratified.tsv')

h2_pabund_norm_unstr_id_2 <- h2 %>% 
  column_to_rownames(var = "# Pathway") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "id") %>% 
  mutate(id = str_replace_all(id, "_S.*", "")) %>%
  # mutate(id = as.double(id)) %>% view()
  semi_join(samples_id_2, by = 'id') %>% 
  arrange(id) %>% 
  column_to_rownames(var = 'id') %>% 
  rownames_to_column(var = 'id_2')

write_delim(h2_pabund_norm_unstr_id_2, file = './Data/chebi_h2_p_abund_norm_unstr_id_2.csv')


# mrc_h2_p_abund_norm_str_id_2 --------------------------------------------
# 
# h2 <- read_delim('./Data/Humann2/humann2_pathabundance_normalised_stratified.tsv')
# 
# h2_pabund_norm_str_id_2 <- h2 %>% 
#   column_to_rownames(var = "# Pathway") %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   rownames_to_column(var = "gut_sample") %>% 
#   mutate(gut_sample = str_replace_all(gut_sample, "_S.*", "")) %>% 
#   mutate(gut_sample = as.double(gut_sample)) %>% 
#   inner_join(samples_id_2, by = 'gut_sample') %>% 
#   arrange(id_2) %>% 
#   column_to_rownames(var = 'id_2') %>% 
#   select(-gut_sample) %>% 
#   rownames_to_column(var = 'id_2') 
# 
# write_delim(h2_pabund_norm_str_id_2, file = './Data/Humann2/mrc_h2_p_abund_norm_str_id_2')

