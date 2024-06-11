data_1 <- read_delim('./Data/chebi_kb_vegan_5_perm_seed123_id_2')
data_g <- read_delim('./Data/chebi_kb_Genus_vegan_5_perm_seed123_id_2')
data_f <- read_delim('./Data/chebi_kb_Family_vegan_5_perm_seed123_id_2')
data_o <- read_delim('./Data/chebi_kb_Order_vegan_5_perm_seed123_id_2')
data_c <- read_delim('./Data/chebi_kb_Class_vegan_5_perm_seed123_id_2')
data_p <- read_delim('./Data/chebi_kb_Phylum_vegan_5_perm_seed123_id_2')

samples <- read.csv2('./Data/meta_chebi.csv')
samples <- samples %>% 
  arrange(participant, trial)
samples$combined <- paste0(samples$group,sep = "-", samples$trial)
vegan <- data_1 %>% column_to_rownames(var = "id_2")
windowsFonts("Times" = windowsFont("Times"))
