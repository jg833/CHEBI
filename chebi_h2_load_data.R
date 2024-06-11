data_1 <- read_delim('./Data/chebi_h2_p_abund_norm_unstr_id_2.csv')
samples <- read.csv2('./Data/meta_chebi.csv')
samples <- samples %>% 
  arrange(id_2)
samples$combined <- paste0(samples$group,sep = "-", samples$trial)
vegan <- data_1 %>% column_to_rownames(var = "id_2")
windowsFonts("Times" = windowsFont("Times"))
data_1$id_2
samples$id_2
rownames(vegan)