library(tidyverse)
acid_extractome <- read_tsv("data-raw/AcidExtractome_Data.txt")
acid_pvals <- read_tsv("data-raw/AcidExtractome_pvals.txt")
gene_expr <- read_tsv("data-raw/RNAseq_norm_count_matrix.txt")


genes_long <- gene_expr %>%
  pivot_longer(-Gene) %>%
  separate(name, into = c("condition", "rep"), sep = " ") %>%
  drop_na() %>%
  rename(Accession = Gene)

saveRDS(genes_long, "data/genes_long.rds")


acid_long <- acid_extractome %>% 
  pivot_longer(c(-Accession, -Description)) %>%
  separate(name, into = c("condition", "rep"), sep = " ") %>%
  drop_na()

fold_change_mean <- acid_long %>%
  filter(condition %in% c("Naive", "Primed")) %>%
  group_by(Accession, condition) %>%
  summarise(mean = mean(value)) %>%
  mutate(mean = if_else(mean == 0, 1, mean)) %>%  # replace 0s with value of 1 - if we go for 0.0000001, we'll get ridiculous fc values
  pivot_wider(names_from = condition, values_from = mean) %>%
  mutate(naive_primed_fc = Naive/Primed) %>%
  mutate(log2fc_naive_primed = log2(naive_primed_fc))
  
acid_pval_fc <- acid_pvals %>%
  select(1:2) %>%
  left_join(fold_change_mean)

acid_meta <- acid_extractome %>%
  select(Accession, Description) %>%
  separate(Description, sep = "GN=", into = c(NA, "Gene"), remove = FALSE) %>%
  separate(Gene, sep = " ", into = c("Gene", NA)) %>%
  select(Accession, Gene, Description) 

# add the extra genes in to table
all_gene_names <- pull(gene_expr, Gene) 
genes_not_in <- tibble::enframe(all_gene_names[!all_gene_names %in% acid_meta$Gene], name = NULL, value = "Gene")

genes_to_add <- genes_not_in %>%
  mutate(Accession = "") %>%
  mutate(Description = "") %>%
  relocate(Accession) 

meta <- bind_rows(acid_meta, genes_to_add)

saveRDS(meta, "data/meta.rds")
saveRDS(acid_long, "data/acid_long.rds")
saveRDS(acid_pval_fc, "data/acid_pval_fc.rds")


#this is a rough check that the values we've got match (ish) the pval dataset that's been provided - 
#not proper way to do stats - don't use these!!
t_test_check <- acid_long %>%
  filter(condition %in% c("Naive", "Primed")) %>%
  #mutate(mean = if_else(value == 0, 0.1, value)) %>% # otherwise t_test won't work for loads of 0 values
  group_by(Accession) %>%
  rstatix::t_test(value~condition) %>%
  ungroup() %>%
  select(Accession, p)
  

pval_check <- acid_pval_fc %>%
  left_join(t_test_check)

saveRDS(pval_check, "data/pval_check_temp.rds")





