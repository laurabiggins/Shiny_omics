library(tidyverse)
acid_extractome <- read_tsv("data-raw/AcidExtractome_Data.txt")
acid_pvals <- read_tsv("data-raw/AcidExtractome_pvals.txt")
gene_expr <- read_tsv("data-raw/RNAseq_norm_count_matrix.txt")


## gene expr data ----

genes_long <- gene_expr %>%
  pivot_longer(-Gene) %>%
  separate(name, into = c("condition", "rep"), sep = " ") %>%
  drop_na() %>%
  rename(Accession = Gene)

saveRDS(genes_long, "data/genes_long.rds")

## chromatin acid extractome data ---- 

acid_long <- acid_extractome %>% 
  pivot_longer(c(-Accession, -Description)) %>%
  separate(name, into = c("condition", "rep"), sep = " ") %>%
  drop_na()

saveRDS(acid_long, "data/acid_long.rds")

## histone data ----

PXGL <- read_tsv("data-raw/hPTM_PXGL.txt")
ENHSM <- read_tsv("data-raw/hPTM_ENHSM.txt")
t2iLGo <- read_tsv("data-raw/hPTM_t2iLGo.txt")

process_histone_data <- function(dataset, medium_name){
  dataset %>% 
    select(c(-Histone, -Sequence, -PTM)) %>%
    pivot_longer(-`Histone mark`) %>%
    mutate(value = as.numeric(value)) %>%
    separate(name, into = c("condition", "rep"), sep = " ") %>%
    drop_na() %>%
    mutate(medium = medium_name)
}

histone_data <- bind_rows(
  process_histone_data(PXGL, "PXGL"),
  process_histone_data(ENHSM, "ENHSM"),
  process_histone_data(t2iLGo, "t2iLGo")
)

histone_data$`Histone mark` <- str_replace(histone_data$`Histone mark`, pattern = "Ac", replacement = "ac")


saveRDS(histone_data, "data/histone_data.rds")


## All data types meta data processing ---- 

acid_meta <- acid_extractome %>%
  select(Accession, Description) %>%
  separate(Description, sep = "GN=", into = c(NA, "Gene"), remove = FALSE) %>%
  separate(Gene, sep = " ", into = c("Gene_expr_id", NA)) %>%
  select(Accession, Gene_expr_id, Description) 

# remove gene ids that don't have expression values
all_gene_names <- pull(gene_expr, Gene) 
acid_meta <- acid_meta %>%
  mutate(Gene_expr_id = if_else(Gene_expr_id %in% all_gene_names, Gene_expr_id, ""))


# add the extra genes in to table
genes_not_in <- tibble::enframe(
  all_gene_names[!all_gene_names %in% acid_meta$Gene_expr_id], 
  name = NULL, 
  value = "Gene_expr_id"
)

genes_to_add <- genes_not_in %>%
  mutate(Accession = "") %>%
  mutate(Description = "") %>%
  relocate(Accession) 

meta <- bind_rows(acid_meta, genes_to_add)


histone_links <- read_tsv("data-raw/Links between histones and chromatin proteins.txt") %>%
  rename(Gene_expr_id = `Edgelist B`) %>%
  rename("histone_mark" = `Edgelist A`)

meta <- meta %>%
  left_join(histone_links)


saveRDS(meta, "data/meta.rds")


# check if the histone matching was good enough
sum(histone_links$histone_mark %in% x$histone_mark)
sum(!histone_links$histone_mark %in% x$histone_mark)

hist_links_matched <- histone_links %>%
  mutate(matched = histone_links$Gene_expr_id %in% x$Gene_expr_id)

# It's mainly where the gene id is not actually a gene id, but a histone where we don't have matches. There are a few genes that don't match though.
# There's extra info where the histones are matched to histone marks which I'm
# ignoring for now.




## p-value processing ----

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





