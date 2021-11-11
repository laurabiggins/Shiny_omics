library(tidyverse)
acid_extractome <- read_tsv("data-raw/AcidExtractome_Data.txt")
gene_expr <- read_tsv("data-raw/RNAseq_norm_count_matrix.txt")


## gene expr data ----

genes_long <- gene_expr %>%
  pivot_longer(-Gene) %>%
  separate(name, into = c("condition", "rep"), sep = " ") %>%
  drop_na() %>%
  rename(Gene_expr_id = Gene)

saveRDS(genes_long, "data/genes_long.rds")

## chromatin acid extractome data ---- 

acid_long <- acid_extractome %>% 
  pivot_longer(c(-Accession, -Description)) %>%
  separate(name, into = c("condition", "rep"), sep = " ") %>%
  drop_na()

#saveRDS(acid_long, "data/acid_long.rds") # this is saved later on

## chep data ----
## this has got p values which I'l deal with afterwards
chep1 <- read_tsv("data-raw/ChEP_naive_primed.txt")
chep2 <- read_tsv("data-raw/ChEP_naive2i.txt")
chep3 <- read_tsv("data-raw/ChEP_primed2i.txt")

chep1 <- chep1 %>%
  select(-1, -`log2 fold change`, -`p-value`) %>%
  pivot_longer(cols = starts_with("LFQ")) %>%
  separate(name, into = c(NA, NA, "condition", "rep"), sep = " ")

chep2 <- chep2 %>%
  select(-1, -`log2 fold change`, -`p-value`) %>%
  pivot_longer(cols = starts_with("LFQ")) %>%
  separate(name, into = c(NA, NA, "condition", "rep"), sep = " ")

chep3 <- chep3 %>%
  select(-1, -`log2 fold change`, -`p-value`) %>%
  pivot_longer(cols = starts_with("LFQ")) %>%
  separate(name, into = c(NA, NA, "condition", "rep"), sep = " ")

# chep_data <- bind_rows(chep1, chep2, chep3) %>%
#   distinct() %>%
#   rename(Gene_id = Gene.names) %>%
#   drop_na(Gene_id)

# extra annotations
chep_curated <- read_delim("data-raw/chep_curation.txt") %>%
  rename(Gene.names = gene_id_chep, Gene_id = `Name in RNA-seq dataset`) 

chep_data <- bind_rows(chep1, chep2, chep3) %>%
  distinct() %>%
  left_join(chep_curated) %>%
  mutate(Gene_id = if_else(is.na(Gene_id), Gene.names, Gene_id))


saveRDS(chep_data, "data/chep_data.rds")

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

histone_data <- dplyr::rename(histone_data, histone_mark = `Histone mark`)
histone_data$histone_mark <- str_replace(histone_data$histone_mark, pattern = "Ac", replacement = "ac")

saveRDS(histone_data, "data/histone_data.rds")


## Meta - all data types ---- 
# we'll join by gene name
all_gene_expr_names <- pull(gene_expr, Gene) 

# extra manual curations to match ids
acid_curated <- read_delim("data-raw/acid_curated.txt") %>%
  rename(Accession = Acid_id)

acid_meta <- acid_extractome %>%
  select(Accession, Description) %>%
  separate(Description, sep = "GN=", into = c(NA, "Gene"), remove = FALSE) %>%
  separate(Gene, sep = " ", into = c("Gene_id", NA)) %>%
  select(Accession, Gene_id, Description) %>%
  left_join(acid_curated) %>%
  mutate(Gene_id = if_else(is.na(gene), Gene_id, gene)) %>%
  select(-gene)

# don't need to keep doing this
#sum(acid_meta$Gene_id %in% all_gene_expr_names)
#sum(!acid_meta$Gene_id %in% all_gene_expr_names)
#View(acid_meta[!acid_meta$Gene_id %in% all_gene_expr_names,])
#
## Functions for splitting up gene ids to check for matches.
## This isn't needed for the acid dataset as there are not multiple gene ids, but it is needed for the #Chep data
#which_match <- function(id){
#  #split the id in case there are multiple ids
#  if(sum(str_detect(all_gene_expr_names, id)) == 0) return (NA)
#  else return(id)
#}
#gene_match <- function(id){
#  print(all_gene_expr_names[str_detect(all_gene_expr_names, id)])
#}
#
#which_multi <- function(id){
#  ids <- str_split(id, ";")[[1]]
#  matches <- unlist(map(ids, which_match))
#  matches[!is.na(matches)][1]
#}
#
#unmatched <- acid_meta$Gene_id[!acid_meta$Gene_id %in% all_gene_expr_names]
#unmatched_acid_ids <- acid_meta[!acid_meta$Gene_id %in% all_gene_expr_names,]
#write_tsv(unmatched_acid_ids, file = "unmatched_acid_ids.txt") # write out the set of unmatched acid gene ids

# No duplicated gene ids in the acid dataset
#sum(duplicated(acid_meta$Gene_id)) 


# The chep data has about 35 duplicated gene names, we can use the gene names 
# to join the tables but we're not using the genes on the volcano plot so we don't need to worry about multiple matches.
# We should add the gene expr data after the protein data has been sorted. 
chep_wide_ids <- chep_data %>%
  select(Majority.protein.IDs, Protein.names, Gene_id) %>%
  distinct() 

#n_distinct(chep_wide_ids$Majority.protein.IDs)
#n_distinct(chep_wide_ids$Protein.names)
#n_distinct(chep_wide_ids$Gene_id)
#
## The majority protein_ids are unique so we'll use those for identifying the ChEP data.
## Have a look at the dup genes
#chep_wide_ids$Gene_id[duplicated(chep_wide_ids$Gene_id)] # quite a few NAs here, that's not much good.
#dups <- chep_wide_ids[duplicated(chep_wide_ids$Gene_id),]

# there are a lot of multiple gene names in the Chep ids so we'll need to separate them and
# see which ones match



# this takes a while so don't keep re-running it, but it does work quite nicely
# matched_gene_id <- map_chr(chep_wide_ids$Gene_id, which_multi)


chep_wide <- chep_wide_ids %>%
  rename(gene_id_chep = Gene_id) %>%
  mutate(Gene_id = if_else(gene_id_chep %in% all_gene_expr_names, gene_id_chep, NA_character_)) 
  
#unmatched_chep_ids <- chep_wide %>%
#  filter(is.na(Gene_id)) %>%
#  select(-Gene_id)
#

# I've looked at the unmatched set and I don't think I can do anymore matching with this
# it'll need some manual input to find any more matches
#write_tsv(unmatched_chep_ids, file = "unmatched_chep_ids.txt")


# join acid to chep, then gene_expr, then histones
meta1 <- full_join(acid_meta, chep_wide, na_matches = "never")
  

# add a column of gene_expr_ids that are in the gene expr dataset. 
# We want to keep the gene_id as a reference but have a separate column that has a missing value if we don't have gene expr data for that one.
gene_set <- gene_expr %>%
  select(Gene) %>%
  rename(Gene_expr_id = Gene) %>%
  mutate(Gene_id = Gene_expr_id)

meta2 <- full_join(meta1, gene_set, na_matches = "never")


## add histones in to meta table

histone_links <- read_tsv("data-raw/Links between histones and chromatin proteins.txt") %>%
  rename(Gene_id = `Edgelist B`) %>%
  rename("histone_mark" = `Edgelist A`)

histone_links$histone_mark <- str_replace(histone_links$histone_mark, pattern = "Ac", replacement = "ac")

# remove any from histone_links that we don't have data for
histone_links <- histone_links %>%
  filter(histone_mark %in% histone_data$histone_mark)


# also add any histones that aren't in this file.
all_histones <- histone_data %>%
  select(histone_mark) %>%
  distinct() %>%
  left_join(histone_links)

#histone_data %>%
#filter(histone_mark == "H3K79ac")

meta <- full_join(meta2, all_histones, na_matches = "never")

table_data <- meta %>%
  relocate(Gene_id) %>%
  arrange(desc(Accession)) %>%
  rowid_to_column()

saveRDS(table_data, "data/table_data.rds")

gene_id_table <- table_data %>%
  mutate(`acid extractome` = if_else(!is.na(Accession), "data available", "no data")) %>%
  mutate(chep = if_else(!is.na(Majority.protein.IDs), "data available", "no data")) %>%
  mutate(`gene expr` = if_else(!is.na(Gene_expr_id), "data available", "no data")) %>%
  mutate(histone = if_else(!is.na(histone_mark), histone_mark, "no data")) %>%
  select(rowid, Gene_id, `acid extractome`:last_col(), Accession, Majority.protein.IDs, gene_id_chep, histone_mark, Description)

saveRDS(gene_id_table, "data/gene_id_table.rds")

# add a gene id column to the acid dataset
gene_ids <- gene_id_table %>%
  select(.data[["Accession"]], Gene_id) %>%
  drop_na() 

acid_long <- acid_long %>%
  left_join(gene_ids)

saveRDS(acid_long, "data/acid_long.rds")


# check if the histone matching was good enough
#sum(histone_links$histone_mark %in% x$histone_mark)
#sum(!histone_links$histone_mark %in% x$histone_mark)

#hist_links_matched <- histone_links %>%
#  mutate(matched = histone_links$Gene_expr_id %in% x$Gene_expr_id)

# It's mainly where the gene id is not actually a gene id, but a histone where we don't have matches. There are a few genes that don't match though.
# There's extra info where the histones are matched to histone marks which I'm
# ignoring for now.




## p-value processing ----
# p-values are provided but we need log2 fc too for the volcano plot
# 
### acid pvalues ----
acid_pvals <- read_tsv("data-raw/AcidExtractome_pvals.txt")

fold_change <- acid_long %>%
  group_by(Accession, condition) %>%
  summarise(mean = mean(value)) %>%
  mutate(mean = if_else(mean == 0, 1, mean)) %>%  # replace 0s with value of 1 - if we go for 0.0000001, we'll get ridiculous fc values
  pivot_wider(names_from = condition, values_from = mean) %>%
  mutate(naive_primed = log2(Naive/Primed)) %>%
  mutate(naive_naive2i = log2(Naive/`Naive+PRC2i`)) %>%
  mutate(primed_primed2i = log2(Primed/`Primed+PRC2i`)) %>%
  ungroup() %>%
  select(-(2:5)) %>%
  pivot_longer(cols = !Accession, values_to = "log2fc", names_to = "condition", values_drop_na = TRUE)

 
acid_pval_fc <- acid_pvals %>%
  pivot_longer(cols = !Accession, names_to = "condition", values_to = "pval") %>%
  left_join(fold_change)

### gene expr pvalues ----
gene_pvals <- read_tsv("data-raw/p_values_RNAseq_naive_primed_PRC2i_overlapping_genes.txt", lazy = FALSE) %>%
  rename(Gene_expr_id = Accession)

gene_fc <- genes_long %>%
  group_by(Gene_expr_id, condition) %>%
  summarise(mean = mean(value)) %>%
  mutate(mean = if_else(mean == 0, 1, mean)) %>%  # replace 0s with value of 1 - if we go for 0.0000001, we'll get ridiculous fc values
  pivot_wider(names_from = condition, values_from = mean) %>%
  mutate(naive_primed = log2(Naive/Primed)) %>%
  mutate(naive_naive2i = log2(Naive/`Naive+PRC2i`)) %>%
  mutate(primed_primed2i = log2(Primed/`Primed+PRC2i`)) %>%
  ungroup() %>%
  select(-(2:5)) %>%
  pivot_longer(cols = -Gene_expr_id, values_to = "log2fc", names_to = "condition", values_drop_na = TRUE)

gene_pval_fc <- gene_pvals %>%
  pivot_longer(cols = -Gene_expr_id, names_to = "condition", values_to = "pval") %>%
  left_join(gene_fc)

### histone_pvals ----
histone_pvals <- read_tsv("data-raw/hPTM_pvals.txt")

histone_fc <- histone_data %>%
  group_by(medium, histone_mark, condition) %>%
  summarise(mean = mean(value)) %>%
  ungroup() %>%
  group_by(medium) %>%
  # mutate(mean = if_else(mean == 0, 1, mean)) %>%  # much smaller values than the acid set
  pivot_wider(names_from = condition, values_from = mean) %>%
  mutate(naive_primed = log2(Naive/Primed)) %>%
  mutate(naive_naive2i = log2(Naive/`Naive+PRC2i`)) %>%
  mutate(primed_primed2i = log2(Primed/`Primed+PRC2i`)) %>%
  ungroup() %>%
  select(-(3:6)) %>%
  pivot_longer(cols = !c(medium, histone_mark), values_to = "log2fc", names_to = "condition", values_drop_na = TRUE)

histone_pval_fc <- histone_pvals %>%
  pivot_longer(cols = !histone_mark, names_to = c("medium", "condition"), names_sep = " ", values_to = "pval") %>%
  mutate(pval = as.numeric(pval)) %>%
  filter(medium == "PXGL") %>%
  left_join(histone_fc) %>%
  drop_na()



### chep pvalues ----
chep1 <- read_tsv("data-raw/ChEP_naive_primed.txt")
chep2 <- read_tsv("data-raw/ChEP_naive2i.txt")
chep3 <- read_tsv("data-raw/ChEP_primed2i.txt")

chep1_p <- chep1 %>%
  rename(log2fc = `log2 fold change`, pval = `p-value`) %>%
  select(Majority.protein.IDs, log2fc, pval) %>%
  mutate (condition = "naive_primed")

chep2_p <- chep2 %>%
  rename(log2fc = `log2 fold change`, pval = `p-value`) %>%
  select(Majority.protein.IDs, log2fc, pval) %>%
  mutate (condition = "naive_naive2i")

chep3_p <- chep3 %>%
  rename(log2fc = `log2 fold change`, pval = `p-value`) %>%
  select(Majority.protein.IDs, log2fc, pval) %>%
  mutate (condition = "primed_primed2i")

chep_pval_fc <- bind_rows(chep1_p, chep2_p, chep3_p)

# these pvalues ahve already been -log10 transformed so I'm changing them back 
# in line with the other data types
chep_pval_fc <- mutate(chep_pval_fc, pval = 10^-pval)

# Add gene id column before saving pvalue datasets
# this doesn't need to be in the reactive - extract it!!  

add_gene_id_column <- function(pval_ds, accession_type, gene_id_table){

  gene_ids <- gene_id_table %>%
    select(.data[[accession_type]], Gene_id) %>%
    drop_na() 

  pval_ds %>%
    drop_na() %>%
    left_join(gene_ids) 
}

histone_pval_fc2 <- add_gene_id_column(histone_pval_fc, "histone_mark", gene_id_table)
acid_pval_fc2 <- add_gene_id_column(acid_pval_fc, "Accession", gene_id_table)
chep_pval_fc2 <- add_gene_id_column(chep_pval_fc, "Majority.protein.IDs", gene_id_table)
gene_pval_fc2 <- mutate(gene_pval_fc, Gene_id = Gene_expr_id)


saveRDS(histone_pval_fc2, "data/histone_pval_fc.rds")
saveRDS(acid_pval_fc2, "data/acid_pval_fc.rds")
saveRDS(chep_pval_fc2, "data/chep_pval_fc.rds")
saveRDS(gene_pval_fc2, "data/gene_pval_fc.rds")







