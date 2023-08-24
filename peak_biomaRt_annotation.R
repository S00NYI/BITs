library(biomaRt)
library(dplyr)
library(stringr)

## NOTE: useast mirror is always more reliable than the main:
GRCh38 = useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "https://useast.ensembl.org")

## List of available attributes and filters:
attributes = listAttributes(GRCh38)
filters = listFilters(GRCh38)

## Choose attributes and filters needed:
selected_attributes = c('chromosome_name', 'start_position', 'end_position', 'external_gene_name', 'gene_biotype', 'strand')
selected_filters = c('chromosome_name')

## from countAnnotations, filter gene_id by i) start with "ENSG" and ii) has length of 15 and get the gene_name and gene_biotype
martOut = getBM(selected_attributes,
                selected_filters,
                c(1:22, 'X', 'Y', 'MT'),
                mart = GRCh38)

martOut = martOut %>% mutate(chromosome_name = paste0('chr', chromosome_name))
martOut = martOut %>% mutate(strand = ifelse(strand == "-1", "-", strand))
martOut = martOut %>% mutate(strand = ifelse(strand == "1", "+", strand))
martOut = martOut %>% mutate(external_gene_name = ifelse(external_gene_name == "", ".", external_gene_name))

martOut = martOut %>% mutate(gene_biotype = ifelse(str_sub(gene_biotype, -10) == "pseudogene", "pseudogene", gene_biotype))

martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "protein_coding", "mRNA", gene_biotype))

# martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "snoRNA", "other ncRNA", gene_biotype))
martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "scaRNA", "ncRNA", gene_biotype))
martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "scRNA", "ncRNA", gene_biotype))
martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "vault_RNA", "ncRNA", gene_biotype))
martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "sRNA", "ncRNA", gene_biotype))

martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "TEC", "other", gene_biotype))
martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "ribozyme", "other", gene_biotype))
martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "TR_V_gene", "other", gene_biotype))
martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "TR_D_gene", "other", gene_biotype))
martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "TR_J_gene", "other", gene_biotype))
martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "TR_C_gene", "other", gene_biotype))
martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "IG_V_gene", "other", gene_biotype))
martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "IG_D_gene", "other", gene_biotype))
martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "IG_J_gene", "other", gene_biotype))
martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "IG_C_gene", "other", gene_biotype))
martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "misc_RNA", "other", gene_biotype))
martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "artifact", "other", gene_biotype))

martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "Mt_tRNA", "mitochondrial", gene_biotype))
martOut = martOut %>% mutate(gene_biotype = ifelse(gene_biotype == "Mt_rRNA", "mitochondrial", gene_biotype))

write.table(martOut, paste0( '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/Annotations.bed'), row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
