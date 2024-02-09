## CoCLIP Analysis:
## snoRNA Analysis
## Written by Soon Yi
## Last edit: 2024-01-31

library(stringr)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsignif)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)


setwd("~/Desktop/Genomics/Annotations")
gtfFile = 'Homo_sapiens.GRCh38.110.gtf'
gtf = import(gtfFile)

gtf_df = as.data.frame(gtf)

## Subset the data_frame by gene/transcript types:
snoRNA = subset(gtf_df, transcript_biotype == "snoRNA" & type == "transcript")
snoRNA_gr = GRanges(seqnames=snoRNA$seqnames, ranges=IRanges(start=snoRNA$start, end=snoRNA$end, names=snoRNA$gene_name), strand=snoRNA$strand)

## Bed Files to Import:

bedBase = "~/Desktop/Genomics/CoCLIP/combined_collapsedBED/collapsedBED"

ID = c('JL0361_Input_G3BP_Input_Mock_1',
       'JL0361_Input_G3BP_Input_Mock_2',
       'JL0361_Input_G3BP_Input_Ars_1',
       'JL0361_Input_G3BP_Input_Ars_2',
       'JL0361_Input_G3BP_Input_Ars_3',
       'JL0361_Enrich_G3BP_Enrich_Mock_1',
       'JL0361_Enrich_G3BP_Enrich_Mock_2',
       'JL0361_Enrich_G3BP_Enrich_Ars_1',
       'JL0361_Enrich_G3BP_Enrich_Ars_2',
       'JL0361_Enrich_G3BP_Enrich_Ars_3',
       'JL0380_Fraction_Cyto_Fraction_Mock_1',
       'JL0380_Fraction_Cyto_Fraction_Mock_2',
       'JL0380_Fraction_Cyto_Fraction_Mock_3',
       'JL0380_Fraction_Cyto_Fraction_Ars_1',
       'JL0380_Fraction_Cyto_Fraction_Ars_2',
       'JL0380_Fraction_Cyto_Fraction_Ars_3',
       'JL0380_Fraction_Nuc_Fraction_Mock_1',
       'JL0380_Fraction_Nuc_Fraction_Mock_2',
       'JL0380_Fraction_Nuc_Fraction_Mock_3',
       'JL0380_Fraction_Nuc_Fraction_Ars_1',
       'JL0380_Fraction_Nuc_Fraction_Ars_2',
       'JL0380_Fraction_Nuc_Fraction_Ars_3',
       'JL0388_Input_NES_Input_Mock_1',
       'JL0388_Input_NES_Input_Ars_1',
       'JL0388_Input_NLS_Input_Mock_1',
       'JL0388_Input_NLS_Input_Ars_1',
       'JL0388_Input_SG_Input_Mock_1',
       'JL0388_Input_SG_Input_Ars_1',
       'JL0388_Enrich_NES_Enrich_Mock_1',
       'JL0388_Enrich_NES_Enrich_Mock_2',
       'JL0388_Enrich_NES_Enrich_Ars_1',
       'JL0388_Enrich_NES_Enrich_Ars_2',
       'JL0388_Enrich_NLS_Enrich_Mock_1',
       'JL0388_Enrich_NLS_Enrich_Mock_2',
       'JL0388_Enrich_NLS_Enrich_Ars_1',
       'JL0388_Enrich_NLS_Enrich_Ars_2',
       'JL0388_Enrich_SG_Enrich_Mock_1',
       'JL0388_Enrich_SG_Enrich_Mock_2',
       'JL0388_Enrich_SG_Enrich_Ars_1',
       'JL0388_Enrich_SG_Enrich_Ars_2',
       'JL1024_Pool_NES_I_Mock_All',
       'JL1024_Pool_NES_I_Ars_All',
       'JL1024_Pool_NLS_I_Mock_All',
       'JL1024_Pool_NLS_I_Ars_All',
       'JL1024_Pool_G3BP_I_Mock_All',
       'JL1024_Pool_G3BP_I_Ars_All',
       'JL1024_Pool_NES_E_Mock_R1',
       'JL1024_Pool_NES_E_Mock_R2',
       'JL1024_Pool_NES_E_Ars_R1',
       'JL1024_Pool_NES_E_Ars_R2',
       'JL1024_Pool_NLS_E_Mock_R1',
       'JL1024_Pool_NLS_E_Mock_R2',
       'JL1024_Pool_NLS_E_Ars_R1',
       'JL1024_Pool_NLS_E_Ars_R2',
       'JL1024_Pool_G3BP_E_Mock_R1',
       'JL1024_Pool_G3BP_E_Mock_R2',
       'JL1024_Pool_G3BP_E_Ars_R1',
       'JL1024_Pool_G3BP_E_Ars_R2')

for (file in ID) {
  bedFile = read_delim(paste0(bedBase, '/', file, '.sorted.collapsed.bed'), show_col_types = F, col_names = F)
  colnames(bedFile) = c('chr', 'start', 'end', 'name', 'score', 'strand')
  depth = nrow(bedFile)
  
  bedFile$chr = gsub("chr", "", bedFile$chr)
  bedFile_gr = GRanges(seqnames=bedFile$chr, ranges=IRanges(start=bedFile$start, end=bedFile$end, names=bedFile$name), strand=bedFile$strand)
  
  read_per_snoRNA = as.data.frame(findOverlaps(query = snoRNA_gr, subject = bedFile_gr, minoverlap = 1, select = 'first'))
  colnames(read_per_snoRNA) = file
  read_per_snoRNA = read_per_snoRNA %>% mutate(!!sym(file) := ifelse(is.na(!!sym(file)), 0, !!sym(file)))
  
  read_per_snoRNA = read_per_snoRNA/depth/1e6
  
  snoRNA = cbind(snoRNA, read_per_snoRNA)
}

# write_csv(snoRNA, '~/Desktop/Genomics/CoCLIP/snoRNA_reads.csv', quote = 'none', col_names = T)

snoRNA = snoRNA %>% mutate(Mock_Nuclear_Fraction = 
                             (JL0380_Fraction_Nuc_Fraction_Mock_1 + 
                             JL0380_Fraction_Nuc_Fraction_Mock_2 +
                             JL0380_Fraction_Nuc_Fraction_Mock_3))

snoRNA = snoRNA %>% mutate(Stress_Nuclear_Fraction = 
                             (JL0380_Fraction_Nuc_Fraction_Ars_1 + 
                             JL0380_Fraction_Nuc_Fraction_Ars_2 +
                             JL0380_Fraction_Nuc_Fraction_Ars_3))

snoRNA = snoRNA %>% mutate(Mock_Cytoplasm_Fraction = 
                             (JL0380_Fraction_Cyto_Fraction_Mock_1 + 
                             JL0380_Fraction_Cyto_Fraction_Mock_2 +
                             JL0380_Fraction_Cyto_Fraction_Mock_3))

snoRNA = snoRNA %>% mutate(Stress_Cytoplasm_Fraction = 
                             (JL0380_Fraction_Cyto_Fraction_Ars_1 + 
                             JL0380_Fraction_Cyto_Fraction_Ars_2 +
                             JL0380_Fraction_Cyto_Fraction_Ars_3))

snoRNA = snoRNA %>% mutate(Mock_Total_Input = 
                             (JL0388_Input_NLS_Input_Mock_1 + 
                             JL0388_Input_NES_Input_Mock_1 +
                             JL0361_Input_G3BP_Input_Mock_1 +
                             JL0361_Input_G3BP_Input_Mock_2 + 
                             JL1024_Pool_NLS_I_Mock_All +
                             JL1024_Pool_NES_I_Mock_All +
                             JL0388_Input_SG_Input_Mock_1 +
                             JL1024_Pool_G3BP_I_Mock_All))

snoRNA = snoRNA %>% mutate(Stress_Total_Input = 
                             (JL0388_Input_NLS_Input_Ars_1 + 
                             JL0388_Input_NES_Input_Ars_1 +
                             JL0361_Input_G3BP_Input_Ars_1 +
                             JL0361_Input_G3BP_Input_Ars_2 +
                             JL0361_Input_G3BP_Input_Ars_3 +
                             JL1024_Pool_NLS_I_Ars_All +
                             JL1024_Pool_NES_I_Ars_All +
                             JL0388_Input_SG_Input_Ars_1 +
                             JL1024_Pool_G3BP_I_Ars_All))

snoRNA = snoRNA %>% mutate(Mock_NLS_Enrich = 
                             (JL0388_Enrich_NLS_Enrich_Mock_1 + 
                             JL0388_Enrich_NLS_Enrich_Mock_2 +
                             JL1024_Pool_NLS_E_Mock_R1 +
                             JL1024_Pool_NLS_E_Mock_R2))

snoRNA = snoRNA %>% mutate(Stress_NLS_Enrich = 
                             (JL0388_Enrich_NLS_Enrich_Ars_1 +
                             JL0388_Enrich_NLS_Enrich_Ars_2 +
                             JL1024_Pool_NLS_E_Ars_R1 +
                             JL1024_Pool_NLS_E_Ars_R2))

snoRNA = snoRNA %>% mutate(Mock_NES_Enrich = 
                             (JL0388_Enrich_NES_Enrich_Mock_1 + 
                             JL0388_Enrich_NES_Enrich_Mock_2 +
                             JL1024_Pool_NES_E_Mock_R1 +
                             JL1024_Pool_NES_E_Mock_R2))

snoRNA = snoRNA %>% mutate(Stress_NES_Enrich = 
                             (JL0388_Enrich_NES_Enrich_Ars_1 +
                             JL0388_Enrich_NES_Enrich_Ars_2 +
                             JL1024_Pool_NES_E_Ars_R1 +
                             JL1024_Pool_NES_E_Ars_R2))

snoRNA = snoRNA %>% mutate(Mock_G3BP_Enrich = 
                             (JL0361_Enrich_G3BP_Enrich_Mock_1 + 
                             JL0361_Enrich_G3BP_Enrich_Mock_2 +
                             JL0388_Enrich_SG_Enrich_Mock_1 +
                             JL0388_Enrich_SG_Enrich_Mock_2 +
                             JL1024_Pool_G3BP_E_Mock_R1 +
                             JL1024_Pool_G3BP_E_Mock_R2))

snoRNA = snoRNA %>% mutate(Stress_G3BP_Enrich = 
                             (JL0361_Enrich_G3BP_Enrich_Ars_1 +
                             JL0361_Enrich_G3BP_Enrich_Ars_2 +
                             JL0361_Enrich_G3BP_Enrich_Ars_3 +
                             JL0388_Enrich_SG_Enrich_Ars_1 +
                             JL0388_Enrich_SG_Enrich_Ars_2 +
                             JL1024_Pool_G3BP_E_Ars_R1 +
                             JL1024_Pool_G3BP_E_Ars_R2))


snoRNA_filtered = snoRNA[, c('seqnames', 'start', 'end', 'width', 'strand', 'gene_id', 'gene_name', 'transcript_name', colnames(snoRNA)[86:97])]



snoRNA_filtered_Mock = snoRNA_filtered[, c('Mock_Total_Input', 'Mock_Nuclear_Fraction', 'Mock_Cytoplasm_Fraction', 'Mock_NLS_Enrich', 'Mock_NES_Enrich', 'Mock_G3BP_Enrich')]
colnames(snoRNA_filtered_Mock) = c("Input", "Nuclear", "Cytoplasm", "NLS", "NES", "G3BP")
rownames(snoRNA_filtered_Mock) = NULL

Group_Order = colnames(snoRNA_filtered_Mock)


snoRNA_filtered_Mock_long = snoRNA_filtered_Mock %>%
  pivot_longer(
    cols = everything(),
    names_to = "Group",
    values_to = "Value" 
  )

snoRNA_filtered_Mock_long$Group = factor(snoRNA_filtered_Mock_long$Group, levels = Group_Order)

GroupColors = c("Input" = "grey", "Nuclear" = "blue", "Cytoplasm" = "green", 
               "NLS" = "skyblue", "NES" = "darkseagreen2", "G3BP" = "salmon")

ggplot(snoRNA_filtered_Mock_long, aes(x = Group, y = Value)) +
  geom_boxplot(notch = T) +
  geom_jitter(aes(color = Group), position = position_jitter(width = 0.1), alpha = 0.2) +
  scale_color_manual(values = GroupColors) +
  labs(title = "Boxplots of Normalized Read Counts in Mock Samples") +
  xlab("Samples") +
  ylab("Depth Normalized Read Counts") +
  scale_y_log10() +
  coord_cartesian(ylim = c(1e-9, 1e-5)) +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


Mock_pValMatrix = matrix(nrow = ncol(snoRNA_filtered_Mock), ncol = ncol(snoRNA_filtered_Mock))
rownames(Mock_pValMatrix) = colnames(Mock_pValMatrix) = names(snoRNA_filtered_Mock)

# Perform pairwise tests and store p-values
for (i in 1:(ncol(Mock_pValMatrix) - 1)) {
  for (j in (i + 1):ncol(Mock_pValMatrix)) {
    i_vector = as.vector(unlist(snoRNA_filtered_Mock[, i]))[!is.na(as.vector(unlist(snoRNA_filtered_Mock[, i])))]
    j_vector = as.vector(unlist(snoRNA_filtered_Mock[, j]))[!is.na(as.vector(unlist(snoRNA_filtered_Mock[, j])))]
    # result = ks.test(i_vector, j_vector)
    result = wilcox.test(i_vector, j_vector)
    
    Mock_pValMatrix[i, j] = result$p.value
    Mock_pValMatrix[j, i] = result$p.value
  }
}

Mock_pValMatrix_FDR = matrix(p.adjust(Mock_pValMatrix, method = 'fdr'), nrow = 6, ncol = 6)
colnames(Mock_pValMatrix_FDR) = colnames(Mock_pValMatrix)
rownames(Mock_pValMatrix_FDR) = rownames(Mock_pValMatrix)


snoRNA_filtered_Stress = snoRNA_filtered[, c('Stress_Total_Input', 'Stress_Nuclear_Fraction', 'Stress_Cytoplasm_Fraction', 'Stress_NLS_Enrich', 'Stress_NES_Enrich', 'Stress_G3BP_Enrich')]
colnames(snoRNA_filtered_Stress) = c("Input", "Nuclear", "Cytoplasm", "NLS", "NES", "G3BP")
rownames(snoRNA_filtered_Stress) = NULL

Group_Order = colnames(snoRNA_filtered_Stress)


snoRNA_filtered_Stress_long = snoRNA_filtered_Stress %>%
  pivot_longer(
    cols = everything(),
    names_to = "Group",
    values_to = "Value" 
  )

snoRNA_filtered_Stress_long$Group = factor(snoRNA_filtered_Stress_long$Group, levels = Group_Order)

GroupColors = c("Input" = "grey", "Nuclear" = "blue", "Cytoplasm" = "green", 
                "NLS" = "skyblue", "NES" = "darkseagreen2", "G3BP" = "salmon")

ggplot(snoRNA_filtered_Stress_long, aes(x = Group, y = Value)) +
  geom_boxplot(notch = T) +
  geom_jitter(aes(color = Group), position = position_jitter(width = 0.1), alpha = 0.2) +
  scale_color_manual(values = GroupColors) +
  labs(title = "Boxplots of Normalized Read Counts in Stress Samples") +
  xlab("Samples") +
  ylab("Depth Normalized Read Counts") +
  scale_y_log10() +
  coord_cartesian(ylim = c(1e-9, 1e-5)) +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


Stress_pValMatrix = matrix(nrow = ncol(snoRNA_filtered_Stress), ncol = ncol(snoRNA_filtered_Stress))
rownames(Stress_pValMatrix) = colnames(Stress_pValMatrix) = names(snoRNA_filtered_Stress)

# Perform pairwise tests and store p-values
for (i in 1:(ncol(Stress_pValMatrix) - 1)) {
  for (j in (i + 1):ncol(Stress_pValMatrix)) {
    i_vector = as.vector(unlist(snoRNA_filtered_Stress[, i]))[!is.na(as.vector(unlist(snoRNA_filtered_Stress[, i])))]
    j_vector = as.vector(unlist(snoRNA_filtered_Stress[, j]))[!is.na(as.vector(unlist(snoRNA_filtered_Stress[, j])))]
    # result = ks.test(i_vector, j_vector)
    result = wilcox.test(i_vector, j_vector)
    
    Stress_pValMatrix[i, j] = result$p.value
    Stress_pValMatrix[j, i] = result$p.value
  }
}

Stress_pValMatrix_FDR = matrix(p.adjust(Stress_pValMatrix, method = 'fdr'), nrow = 6, ncol = 6)
colnames(Stress_pValMatrix_FDR) = colnames(Stress_pValMatrix)
rownames(Stress_pValMatrix_FDR) = rownames(Stress_pValMatrix)



