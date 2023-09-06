
library(stringr)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(corrplot)
library(pheatmap)

## Load peak matrix and clean up:
####################################################################################################################
peaksMatrix_PATH = 'L:/.shortcut-targets-by-id/13hY9t_p6eUdnvP-c2OClHbzyeWMOruD_/CoCLIP_HuR_Paper/Data/'    ## Use this for windows machine
# peaksMatrix_PATH = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/'
peaksMatrix_FILE = 'Combined_peakCoverage_groomed_normalized_annotated.txt'

peaksMatrix = read_delim(paste0(peaksMatrix_PATH, peaksMatrix_FILE), show_col_types = FALSE)
peaksMatrix = peaksMatrix %>% mutate_at('TOTAL_BC', as.numeric)
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'unannotated', 'UnAn', finalized_annotation))
peaksMatrix = peaksMatrix %>% mutate(grouped_annotation = ifelse(grouped_annotation == 'unannotated', 'UnAn', grouped_annotation))
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'downstream 10K', 'DS10K', finalized_annotation))
peaksMatrix = peaksMatrix %>% mutate(grouped_annotation = ifelse(grouped_annotation == 'downstream 10K', 'DS10K', grouped_annotation))
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'ncRNA_Retained_intron', 'nC_RI', finalized_annotation))
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'CDS_Retained_intron', 'CDS_RI', finalized_annotation))

## Column organization:
inert_columns = c('chrom', 'start', 'end', 'peak_names', 'score', 'strand', 
                  "gene", "external_gene_name", "annotation", "finalized_annotation", "grouped_annotation", "annotation_count", "TOTAL_TagCount")
BC_columns = c("Nuc_F_M_BC", "Nuc_F_S_BC", "Cyto_F_M_BC", "Cyto_F_S_BC", 
               "NLS_I_M_BC", "NLS_I_S_BC", "NES_I_M_BC", "NES_I_S_BC", "G3BP_I_M_BC", "G3BP_I_S_BC",
               "NLS_E_M_BC", "NLS_E_S_BC", "NES_E_M_BC", "NES_E_S_BC", "G3BP_E_M_BC", "G3BP_E_S_BC",
               "TOTAL_BC")

Nuc_F_M = c('Nuc_F_M_1', 'Nuc_F_M_2', 'Nuc_F_M_3')
Nuc_F_S = c('Nuc_F_S_1', 'Nuc_F_S_2', 'Nuc_F_S_3')
Cyto_F_M = c('Cyto_F_M_1', 'Cyto_F_M_2', 'Cyto_F_M_3')
Cyto_F_S = c('Cyto_F_S_1', 'Cyto_F_S_2', 'Cyto_F_S_3')

NLS_I_M = c('NLS_I_M_1', 'NLS_I_M_2')
NLS_I_S = c('NLS_I_S_1', 'NLS_I_S_2')
NES_I_M = c('NES_I_M_1', 'NES_I_M_2')
NES_I_S = c('NES_I_S_1', 'NES_I_S_2')
G3BP_I_M = c('G3BP_I_M_1', 'G3BP_I_M_2', 'G3BP_I_M_3', 'G3BP_I_M_4')
G3BP_I_S = c('G3BP_I_S_1', 'G3BP_I_S_2', 'G3BP_I_S_3', 'G3BP_I_S_4', 'G3BP_I_S_5')

NLS_E_M = c('NLS_E_M_1', 'NLS_E_M_2', 'NLS_E_M_3', 'NLS_E_M_4')
NLS_E_S = c('NLS_E_S_1', 'NLS_E_S_2', 'NLS_E_S_3', 'NLS_E_S_4')
NES_E_M = c('NES_E_M_1', 'NES_E_M_2', 'NES_E_M_3', 'NES_E_M_4')
NES_E_S = c('NES_E_S_1', 'NES_E_S_2', 'NES_E_S_3', 'NES_E_S_4')
G3BP_E_M = c('G3BP_E_M_1', 'G3BP_E_M_2', 'G3BP_E_M_3', 'G3BP_E_M_4', 'G3BP_E_M_5', 'G3BP_E_M_6')
G3BP_E_S = c('G3BP_E_S_1', 'G3BP_E_S_2', 'G3BP_E_S_3', 'G3BP_E_S_4', 'G3BP_E_S_5', 'G3BP_E_S_6', 'G3BP_E_S_7')

## Add row sum columns for further filtering:
peaksMatrix$F_rowSum = rowSums(peaksMatrix[, c(Nuc_F_M, Nuc_F_S, Cyto_F_M, Cyto_F_S)])
peaksMatrix$I_rowSum = rowSums(peaksMatrix[, c(NLS_I_M, NLS_I_S, NES_I_M, NES_I_S, G3BP_I_M, G3BP_I_S)])
peaksMatrix$E_rowSum = rowSums(peaksMatrix[, c(NLS_E_M, NLS_E_S, NES_E_M, NES_E_S, G3BP_E_M, G3BP_E_S)])
rowSum_columns = c('F_rowSum', 'I_rowSum', 'E_rowSum')

## Add pseudocount instead of filtering by median. This might be better because we will be filtering by BC down the road anyway.
pseudoCount = min(peaksMatrix[, colnames(peaksMatrix)[6:63]][peaksMatrix[, colnames(peaksMatrix)[6:63]] != 0], na.rm = TRUE)
peaksMatrix[, colnames(peaksMatrix)[6:63]] = peaksMatrix[, colnames(peaksMatrix)[6:63]] + pseudoCount
####################################################################################################################

## Filter Criteria:
####################################################################################################################
BC_Threshold_F = 2
BC_Threshold_I = 1
BC_Threshold_I_SG = 2
BC_Threshold_E = 2
BC_Threshold_E_SG = 3

rowSum_Multiplier_F = 4
rowSum_Multiplier_I = 3
rowSum_Multiplier_E = 2
####################################################################################################################

## Subset of Peaks for downstream analysis:
####################################################################################################################
Peak_F_Nuc_M = (peaksMatrix[, c(inert_columns, Nuc_F_M, BC_columns, rowSum_columns)] 
                %>% filter(Nuc_F_M_BC >= BC_Threshold_F & rowSums(peaksMatrix[, c(Nuc_F_M)]) / length(Nuc_F_M) > median(rowSums(peaksMatrix[, c(Nuc_F_M)]) / length(Nuc_F_M)) * rowSum_Multiplier_F))

Peak_F_Nuc_S = (peaksMatrix[, c(inert_columns, Nuc_F_S, BC_columns, rowSum_columns)] 
                %>% filter(Nuc_F_S_BC >= BC_Threshold_F & rowSums(peaksMatrix[, c(Nuc_F_S)]) / length(Nuc_F_S) > median(rowSums(peaksMatrix[, c(Nuc_F_S)]) / length(Nuc_F_S)) * rowSum_Multiplier_F))

Peak_F_Cyt_M = (peaksMatrix[, c(inert_columns, Cyto_F_M, BC_columns, rowSum_columns)] 
                %>% filter(Cyto_F_M_BC >= BC_Threshold_F & rowSums(peaksMatrix[, c(Cyto_F_M)]) / length(Cyto_F_M) > median(rowSums(peaksMatrix[, c(Cyto_F_M)]) / length(Cyto_F_M)) * rowSum_Multiplier_F))

Peak_F_Cyt_S = (peaksMatrix[, c(inert_columns, Cyto_F_S, BC_columns, rowSum_columns)] 
                %>% filter(Cyto_F_S_BC >= BC_Threshold_F & rowSums(peaksMatrix[, c(Cyto_F_S)]) / length(Cyto_F_S) > median(rowSums(peaksMatrix[, c(Cyto_F_S)]) / length(Cyto_F_S)) * rowSum_Multiplier_F))

Peak_Co_Input_M = (peaksMatrix[, c(inert_columns, NLS_I_M, NES_I_M, G3BP_I_M, BC_columns, rowSum_columns)] 
                   %>% filter(NLS_I_M_BC >= BC_Threshold_I & 
                                NES_I_M_BC >= BC_Threshold_I & 
                                G3BP_I_M_BC >= BC_Threshold_I_SG &
                                rowSums(peaksMatrix[, c(NLS_I_M, NES_I_M, G3BP_I_M)]) / length(c(NLS_I_M, NES_I_M, G3BP_I_M)) > median(rowSums(peaksMatrix[, c(NLS_I_M, NES_I_M, G3BP_I_M)]) / length(c(NLS_I_M, NES_I_M, G3BP_I_M))) * rowSum_Multiplier_F))

Peak_Co_Input_S = (peaksMatrix[, c(inert_columns, NLS_I_S, NES_I_S, G3BP_I_S, BC_columns, rowSum_columns)] 
                   %>% filter(NLS_I_S_BC >= BC_Threshold_I & 
                                NES_I_S_BC >= BC_Threshold_I & 
                                G3BP_I_S_BC >= BC_Threshold_I_SG &
                                rowSums(peaksMatrix[, c(NLS_I_S, NES_I_S, G3BP_I_S)]) / length(c(NLS_I_S, NES_I_S, G3BP_I_S)) > median(rowSums(peaksMatrix[, c(NLS_I_S, NES_I_S, G3BP_I_S)]) / length(c(NLS_I_S, NES_I_S, G3BP_I_S))) * rowSum_Multiplier_F))

Peak_Co_NLS_M = (peaksMatrix[, c(inert_columns, NLS_E_M, BC_columns, rowSum_columns)] 
                 %>% filter(NLS_E_M_BC >= BC_Threshold_E &
                              rowSums(peaksMatrix[, c(NLS_E_M)]) / length(c(NLS_E_M)) > median(rowSums(peaksMatrix[, c(NLS_E_M)]) / length(c(NLS_E_M))) * rowSum_Multiplier_F))

Peak_Co_NLS_S = (peaksMatrix[, c(inert_columns, NLS_E_S, BC_columns, rowSum_columns)] 
                 %>% filter(NLS_E_S_BC >= BC_Threshold_E &
                              rowSums(peaksMatrix[, c(NLS_E_S)]) / length(c(NLS_E_S)) > median(rowSums(peaksMatrix[, c(NLS_E_S)]) / length(c(NLS_E_S))) * rowSum_Multiplier_F))

Peak_Co_NES_M = (peaksMatrix[, c(inert_columns, NES_E_M, BC_columns, rowSum_columns)] 
                 %>% filter(NES_E_M_BC >= BC_Threshold_E &
                              rowSums(peaksMatrix[, c(NES_E_M)]) / length(c(NES_E_M)) > median(rowSums(peaksMatrix[, c(NES_E_M)]) / length(c(NES_E_M))) * rowSum_Multiplier_F))

Peak_Co_NES_S = (peaksMatrix[, c(inert_columns, NES_E_S, BC_columns, rowSum_columns)] 
                 %>% filter(NES_E_S_BC >= BC_Threshold_E &
                              rowSums(peaksMatrix[, c(NES_E_S)]) / length(c(NES_E_S)) > median(rowSums(peaksMatrix[, c(NES_E_S)]) / length(c(NES_E_S))) * rowSum_Multiplier_F))

Peak_Co_G3BP_M = (peaksMatrix[, c(inert_columns, G3BP_E_M, BC_columns, rowSum_columns)] 
                  %>% filter(G3BP_E_M_BC >= BC_Threshold_E_SG &
                               rowSums(peaksMatrix[, c(G3BP_E_M)]) / length(c(G3BP_E_M)) > median(rowSums(peaksMatrix[, c(G3BP_E_M)]) / length(c(G3BP_E_M))) * rowSum_Multiplier_F))

Peak_Co_G3BP_S = (peaksMatrix[, c(inert_columns, G3BP_E_S, BC_columns, rowSum_columns)] 
                  %>% filter(G3BP_E_S_BC >= BC_Threshold_E_SG &
                               rowSums(peaksMatrix[, c(G3BP_E_S)]) / length(c(G3BP_E_S)) > median(rowSums(peaksMatrix[, c(G3BP_E_S)]) / length(c(G3BP_E_S))) * rowSum_Multiplier_F))
# 
# ## Nuclear Fraction Mock tags
# Counts_perGene_F_N_M = Peak_F_Nuc_M %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(Nuc_F_M_1 + Nuc_F_M_2 + Nuc_F_M_3), peakCounts = n())
# Counts_perGene_F_N_M = Counts_perGene_F_N_M %>% mutate(tagDensity = tagCounts/peakCounts)
# 
# ## Nuclear Fraction Stress tags
# Counts_perGene_F_N_S = Peak_F_Nuc_S %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(Nuc_F_S_1 + Nuc_F_S_2 + Nuc_F_S_3), peakCounts = n())
# Counts_perGene_F_N_S = Counts_perGene_F_N_S %>% mutate(tagDensity = tagCounts/peakCounts)
# 
# ## Cytoplasm Fraction Mock tags
# Counts_perGene_F_C_M = Peak_F_Cyt_M %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(Cyto_F_M_1 + Cyto_F_M_2 + Cyto_F_M_3), peakCounts = n())
# Counts_perGene_F_C_M = Counts_perGene_F_C_M %>% mutate(tagDensity = tagCounts/peakCounts)
# 
# ## Cytoplasm Fraction Stress tags
# Counts_perGene_F_C_S = Peak_F_Cyt_S %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(Cyto_F_S_1 + Cyto_F_S_2 + Cyto_F_S_3), peakCounts = n())
# Counts_perGene_F_C_S = Counts_perGene_F_C_S %>% mutate(tagDensity = tagCounts/peakCounts)
# 
# ## CoCLIP Input Mock tags
# Counts_perGene_I_M = Peak_Co_Input_M %>% filter(NLS_I_M_BC >= 2) %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(NLS_I_M_1 + NLS_I_M_2 + NES_I_M_1 + NES_I_M_2 + G3BP_I_M_1 + G3BP_I_M_2 + G3BP_I_M_3 + G3BP_I_M_4), peakCounts = n())
# Counts_perGene_I_M = Counts_perGene_I_M %>% mutate(tagDensity = tagCounts/peakCounts)
# 
# ## CoCLIP Input Stress tags
# Counts_perGene_I_S = Peak_Co_Input_S %>% filter(NLS_I_S_BC >= 2) %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(NLS_I_S_1 + NLS_I_S_2 + NES_I_S_1 + NES_I_S_2 + G3BP_I_S_1 + G3BP_I_S_2 + G3BP_I_S_3 + G3BP_I_S_4), peakCounts = n())
# Counts_perGene_I_S = Counts_perGene_I_S %>% mutate(tagDensity = tagCounts/peakCounts)
# 
# ## CoCLIP NLS Mock tags
# Counts_perGene_E_NLS_M = Peak_Co_NLS_M %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(NLS_E_M_1 + NLS_E_M_2 + NLS_E_M_3 + NLS_E_M_4), peakCounts = n())
# Counts_perGene_E_NLS_M = Counts_perGene_E_NLS_M %>% mutate(tagDensity = tagCounts/peakCounts)
# 
# ## CoCLIP NLS Stress tags
# Counts_perGene_E_NLS_S = Peak_Co_NLS_S %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(NLS_E_S_1 + NLS_E_S_2 + NLS_E_S_3 + NLS_E_S_4), peakCounts = n())
# Counts_perGene_E_NLS_S = Counts_perGene_E_NLS_S %>% mutate(tagDensity = tagCounts/peakCounts)
# 
# ## CoCLIP NES Mock tags
# Counts_perGene_E_NES_M = Peak_Co_NES_M %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(NES_E_M_1 + NES_E_M_2 + NES_E_M_3 + NES_E_M_4), peakCounts = n())
# Counts_perGene_E_NES_M = Counts_perGene_E_NES_M %>% mutate(tagDensity = tagCounts/peakCounts)
# 
# ## CoCLIP NES Stress tags
# Counts_perGene_E_NES_S = Peak_Co_NES_S %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(NES_E_S_1 + NES_E_S_2 + NES_E_S_3 + NES_E_S_4), peakCounts = n())
# Counts_perGene_E_NES_S = Counts_perGene_E_NES_S %>% mutate(tagDensity = tagCounts/peakCounts)
# 
# ## CoCLIP G3BP Mock tags
# Counts_perGene_E_G3BP_M = Peak_Co_G3BP_M %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(G3BP_E_M_1 + G3BP_E_M_2 + G3BP_E_M_3 + G3BP_E_M_4 + G3BP_E_M_5 + G3BP_E_M_6), peakCounts = n())
# Counts_perGene_E_G3BP_M = Counts_perGene_E_G3BP_M %>% mutate(tagDensity = tagCounts/peakCounts)
# 
# ## CoCLIP G3BP Stress tags
# Counts_perGene_E_G3BP_S = Peak_Co_G3BP_S %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(G3BP_E_S_1 + G3BP_E_S_2 + G3BP_E_S_3 + G3BP_E_S_4 + G3BP_E_S_5 + G3BP_E_S_6 + G3BP_E_S_7), peakCounts = n())
# Counts_perGene_E_G3BP_S = Counts_perGene_E_G3BP_S %>% mutate(tagDensity = tagCounts/peakCounts)

####################################################################################################################

## Start Building Enrichment Table:
####################################################################################################################
peakRowSum = peaksMatrix[, c(inert_columns, rowSum_columns)]

peakRowSum = peakRowSum %>% mutate(Nuc_F_M = rowSums(peaksMatrix[, Nuc_F_M])/length(Nuc_F_M) * 1e6)
peakRowSum = peakRowSum %>% mutate(Nuc_F_S = rowSums(peaksMatrix[, Nuc_F_S])/length(Nuc_F_S) * 1e6)
peakRowSum = peakRowSum %>% mutate(Cyto_F_M = rowSums(peaksMatrix[, Cyto_F_M])/length(Cyto_F_M) * 1e6)
peakRowSum = peakRowSum %>% mutate(Cyto_F_S = rowSums(peaksMatrix[, Cyto_F_S])/length(Cyto_F_S) * 1e6)

peakRowSum = peakRowSum %>% mutate(NLS_I_M = rowSums(peaksMatrix[, NLS_I_M])/length(NLS_I_M) * 1e6)
peakRowSum = peakRowSum %>% mutate(NES_I_M = rowSums(peaksMatrix[, NES_I_M])/length(NES_I_M) * 1e6)
peakRowSum = peakRowSum %>% mutate(G3BP_I_M = rowSums(peaksMatrix[, G3BP_I_M])/length(G3BP_I_M) * 1e6)
temp = rowSums(peaksMatrix[, c(NLS_I_M, NES_I_M, G3BP_I_M)])/length(c(NLS_I_M, NES_I_M, G3BP_I_M)) * 1e6
peakRowSum = peakRowSum %>% mutate(I_M = temp)

peakRowSum = peakRowSum %>% mutate(NLS_I_S = rowSums(peaksMatrix[, NLS_I_S])/length(NLS_I_S) * 1e6)
peakRowSum = peakRowSum %>% mutate(NES_I_S = rowSums(peaksMatrix[, NES_I_S])/length(NES_I_S) * 1e6)
peakRowSum = peakRowSum %>% mutate(G3BP_I_S = rowSums(peaksMatrix[, G3BP_I_S])/length(G3BP_I_S) * 1e6)
temp = rowSums(peaksMatrix[, c(NLS_I_S, NES_I_S, G3BP_I_S)])/length(c(NLS_I_S, NES_I_S, G3BP_I_S)) * 1e6
peakRowSum = peakRowSum %>% mutate(I_S = temp)

peakRowSum = peakRowSum %>% mutate(NLS_E_M = rowSums(peaksMatrix[, NLS_E_M])/length(NLS_E_M) * 1e6)
peakRowSum = peakRowSum %>% mutate(NLS_E_S = rowSums(peaksMatrix[, NLS_E_S])/length(NLS_E_S) * 1e6)
peakRowSum = peakRowSum %>% mutate(NES_E_M = rowSums(peaksMatrix[, NES_E_M])/length(NES_E_M) * 1e6)
peakRowSum = peakRowSum %>% mutate(NES_E_S = rowSums(peaksMatrix[, NES_E_S])/length(NES_E_S) * 1e6)
peakRowSum = peakRowSum %>% mutate(G3BP_E_M = rowSums(peaksMatrix[, G3BP_E_M])/length(G3BP_E_M) * 1e6)
peakRowSum = peakRowSum %>% mutate(G3BP_E_S = rowSums(peaksMatrix[, G3BP_E_S])/length(G3BP_E_S) * 1e6)

peakRowSum = cbind(peakRowSum, peaksMatrix[, BC_columns])

# write.table(peakRowSum, paste0(peaksMatrix_PATH, str_replace(peaksMatrix_FILE, ".txt", "_rowSum.txt")), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t', na = "")

peakEnrichment = peaksMatrix[, c(inert_columns, BC_columns, rowSum_columns)]

peakEnrichment = peakEnrichment %>% mutate(NLS_EvI_M = peakRowSum$NLS_E_M / peakRowSum$I_M)
peakEnrichment = peakEnrichment %>% mutate(NES_EvI_M = peakRowSum$NES_E_M / peakRowSum$I_M)
peakEnrichment = peakEnrichment %>% mutate(G3BP_EvI_M = peakRowSum$G3BP_E_M / peakRowSum$I_M)
peakEnrichment = peakEnrichment %>% mutate(NLS_EvI_S = peakRowSum$NLS_E_S / peakRowSum$I_S)
peakEnrichment = peakEnrichment %>% mutate(NES_EvI_S = peakRowSum$NES_E_S / peakRowSum$I_S)
peakEnrichment = peakEnrichment %>% mutate(G3BP_EvI_S = peakRowSum$G3BP_E_S / peakRowSum$I_S)

peakEnrichment = peakEnrichment %>% mutate(Nuc_EvI_M = peakRowSum$Nuc_F_M / peakRowSum$I_M)
peakEnrichment = peakEnrichment %>% mutate(Cyto_EvI_M = peakRowSum$Cyto_F_M / peakRowSum$I_M)
peakEnrichment = peakEnrichment %>% mutate(Nuc_EvI_S = peakRowSum$Nuc_F_S / peakRowSum$I_S)
peakEnrichment = peakEnrichment %>% mutate(Cyto_EvI_S = peakRowSum$Cyto_F_S / peakRowSum$I_S)

peakEnrichment = peakEnrichment %>% mutate(E_NvC_M = peakRowSum$NLS_E_M / peakRowSum$NES_E_M)
peakEnrichment = peakEnrichment %>% mutate(F_NvC_M = peakRowSum$Nuc_F_M / peakRowSum$Cyto_F_M)
peakEnrichment = peakEnrichment %>% mutate(E_NvC_S = peakRowSum$NLS_E_S / peakRowSum$NES_E_S)
peakEnrichment = peakEnrichment %>% mutate(F_NvC_S = peakRowSum$Nuc_F_S / peakRowSum$Cyto_F_S)

peakEnrichment = cbind(peakEnrichment, peakRowSum[, colnames(peakRowSum)[17:34]])

# write.table(peakEnrichment, paste0(peaksMatrix_PATH, str_replace(peaksMatrix_FILE, ".txt", "_Enrichment.txt")), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t', na = "")
####################################################################################################################


## Compartment Comparison
####################################################################################################################
# Mock NLS vs NES:
M_Peaks = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                      ((NLS_I_M_BC >= BC_Threshold_I & NLS_E_M_BC >= BC_Threshold_E) | 
                                         (NES_I_M_BC >= BC_Threshold_I & NES_E_M_BC >= BC_Threshold_E)) &
                                      ((I_rowSum >= median(peakEnrichment$I_rowSum)*rowSum_Multiplier_I) &
                                         (E_rowSum >= median(peakEnrichment$E_rowSum)*rowSum_Multiplier_E)))

plotData = data.frame(NLS_EvI_M = M_Peaks$NLS_EvI_M, NES_EvI_M = M_Peaks$NES_EvI_M)
plotData$annotation = factor(M_Peaks$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(plotData, aes(x = log2(NLS_EvI_M), y = log2(NES_EvI_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Cytoplasm Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "red") +
  ggtitle(paste0('Mock HuR Peaks: Nuclear vs Cytoplasm of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Per Gene Peaks so that we can filter:
M_Peaks_per_gene = M_Peaks %>% group_by(gene, external_gene_name) %>% summarise(peakCounts = n())

## Peaks per location:
M_Peaks_NLS = M_Peaks %>% filter((NLS_EvI_M >= 2 & NES_EvI_M < 2) & 
                                   (NLS_E_M_BC >= BC_Threshold_E))

M_Peaks_NES = M_Peaks %>% filter((NES_EvI_M >= 2 & NLS_EvI_M < 2) &
                                   (NES_E_M_BC >= BC_Threshold_E))

M_Peaks_Shared = M_Peaks %>% filter((NLS_EvI_M >= 2 & NES_EvI_M >= 2) &
                                      (NLS_E_M_BC >= BC_Threshold_E & NES_E_M_BC >= BC_Threshold_E))

## Gene Level Grouping:
M_Genes_NLS = unique(M_Peaks_NLS$gene)
M_Genes_NES = unique(M_Peaks_NES$gene)

M_Peaks_NLS_UniqueGenes = M_Peaks_NLS %>% filter(!gene %in% c(M_Genes_NES))
M_Peaks_NES_UniqueGenes = M_Peaks_NES %>% filter(!gene %in% c(M_Genes_NLS))

M_Peaks_SharedGenes = rbind(M_Peaks_NLS %>% filter(gene %in% c(M_Genes_NES)),
                            M_Peaks_NES %>% filter(gene %in% c(M_Genes_NLS)),
                            M_Peaks_Shared)

M_Peaks_NLS_UniqueGenes = M_Peaks_NLS_UniqueGenes[, c(inert_columns[1:length(inert_columns)-1], "NLS_EvI_M", "NES_EvI_M", "NLS_E_M", "NLS_I_M", "NES_E_M", "NES_I_M", "NLS_E_M_BC", "NLS_I_M_BC", "NES_E_M_BC", "NES_I_M_BC")]
M_Peaks_NES_UniqueGenes = M_Peaks_NES_UniqueGenes[, c(inert_columns[1:length(inert_columns)-1], "NLS_EvI_M", "NES_EvI_M", "NLS_E_M", "NLS_I_M", "NES_E_M", "NES_I_M", "NLS_E_M_BC", "NLS_I_M_BC", "NES_E_M_BC", "NES_I_M_BC")]
M_Peaks_SharedGenes = M_Peaks_SharedGenes[, c(inert_columns[1:length(inert_columns)-1], "NLS_EvI_M", "NES_EvI_M", "NLS_E_M", "NLS_I_M", "NES_E_M", "NES_I_M", "NLS_E_M_BC", "NLS_I_M_BC", "NES_E_M_BC", "NES_I_M_BC")]

## These can be filtered again by the Peaks per gene table


# Stress NLS vs NES:
S_Peaks = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                      ((NLS_I_S_BC >= BC_Threshold_I & NLS_E_S_BC >= BC_Threshold_E) | 
                                         (NES_I_S_BC >= BC_Threshold_I & NES_E_S_BC >= BC_Threshold_E)) &
                                      ((I_rowSum >= median(peakEnrichment$I_rowSum)*rowSum_Multiplier_I) &
                                         (E_rowSum >= median(peakEnrichment$E_rowSum)*rowSum_Multiplier_E)))

plotData = data.frame(NLS_EvI_S = S_Peaks$NLS_EvI_S, NES_EvI_S = S_Peaks$NES_EvI_S)
plotData$annotation = factor(S_Peaks$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(plotData, aes(x = log2(NLS_EvI_S), y = log2(NES_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Cytoplasm Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "red") +
  ggtitle(paste0('Stress HuR Peaks: Nuclear vs Cytoplasm of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Per Gene Peaks so that we can filter:
S_Peaks_per_gene = S_Peaks %>% group_by(gene, external_gene_name) %>% summarise(peakCounts = n())


## Peaks per location:
S_Peaks_NLS = S_Peaks %>% filter((NLS_EvI_S >= 2 & NES_EvI_S < 2) &
                                   (NLS_E_S_BC >= BC_Threshold_E))

S_Peaks_NES = S_Peaks %>% filter((NES_EvI_S >= 2 & NLS_EvI_S < 2) &
                                   (NES_E_S_BC >= BC_Threshold_E))

S_Peaks_Shared = S_Peaks %>% filter((NLS_EvI_S >= 2 & NES_EvI_S >= 2) &
                                      (NLS_E_S_BC >= BC_Threshold_E & NES_E_S_BC >= BC_Threshold_E))

## Gene Level Grouping:
S_Genes_NLS = unique(S_Peaks_NLS$gene)
S_Genes_NES = unique(S_Peaks_NES$gene)

S_Peaks_NLS_UniqueGenes = S_Peaks_NLS %>% filter(!gene %in% c(S_Genes_NES))
S_Peaks_NES_UniqueGenes = S_Peaks_NES %>% filter(!gene %in% c(S_Genes_NLS))

S_Peaks_SharedGenes = rbind(S_Peaks_NLS %>% filter(gene %in% c(S_Genes_NES)),
                            S_Peaks_NES %>% filter(gene %in% c(S_Genes_NLS)),
                            S_Peaks_Shared)

S_Peaks_NLS_UniqueGenes = S_Peaks_NLS_UniqueGenes[, c(inert_columns[1:length(inert_columns)-1], "NLS_EvI_S", "NES_EvI_S", "NLS_E_S", "NLS_I_S", "NES_E_S", "NES_I_S", "NLS_E_S_BC", "NLS_I_S_BC", "NES_E_S_BC", "NES_I_S_BC")]
S_Peaks_NES_UniqueGenes = S_Peaks_NES_UniqueGenes[, c(inert_columns[1:length(inert_columns)-1], "NLS_EvI_S", "NES_EvI_S", "NLS_E_S", "NLS_I_S", "NES_E_S", "NES_I_S", "NLS_E_S_BC", "NLS_I_S_BC", "NES_E_S_BC", "NES_I_S_BC")]
S_Peaks_SharedGenes = S_Peaks_SharedGenes[, c(inert_columns[1:length(inert_columns)-1], "NLS_EvI_S", "NES_EvI_S", "NLS_E_S", "NLS_I_S", "NES_E_S", "NES_I_S", "NLS_E_S_BC", "NLS_I_S_BC", "NES_E_S_BC", "NES_I_S_BC")]

## These can be filtered again by the Peaks per gene table


## Mock to Stress Change at Gene Level:
Genes_M2S_NLS_NLS = intersect(M_Peaks_NLS_UniqueGenes$gene, S_Peaks_NLS_UniqueGenes$gene)
Genes_M2S_NLS_NES = intersect(M_Peaks_NLS_UniqueGenes$gene, S_Peaks_NES_UniqueGenes$gene)

Genes_M2S_NES_NES = intersect(M_Peaks_NES_UniqueGenes$gene, S_Peaks_NES_UniqueGenes$gene)
Genes_M2S_NES_NLS = intersect(M_Peaks_NES_UniqueGenes$gene, S_Peaks_NLS_UniqueGenes$gene)






# Mock NLS vs G3BP:
NLS_G3BP_Peaks_Filtered = peakEnrichment %>% filter(grouped_annotation != 'UnAn' & ((NLS_I_M_BC >= BC_Threshold_I & 
                                                                                       NLS_E_M_BC >= BC_Threshold_E) | 
                                                                                      (G3BP_I_M_BC >= BC_Threshold_I_SG & 
                                                                                         G3BP_E_M_BC >= BC_Threshold_E_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_G3BP_Peaks_Filtered$NLS_EvI_M, G3BP_EvI_M = NLS_G3BP_Peaks_Filtered$G3BP_EvI_M)
data$annotation = factor(NLS_G3BP_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(NLS_EvI_M), y = log2(G3BP_EvI_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Stress Granule Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: Nuclear vs Stress Granule of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Mock NLS vs G3BP - mRNA specific:
NLS_G3BP_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                       finalized_annotation == "3'UTR" | 
                                                       finalized_annotation == "CDS" | 
                                                       finalized_annotation == "intron" | 
                                                       finalized_annotation == "CDS_RI" |
                                                       finalized_annotation == "DS10K") &
                                                      ((NLS_I_M_BC >= BC_Threshold_I & 
                                                          NLS_E_M_BC >= BC_Threshold_E) | 
                                                         (G3BP_I_M_BC >= BC_Threshold_I_SG & 
                                                            G3BP_E_M_BC >= BC_Threshold_E_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_G3BP_Peaks_Filtered$NLS_EvI_M, G3BP_EvI_M = NLS_G3BP_Peaks_Filtered$G3BP_EvI_M)
data$annotation = factor(NLS_G3BP_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(NLS_EvI_M), y = log2(G3BP_EvI_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Stress Granule Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: Nuclear vs Stress Granule of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Mock NLS vs G3BP - non-mRNA specific:
NLS_G3BP_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                       finalized_annotation != "3'UTR" & 
                                                       finalized_annotation != "CDS" & 
                                                       finalized_annotation != "intron" &  
                                                       finalized_annotation != "CDS_RI" &
                                                       finalized_annotation != "DS10K" &
                                                       finalized_annotation != "UnAn") & 
                                                      ((NLS_I_M_BC >= BC_Threshold_I & 
                                                          NLS_E_M_BC >= BC_Threshold_E) | 
                                                         (G3BP_I_M_BC >= BC_Threshold_I_SG & 
                                                            G3BP_E_M_BC >= BC_Threshold_E_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_G3BP_Peaks_Filtered$NLS_EvI_M, G3BP_EvI_M = NLS_G3BP_Peaks_Filtered$G3BP_EvI_M)
data$annotation = factor(NLS_G3BP_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(NLS_EvI_M), y = log2(G3BP_EvI_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Stress Granule Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: Nuclear vs Stress Granule of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


# Mock NES vs G3BP:
NES_G3BP_Peaks_Filtered = peakEnrichment %>% filter(grouped_annotation != 'UnAn' & ((NES_I_M_BC >= BC_Threshold_I & 
                                                                                       NES_E_M_BC >= BC_Threshold_E) | 
                                                                                      (G3BP_I_M_BC >= BC_Threshold_I_SG & 
                                                                                         G3BP_E_M_BC >= BC_Threshold_E_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_M = NES_G3BP_Peaks_Filtered$NES_EvI_M, G3BP_EvI_M = NES_G3BP_Peaks_Filtered$G3BP_EvI_M)
data$annotation = factor(NES_G3BP_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(NES_EvI_M), y = log2(G3BP_EvI_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Cytoplasm Enriched/Input)', y = 'log2(Stress Granule Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: Cytoplasm vs Stress Granule of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Mock NES vs G3BP - mRNA specific:
NES_G3BP_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                       finalized_annotation == "3'UTR" | 
                                                       finalized_annotation == "CDS" | 
                                                       finalized_annotation == "intron" | 
                                                       finalized_annotation == "CDS_RI" |
                                                       finalized_annotation == "DS10K") &
                                                      ((NES_I_M_BC >= BC_Threshold_I & 
                                                          NES_E_M_BC >= BC_Threshold_E) | 
                                                         (G3BP_I_M_BC >= BC_Threshold_I_SG & 
                                                            G3BP_E_M_BC >= BC_Threshold_E_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_M = NES_G3BP_Peaks_Filtered$NES_EvI_M, G3BP_EvI_M = NES_G3BP_Peaks_Filtered$G3BP_EvI_M)
data$annotation = factor(NES_G3BP_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(NES_EvI_M), y = log2(G3BP_EvI_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Cytoplasm Enriched/Input)', y = 'log2(Stress Granule Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: Cytoplasm vs Stress Granule of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Mock NES vs G3BP - non-mRNA specific:
NES_G3BP_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                       finalized_annotation != "3'UTR" & 
                                                       finalized_annotation != "CDS" & 
                                                       finalized_annotation != "intron" &  
                                                       finalized_annotation != "CDS_RI" &
                                                       finalized_annotation != "DS10K" &
                                                       finalized_annotation != "UnAn") & 
                                                      ((NES_I_M_BC >= BC_Threshold_I & 
                                                          NES_E_M_BC >= BC_Threshold_E) | 
                                                         (G3BP_I_M_BC >= BC_Threshold_I_SG & 
                                                            G3BP_E_M_BC >= BC_Threshold_E_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_M = NES_G3BP_Peaks_Filtered$NES_EvI_M, G3BP_EvI_M = NES_G3BP_Peaks_Filtered$G3BP_EvI_M)
data$annotation = factor(NES_G3BP_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(NES_EvI_M), y = log2(G3BP_EvI_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Cytoplasm Enriched/Input)', y = 'log2(Stress Granule Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: Cytoplasm vs Stress Granule of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress NLS vs NES:
NLS_NES_Peaks_Filtered = peakEnrichment %>% filter(grouped_annotation != 'UnAn' & ((NLS_I_S_BC >= BC_Threshold_I & 
                                                                                      NLS_E_S_BC >= BC_Threshold_E) | 
                                                                                     (NES_I_S_BC >= BC_Threshold_I & 
                                                                                        NES_E_S_BC >= BC_Threshold_E)) &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                     E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_S = NLS_NES_Peaks_Filtered$NLS_EvI_S, NES_EvI_S = NLS_NES_Peaks_Filtered$NES_EvI_S)
data$annotation = factor(NLS_NES_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(NLS_EvI_S), y = log2(NES_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Cytoplasm Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "red") +
  ggtitle(paste0('Stress HuR Peaks: Nuclear vs Cytoplasm of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress NLS vs NES - mRNA specific:
NLS_NES_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                      finalized_annotation == "3'UTR" | 
                                                      finalized_annotation == "CDS" | 
                                                      finalized_annotation == "intron" | 
                                                      finalized_annotation == "CDS_RI" |
                                                      finalized_annotation == "DS10K") &
                                                     ((NLS_I_S_BC >= BC_Threshold_I & 
                                                         NLS_E_S_BC >= BC_Threshold_E) | 
                                                        (NES_I_S_BC >= BC_Threshold_I & 
                                                           NES_E_S_BC >= BC_Threshold_E)) &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                     E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_S = NLS_NES_Peaks_Filtered$NLS_EvI_S, NES_EvI_S = NLS_NES_Peaks_Filtered$NES_EvI_S)
data$annotation = factor(NLS_NES_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(NLS_EvI_S), y = log2(NES_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Cytoplasm Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Stress HuR Peaks: Nuclear vs Cytoplasm of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress NLS vs NES - non-mRNA specific:
NLS_NES_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                      finalized_annotation != "3'UTR" & 
                                                      finalized_annotation != "CDS" & 
                                                      finalized_annotation != "intron" &  
                                                      finalized_annotation != "CDS_RI" &
                                                      finalized_annotation != "DS10K" &
                                                      finalized_annotation != "UnAn") & 
                                                     ((NLS_I_S_BC >= BC_Threshold_I & 
                                                         NLS_E_S_BC >= BC_Threshold_E) | 
                                                        (NES_I_S_BC >= BC_Threshold_I & 
                                                           NES_E_S_BC >= BC_Threshold_E)) &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                     E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_S = NLS_NES_Peaks_Filtered$NLS_EvI_S, NES_EvI_S = NLS_NES_Peaks_Filtered$NES_EvI_S)
data$annotation = factor(NLS_NES_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(NLS_EvI_S), y = log2(NES_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Cytoplasm Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Stress HuR Peaks: Nuclear vs Cytoplasm of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress NLS vs G3BP:
NLS_G3BP_Peaks_Filtered = peakEnrichment %>% filter(grouped_annotation != 'UnAn' & ((NLS_I_S_BC >= BC_Threshold_I & 
                                                                                       NLS_E_S_BC >= BC_Threshold_E) | 
                                                                                      (G3BP_I_S_BC >= BC_Threshold_I_SG & 
                                                                                         G3BP_E_S_BC >= BC_Threshold_E_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_S = NLS_G3BP_Peaks_Filtered$NLS_EvI_S, G3BP_EvI_S = NLS_G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(NLS_G3BP_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(NLS_EvI_S), y = log2(G3BP_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Stress Granule Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Stress HuR Peaks: Nuclear vs Stress Granule of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress NLS vs G3BP - mRNA specific:
NLS_G3BP_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                       finalized_annotation == "3'UTR" | 
                                                       finalized_annotation == "CDS" | 
                                                       finalized_annotation == "intron" | 
                                                       finalized_annotation == "CDS_RI" |
                                                       finalized_annotation == "DS10K") &
                                                      ((NLS_I_S_BC >= BC_Threshold_I & 
                                                          NLS_E_S_BC >= BC_Threshold_E) | 
                                                         (G3BP_I_S_BC >= BC_Threshold_I_SG & 
                                                            G3BP_E_S_BC >= BC_Threshold_E_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_S = NLS_G3BP_Peaks_Filtered$NLS_EvI_S, G3BP_EvI_S = NLS_G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(NLS_G3BP_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(NLS_EvI_S), y = log2(G3BP_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Stress Granule Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Stress HuR Peaks: Nuclear vs Stress Granule of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress NLS vs G3BP - non-mRNA specific:
NLS_G3BP_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                       finalized_annotation != "3'UTR" & 
                                                       finalized_annotation != "CDS" & 
                                                       finalized_annotation != "intron" &  
                                                       finalized_annotation != "CDS_RI" &
                                                       finalized_annotation != "DS10K" &
                                                       finalized_annotation != "UnAn") & 
                                                      ((NLS_I_S_BC >= BC_Threshold_I & 
                                                          NLS_E_S_BC >= BC_Threshold_E) | 
                                                         (G3BP_I_S_BC >= BC_Threshold_I_SG & 
                                                            G3BP_E_S_BC >= BC_Threshold_E_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_S = NLS_G3BP_Peaks_Filtered$NLS_EvI_S, G3BP_EvI_S = NLS_G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(NLS_G3BP_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(NLS_EvI_S), y = log2(G3BP_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Stress Granule Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Stress HuR Peaks: Nuclear vs Stress Granule of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress NES vs G3BP:
NES_G3BP_Peaks_Filtered = peakEnrichment %>% filter(grouped_annotation != 'UnAn' & ((NES_I_S_BC >= BC_Threshold_I & 
                                                                                       NES_E_S_BC >= BC_Threshold_E) | 
                                                                                      (G3BP_I_S_BC >= BC_Threshold_I_SG & 
                                                                                         G3BP_E_S_BC >= BC_Threshold_E_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_S = NES_G3BP_Peaks_Filtered$NES_EvI_S, G3BP_EvI_S = NES_G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(NES_G3BP_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(NES_EvI_S), y = log2(G3BP_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Cytoplasm Enriched/Input)', y = 'log2(Stress Granule Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Stress HuR Peaks: Cytoplasm vs Stress Granule of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress NES vs G3BP - mRNA specific:
NES_G3BP_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                       finalized_annotation == "3'UTR" | 
                                                       finalized_annotation == "CDS" | 
                                                       finalized_annotation == "intron" | 
                                                       finalized_annotation == "CDS_RI" |
                                                       finalized_annotation == "DS10K") &
                                                      ((NES_I_S_BC >= BC_Threshold_I & 
                                                          NES_E_S_BC >= BC_Threshold_E) | 
                                                         (G3BP_I_S_BC >= BC_Threshold_I_SG & 
                                                            G3BP_E_S_BC >= BC_Threshold_E_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_S = NES_G3BP_Peaks_Filtered$NES_EvI_S, G3BP_EvI_S = NES_G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(NES_G3BP_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(NES_EvI_S), y = log2(G3BP_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Cytoplasm Enriched/Input)', y = 'log2(Stress Granule Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Stress HuR Peaks: Cytoplasm vs Stress Granule of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress NES vs G3BP - non-mRNA specific:
NES_G3BP_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                       finalized_annotation != "3'UTR" & 
                                                       finalized_annotation != "CDS" & 
                                                       finalized_annotation != "intron" &  
                                                       finalized_annotation != "CDS_RI" &
                                                       finalized_annotation != "DS10K" &
                                                       finalized_annotation != "UnAn") & 
                                                      ((NES_I_S_BC >= BC_Threshold_I & 
                                                          NES_E_S_BC >= BC_Threshold_E) | 
                                                         (G3BP_I_S_BC >= BC_Threshold_I_SG & 
                                                            G3BP_E_S_BC >= BC_Threshold_E_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_S = NES_G3BP_Peaks_Filtered$NES_EvI_S, G3BP_EvI_S = NES_G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(NES_G3BP_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(NES_EvI_S), y = log2(G3BP_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Cytoplasm Enriched/Input)', y = 'log2(Stress Granule Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Stress HuR Peaks: Cytoplasm vs Stress Granule of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

####################################################################################################################

## Fractionation Compartment Comparison
####################################################################################################################
# Mock Nuclear vs Cytoplasm:
Mock_Peaks_Filtered = peakEnrichment %>% filter(grouped_annotation != 'UnAn' & ((Cyto_F_M_BC >= BC_Threshold_F & 
                                                                                   NES_I_M_BC >= BC_Threshold_I) | 
                                                                                  (Nuc_F_M_BC >= BC_Threshold_F & 
                                                                                     NLS_I_M_BC >= BC_Threshold_I)) &
                                                  (F_rowSum >= median(peakEnrichment$F_rowSum)*rowSum_Multiplier_F &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*rowSum_Multiplier_F))

data = data.frame(Nuc_EvI_M = Mock_Peaks_Filtered$Nuc_EvI_M, Cyto_EvI_M = Mock_Peaks_Filtered$Cyto_EvI_M)
data$annotation = factor(Mock_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(Nuc_EvI_M), y = log2(Cyto_EvI_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Cytoplasm Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: Nuclear vs Cytoplasm of Fraction/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Mock Nuclear vs Cytoplasm - mRNA specific:
Mock_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                   finalized_annotation == "3'UTR" | 
                                                   finalized_annotation == "CDS" | 
                                                   finalized_annotation == "intron" | 
                                                   finalized_annotation == "CDS_RI" |
                                                   finalized_annotation == "DS10K") &
                                                  ((Cyto_F_M_BC >= BC_Threshold_F & 
                                                      NES_I_M_BC >= BC_Threshold_I) | 
                                                     (Nuc_F_M_BC >= BC_Threshold_F & 
                                                        NLS_I_M_BC >= BC_Threshold_I)) &
                                                  (F_rowSum >= median(peakEnrichment$F_rowSum)*rowSum_Multiplier_F &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*rowSum_Multiplier_F))

data = data.frame(Nuc_EvI_M = Mock_Peaks_Filtered$Nuc_EvI_M, Cyto_EvI_M = Mock_Peaks_Filtered$Cyto_EvI_M)
data$annotation = factor(Mock_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(Nuc_EvI_M), y = log2(Cyto_EvI_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Cytoplasm Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: Nuclear vs Cytoplasm of Fraction/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Mock Nuclear vs Cytoplasm - non-mRNA specific:
Mock_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                   finalized_annotation != "3'UTR" & 
                                                   finalized_annotation != "CDS" & 
                                                   finalized_annotation != "intron" &  
                                                   finalized_annotation != "CDS_RI" &
                                                   finalized_annotation != "DS10K" &
                                                   finalized_annotation != "UnAn") & 
                                                  ((Cyto_F_M_BC >= BC_Threshold_F & 
                                                      NES_I_M_BC >= BC_Threshold_I) | 
                                                     (Nuc_F_M_BC >= BC_Threshold_F & 
                                                        NLS_I_M_BC >= BC_Threshold_I)) &
                                                  (F_rowSum >= median(peakEnrichment$F_rowSum)*rowSum_Multiplier_F &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*rowSum_Multiplier_F))

data = data.frame(Nuc_EvI_M = Mock_Peaks_Filtered$Nuc_EvI_M, Cyto_EvI_M = Mock_Peaks_Filtered$Cyto_EvI_M)
data$annotation = factor(Mock_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(Nuc_EvI_M), y = log2(Cyto_EvI_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Cytoplasm Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: Nuclear vs Cytoplasm of Fraction/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress Nuclear vs Cytoplasm:
Stress_Peaks_Filtered = peakEnrichment %>% filter(grouped_annotation != 'UnAn' & ((Cyto_F_S_BC >= BC_Threshold_F & 
                                                                                     NES_I_S_BC >= BC_Threshold_I) | 
                                                                                    (Nuc_F_S_BC >= BC_Threshold_F & 
                                                                                       NLS_I_S_BC >= BC_Threshold_I)) &
                                                    (F_rowSum >= median(peakEnrichment$F_rowSum)*rowSum_Multiplier_F &
                                                       I_rowSum >= median(peakEnrichment$I_rowSum)*rowSum_Multiplier_F))

data = data.frame(Nuc_EvI_S = Stress_Peaks_Filtered$Nuc_EvI_S, Cyto_EvI_S = Stress_Peaks_Filtered$Cyto_EvI_S)
data$annotation = factor(Stress_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(Nuc_EvI_S), y = log2(Cyto_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Cytoplasm Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Stress HuR Peaks: Nuclear vs Cytoplasm of Fraction/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress Nuclear vs Cytoplasm - mRNA specific:
Stress_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                     finalized_annotation == "3'UTR" | 
                                                     finalized_annotation == "CDS" | 
                                                     finalized_annotation == "intron" | 
                                                     finalized_annotation == "CDS_RI" |
                                                     finalized_annotation == "DS10K") &
                                                    ((Cyto_F_S_BC >= BC_Threshold_F & 
                                                        NES_I_S_BC >= BC_Threshold_I) | 
                                                       (Nuc_F_S_BC >= BC_Threshold_F & 
                                                          NLS_I_S_BC >= BC_Threshold_I)) &
                                                    (F_rowSum >= median(peakEnrichment$F_rowSum)*rowSum_Multiplier_F &
                                                       I_rowSum >= median(peakEnrichment$I_rowSum)*rowSum_Multiplier_F))

data = data.frame(Nuc_EvI_S = Stress_Peaks_Filtered$Nuc_EvI_S, Cyto_EvI_S = Stress_Peaks_Filtered$Cyto_EvI_S)
data$annotation = factor(Stress_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(Nuc_EvI_S), y = log2(Cyto_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Cytoplasm Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Stress HuR Peaks: Nuclear vs Cytoplasm of Fraction/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress Nuclear vs Cytoplasm - non-mRNA specific:
Stress_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                     finalized_annotation != "3'UTR" & 
                                                     finalized_annotation != "CDS" & 
                                                     finalized_annotation != "intron" &  
                                                     finalized_annotation != "CDS_RI" &
                                                     finalized_annotation != "DS10K" &
                                                     finalized_annotation != "UnAn") & 
                                                    ((Cyto_F_S_BC >= BC_Threshold_F & 
                                                        NES_I_S_BC >= BC_Threshold_I) | 
                                                       (Nuc_F_S_BC >= BC_Threshold_F & 
                                                          NLS_I_S_BC >= BC_Threshold_I)) &
                                                    (F_rowSum >= median(peakEnrichment$F_rowSum)*rowSum_Multiplier_F &
                                                       I_rowSum >= median(peakEnrichment$I_rowSum)*rowSum_Multiplier_F))

data = data.frame(Nuc_EvI_S = Stress_Peaks_Filtered$Nuc_EvI_S, Cyto_EvI_S = Stress_Peaks_Filtered$Cyto_EvI_S)
data$annotation = factor(Stress_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(Nuc_EvI_S), y = log2(Cyto_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Cytoplasm Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Stress HuR Peaks: Nuclear vs Cytoplasm of Fraction/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

####################################################################################################################