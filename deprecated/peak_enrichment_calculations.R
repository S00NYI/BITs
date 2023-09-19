## CoCLIP Analysis: 
## Peak Enrichment Calculation
## Written by Soon Yi
## Last Edit: 2023-09-04

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
# peaksMatrix_PATH = 'L:/My Drive/CWRU/PhD/Luna Lab/1. coCLIP/Analysis/peaks/'    ## Use this for windows machine
peaksMatrix_PATH = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/'
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

## Make an peak density table:
## Actually this is just peak row sums by sample, not "density"
avgTagCounts = peaksMatrix[, c(inert_columns, BC_columns, rowSum_columns)]
avgTagCounts = avgTagCounts %>% mutate(Nuc_F_M = rowSums(peaksMatrix[, c(Nuc_F_M)]) / length(Nuc_F_M))
avgTagCounts = avgTagCounts %>% mutate(Nuc_F_S = rowSums(peaksMatrix[, c(Nuc_F_S)]) / length(Nuc_F_S))
avgTagCounts = avgTagCounts %>% mutate(Cyto_F_M = rowSums(peaksMatrix[, c(Cyto_F_M)]) / length(Cyto_F_M))
avgTagCounts = avgTagCounts %>% mutate(Cyto_F_S = rowSums(peaksMatrix[, c(Cyto_F_S)]) / length(Cyto_F_S))

avgTagCounts = avgTagCounts %>% mutate(NLS_I_M = rowSums(peaksMatrix[, c(NLS_I_M)]) / length(NLS_I_M))
avgTagCounts = avgTagCounts %>% mutate(NLS_I_S = rowSums(peaksMatrix[, c(NLS_I_S)]) / length(NLS_I_S))
avgTagCounts = avgTagCounts %>% mutate(NES_I_M = rowSums(peaksMatrix[, c(NES_I_M)]) / length(NES_I_M))
avgTagCounts = avgTagCounts %>% mutate(NES_I_S = rowSums(peaksMatrix[, c(NES_I_S)]) / length(NES_I_S))
avgTagCounts = avgTagCounts %>% mutate(G3BP_I_M = rowSums(peaksMatrix[, c(G3BP_I_M)]) / length(G3BP_I_M))
avgTagCounts = avgTagCounts %>% mutate(G3BP_I_S = rowSums(peaksMatrix[, c(G3BP_I_S)]) / length(G3BP_I_S))

avgTagCounts = avgTagCounts %>% mutate(NLS_E_M = rowSums(peaksMatrix[, c(NLS_E_M)]) / length(NLS_E_M))
avgTagCounts = avgTagCounts %>% mutate(NLS_E_S = rowSums(peaksMatrix[, c(NLS_E_S)]) / length(NLS_E_S))
avgTagCounts = avgTagCounts %>% mutate(NES_E_M = rowSums(peaksMatrix[, c(NES_E_M)]) / length(NES_E_M))
avgTagCounts = avgTagCounts %>% mutate(NES_E_S = rowSums(peaksMatrix[, c(NES_E_S)]) / length(NES_E_S))
avgTagCounts = avgTagCounts %>% mutate(G3BP_E_M = rowSums(peaksMatrix[, c(G3BP_E_M)]) / length(G3BP_E_M))
avgTagCounts = avgTagCounts %>% mutate(G3BP_E_S = rowSums(peaksMatrix[, c(G3BP_E_S)]) / length(G3BP_E_S))

## Add pseudocount instead of filtering by median. This might be better because we will be filtering by BC down the road anyway.
# peakCount_median = median(peaksMatrix[, colnames(peaksMatrix)[6:63]][peaksMatrix[, colnames(peaksMatrix)[6:63]] != 0], na.rm = TRUE)
# peaksMatrix = peaksMatrix %>% filter(if_any(all_of(comparison_vector), ~ . > peakCount_median ))
pseudoCount = min(peaksMatrix[, colnames(peaksMatrix)[6:63]][peaksMatrix[, colnames(peaksMatrix)[6:63]] != 0], na.rm = TRUE)
peaksMatrix[, colnames(peaksMatrix)[6:63]] = peaksMatrix[, colnames(peaksMatrix)[6:63]] + pseudoCount

####################################################################################################################

## PCA of Normalized Peaks:
####################################################################################################################
PCA_data = peaksMatrix[, c(Nuc_F_M, Nuc_F_S, Cyto_F_M, Cyto_F_S, NLS_I_M, NLS_I_S, NES_I_M, NES_I_S, G3BP_I_M, G3BP_I_S, NLS_E_M, NLS_E_S, NES_E_M, NES_E_S, G3BP_E_M, G3BP_E_S)]
# PCA_data = scale(PCA_data)

PCA_result = prcomp(t(PCA_data))

num_PCs = 10
PCs = PCA_result$x[, 1:num_PCs]

dataLabel= data.frame(Type = c('F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 
                               'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I',
                               'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 
                               'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E'), 
                      Localization = c('N', 'N', 'N', 'N', 'N', 'N', 'C', 'C', 'C', 'C', 'C', 'C',
                                       'N', 'N', 'N', 'N', 'C', 'C', 'C', 'C', 'SG', 'SG', 'SG', 'SG', 'SG', 'SG', 'SG', 'SG', 'SG',
                                       'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 
                                       'SG', 'SG', 'SG', 'SG', 'SG', 'SG', 'SG', 'SG', 'SG', 'SG', 'SG', 'SG', 'SG'), 
                      Condition = c('M', 'M', 'M', 'S', 'S', 'S', 'M', 'M', 'M', 'S', 'S', 'S',
                                    'M', 'M', 'S', 'S', 'M', 'M', 'S', 'S', 'M', 'M', 'M', 'M', 'S', 'S', 'S', 'S', 'S',
                                    'M', 'M', 'M', 'M', 'S', 'S', 'S', 'S', 'M', 'M', 'M', 'M', 'S', 'S', 'S', 'S',
                                    'M', 'M', 'M', 'M', 'M', 'M', 'S', 'S', 'S', 'S', 'S', 'S', 'S'),
                      Replicate = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 
                                    1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 3, 4, 1, 2, 3, 4, 5,
                                    1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 
                                    1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 7))
rownames(dataLabel) = rownames(PCs)
PCs = cbind(PCs, dataLabel)

color_palette = c("M" = 'green', "S" = 'red')
shape_palette = c("E" = 16, "I" = 17)

# Create the plot
ggplot(PCs %>% filter(Type != 'F'), aes(x = PC1, y = PC2, color = Condition, shape = Type)) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  labs(title = "PCA Plot of PC1 vs PC2", x = "PC1", y = "PC2",
       color = "Localization",
       shape = "Type") +
  theme_minimal()

####################################################################################################################

## Correlation plot of Normalized Peaks:
####################################################################################################################
CorrMatrix = cor(PCA_data)

CorrMatrix = matrix(round(CorrMatrix,2), nrow = 58)
colnames(CorrMatrix) = colnames(PCA_data)
rownames(CorrMatrix) = colnames(PCA_data)
pheatmap(CorrMatrix)
####################################################################################################################

## Filter Criteria:
####################################################################################################################
# peaksMatrix = peaksMatrix %>% filter(annotation_count == 1)

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
# Peak_all_M = (peaksMatrix[, c(inert_columns, Nuc_F_M, Cyto_F_M, NLS_I_M, NES_I_M, G3BP_I_M, NLS_E_M, NES_E_M, G3BP_E_M, BC_columns)] 
#               %>% filter((Nuc_F_M_BC >= BC_Threshold_F & Cyto_F_M_BC >= BC_Threshold_F) | 
#                            (NLS_I_M_BC >= BC_Threshold_I & NES_I_M_BC >= BC_Threshold_I & G3BP_I_M_BC >= BC_Threshold_I_SG) |
#                            (NLS_E_M_BC >= BC_Threshold_E & NES_E_M_BC >= BC_Threshold_E & G3BP_E_M_BC >= BC_Threshold_E_SG)))
# 
# Peak_all_S = (peaksMatrix[, c(inert_columns, Nuc_F_S, Cyto_F_S, NLS_I_S, NES_I_S, G3BP_I_S, NLS_E_S, NES_E_S, G3BP_E_S, BC_columns)] 
#               %>% filter((Nuc_F_S_BC >= BC_Threshold_F & Cyto_F_S_BC >= BC_Threshold_F) |
#                            (NLS_I_S_BC >= BC_Threshold_I & NES_I_S_BC >= BC_Threshold_I & G3BP_I_S_BC >= BC_Threshold_E_SG) | 
#                            (NLS_E_S_BC >= BC_Threshold_E & NES_E_S_BC >= BC_Threshold_E & G3BP_E_S_BC >= BC_Threshold_E_SG)))
# 
# Peak_all = (peaksMatrix[, c(inert_columns, Nuc_F_S, Cyto_F_S, NLS_I_S, NES_I_S, G3BP_I_S, NLS_E_S, NES_E_S, G3BP_E_S, BC_columns)] 
#             %>% filter(((Nuc_F_M_BC >= BC_Threshold_F & Cyto_F_M_BC >= BC_Threshold_F) | 
#                           (NLS_I_M_BC >= BC_Threshold_I & NES_I_M_BC >= BC_Threshold_I & G3BP_I_M_BC >= BC_Threshold_E_SG) | 
#                           (NLS_E_M_BC >= BC_Threshold_E & NES_E_M_BC >= BC_Threshold_E & G3BP_E_M_BC >= BC_Threshold_E_SG)) | 
#                          ((Nuc_F_S_BC >= BC_Threshold_F & Cyto_F_S_BC >= BC_Threshold_F) | 
#                             (NLS_I_S_BC >= BC_Threshold_I & NES_I_S_BC >= BC_Threshold_I & G3BP_I_S_BC >= BC_Threshold_E_SG) | 
#                             NLS_E_S_BC >= BC_Threshold_E & NES_E_S_BC >= BC_Threshold_E & G3BP_E_S_BC >= BC_Threshold_E_SG)))

Peak_F_Nuc_M = (peaksMatrix[, c(inert_columns, Nuc_F_M, BC_columns, rowSum_columns)] 
                %>% filter(Nuc_F_M_BC >= BC_Threshold_F & rowSums(peaksMatrix[, c(Nuc_F_M)]) / length(Nuc_F_M) > median(rowSums(peaksMatrix[, c(Nuc_F_M)]) / length(Nuc_F_M)) * rowSum_Multiplier_F))

Peak_F_Nuc_S = (peaksMatrix[, c(inert_columns, Nuc_F_S, BC_columns, rowSum_columns)] 
                %>% filter(Nuc_F_S_BC >= BC_Threshold_F & rowSums(peaksMatrix[, c(Nuc_F_S)]) / length(Nuc_F_S) > median(rowSums(peaksMatrix[, c(Nuc_F_S)]) / length(Nuc_F_S)) * rowSum_Multiplier_F))

Peak_F_Cyt_M = (peaksMatrix[, c(inert_columns, Cyto_F_M, BC_columns, rowSum_columns)] 
                %>% filter(Cyto_F_M_BC >= BC_Threshold_F & rowSums(peaksMatrix[, c(Cyto_F_M)]) / length(Cyto_F_M) > median(rowSums(peaksMatrix[, c(Cyto_F_M)]) / length(Cyto_F_M)) * rowSum_Multiplier_F))

Peak_F_Cyt_S = (peaksMatrix[, c(inert_columns, Cyto_F_S, BC_columns, rowSum_columns)] 
                %>% filter(Cyto_F_S_BC >= BC_Threshold_F & rowSums(peaksMatrix[, c(Cyto_F_S)]) / length(Cyto_F_S) > median(rowSums(peaksMatrix[, c(Cyto_F_S)]) / length(Cyto_F_S)) * rowSum_Multiplier_F))

# Peak_F_all_M = (peaksMatrix[, c(inert_columns, Nuc_F_M, Cyto_F_M, BC_columns, rowSum_columns)] 
#                 %>% filter(Nuc_F_M_BC >= BC_Threshold_F &
#                            Cyto_F_M_BC >= BC_Threshold_F & F_rowSum > median(F_rowSum) * rowSum_Multiplier_F))
# 
# Peak_F_all_S = (peaksMatrix[, c(inert_columns, Nuc_F_S, Cyto_F_S, BC_columns, rowSum_columns)] 
#                 %>% filter(Nuc_F_S_BC >= BC_Threshold_F & 
#                            Cyto_F_S_BC >= BC_Threshold_F & F_rowSum > median(F_rowSum) * rowSum_Multiplier_F))

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

# These outputs are in peaks format suitable for Homer motif analysis:
# write.table(Peak_all_M[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_all_M.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
# write.table(Peak_all_S[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_all_S.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
# write.table(Peak_all[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_all.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
# 
# write.table(Peak_F_Nuc_M[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_F_Nuc_M.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
# write.table(Peak_F_Nuc_S[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_F_Nuc_S.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
# write.table(Peak_F_Cyt_M[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_F_Cyt_M.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
# write.table(Peak_F_Cyt_S[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_F_Cyt_S.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
# write.table(Peak_F_all_M[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_F_all_M.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
# write.table(Peak_F_all_S[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_F_all_S.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
# 
# write.table(Peak_Co_Input_M[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_Co_Input_M.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
# write.table(Peak_Co_Input_S[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_Co_Input_S.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
# write.table(Peak_Co_NLS_M[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_Co_NLS_M.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
# write.table(Peak_Co_NLS_S[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_Co_NLS_S.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
# write.table(Peak_Co_NES_M[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_Co_NES_M.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
# write.table(Peak_Co_NES_S[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_Co_NES_S.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
# write.table(Peak_Co_G3BP_M[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_Co_G3BP_M.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
# write.table(Peak_Co_G3BP_S[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_Co_G3BP_S.txt'), quote = F, sep = '\t', row.names = F, col.names = F)

## per Gene information:
# Counts_perGene = Peak_all %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(TOTAL_TagCount), peakCounts = n())
# Counts_perGene = Counts_perGene %>% mutate(tagDensity = tagCounts/peakCounts)
# 
# ## All Mock tags
# Counts_perGene_A_M= Peak_all_M %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(Nuc_F_M_1 + Nuc_F_M_2 + Nuc_F_M_3 + 
#                                                                                                       Cyto_F_M_1 + Cyto_F_M_2 + Cyto_F_M_3 + 
#                                                                                                       NLS_I_M_1 + NLS_I_M_2 + NES_I_M_1 + NES_I_M_2 + G3BP_I_M_1 + G3BP_I_M_2 + G3BP_I_M_3 + G3BP_I_M_4 +
#                                                                                                       NLS_E_M_1 + NLS_E_M_2 + NLS_E_M_3 + NLS_E_M_4 + 
#                                                                                                       NES_E_M_1 + NES_E_M_2 + NES_E_M_3 + NES_E_M_4 +
#                                                                                                       G3BP_E_M_1 + G3BP_E_M_2 + G3BP_E_M_3 + G3BP_E_M_4 + G3BP_E_M_5 + G3BP_E_M_6), peakCounts = n())
# Counts_perGene_A_M = Counts_perGene_A_M %>% mutate(tagDensity = tagCounts/peakCounts)
# 
# ## All Stress tags
# Counts_perGene_A_S= Peak_all_S %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(Nuc_F_S_1 + Nuc_F_S_2 + Nuc_F_S_3 +
#                                                                                                       Cyto_F_S_1 + Cyto_F_S_2 + Cyto_F_S_3 +
#                                                                                                       NLS_I_S_1 + NLS_I_S_2 + NES_I_S_1 + NES_I_S_2 + G3BP_I_S_1 + G3BP_I_S_2 + G3BP_I_S_3 + G3BP_I_S_4 +
#                                                                                                       NLS_E_S_1 + NLS_E_S_2 + NLS_E_S_3 + NLS_E_S_4 +
#                                                                                                       NES_E_S_1 + NES_E_S_2 + NES_E_S_3 + NES_E_S_4 +
#                                                                                                       G3BP_E_S_1 + G3BP_E_S_2 + G3BP_E_S_3 + G3BP_E_S_4 + G3BP_E_S_5 + G3BP_E_S_6 + G3BP_E_S_7), peakCounts = n())
# Counts_perGene_A_S = Counts_perGene_A_S %>% mutate(tagDensity = tagCounts/peakCounts)

## Nuclear Fraction Mock tags
Counts_perGene_F_N_M = Peak_F_Nuc_M %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(Nuc_F_M_1 + Nuc_F_M_2 + Nuc_F_M_3), peakCounts = n())
Counts_perGene_F_N_M = Counts_perGene_F_N_M %>% mutate(tagDensity = tagCounts/peakCounts)

## Nuclear Fraction Stress tags
Counts_perGene_F_N_S = Peak_F_Nuc_S %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(Nuc_F_S_1 + Nuc_F_S_2 + Nuc_F_S_3), peakCounts = n())
Counts_perGene_F_N_S = Counts_perGene_F_N_S %>% mutate(tagDensity = tagCounts/peakCounts)

## Cytoplasm Fraction Mock tags
Counts_perGene_F_C_M = Peak_F_Cyt_M %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(Cyto_F_M_1 + Cyto_F_M_2 + Cyto_F_M_3), peakCounts = n())
Counts_perGene_F_C_M = Counts_perGene_F_C_M %>% mutate(tagDensity = tagCounts/peakCounts)

## Cytoplasm Fraction Stress tags
Counts_perGene_F_C_S = Peak_F_Cyt_S %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(Cyto_F_S_1 + Cyto_F_S_2 + Cyto_F_S_3), peakCounts = n())
Counts_perGene_F_C_S = Counts_perGene_F_C_S %>% mutate(tagDensity = tagCounts/peakCounts)

## CoCLIP Input Mock tags
Counts_perGene_I_M = Peak_Co_Input_M %>% filter(NLS_I_M_BC >= 2) %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(NLS_I_M_1 + NLS_I_M_2 + NES_I_M_1 + NES_I_M_2 + G3BP_I_M_1 + G3BP_I_M_2 + G3BP_I_M_3 + G3BP_I_M_4), peakCounts = n())
Counts_perGene_I_M = Counts_perGene_I_M %>% mutate(tagDensity = tagCounts/peakCounts)

## CoCLIP Input Stress tags
Counts_perGene_I_S = Peak_Co_Input_S %>% filter(NLS_I_S_BC >= 2) %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(NLS_I_S_1 + NLS_I_S_2 + NES_I_S_1 + NES_I_S_2 + G3BP_I_S_1 + G3BP_I_S_2 + G3BP_I_S_3 + G3BP_I_S_4), peakCounts = n())
Counts_perGene_I_S = Counts_perGene_I_S %>% mutate(tagDensity = tagCounts/peakCounts)

## CoCLIP NLS Mock tags
Counts_perGene_E_NLS_M = Peak_Co_NLS_M %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(NLS_E_M_1 + NLS_E_M_2 + NLS_E_M_3 + NLS_E_M_4), peakCounts = n())
Counts_perGene_E_NLS_M = Counts_perGene_E_NLS_M %>% mutate(tagDensity = tagCounts/peakCounts)

## CoCLIP NLS Stress tags
Counts_perGene_E_NLS_S = Peak_Co_NLS_S %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(NLS_E_S_1 + NLS_E_S_2 + NLS_E_S_3 + NLS_E_S_4), peakCounts = n())
Counts_perGene_E_NLS_S = Counts_perGene_E_NLS_S %>% mutate(tagDensity = tagCounts/peakCounts)

## CoCLIP NES Mock tags
Counts_perGene_E_NES_M = Peak_Co_NES_M %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(NES_E_M_1 + NES_E_M_2 + NES_E_M_3 + NES_E_M_4), peakCounts = n())
Counts_perGene_E_NES_M = Counts_perGene_E_NES_M %>% mutate(tagDensity = tagCounts/peakCounts)

## CoCLIP NES Stress tags
Counts_perGene_E_NES_S = Peak_Co_NES_S %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(NES_E_S_1 + NES_E_S_2 + NES_E_S_3 + NES_E_S_4), peakCounts = n())
Counts_perGene_E_NES_S = Counts_perGene_E_NES_S %>% mutate(tagDensity = tagCounts/peakCounts)

## CoCLIP G3BP Mock tags
Counts_perGene_E_G3BP_M = Peak_Co_G3BP_M %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(G3BP_E_M_1 + G3BP_E_M_2 + G3BP_E_M_3 + G3BP_E_M_4 + G3BP_E_M_5 + G3BP_E_M_6), peakCounts = n())
Counts_perGene_E_G3BP_M = Counts_perGene_E_G3BP_M %>% mutate(tagDensity = tagCounts/peakCounts)

## CoCLIP G3BP Stress tags
Counts_perGene_E_G3BP_S = Peak_Co_G3BP_S %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(G3BP_E_S_1 + G3BP_E_S_2 + G3BP_E_S_3 + G3BP_E_S_4 + G3BP_E_S_5 + G3BP_E_S_6 + G3BP_E_S_7), peakCounts = n())
Counts_perGene_E_G3BP_S = Counts_perGene_E_G3BP_S %>% mutate(tagDensity = tagCounts/peakCounts)

####################################################################################################################

## Exploratory Stacked Bar Plots Across All Peaks:
####################################################################################################################
# ## Mock vs Stress:
# PeakDistribution_all_M = data.frame(table(Peak_all_M$grouped_annotation), row.names = 1)
# PeakDistribution_all_S = data.frame(table(Peak_all_S$grouped_annotation), row.names = 1)
# PeakDistribution_all = data.frame(table(Peak_all$grouped_annotation), row.names = 1)
# 
# colnames(PeakDistribution_all_M) = c('All_M')
# colnames(PeakDistribution_all_S) = c('All_S')
# colnames(PeakDistribution_all) = c('All_M+S')
# 
# ## Peak Counts Distribution Stacked Bar Graph
# PeakDistribution_all_combined = cbind(PeakDistribution_all_M, 
#                                       PeakDistribution_all_S, 
#                                       PeakDistribution_all)
# PeakDistribution_all_combined$Annotation = rownames(PeakDistribution_all_combined)
# 
# PeakDistribution_all_combined = PeakDistribution_all_combined %>%
#   gather(key = "Source", value = "Freq", 'All_M', 'All_S', 'All_M+S') %>%
#   select(Source, Freq, Annotation)
# 
# PeakDistribution_all_combined$Source = factor(PeakDistribution_all_combined$Source, levels = c('All_M+S', 'All_M', 'All_S'))
# PeakDistribution_all_combined$Annotation = factor(PeakDistribution_all_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "ncRNA", "TE", "Other", "DS10K", 'UnAn'))
# 
# ggplot(PeakDistribution_all_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
#   geom_bar(position='stack', stat='identity') +
#   scale_x_discrete(labels= c('Mock + Stress', 'Mock Only', 'Stress Only')) +
#   ggtitle('All Raw Peak Counts Distributions') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_brewer(palette = "Set3")
# 
# ## Fraction Distribution Stacked Bar Graph
# PeakDistribution_all_combined = cbind(PeakDistribution_all_M/sum(PeakDistribution_all_M), 
#                                       PeakDistribution_all_S/sum(PeakDistribution_all_S), 
#                                       PeakDistribution_all/sum(PeakDistribution_all))
# PeakDistribution_all_combined$Annotation = rownames(PeakDistribution_all_combined)
# 
# PeakDistribution_all_combined = PeakDistribution_all_combined %>%
#   gather(key = "Source", value = "Freq", 'All_M', 'All_S', 'All_M+S') %>%
#   select(Source, Freq, Annotation)
# 
# PeakDistribution_all_combined$Source = factor(PeakDistribution_all_combined$Source, levels = c('All_M+S', 'All_M', 'All_S'))
# PeakDistribution_all_combined$Annotation = factor(PeakDistribution_all_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "ncRNA", "TE", "Other", "DS10K", 'UnAn'))
# 
# ggplot(PeakDistribution_all_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
#   geom_bar(position='stack', stat='identity') +
#   scale_x_discrete(labels= c('Mock + Stress', 'Mock Only', 'Stress Only')) +
#   ggtitle('All Peak Fractions Distributions') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_brewer(palette = "Set3")
####################################################################################################################

## Exploratory Stacked Bar Plots For Fractionation CLIP Peaks:
####################################################################################################################
## Mock vs Stress for Each Fraction
PeakDistribution_F_Nuc_M = data.frame(table(Peak_F_Nuc_M$grouped_annotation), row.names = 1)
PeakDistribution_F_Nuc_S = data.frame(table(Peak_F_Nuc_S$grouped_annotation), row.names = 1)

PeakDistribution_F_Cyt_M = data.frame(table(Peak_F_Cyt_M$grouped_annotation), row.names = 1)
PeakDistribution_F_Cyt_S = data.frame(table(Peak_F_Cyt_S$grouped_annotation), row.names = 1)

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Nuc_S, 
                                    PeakDistribution_F_Cyt_M, PeakDistribution_F_Cyt_S)
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K", 'UnAn'))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
  ggtitle('Raw Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## Fraction Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), 
                                    PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S))
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K", 'UnAn'))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
  ggtitle('Raw Peak Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")
####################################################################################################################

## Exploratory Stacked Bar Plots For Fractionation CLIP Peaks: Specific RNA biotypes ONLY
####################################################################################################################
## mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_F_Nuc_M_mRNA = Peak_F_Nuc_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_F_Nuc_S_mRNA = Peak_F_Nuc_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

Peak_F_Cyt_M_mRNA = Peak_F_Cyt_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_F_Cyt_S_mRNA = Peak_F_Cyt_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

PeakDistribution_F_Nuc_M = data.frame(table(Peak_F_Nuc_M_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_Nuc_S = data.frame(table(Peak_F_Nuc_S_mRNA$finalized_annotation), row.names = 1)

PeakDistribution_F_Cyt_M = data.frame(table(Peak_F_Cyt_M_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_Cyt_S = data.frame(table(Peak_F_Cyt_S_mRNA$finalized_annotation), row.names = 1)

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Nuc_S, 
                                    PeakDistribution_F_Cyt_M, PeakDistribution_F_Cyt_S)
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
  ggtitle('mRNA Raw Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## Fraction Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), 
                                    PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S))
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
  ggtitle('mRNA Raw Peak Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## non-mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_F_Nuc_M_Not_mRNA = Peak_F_Nuc_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_F_Nuc_S_Not_mRNA = Peak_F_Nuc_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

Peak_F_Cyt_M_Not_mRNA = Peak_F_Cyt_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_F_Cyt_S_Not_mRNA = Peak_F_Cyt_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

PeakDistribution_F_Nuc_M = data.frame(table(Peak_F_Nuc_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_Nuc_S = data.frame(table(Peak_F_Nuc_S_Not_mRNA$finalized_annotation), row.names = 1)

PeakDistribution_F_Cyt_M = data.frame(table(Peak_F_Cyt_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_Cyt_S = data.frame(table(Peak_F_Cyt_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_Cyt_S = rbind(PeakDistribution_F_Cyt_S, c(0))
row.names(PeakDistribution_F_Cyt_S) = c(row.names(PeakDistribution_F_Cyt_S)[1:9], 'miRNA')
PeakDistribution_F_Cyt_S = PeakDistribution_F_Cyt_S[row.names(PeakDistribution_F_Cyt_M), 'Freq', drop = FALSE]

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Nuc_S, 
                                    PeakDistribution_F_Cyt_M, PeakDistribution_F_Cyt_S)
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
  ggtitle('non-mRNA Raw Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## Fraction Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), 
                                    PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S))
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
  ggtitle('non-mRNA Raw Peak Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")
####################################################################################################################

## Exploratory Stacked Bar Plots for Fractionation CLIP Tags:
####################################################################################################################
## Mock vs Stress for Each Fraction
PeakDistribution_F_Nuc_M = data.frame(data.frame(annotation = Peak_F_Nuc_M$grouped_annotation, 
                                                 tagSum = rowSums(Peak_F_Nuc_M[, Nuc_F_M])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_Nuc_S = data.frame(data.frame(annotation = Peak_F_Nuc_S$grouped_annotation, 
                                                 tagSum = rowSums(Peak_F_Nuc_S[, Nuc_F_S])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_Cyt_M = data.frame(data.frame(annotation = Peak_F_Cyt_M$grouped_annotation, 
                                                 tagSum = rowSums(Peak_F_Cyt_M[, Cyto_F_M])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_Cyt_S = data.frame(data.frame(annotation = Peak_F_Cyt_S$grouped_annotation, 
                                                 tagSum = rowSums(Peak_F_Cyt_S[, Cyto_F_S])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')

## Tag Counts Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Nuc_S, 
                                    PeakDistribution_F_Cyt_M, PeakDistribution_F_Cyt_S)
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K", 'UnAn'))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
  ggtitle('Normalized Tag Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## Fraction Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), 
                                    PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S))
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K", 'UnAn'))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
  ggtitle('Normalized Tag Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

####################################################################################################################

## Exploratory Stacked Bar Plots For Fractionation CLIP Tags: Specific RNA biotypes ONLY
####################################################################################################################
## mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_F_Nuc_M_mRNA = Peak_F_Nuc_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_F_Nuc_S_mRNA = Peak_F_Nuc_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

Peak_F_Cyt_M_mRNA = Peak_F_Cyt_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_F_Cyt_S_mRNA = Peak_F_Cyt_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

PeakDistribution_F_Nuc_M = data.frame(data.frame(annotation = Peak_F_Nuc_M_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_Nuc_M_mRNA[, Nuc_F_M])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_Nuc_S = data.frame(data.frame(annotation = Peak_F_Nuc_S_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_Nuc_S_mRNA[, Nuc_F_S])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_Cyt_M = data.frame(data.frame(annotation = Peak_F_Cyt_M_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_Cyt_M_mRNA[, Cyto_F_M])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_Cyt_S = data.frame(data.frame(annotation = Peak_F_Cyt_S_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_Cyt_S_mRNA[, Cyto_F_S])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Nuc_S, 
                                    PeakDistribution_F_Cyt_M, PeakDistribution_F_Cyt_S)
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
  ggtitle('mRNA Normalized Tag Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## Fraction Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), 
                                    PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S))
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
  ggtitle('Normalized Tag Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## non-mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_F_Nuc_M_Not_mRNA = Peak_F_Nuc_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_F_Nuc_S_Not_mRNA = Peak_F_Nuc_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

Peak_F_Cyt_M_Not_mRNA = Peak_F_Cyt_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_F_Cyt_S_Not_mRNA = Peak_F_Cyt_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

PeakDistribution_F_Nuc_M = data.frame(data.frame(annotation = Peak_F_Nuc_M_Not_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_Nuc_M_Not_mRNA[, Nuc_F_M])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_Nuc_S = data.frame(data.frame(annotation = Peak_F_Nuc_S_Not_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_Nuc_S_Not_mRNA[, Nuc_F_S])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_Cyt_M = data.frame(data.frame(annotation = Peak_F_Cyt_M_Not_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_Cyt_M_Not_mRNA[, Cyto_F_M])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_Cyt_S = data.frame(data.frame(annotation = Peak_F_Cyt_S_Not_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_Cyt_S_Not_mRNA[, Cyto_F_S])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_F_Cyt_S = rbind(PeakDistribution_F_Cyt_S, c(0))
row.names(PeakDistribution_F_Cyt_S) = c(row.names(PeakDistribution_F_Cyt_S)[1:9], 'miRNA')
PeakDistribution_F_Cyt_S = PeakDistribution_F_Cyt_S[row.names(PeakDistribution_F_Cyt_M), 'tags', drop = FALSE]

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Nuc_S, 
                                    PeakDistribution_F_Cyt_M, PeakDistribution_F_Cyt_S)
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
  ggtitle('non-mRNA Normalized Tag Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## Fraction Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), 
                                    PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S))
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
  ggtitle('non-mRNA Normalized Tag Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")
####################################################################################################################

## Exploratory Stacked Bar Plots for Fractionation CLIP average Tags:
####################################################################################################################
# ## Mock vs Stress for Each Fraction
# PeakDistribution_F_Nuc_M = data.frame(data.frame(annotation = avgTagCounts$grouped_annotation, 
#                                                  AvgTagCounts = avgTagCounts$Nuc_F_M, 
#                                                  BC = avgTagCounts$Nuc_F_M_BC)
#                                       %>% filter(BC >= BC_Threshold_F & AvgTagCounts >= median(avgTagCounts$Nuc_F_M)*rowSum_Multiplier_F)
#                                       %>% group_by(annotation) 
#                                       %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_F_Nuc_S = data.frame(data.frame(annotation = avgTagCounts$grouped_annotation, 
#                                                  AvgTagCounts = avgTagCounts$Nuc_F_S, 
#                                                  BC = avgTagCounts$Nuc_F_S_BC)
#                                       %>% filter(BC >= BC_Threshold_F & AvgTagCounts >= median(avgTagCounts$Nuc_F_S)*rowSum_Multiplier_F)
#                                       %>% group_by(annotation) 
#                                       %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_F_Cyt_M = data.frame(data.frame(annotation = avgTagCounts$grouped_annotation, 
#                                                  AvgTagCounts = avgTagCounts$Cyto_F_M, 
#                                                  BC = avgTagCounts$Cyto_F_M_BC)
#                                       %>% filter(BC >= BC_Threshold_F & AvgTagCounts >= median(avgTagCounts$Cyto_F_M)*rowSum_Multiplier_F)
#                                       %>% group_by(annotation) 
#                                       %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_F_Cyt_S = data.frame(data.frame(annotation = avgTagCounts$grouped_annotation, 
#                                                  AvgTagCounts = avgTagCounts$Cyto_F_S, 
#                                                  BC = avgTagCounts$Cyto_F_S_BC)
#                                       %>% filter(BC >= BC_Threshold_F & AvgTagCounts >= median(avgTagCounts$Cyto_F_S)*rowSum_Multiplier_F)
#                                       %>% group_by(annotation) 
#                                       %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
# colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
# colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
# colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')
# 
# ## Tag Counts Distribution Stacked Bar Graph
# PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Nuc_S,
#                                     PeakDistribution_F_Cyt_M, PeakDistribution_F_Cyt_S)
# PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)
# 
# PeakDistribution_F_combined = PeakDistribution_F_combined %>%
#   gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
#   select(Source, Freq, Annotation)
# 
# PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
# PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn"))
# 
# ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
#   geom_bar(position='stack', stat='identity') +
#   scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
#   ggtitle('Peak Density Distributions') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_brewer(palette = "Set3")
# 
# ## Fraction Distribution Stacked Bar Graph
# PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), 
#                                     PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S))
# PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)
# 
# PeakDistribution_F_combined = PeakDistribution_F_combined %>%
#   gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
#   select(Source, Freq, Annotation)
# 
# PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
# PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn"))
# 
# ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
#   geom_bar(position='stack', stat='identity') +
#   scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
#   ggtitle('Peak Density Distributions by Fraction') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_brewer(palette = "Set3")

####################################################################################################################

## Exploratory Stacked Bar Plots for Fractionation CLIP average Tags: Specific RNA biotypes ONLY
####################################################################################################################
# ## mRNA Features Only
# ## Mock vs Stress for Each Fraction
# avgTagCounts_mRNA = avgTagCounts %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == "CDS_RI" | finalized_annotation == 'DS10K')
# 
# ## Mock vs Stress for Each Fraction
# PeakDistribution_F_Nuc_M = data.frame(data.frame(annotation = avgTagCounts_mRNA$finalized_annotation, 
#                                                  AvgTagCounts = avgTagCounts_mRNA$Nuc_F_M, 
#                                                  BC = avgTagCounts_mRNA$Nuc_F_M_BC)
#                                       %>% filter(BC >= BC_Threshold_F & AvgTagCounts >= median(avgTagCounts_mRNA$Nuc_F_M)*rowSum_Multiplier_F)
#                                       %>% group_by(annotation) 
#                                       %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_F_Nuc_S = data.frame(data.frame(annotation = avgTagCounts_mRNA$finalized_annotation, 
#                                                  AvgTagCounts = avgTagCounts_mRNA$Nuc_F_S, 
#                                                  BC = avgTagCounts_mRNA$Nuc_F_S_BC)
#                                       %>% filter(BC >= BC_Threshold_F & AvgTagCounts >= median(avgTagCounts_mRNA$Nuc_F_S)*rowSum_Multiplier_F)
#                                       %>% group_by(annotation) 
#                                       %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_F_Cyt_M = data.frame(data.frame(annotation = avgTagCounts_mRNA$finalized_annotation, 
#                                                  AvgTagCounts = avgTagCounts_mRNA$Cyto_F_M, 
#                                                  BC = avgTagCounts_mRNA$Cyto_F_M_BC)
#                                       %>% filter(BC >= BC_Threshold_F & AvgTagCounts >= median(avgTagCounts_mRNA$Cyto_F_M)*rowSum_Multiplier_F)
#                                       %>% group_by(annotation) 
#                                       %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_F_Cyt_S = data.frame(data.frame(annotation = avgTagCounts_mRNA$finalized_annotation, 
#                                                  AvgTagCounts = avgTagCounts_mRNA$Cyto_F_S, 
#                                                  BC = avgTagCounts_mRNA$Cyto_F_S_BC)
#                                       %>% filter(BC >= BC_Threshold_F & AvgTagCounts >= median(avgTagCounts_mRNA$Cyto_F_S)*rowSum_Multiplier_F)
#                                       %>% group_by(annotation) 
#                                       %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
# colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
# colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
# colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')
# 
# ## Tag Counts Distribution Stacked Bar Graph
# PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Nuc_S,
#                                     PeakDistribution_F_Cyt_M, PeakDistribution_F_Cyt_S)
# PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)
# 
# PeakDistribution_F_combined = PeakDistribution_F_combined %>%
#   gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
#   select(Source, Freq, Annotation)
# 
# PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
# PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))
# 
# ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
#   geom_bar(position='stack', stat='identity') +
#   scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
#   ggtitle('mRNA Peak Density Distributions') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_brewer(palette = "Set3")
# 
# ## Fraction Distribution Stacked Bar Graph
# PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), 
#                                     PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S))
# PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)
# 
# PeakDistribution_F_combined = PeakDistribution_F_combined %>%
#   gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
#   select(Source, Freq, Annotation)
# 
# PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
# PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))
# 
# ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
#   geom_bar(position='stack', stat='identity') +
#   scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
#   ggtitle('mRNA Peak Density Distributions by Fraction') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_brewer(palette = "Set3")
# 
# ## nonmRNA Features Only
# ## Mock vs Stress for Each Fraction
# avgTagCounts_Not_mRNA = avgTagCounts %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
# 
# ## Mock vs Stress for Each Fraction
# PeakDistribution_F_Nuc_M = data.frame(data.frame(annotation = avgTagCounts_Not_mRNA$finalized_annotation, 
#                                                  AvgTagCounts = avgTagCounts_Not_mRNA$Nuc_F_M, 
#                                                  BC = avgTagCounts_Not_mRNA$Nuc_F_M_BC)
#                                       %>% filter(BC >= BC_Threshold_F & AvgTagCounts >= median(avgTagCounts_Not_mRNA$Nuc_F_M)*rowSum_Multiplier_F)
#                                       %>% group_by(annotation) 
#                                       %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_F_Nuc_S = data.frame(data.frame(annotation = avgTagCounts_Not_mRNA$finalized_annotation, 
#                                                  AvgTagCounts = avgTagCounts_Not_mRNA$Nuc_F_S, 
#                                                  BC = avgTagCounts_Not_mRNA$Nuc_F_S_BC)
#                                       %>% filter(BC >= BC_Threshold_F & AvgTagCounts >= median(avgTagCounts_Not_mRNA$Nuc_F_S)*rowSum_Multiplier_F)
#                                       %>% group_by(annotation) 
#                                       %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_F_Cyt_M = data.frame(data.frame(annotation = avgTagCounts_Not_mRNA$finalized_annotation, 
#                                                  AvgTagCounts = avgTagCounts_Not_mRNA$Cyto_F_M, 
#                                                  BC = avgTagCounts_Not_mRNA$Cyto_F_M_BC)
#                                       %>% filter(BC >= BC_Threshold_F & AvgTagCounts >= median(avgTagCounts_Not_mRNA$Cyto_F_M)*rowSum_Multiplier_F)
#                                       %>% group_by(annotation) 
#                                       %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# PeakDistribution_F_Cyt_M = rbind(PeakDistribution_F_Cyt_M, c(0))
# row.names(PeakDistribution_F_Cyt_M) = c(row.names(PeakDistribution_F_Cyt_M)[1:10], 'scaRNA')
# PeakDistribution_F_Cyt_M = PeakDistribution_F_Cyt_M[row.names(PeakDistribution_F_Nuc_S), 'tags', drop = FALSE]
# 
# PeakDistribution_F_Cyt_S = data.frame(data.frame(annotation = avgTagCounts_Not_mRNA$finalized_annotation, 
#                                                  AvgTagCounts = avgTagCounts_Not_mRNA$Cyto_F_S, 
#                                                  BC = avgTagCounts_Not_mRNA$Cyto_F_S_BC)
#                                       %>% filter(BC >= BC_Threshold_F & AvgTagCounts >= median(avgTagCounts_Not_mRNA$Cyto_F_S)*rowSum_Multiplier_F)
#                                       %>% group_by(annotation) 
#                                       %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# PeakDistribution_F_Cyt_S = rbind(PeakDistribution_F_Cyt_S, c(0), c(0))
# row.names(PeakDistribution_F_Cyt_S) = c(row.names(PeakDistribution_F_Cyt_S)[1:9], 'scaRNA', 'miRNA')
# PeakDistribution_F_Cyt_S = PeakDistribution_F_Cyt_S[row.names(PeakDistribution_F_Nuc_S), 'tags', drop = FALSE]
# 
# colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
# colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
# colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
# colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')
# 
# ## Tag Counts Distribution Stacked Bar Graph
# PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Nuc_S,
#                                     PeakDistribution_F_Cyt_M, PeakDistribution_F_Cyt_S)
# PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)
# 
# PeakDistribution_F_combined = PeakDistribution_F_combined %>%
#   gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
#   select(Source, Freq, Annotation)
# 
# PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
# PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))
# 
# ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
#   geom_bar(position='stack', stat='identity') +
#   scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
#   ggtitle('non-mRNA Peak Density Distributions') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_brewer(palette = "Set3")
# 
# ## Fraction Distribution Stacked Bar Graph
# PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), 
#                                     PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S))
# PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)
# 
# PeakDistribution_F_combined = PeakDistribution_F_combined %>%
#   gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt') %>%
#   select(Source, Freq, Annotation)
# 
# PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt'))
# PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))
# 
# ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
#   geom_bar(position='stack', stat='identity') +
#   scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress')) +
#   ggtitle('non-mRNA Peak Density Distributions by Fraction') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_brewer(palette = "Set3")
####################################################################################################################

## FIGURE1 Exploratory Stacked Bar Plots For CoCLIP Peaks:
####################################################################################################################
## Mock vs Stress for Each Fraction
PeakDistribution_Co_Input_M = data.frame(table(Peak_Co_Input_M$grouped_annotation), row.names = 1)
PeakDistribution_Co_Input_S = data.frame(table(Peak_Co_Input_S$grouped_annotation), row.names = 1)

PeakDistribution_Co_NLS_M = data.frame(table(Peak_Co_NLS_M$grouped_annotation), row.names = 1)
PeakDistribution_Co_NLS_S = data.frame(table(Peak_Co_NLS_S$grouped_annotation), row.names = 1)

PeakDistribution_Co_NES_M = data.frame(table(Peak_Co_NES_M$grouped_annotation), row.names = 1)
PeakDistribution_Co_NES_S = data.frame(table(Peak_Co_NES_S$grouped_annotation), row.names = 1)

PeakDistribution_Co_G3BP_M = data.frame(table(Peak_Co_G3BP_M$grouped_annotation), row.names = 1)
PeakDistribution_Co_G3BP_S = data.frame(table(Peak_Co_G3BP_S$grouped_annotation), row.names = 1)

colnames(PeakDistribution_Co_Input_M) = c('Co_M_Input')
colnames(PeakDistribution_Co_Input_S) = c('Co_S_Input')
colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')
colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')
colnames(PeakDistribution_Co_G3BP_M) = c('Co_M_G3BP')
colnames(PeakDistribution_Co_G3BP_S) = c('Co_S_G3BP')

PeakDistribution_Co_Input_M = PeakDistribution_Co_Input_M %>% slice(-which(rownames(PeakDistribution_Co_Input_M) == 'UnAn'))
PeakDistribution_Co_Input_S = PeakDistribution_Co_Input_S %>% slice(-which(rownames(PeakDistribution_Co_Input_S) == 'UnAn'))
PeakDistribution_Co_NLS_M = PeakDistribution_Co_NLS_M %>% slice(-which(rownames(PeakDistribution_Co_NLS_M) == 'UnAn'))
PeakDistribution_Co_NLS_S = PeakDistribution_Co_NLS_S %>% slice(-which(rownames(PeakDistribution_Co_NLS_S) == 'UnAn'))
PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M %>% slice(-which(rownames(PeakDistribution_Co_NES_M) == 'UnAn'))
PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S %>% slice(-which(rownames(PeakDistribution_Co_NES_S) == 'UnAn'))
PeakDistribution_Co_G3BP_M = PeakDistribution_Co_G3BP_M %>% slice(-which(rownames(PeakDistribution_Co_G3BP_M) == 'UnAn'))
PeakDistribution_Co_G3BP_S = PeakDistribution_Co_G3BP_S %>% slice(-which(rownames(PeakDistribution_Co_G3BP_S) == 'UnAn'))

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
                                     PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
                                     PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
                                     PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", 'ncRNA', "TE", "Other", "DS10K", 'UnAn'))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_M_Input' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels= c('Input', 'NLS', 'NES')) +
  ggtitle('CoCLIP Mock Peak Distribution') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Fraction Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M/sum(PeakDistribution_Co_Input_M), PeakDistribution_Co_Input_S/sum(PeakDistribution_Co_Input_S), 
                                     PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), 
                                     PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S),
                                     PeakDistribution_Co_G3BP_M/sum(PeakDistribution_Co_G3BP_M), PeakDistribution_Co_G3BP_S/sum(PeakDistribution_Co_G3BP_S))
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", 'ncRNA', "TE", "Other", "DS10K", 'UnAn'))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_M_Input' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels= c('Input', 'NLS', 'NES')) +
  ggtitle('CoCLIP Mock Peak Distribution By Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
####################################################################################################################

## FIGURE1 Exploratory Stacked Bar Plots For CoCLIP Peaks: Specific RNA biotypes ONLY
####################################################################################################################
## mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_Co_Input_M_mRNA = Peak_Co_Input_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_Input_S_mRNA = Peak_Co_Input_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

Peak_Co_NLS_M_mRNA = Peak_Co_NLS_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_NLS_S_mRNA = Peak_Co_NLS_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

Peak_Co_NES_M_mRNA = Peak_Co_NES_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_NES_S_mRNA = Peak_Co_NES_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

Peak_Co_G3BP_M_mRNA = Peak_Co_G3BP_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_G3BP_S_mRNA = Peak_Co_G3BP_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

PeakDistribution_Co_Input_M = data.frame(table(Peak_Co_Input_M_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_Input_S = data.frame(table(Peak_Co_Input_S_mRNA$finalized_annotation), row.names = 1)

PeakDistribution_Co_NLS_M = data.frame(table(Peak_Co_NLS_M_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NLS_S = data.frame(table(Peak_Co_NLS_S_mRNA$finalized_annotation), row.names = 1)

PeakDistribution_Co_NES_M = data.frame(table(Peak_Co_NES_M_mRNA$finalized_annotation), row.names = 1)
# PeakDistribution_Co_NES_M = rbind(PeakDistribution_Co_NES_M, c(0))
# row.names(PeakDistribution_Co_NES_M) = c(row.names(PeakDistribution_Co_NES_M)[1:6], 'nC_RI')
# PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

PeakDistribution_Co_NES_S = data.frame(table(Peak_Co_NES_S_mRNA$finalized_annotation), row.names = 1)
# PeakDistribution_Co_NES_S = rbind(PeakDistribution_Co_NES_S, c(0))
# row.names(PeakDistribution_Co_NES_S) = c(row.names(PeakDistribution_Co_NES_S)[1:6], 'nC_RI')
# PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

PeakDistribution_Co_G3BP_M = data.frame(table(Peak_Co_G3BP_M_mRNA$finalized_annotation), row.names = 1)
# PeakDistribution_Co_G3BP_M = rbind(PeakDistribution_Co_G3BP_M, c(0))
# row.names(PeakDistribution_Co_G3BP_M) = c(row.names(PeakDistribution_Co_G3BP_M)[1:6], 'nC_RI')
# PeakDistribution_Co_G3BP_M = PeakDistribution_Co_G3BP_M[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

PeakDistribution_Co_G3BP_S = data.frame(table(Peak_Co_G3BP_S_mRNA$finalized_annotation), row.names = 1)
# PeakDistribution_Co_G3BP_S = rbind(PeakDistribution_Co_G3BP_S, c(0))
# row.names(PeakDistribution_Co_G3BP_S) = c(row.names(PeakDistribution_Co_G3BP_S)[1:6], 'nC_RI')
# PeakDistribution_Co_G3BP_S = PeakDistribution_Co_G3BP_S[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

colnames(PeakDistribution_Co_Input_M) = c('Co_M_Input')
colnames(PeakDistribution_Co_Input_S) = c('Co_S_Input')
colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')
colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')
colnames(PeakDistribution_Co_G3BP_M) = c('Co_M_G3BP')
colnames(PeakDistribution_Co_G3BP_S) = c('Co_S_G3BP')

PeakDistribution_Co_Input_M = PeakDistribution_Co_Input_M %>% slice(-which(rownames(PeakDistribution_Co_Input_M) == 'UnAn'))
PeakDistribution_Co_Input_S = PeakDistribution_Co_Input_S %>% slice(-which(rownames(PeakDistribution_Co_Input_S) == 'UnAn'))
PeakDistribution_Co_NLS_M = PeakDistribution_Co_NLS_M %>% slice(-which(rownames(PeakDistribution_Co_NLS_M) == 'UnAn'))
PeakDistribution_Co_NLS_S = PeakDistribution_Co_NLS_S %>% slice(-which(rownames(PeakDistribution_Co_NLS_S) == 'UnAn'))
PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M %>% slice(-which(rownames(PeakDistribution_Co_NES_M) == 'UnAn'))
PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S %>% slice(-which(rownames(PeakDistribution_Co_NES_S) == 'UnAn'))
PeakDistribution_Co_G3BP_M = PeakDistribution_Co_G3BP_M %>% slice(-which(rownames(PeakDistribution_Co_G3BP_M) == 'UnAn'))
PeakDistribution_Co_G3BP_S = PeakDistribution_Co_G3BP_S %>% slice(-which(rownames(PeakDistribution_Co_G3BP_S) == 'UnAn'))

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
                                     PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
                                     PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
                                     PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_M_Input' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels= c('Input', 'NLS', 'NES')) +
  ggtitle('CoCLIP Mock mRNA Peak Distribution') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Fraction Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M/sum(PeakDistribution_Co_Input_M), PeakDistribution_Co_Input_S/sum(PeakDistribution_Co_Input_S), 
                                     PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), 
                                     PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S),
                                     PeakDistribution_Co_G3BP_M/sum(PeakDistribution_Co_G3BP_M), PeakDistribution_Co_G3BP_S/sum(PeakDistribution_Co_G3BP_S))
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_M_Input' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels= c('Input', 'NLS', 'NES')) +
  ggtitle('CoCLIP Mock mRNA Peak Distribution by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## non-mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_Co_Input_M_Not_mRNA = Peak_Co_Input_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_Input_S_Not_mRNA = Peak_Co_Input_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

Peak_Co_NLS_M_Not_mRNA = Peak_Co_NLS_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_NLS_S_Not_mRNA = Peak_Co_NLS_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

Peak_Co_NES_M_Not_mRNA = Peak_Co_NES_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_NES_S_Not_mRNA = Peak_Co_NES_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

Peak_Co_G3BP_M_Not_mRNA = Peak_Co_G3BP_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_G3BP_S_Not_mRNA = Peak_Co_G3BP_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

PeakDistribution_Co_NLS_M = data.frame(table(Peak_Co_NLS_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NLS_M = rbind(PeakDistribution_Co_NLS_M, c(0))
row.names(PeakDistribution_Co_NLS_M) = c(row.names(PeakDistribution_Co_NLS_M)[1:10], 'scaRNA')

PeakDistribution_Co_NLS_S = data.frame(table(Peak_Co_NLS_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NLS_S = rbind(PeakDistribution_Co_NLS_S, c(0))
row.names(PeakDistribution_Co_NLS_S) = c(row.names(PeakDistribution_Co_NLS_S)[1:10], 'scaRNA')
PeakDistribution_Co_NLS_S = PeakDistribution_Co_NLS_S[row.names(PeakDistribution_Co_NLS_M), 'Freq', drop = FALSE]

PeakDistribution_Co_Input_M = data.frame(table(Peak_Co_Input_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_Input_M = rbind(PeakDistribution_Co_Input_M, c(0), c(0))
row.names(PeakDistribution_Co_Input_M) = c(row.names(PeakDistribution_Co_Input_M)[1:9], 'scaRNA', 'nC_RI')
PeakDistribution_Co_Input_M = PeakDistribution_Co_Input_M[row.names(PeakDistribution_Co_NLS_M), 'Freq', drop = FALSE]

PeakDistribution_Co_Input_S = data.frame(table(Peak_Co_Input_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_Input_S = rbind(PeakDistribution_Co_Input_S, c(0), c(0))
row.names(PeakDistribution_Co_Input_S) = c(row.names(PeakDistribution_Co_Input_S)[1:9], 'scaRNA', 'nC_RI')
PeakDistribution_Co_Input_S = PeakDistribution_Co_Input_S[row.names(PeakDistribution_Co_NLS_M), 'Freq', drop = FALSE]

PeakDistribution_Co_NES_M = data.frame(table(Peak_Co_NES_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NES_M = rbind(PeakDistribution_Co_NES_M, c(0), c(0), c(0))
row.names(PeakDistribution_Co_NES_M) = c(row.names(PeakDistribution_Co_NES_M)[1:8], 'miRNA', 'nC_RI', 'scaRNA')
PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

PeakDistribution_Co_NES_S = data.frame(table(Peak_Co_NES_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NES_S = rbind(PeakDistribution_Co_NES_S, c(0), c(0))
row.names(PeakDistribution_Co_NES_S) = c(row.names(PeakDistribution_Co_NES_S)[1:9], 'scaRNA', 'nC_RI')
PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S[row.names(PeakDistribution_Co_NLS_M), 'Freq', drop = FALSE]

PeakDistribution_Co_G3BP_M = data.frame(table(Peak_Co_G3BP_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_G3BP_M = rbind(PeakDistribution_Co_G3BP_M, c(0), c(0), c(0))
row.names(PeakDistribution_Co_G3BP_M) = c(row.names(PeakDistribution_Co_G3BP_M)[1:8], 'miRNA', 'nC_RI', 'scaRNA')
PeakDistribution_Co_G3BP_M = PeakDistribution_Co_G3BP_M[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

PeakDistribution_Co_G3BP_S = data.frame(table(Peak_Co_G3BP_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_G3BP_S = rbind(PeakDistribution_Co_G3BP_S, c(0), c(0), c(0))
row.names(PeakDistribution_Co_G3BP_S) = c(row.names(PeakDistribution_Co_G3BP_S)[1:8], 'miRNA', 'nC_RI', 'scaRNA')
PeakDistribution_Co_G3BP_S = PeakDistribution_Co_G3BP_S[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

colnames(PeakDistribution_Co_Input_M) = c('Co_M_Input')
colnames(PeakDistribution_Co_Input_S) = c('Co_S_Input')
colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')
colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')
colnames(PeakDistribution_Co_G3BP_M) = c('Co_M_G3BP')
colnames(PeakDistribution_Co_G3BP_S) = c('Co_S_G3BP')

PeakDistribution_Co_Input_M = PeakDistribution_Co_Input_M %>% slice(-which(rownames(PeakDistribution_Co_Input_M) == 'UnAn'))
PeakDistribution_Co_Input_S = PeakDistribution_Co_Input_S %>% slice(-which(rownames(PeakDistribution_Co_Input_S) == 'UnAn'))
PeakDistribution_Co_NLS_M = PeakDistribution_Co_NLS_M %>% slice(-which(rownames(PeakDistribution_Co_NLS_M) == 'UnAn'))
PeakDistribution_Co_NLS_S = PeakDistribution_Co_NLS_S %>% slice(-which(rownames(PeakDistribution_Co_NLS_S) == 'UnAn'))
PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M %>% slice(-which(rownames(PeakDistribution_Co_NES_M) == 'UnAn'))
PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S %>% slice(-which(rownames(PeakDistribution_Co_NES_S) == 'UnAn'))
PeakDistribution_Co_G3BP_M = PeakDistribution_Co_G3BP_M %>% slice(-which(rownames(PeakDistribution_Co_G3BP_M) == 'UnAn'))
PeakDistribution_Co_G3BP_S = PeakDistribution_Co_G3BP_S %>% slice(-which(rownames(PeakDistribution_Co_G3BP_S) == 'UnAn'))

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
                                     PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
                                     PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
                                     PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_M_Input' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels = c('Input', 'NLS', 'NES')) +
  ggtitle('CoCLIP Mock ncRNA Peak Distribution') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Fraction Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M/sum(PeakDistribution_Co_Input_M), PeakDistribution_Co_Input_S/sum(PeakDistribution_Co_Input_S), 
                                     PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), 
                                     PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S),
                                     PeakDistribution_Co_G3BP_M/sum(PeakDistribution_Co_G3BP_M), PeakDistribution_Co_G3BP_S/sum(PeakDistribution_Co_G3BP_S))
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_M_Input' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels = c('Input', 'NLS', 'NES')) +
  ggtitle('CoCLIP Mock ncRNA Peak Distribution by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
####################################################################################################################

## FIGURE5 Exploratory Stacked Bar Plots For G3BP CoCLIP Peaks: Specific RNA biotypes ONLY
####################################################################################################################
## mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_Co_Input_M_mRNA = Peak_Co_Input_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_Input_S_mRNA = Peak_Co_Input_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

Peak_Co_NLS_M_mRNA = Peak_Co_NLS_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_NLS_S_mRNA = Peak_Co_NLS_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

Peak_Co_NES_M_mRNA = Peak_Co_NES_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_NES_S_mRNA = Peak_Co_NES_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

Peak_Co_G3BP_M_mRNA = Peak_Co_G3BP_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_G3BP_S_mRNA = Peak_Co_G3BP_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

PeakDistribution_Co_Input_M = data.frame(table(Peak_Co_Input_M_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_Input_S = data.frame(table(Peak_Co_Input_S_mRNA$finalized_annotation), row.names = 1)

PeakDistribution_Co_NLS_M = data.frame(table(Peak_Co_NLS_M_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NLS_S = data.frame(table(Peak_Co_NLS_S_mRNA$finalized_annotation), row.names = 1)

PeakDistribution_Co_NES_M = data.frame(table(Peak_Co_NES_M_mRNA$finalized_annotation), row.names = 1)
# PeakDistribution_Co_NES_M = rbind(PeakDistribution_Co_NES_M, c(0))
# row.names(PeakDistribution_Co_NES_M) = c(row.names(PeakDistribution_Co_NES_M)[1:6], 'nC_RI')
# PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

PeakDistribution_Co_NES_S = data.frame(table(Peak_Co_NES_S_mRNA$finalized_annotation), row.names = 1)
# PeakDistribution_Co_NES_S = rbind(PeakDistribution_Co_NES_S, c(0))
# row.names(PeakDistribution_Co_NES_S) = c(row.names(PeakDistribution_Co_NES_S)[1:6], 'nC_RI')
# PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

PeakDistribution_Co_G3BP_M = data.frame(table(Peak_Co_G3BP_M_mRNA$finalized_annotation), row.names = 1)
# PeakDistribution_Co_G3BP_M = rbind(PeakDistribution_Co_G3BP_M, c(0))
# row.names(PeakDistribution_Co_G3BP_M) = c(row.names(PeakDistribution_Co_G3BP_M)[1:6], 'nC_RI')
# PeakDistribution_Co_G3BP_M = PeakDistribution_Co_G3BP_M[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

PeakDistribution_Co_G3BP_S = data.frame(table(Peak_Co_G3BP_S_mRNA$finalized_annotation), row.names = 1)
# PeakDistribution_Co_G3BP_S = rbind(PeakDistribution_Co_G3BP_S, c(0))
# row.names(PeakDistribution_Co_G3BP_S) = c(row.names(PeakDistribution_Co_G3BP_S)[1:6], 'nC_RI')
# PeakDistribution_Co_G3BP_S = PeakDistribution_Co_G3BP_S[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

colnames(PeakDistribution_Co_Input_M) = c('Co_M_Input')
colnames(PeakDistribution_Co_Input_S) = c('Co_S_Input')
colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')
colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')
colnames(PeakDistribution_Co_G3BP_M) = c('Co_M_G3BP')
colnames(PeakDistribution_Co_G3BP_S) = c('Co_S_G3BP')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
                                     PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
                                     PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
                                     PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_M_Input' | Source == 'Co_M_G3BP'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels= c('Input', 'G3BP')) +
  ggtitle('CoCLIP Mock mRNA Peak Distribution') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_S_Input' | Source == 'Co_S_G3BP'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels= c('Input', 'G3BP')) +
  ggtitle('CoCLIP Arsenite mRNA Peak Distribution') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Fraction Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M/sum(PeakDistribution_Co_Input_M), PeakDistribution_Co_Input_S/sum(PeakDistribution_Co_Input_S), 
                                     PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), 
                                     PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S),
                                     PeakDistribution_Co_G3BP_M/sum(PeakDistribution_Co_G3BP_M), PeakDistribution_Co_G3BP_S/sum(PeakDistribution_Co_G3BP_S))
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_M_Input' | Source == 'Co_M_G3BP'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels = c('Input', 'G3BP'))+
  ggtitle('CoCLIP Mock mRNA Peak Distribution by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_S_Input' | Source == 'Co_S_G3BP'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels = c('Input', 'G3BP')) +
  ggtitle('CoCLIP Arsenite mRNA Peak Distribution by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## non-mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_Co_Input_M_Not_mRNA = Peak_Co_Input_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_Input_S_Not_mRNA = Peak_Co_Input_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

Peak_Co_NLS_M_Not_mRNA = Peak_Co_NLS_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_NLS_S_Not_mRNA = Peak_Co_NLS_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

Peak_Co_NES_M_Not_mRNA = Peak_Co_NES_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_NES_S_Not_mRNA = Peak_Co_NES_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

Peak_Co_G3BP_M_Not_mRNA = Peak_Co_G3BP_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_G3BP_S_Not_mRNA = Peak_Co_G3BP_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

PeakDistribution_Co_NLS_M = data.frame(table(Peak_Co_NLS_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NLS_M = rbind(PeakDistribution_Co_NLS_M, c(0))
row.names(PeakDistribution_Co_NLS_M) = c(row.names(PeakDistribution_Co_NLS_M)[1:10], 'scaRNA')

PeakDistribution_Co_NLS_S = data.frame(table(Peak_Co_NLS_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NLS_S = rbind(PeakDistribution_Co_NLS_S, c(0))
row.names(PeakDistribution_Co_NLS_S) = c(row.names(PeakDistribution_Co_NLS_S)[1:10], 'scaRNA')
PeakDistribution_Co_NLS_S = PeakDistribution_Co_NLS_S[row.names(PeakDistribution_Co_NLS_M), 'Freq', drop = FALSE]

PeakDistribution_Co_Input_M = data.frame(table(Peak_Co_Input_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_Input_M = rbind(PeakDistribution_Co_Input_M, c(0), c(0))
row.names(PeakDistribution_Co_Input_M) = c(row.names(PeakDistribution_Co_Input_M)[1:9], 'scaRNA', 'nC_RI')
PeakDistribution_Co_Input_M = PeakDistribution_Co_Input_M[row.names(PeakDistribution_Co_NLS_M), 'Freq', drop = FALSE]

PeakDistribution_Co_Input_S = data.frame(table(Peak_Co_Input_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_Input_S = rbind(PeakDistribution_Co_Input_S, c(0), c(0))
row.names(PeakDistribution_Co_Input_S) = c(row.names(PeakDistribution_Co_Input_S)[1:9], 'scaRNA', 'nC_RI')
PeakDistribution_Co_Input_S = PeakDistribution_Co_Input_S[row.names(PeakDistribution_Co_NLS_M), 'Freq', drop = FALSE]

PeakDistribution_Co_NES_M = data.frame(table(Peak_Co_NES_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NES_M = rbind(PeakDistribution_Co_NES_M, c(0), c(0), c(0))
row.names(PeakDistribution_Co_NES_M) = c(row.names(PeakDistribution_Co_NES_M)[1:8], 'miRNA', 'nC_RI', 'scaRNA')
PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

PeakDistribution_Co_NES_S = data.frame(table(Peak_Co_NES_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NES_S = rbind(PeakDistribution_Co_NES_S, c(0), c(0))
row.names(PeakDistribution_Co_NES_S) = c(row.names(PeakDistribution_Co_NES_S)[1:9], 'scaRNA', 'nC_RI')
PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S[row.names(PeakDistribution_Co_NLS_M), 'Freq', drop = FALSE]

PeakDistribution_Co_G3BP_M = data.frame(table(Peak_Co_G3BP_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_G3BP_M = rbind(PeakDistribution_Co_G3BP_M, c(0), c(0), c(0))
row.names(PeakDistribution_Co_G3BP_M) = c(row.names(PeakDistribution_Co_G3BP_M)[1:8], 'miRNA', 'nC_RI', 'scaRNA')
PeakDistribution_Co_G3BP_M = PeakDistribution_Co_G3BP_M[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

PeakDistribution_Co_G3BP_S = data.frame(table(Peak_Co_G3BP_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_G3BP_S = rbind(PeakDistribution_Co_G3BP_S, c(0), c(0), c(0))
row.names(PeakDistribution_Co_G3BP_S) = c(row.names(PeakDistribution_Co_G3BP_S)[1:8], 'miRNA', 'nC_RI', 'scaRNA')
PeakDistribution_Co_G3BP_S = PeakDistribution_Co_G3BP_S[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

colnames(PeakDistribution_Co_Input_M) = c('Co_M_Input')
colnames(PeakDistribution_Co_Input_S) = c('Co_S_Input')
colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')
colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')
colnames(PeakDistribution_Co_G3BP_M) = c('Co_M_G3BP')
colnames(PeakDistribution_Co_G3BP_S) = c('Co_S_G3BP')

PeakDistribution_Co_Input_M = PeakDistribution_Co_Input_M %>% slice(-which(rownames(PeakDistribution_Co_Input_M) == 'UnAn'))
PeakDistribution_Co_Input_S = PeakDistribution_Co_Input_S %>% slice(-which(rownames(PeakDistribution_Co_Input_S) == 'UnAn'))
PeakDistribution_Co_NLS_M = PeakDistribution_Co_NLS_M %>% slice(-which(rownames(PeakDistribution_Co_NLS_M) == 'UnAn'))
PeakDistribution_Co_NLS_S = PeakDistribution_Co_NLS_S %>% slice(-which(rownames(PeakDistribution_Co_NLS_S) == 'UnAn'))
PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M %>% slice(-which(rownames(PeakDistribution_Co_NES_M) == 'UnAn'))
PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S %>% slice(-which(rownames(PeakDistribution_Co_NES_S) == 'UnAn'))
PeakDistribution_Co_G3BP_M = PeakDistribution_Co_G3BP_M %>% slice(-which(rownames(PeakDistribution_Co_G3BP_M) == 'UnAn'))
PeakDistribution_Co_G3BP_S = PeakDistribution_Co_G3BP_S %>% slice(-which(rownames(PeakDistribution_Co_G3BP_S) == 'UnAn'))

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
                                     PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
                                     PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
                                     PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_M_Input' | Source == 'Co_M_G3BP'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels = c('Input', 'G3BP')) +
  ggtitle('CoCLIP Mock ncRNA Peak Distribution') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_S_Input' | Source == 'Co_S_G3BP'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels = c('Input', 'G3BP')) +
  ggtitle('CoCLIP Arsenite ncRNA Peak Distribution') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Fraction Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M/sum(PeakDistribution_Co_Input_M), PeakDistribution_Co_Input_S/sum(PeakDistribution_Co_Input_S), 
                                     PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), 
                                     PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S),
                                     PeakDistribution_Co_G3BP_M/sum(PeakDistribution_Co_G3BP_M), PeakDistribution_Co_G3BP_S/sum(PeakDistribution_Co_G3BP_S))
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_M_Input' | Source == 'Co_M_G3BP'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels = c('Input', 'G3BP')) +
  ggtitle('CoCLIP Mock ncRNA Peak Distribution by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_S_Input' | Source == 'Co_S_G3BP'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels = c('Input', 'G3BP')) +
  ggtitle('CoCLIP Arsenite ncRNA Peak Distribution by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

####################################################################################################################

## Exploratory Stacked Bar Plots For CoCLIP Tags:
####################################################################################################################
## Mock vs Stress for Each Fraction
PeakDistribution_Co_Input_M = data.frame(data.frame(annotation = Peak_Co_Input_M$grouped_annotation, 
                                                    tagSum = rowSums(Peak_Co_Input_M[, c(NLS_I_M, NES_I_M, G3BP_I_M)])) 
                                         %>% group_by(annotation) 
                                         %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_Input_S = data.frame(data.frame(annotation = Peak_Co_Input_S$grouped_annotation, 
                                                    tagSum = rowSums(Peak_Co_Input_S[, c(NLS_I_S, NES_I_S, G3BP_I_S)])) 
                                         %>% group_by(annotation) 
                                         %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_NLS_M = data.frame(data.frame(annotation = Peak_Co_NLS_M$grouped_annotation, 
                                                  tagSum = rowSums(Peak_Co_NLS_M[, c(NLS_E_M)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_NLS_S = data.frame(data.frame(annotation = Peak_Co_NLS_S$grouped_annotation, 
                                                  tagSum = rowSums(Peak_Co_NLS_S[, c(NLS_E_S)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_NES_M = data.frame(data.frame(annotation = Peak_Co_NES_M$grouped_annotation, 
                                                  tagSum = rowSums(Peak_Co_NES_M[, c(NES_E_M)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_NES_S = data.frame(data.frame(annotation = Peak_Co_NES_S$grouped_annotation, 
                                                  tagSum = rowSums(Peak_Co_NES_S[, c(NES_E_S)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_G3BP_M = data.frame(data.frame(annotation = Peak_Co_G3BP_M$grouped_annotation, 
                                                   tagSum = rowSums(Peak_Co_G3BP_M[, c(G3BP_E_M)])) 
                                        %>% group_by(annotation) 
                                        %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_G3BP_S = data.frame(data.frame(annotation = Peak_Co_G3BP_S$grouped_annotation, 
                                                   tagSum = rowSums(Peak_Co_G3BP_S[, c(G3BP_E_S)])) 
                                        %>% group_by(annotation) 
                                        %>% summarise(tags = sum(tagSum)), row.names = 1)

colnames(PeakDistribution_Co_Input_M) = c('Co_M_Input')
colnames(PeakDistribution_Co_Input_S) = c('Co_S_Input')
colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')
colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')
colnames(PeakDistribution_Co_G3BP_M) = c('Co_M_G3BP')
colnames(PeakDistribution_Co_G3BP_S) = c('Co_S_G3BP')

## Tag Counts Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
                                     PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
                                     PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
                                     PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", 'snoRNA', "ncRNA", "TE", "Other", "DS10K", 'UnAn'))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_M_Input' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels= c('Input', 'NLS', 'NES')) +
  ggtitle('CoCLIP Mock Normalized Tag Counts Distribution') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Fraction Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M/sum(PeakDistribution_Co_Input_M), PeakDistribution_Co_Input_S/sum(PeakDistribution_Co_Input_S), 
                                     PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), 
                                     PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S),
                                     PeakDistribution_Co_G3BP_M/sum(PeakDistribution_Co_G3BP_M), PeakDistribution_Co_G3BP_S/sum(PeakDistribution_Co_G3BP_S))
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", 'snoRNA', "ncRNA", "TE", "Other", "DS10K", 'UnAn'))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_M_Input' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels= c('Input', 'NLS', 'NES')) +
  ggtitle('CoCLIP Mock Normalized Tag Counts Distribution by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
####################################################################################################################

## Exploratory Stacked Bar Plots For CoCLIP Tags: Specific RNA biotypes ONLY
####################################################################################################################
## mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_Co_Input_M_mRNA = Peak_Co_Input_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_Input_S_mRNA = Peak_Co_Input_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

Peak_Co_NLS_M_mRNA = Peak_Co_NLS_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_NLS_S_mRNA = Peak_Co_NLS_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

Peak_Co_NES_M_mRNA = Peak_Co_NES_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_NES_S_mRNA = Peak_Co_NES_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

Peak_Co_G3BP_M_mRNA = Peak_Co_G3BP_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_G3BP_S_mRNA = Peak_Co_G3BP_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

PeakDistribution_Co_Input_M = data.frame(data.frame(annotation = Peak_Co_Input_M_mRNA$finalized_annotation, 
                                                    tagSum = rowSums(Peak_Co_Input_M_mRNA[, c(NLS_I_M, NES_I_M, G3BP_I_M)])) 
                                         %>% group_by(annotation) 
                                         %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_Input_S = data.frame(data.frame(annotation = Peak_Co_Input_S_mRNA$finalized_annotation, 
                                                    tagSum = rowSums(Peak_Co_Input_S_mRNA[, c(NLS_I_S, NES_I_S, G3BP_I_S)])) 
                                         %>% group_by(annotation) 
                                         %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_NLS_M = data.frame(data.frame(annotation = Peak_Co_NLS_M_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NLS_M_mRNA[, c(NLS_E_M)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_NLS_S = data.frame(data.frame(annotation = Peak_Co_NLS_S_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NLS_S_mRNA[, c(NLS_E_S)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_NES_M = data.frame(data.frame(annotation = Peak_Co_NES_M_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NES_M_mRNA[, c(NES_E_M)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)
# PeakDistribution_Co_NES_M = rbind(PeakDistribution_Co_NES_M, c(0))
# row.names(PeakDistribution_Co_NES_M) = c(row.names(PeakDistribution_Co_NES_M)[1:6], 'nC_RI')
# PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M[row.names(PeakDistribution_Co_NLS_S), 'tags', drop = FALSE]

PeakDistribution_Co_NES_S = data.frame(data.frame(annotation = Peak_Co_NES_S_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NES_S_mRNA[, c(NES_E_S)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)
# PeakDistribution_Co_NES_S = rbind(PeakDistribution_Co_NES_S, c(0))
# row.names(PeakDistribution_Co_NES_S) = c(row.names(PeakDistribution_Co_NES_S)[1:6], 'nC_RI')
# PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S[row.names(PeakDistribution_Co_NLS_S), 'tags', drop = FALSE]

PeakDistribution_Co_G3BP_M = data.frame(data.frame(annotation = Peak_Co_G3BP_M_mRNA$finalized_annotation, 
                                                   tagSum = rowSums(Peak_Co_G3BP_M_mRNA[, c(G3BP_E_M)])) 
                                        %>% group_by(annotation) 
                                        %>% summarise(tags = sum(tagSum)), row.names = 1)
# PeakDistribution_Co_G3BP_M = rbind(PeakDistribution_Co_G3BP_M, c(0))
# row.names(PeakDistribution_Co_G3BP_M) = c(row.names(PeakDistribution_Co_G3BP_M)[1:6], 'nC_RI')
# PeakDistribution_Co_G3BP_M = PeakDistribution_Co_G3BP_M[row.names(PeakDistribution_Co_NLS_S), 'tags', drop = FALSE]

PeakDistribution_Co_G3BP_S = data.frame(data.frame(annotation = Peak_Co_G3BP_S_mRNA$finalized_annotation, 
                                                   tagSum = rowSums(Peak_Co_G3BP_S_mRNA[, c(G3BP_E_S)])) 
                                        %>% group_by(annotation) 
                                        %>% summarise(tags = sum(tagSum)), row.names = 1)
# PeakDistribution_Co_G3BP_S = rbind(PeakDistribution_Co_G3BP_S, c(0))
# row.names(PeakDistribution_Co_G3BP_S) = c(row.names(PeakDistribution_Co_G3BP_S)[1:6], 'nC_RI')
# PeakDistribution_Co_G3BP_S = PeakDistribution_Co_G3BP_S[row.names(PeakDistribution_Co_NLS_S), 'tags', drop = FALSE]

colnames(PeakDistribution_Co_Input_M) = c('Co_M_Input')
colnames(PeakDistribution_Co_Input_S) = c('Co_S_Input')
colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')
colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')
colnames(PeakDistribution_Co_G3BP_M) = c('Co_M_G3BP')
colnames(PeakDistribution_Co_G3BP_S) = c('Co_S_G3BP')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
                                     PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
                                     PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
                                     PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_M_Input' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels= c('Input', 'NLS', 'NES')) +
  ggtitle('CoCLIP Mock mRNA Normalized Tag Counts Distribution') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Fraction Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M/sum(PeakDistribution_Co_Input_M), PeakDistribution_Co_Input_S/sum(PeakDistribution_Co_Input_S), 
                                     PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), 
                                     PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S),
                                     PeakDistribution_Co_G3BP_M/sum(PeakDistribution_Co_G3BP_M), PeakDistribution_Co_G3BP_S/sum(PeakDistribution_Co_G3BP_S))
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_M_Input' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels= c('Input', 'NLS', 'NES')) +
  ggtitle('CoCLIP Mock mRNA Normalized Tag Counts Distribution by') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## non-mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_Co_Input_M_Not_mRNA = Peak_Co_Input_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_Input_S_Not_mRNA = Peak_Co_Input_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

Peak_Co_NLS_M_Not_mRNA = Peak_Co_NLS_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_NLS_S_Not_mRNA = Peak_Co_NLS_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

Peak_Co_NES_M_Not_mRNA = Peak_Co_NES_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_NES_S_Not_mRNA = Peak_Co_NES_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

Peak_Co_G3BP_M_Not_mRNA = Peak_Co_G3BP_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_G3BP_S_Not_mRNA = Peak_Co_G3BP_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

PeakDistribution_Co_NLS_M = data.frame(data.frame(annotation = Peak_Co_NLS_M_Not_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NLS_M_Not_mRNA[, c(NLS_E_M)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_NLS_M = rbind(PeakDistribution_Co_NLS_M, c(0))
row.names(PeakDistribution_Co_NLS_M) = c(row.names(PeakDistribution_Co_NLS_M)[1:10], 'scaRNA')

PeakDistribution_Co_NLS_S = data.frame(data.frame(annotation = Peak_Co_NLS_S_Not_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NLS_S_Not_mRNA[, c(NLS_E_S)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_NLS_S = rbind(PeakDistribution_Co_NLS_S, c(0))
row.names(PeakDistribution_Co_NLS_S) = c(row.names(PeakDistribution_Co_NLS_S)[1:10], 'scaRNA')
PeakDistribution_Co_NLS_S = PeakDistribution_Co_NLS_S[row.names(PeakDistribution_Co_NLS_M), 'tags', drop = FALSE]

PeakDistribution_Co_Input_M = data.frame(data.frame(annotation = Peak_Co_Input_M_Not_mRNA$finalized_annotation, 
                                                    tagSum = rowSums(Peak_Co_Input_M_Not_mRNA[, c(NLS_I_M, NES_I_M, G3BP_I_M)])) 
                                         %>% group_by(annotation) 
                                         %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_Input_M = rbind(PeakDistribution_Co_Input_M, c(0), c(0))
row.names(PeakDistribution_Co_Input_M) = c(row.names(PeakDistribution_Co_Input_M)[1:9], 'scaRNA', 'nC_RI')
PeakDistribution_Co_Input_M = PeakDistribution_Co_Input_M[row.names(PeakDistribution_Co_NLS_M), 'tags', drop = FALSE]

PeakDistribution_Co_Input_S = data.frame(data.frame(annotation = Peak_Co_Input_S_Not_mRNA$finalized_annotation, 
                                                    tagSum = rowSums(Peak_Co_Input_S_Not_mRNA[, c(NLS_I_S, NES_I_S, G3BP_I_S)])) 
                                         %>% group_by(annotation) 
                                         %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_Input_S = rbind(PeakDistribution_Co_Input_S, c(0), c(0))
row.names(PeakDistribution_Co_Input_S) = c(row.names(PeakDistribution_Co_Input_S)[1:9], 'scaRNA', 'nC_RI')
PeakDistribution_Co_Input_S = PeakDistribution_Co_Input_S[row.names(PeakDistribution_Co_NLS_M), 'tags', drop = FALSE]

PeakDistribution_Co_NES_M = data.frame(data.frame(annotation = Peak_Co_NES_M_Not_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NES_M_Not_mRNA[, c(NES_E_M)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_NES_M = rbind(PeakDistribution_Co_NES_M, c(0), c(0), c(0))
row.names(PeakDistribution_Co_NES_M) = c(row.names(PeakDistribution_Co_NES_M)[1:8], 'miRNA', 'nC_RI', 'scaRNA')
PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M[row.names(PeakDistribution_Co_NLS_S), 'tags', drop = FALSE]

PeakDistribution_Co_NES_S = data.frame(data.frame(annotation = Peak_Co_NES_S_Not_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NES_S_Not_mRNA[, c(NES_E_S)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_NES_S = rbind(PeakDistribution_Co_NES_S, c(0), c(0))
row.names(PeakDistribution_Co_NES_S) = c(row.names(PeakDistribution_Co_NES_S)[1:9], 'scaRNA', 'nC_RI')
PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S[row.names(PeakDistribution_Co_NLS_M), 'tags', drop = FALSE]

PeakDistribution_Co_G3BP_M = data.frame(data.frame(annotation = Peak_Co_G3BP_M_Not_mRNA$finalized_annotation, 
                                                   tagSum = rowSums(Peak_Co_G3BP_M_Not_mRNA[, c(G3BP_E_M)])) 
                                        %>% group_by(annotation) 
                                        %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_G3BP_M = rbind(PeakDistribution_Co_G3BP_M, c(0), c(0), c(0))
row.names(PeakDistribution_Co_G3BP_M) = c(row.names(PeakDistribution_Co_G3BP_M)[1:8], 'miRNA', 'nC_RI', 'scaRNA')
PeakDistribution_Co_G3BP_M = PeakDistribution_Co_G3BP_M[row.names(PeakDistribution_Co_NLS_S), 'tags', drop = FALSE]

PeakDistribution_Co_G3BP_S = data.frame(data.frame(annotation = Peak_Co_G3BP_S_Not_mRNA$finalized_annotation, 
                                                   tagSum = rowSums(Peak_Co_G3BP_S_Not_mRNA[, c(G3BP_E_S)])) 
                                        %>% group_by(annotation) 
                                        %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_G3BP_S = rbind(PeakDistribution_Co_G3BP_S, c(0), c(0), c(0))
row.names(PeakDistribution_Co_G3BP_S) = c(row.names(PeakDistribution_Co_G3BP_S)[1:8], 'miRNA', 'nC_RI', 'scaRNA')
PeakDistribution_Co_G3BP_S = PeakDistribution_Co_G3BP_S[row.names(PeakDistribution_Co_NLS_S), 'tags', drop = FALSE]


colnames(PeakDistribution_Co_Input_M) = c('Co_M_Input')
colnames(PeakDistribution_Co_Input_S) = c('Co_S_Input')
colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')
colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')
colnames(PeakDistribution_Co_G3BP_M) = c('Co_M_G3BP')
colnames(PeakDistribution_Co_G3BP_S) = c('Co_S_G3BP')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
                                     PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
                                     PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
                                     PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_M_Input' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels= c('Input', 'NLS', 'NES')) +
  ggtitle('CoCLIP Mock ncRNA Normalized Tag Counts Distribution') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Fraction Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M/sum(PeakDistribution_Co_Input_M), PeakDistribution_Co_Input_S/sum(PeakDistribution_Co_Input_S), 
                                     PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), 
                                     PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S),
                                     PeakDistribution_Co_G3BP_M/sum(PeakDistribution_Co_G3BP_M), PeakDistribution_Co_G3BP_S/sum(PeakDistribution_Co_G3BP_S))
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))

ggplot(PeakDistribution_Co_combined %>% filter(Source == 'Co_M_Input' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  # scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  scale_x_discrete(labels= c('Input', 'NLS', 'NES')) +
  ggtitle('CoCLIP Mock ncRNA Normalized Tag Counts Distribution by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
####################################################################################################################

## Exploratory Stacked Bar Plots For CoCLIP average Tags:
####################################################################################################################
# ## Mock vs Stress for Each Fraction
# PeakDistribution_Co_Input_M = data.frame(data.frame(annotation = avgTagCounts$grouped_annotation, 
#                                                     AvgTagCounts = rowSums(avgTagCounts[, c('NLS_I_M', 'NES_I_M', 'G3BP_I_M')])/3,
#                                                     BC = rowSums(avgTagCounts[, c('NLS_I_M_BC', 'NES_I_M_BC', 'G3BP_I_M_BC')]))
#                                          %>% filter(BC >= BC_Threshold_I*3)
#                                          %>% group_by(annotation) 
#                                          %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_Input_S = data.frame(data.frame(annotation = avgTagCounts$grouped_annotation, 
#                                                     AvgTagCounts = rowSums(avgTagCounts[, c('NLS_I_S', 'NES_I_S', 'G3BP_I_S')])/3,
#                                                     BC = rowSums(avgTagCounts[, c('NLS_I_S_BC', 'NES_I_S_BC', 'G3BP_I_S_BC')]))
#                                          %>% filter(BC >= BC_Threshold_I*3)
#                                          %>% group_by(annotation) 
#                                          %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_NLS_M = data.frame(data.frame(annotation = avgTagCounts$grouped_annotation, 
#                                                   AvgTagCounts = avgTagCounts$NLS_E_M,
#                                                   BC = avgTagCounts$NLS_E_M_BC)
#                                        %>% filter(BC >= BC_Threshold_E)
#                                        %>% group_by(annotation) 
#                                        %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_NLS_S = data.frame(data.frame(annotation = avgTagCounts$grouped_annotation, 
#                                                   AvgTagCounts = avgTagCounts$NLS_E_S,
#                                                   BC = avgTagCounts$NLS_E_S_BC)
#                                        %>% filter(BC >= BC_Threshold_E)
#                                        %>% group_by(annotation) 
#                                        %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_NES_M = data.frame(data.frame(annotation = avgTagCounts$grouped_annotation, 
#                                                   AvgTagCounts = avgTagCounts$NES_E_M,
#                                                   BC = avgTagCounts$NES_E_M_BC)
#                                        %>% filter(BC >= BC_Threshold_E)
#                                        %>% group_by(annotation) 
#                                        %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_NES_S = data.frame(data.frame(annotation = avgTagCounts$grouped_annotation, 
#                                                   AvgTagCounts = avgTagCounts$NES_E_S,
#                                                   BC = avgTagCounts$NES_E_S_BC)
#                                        %>% filter(BC >= BC_Threshold_E)
#                                        %>% group_by(annotation) 
#                                        %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_G3BP_M = data.frame(data.frame(annotation = avgTagCounts$grouped_annotation, 
#                                                   AvgTagCounts = avgTagCounts$G3BP_E_M,
#                                                   BC = avgTagCounts$G3BP_E_M_BC)
#                                        %>% filter(BC >= BC_Threshold_E)
#                                        %>% group_by(annotation) 
#                                        %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_G3BP_S = data.frame(data.frame(annotation = avgTagCounts$grouped_annotation, 
#                                                   AvgTagCounts = avgTagCounts$G3BP_E_S,
#                                                   BC = avgTagCounts$G3BP_E_S_BC)
#                                        %>% filter(BC >= BC_Threshold_E)
#                                        %>% group_by(annotation) 
#                                        %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# colnames(PeakDistribution_Co_Input_M) = c('Co_M_Input')
# colnames(PeakDistribution_Co_Input_S) = c('Co_S_Input')
# colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
# colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
# colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')
# colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')
# colnames(PeakDistribution_Co_G3BP_M) = c('Co_M_G3BP')
# colnames(PeakDistribution_Co_G3BP_S) = c('Co_S_G3BP')
# 
# ## Tag Counts Distribution Stacked Bar Graph
# PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
#                                      PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
#                                      PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
#                                      PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
# PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)
# 
# PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
#   gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
#   select(Source, Freq, Annotation)
# 
# PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
# PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", 'snoRNA', "ncRNA", "TE", "Other", "DS10K", 'UnAn'))
# 
# ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
#   geom_bar(position='stack', stat='identity') +
#   scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
#   ggtitle('Peak Density Distributions') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_brewer(palette = "Set3")
# 
# ## Fraction Distribution Stacked Bar Graph
# PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M/sum(PeakDistribution_Co_Input_M), PeakDistribution_Co_Input_S/sum(PeakDistribution_Co_Input_S), 
#                                      PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), 
#                                      PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S),
#                                      PeakDistribution_Co_G3BP_M/sum(PeakDistribution_Co_G3BP_M), PeakDistribution_Co_G3BP_S/sum(PeakDistribution_Co_G3BP_S))
# PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)
# 
# PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
#   gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
#   select(Source, Freq, Annotation)
# 
# PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
# PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", 'snoRNA', "ncRNA", "TE", "Other", "DS10K", 'UnAn'))
# 
# ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
#   geom_bar(position='stack', stat='identity') +
#   scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
#   ggtitle('Peak Density Distributions by Fraction') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_brewer(palette = "Set3")
####################################################################################################################

## Exploratory Stacked Bar Plots For CoCLIP average Tags: Specific RNA biotypes ONLY
####################################################################################################################
# ## mRNA Features Only
# ## Mock vs Stress for Each Fraction
# avgTagCounts_mRNA = avgTagCounts %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == "CDS_RI" | finalized_annotation == 'DS10K')
# 
# PeakDistribution_Co_Input_M = data.frame(data.frame(annotation = avgTagCounts_mRNA$finalized_annotation, 
#                                                     AvgTagCounts = rowSums(avgTagCounts_mRNA[, c('NLS_I_M', 'NES_I_M', 'G3BP_I_M')])/3,
#                                                     BC = rowSums(avgTagCounts_mRNA[, c('NLS_I_M_BC', 'NES_I_M_BC', 'G3BP_I_M_BC')]))
#                                          %>% filter(BC >= BC_Threshold_I*3)
#                                          %>% group_by(annotation) 
#                                          %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_Input_S = data.frame(data.frame(annotation = avgTagCounts_mRNA$finalized_annotation, 
#                                                     AvgTagCounts = rowSums(avgTagCounts_mRNA[, c('NLS_I_S', 'NES_I_S', 'G3BP_I_S')])/3,
#                                                     BC = rowSums(avgTagCounts_mRNA[, c('NLS_I_S_BC', 'NES_I_S_BC', 'G3BP_I_S_BC')]))
#                                          %>% filter(BC >= BC_Threshold_I*3)
#                                          %>% group_by(annotation) 
#                                          %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_NLS_M = data.frame(data.frame(annotation = avgTagCounts_mRNA$finalized_annotation, 
#                                                   AvgTagCounts = avgTagCounts_mRNA$NLS_E_M,
#                                                   BC = avgTagCounts_mRNA$NLS_E_M_BC)
#                                        %>% filter(BC >= BC_Threshold_E)
#                                        %>% group_by(annotation) 
#                                        %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_NLS_S = data.frame(data.frame(annotation = avgTagCounts_mRNA$finalized_annotation, 
#                                                   AvgTagCounts = avgTagCounts_mRNA$NLS_E_S,
#                                                   BC = avgTagCounts_mRNA$NLS_E_S_BC)
#                                        %>% filter(BC >= BC_Threshold_E)
#                                        %>% group_by(annotation) 
#                                        %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_NES_M = data.frame(data.frame(annotation = avgTagCounts_mRNA$finalized_annotation, 
#                                                   AvgTagCounts = avgTagCounts_mRNA$NES_E_M,
#                                                   BC = avgTagCounts_mRNA$NES_E_M_BC)
#                                        %>% filter(BC >= BC_Threshold_E)
#                                        %>% group_by(annotation) 
#                                        %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_NES_S = data.frame(data.frame(annotation = avgTagCounts_mRNA$finalized_annotation, 
#                                                   AvgTagCounts = avgTagCounts_mRNA$NES_E_S,
#                                                   BC = avgTagCounts_mRNA$NES_E_S_BC)
#                                        %>% filter(BC >= BC_Threshold_E)
#                                        %>% group_by(annotation) 
#                                        %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_G3BP_M = data.frame(data.frame(annotation = avgTagCounts_mRNA$finalized_annotation, 
#                                                    AvgTagCounts = avgTagCounts_mRNA$G3BP_E_M,
#                                                    BC = avgTagCounts_mRNA$G3BP_E_M_BC)
#                                         %>% filter(BC >= BC_Threshold_E)
#                                         %>% group_by(annotation) 
#                                         %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_G3BP_S = data.frame(data.frame(annotation = avgTagCounts_mRNA$finalized_annotation, 
#                                                    AvgTagCounts = avgTagCounts_mRNA$G3BP_E_S,
#                                                    BC = avgTagCounts_mRNA$G3BP_E_S_BC)
#                                         %>% filter(BC >= BC_Threshold_E)
#                                         %>% group_by(annotation) 
#                                         %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# colnames(PeakDistribution_Co_Input_M) = c('Co_M_Input')
# colnames(PeakDistribution_Co_Input_S) = c('Co_S_Input')
# colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
# colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
# colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')
# colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')
# colnames(PeakDistribution_Co_G3BP_M) = c('Co_M_G3BP')
# colnames(PeakDistribution_Co_G3BP_S) = c('Co_S_G3BP')
# 
# ## Tag Counts Distribution Stacked Bar Graph
# PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
#                                      PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
#                                      PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
#                                      PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
# PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)
# 
# PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
#   gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
#   select(Source, Freq, Annotation)
# 
# PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
# PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))
# 
# ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
#   geom_bar(position='stack', stat='identity') +
#   scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
#   ggtitle('mRNA Peak Density Distributions') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_brewer(palette = "Set3")
# 
# ## Fraction Distribution Stacked Bar Graph
# PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M/sum(PeakDistribution_Co_Input_M), PeakDistribution_Co_Input_S/sum(PeakDistribution_Co_Input_S), 
#                                      PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), 
#                                      PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S),
#                                      PeakDistribution_Co_G3BP_M/sum(PeakDistribution_Co_G3BP_M), PeakDistribution_Co_G3BP_S/sum(PeakDistribution_Co_G3BP_S))
# PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)
# 
# PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
#   gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
#   select(Source, Freq, Annotation)
# 
# PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
# PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))
# 
# ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
#   geom_bar(position='stack', stat='identity') +
#   scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
#   ggtitle('mRNA Peak Density Distributions by Fraction') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_brewer(palette = "Set3")
# 
# ## nonmRNA Features Only
# ## Mock vs Stress for Each Fraction
# avgTagCounts_Not_mRNA = avgTagCounts %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
# 
# PeakDistribution_Co_Input_M = data.frame(data.frame(annotation = avgTagCounts_Not_mRNA$finalized_annotation, 
#                                                     AvgTagCounts = rowSums(avgTagCounts_Not_mRNA[, c('NLS_I_M', 'NES_I_M', 'G3BP_I_M')])/3,
#                                                     BC = rowSums(avgTagCounts_Not_mRNA[, c('NLS_I_M_BC', 'NES_I_M_BC', 'G3BP_I_M_BC')]))
#                                          %>% filter(BC >= BC_Threshold_I*3)
#                                          %>% group_by(annotation) 
#                                          %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_Input_S = data.frame(data.frame(annotation = avgTagCounts_Not_mRNA$finalized_annotation, 
#                                                     AvgTagCounts = rowSums(avgTagCounts_Not_mRNA[, c('NLS_I_S', 'NES_I_S', 'G3BP_I_S')])/3,
#                                                     BC = rowSums(avgTagCounts_Not_mRNA[, c('NLS_I_S_BC', 'NES_I_S_BC', 'G3BP_I_S_BC')]))
#                                          %>% filter(BC >= BC_Threshold_I*3)
#                                          %>% group_by(annotation) 
#                                          %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_NLS_M = data.frame(data.frame(annotation = avgTagCounts_Not_mRNA$finalized_annotation, 
#                                                   AvgTagCounts = avgTagCounts_Not_mRNA$NLS_E_M,
#                                                   BC = avgTagCounts_Not_mRNA$NLS_E_M_BC)
#                                        %>% filter(BC >= BC_Threshold_E)
#                                        %>% group_by(annotation) 
#                                        %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# 
# PeakDistribution_Co_NLS_S = data.frame(data.frame(annotation = avgTagCounts_Not_mRNA$finalized_annotation, 
#                                                   AvgTagCounts = avgTagCounts_Not_mRNA$NLS_E_S,
#                                                   BC = avgTagCounts_Not_mRNA$NLS_E_S_BC)
#                                        %>% filter(BC >= BC_Threshold_E)
#                                        %>% group_by(annotation) 
#                                        %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# PeakDistribution_Co_NLS_S = rbind(PeakDistribution_Co_NLS_S, c(0))
# row.names(PeakDistribution_Co_NLS_S) = c(row.names(PeakDistribution_Co_NLS_S)[1:10], 'scaRNA')
# PeakDistribution_Co_NLS_S = PeakDistribution_Co_NLS_S[row.names(PeakDistribution_Co_NLS_M), 'tags', drop = FALSE]
# 
# 
# PeakDistribution_Co_NES_M = data.frame(data.frame(annotation = avgTagCounts_Not_mRNA$finalized_annotation, 
#                                                   AvgTagCounts = avgTagCounts_Not_mRNA$NES_E_M,
#                                                   BC = avgTagCounts_Not_mRNA$NES_E_M_BC)
#                                        %>% filter(BC >= BC_Threshold_E)
#                                        %>% group_by(annotation) 
#                                        %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# PeakDistribution_Co_NES_M = rbind(PeakDistribution_Co_NES_M, c(0), c(0), c(0))
# row.names(PeakDistribution_Co_NES_M) = c(row.names(PeakDistribution_Co_NES_M)[1:8], 'scaRNA', 'nC_RI', 'miRNA')
# PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M[row.names(PeakDistribution_Co_NLS_M), 'tags', drop = FALSE]
# 
# PeakDistribution_Co_NES_S = data.frame(data.frame(annotation = avgTagCounts_Not_mRNA$finalized_annotation, 
#                                                   AvgTagCounts = avgTagCounts_Not_mRNA$NES_E_S,
#                                                   BC = avgTagCounts_Not_mRNA$NES_E_S_BC)
#                                        %>% filter(BC >= BC_Threshold_E)
#                                        %>% group_by(annotation) 
#                                        %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# PeakDistribution_Co_NES_S = rbind(PeakDistribution_Co_NES_S, c(0), c(0))
# row.names(PeakDistribution_Co_NES_S) = c(row.names(PeakDistribution_Co_NES_S)[1:9], 'scaRNA', 'nC_RI')
# PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S[row.names(PeakDistribution_Co_NLS_M), 'tags', drop = FALSE]
# 
# PeakDistribution_Co_G3BP_M = data.frame(data.frame(annotation = avgTagCounts_Not_mRNA$finalized_annotation, 
#                                                    AvgTagCounts = avgTagCounts_Not_mRNA$G3BP_E_M,
#                                                    BC = avgTagCounts_Not_mRNA$G3BP_E_M_BC)
#                                         %>% filter(BC >= BC_Threshold_E)
#                                         %>% group_by(annotation) 
#                                         %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# PeakDistribution_Co_G3BP_M = rbind(PeakDistribution_Co_G3BP_M, c(0), c(0))
# row.names(PeakDistribution_Co_G3BP_M) = c(row.names(PeakDistribution_Co_G3BP_M)[1:9], 'scaRNA', 'miRNA')
# PeakDistribution_Co_G3BP_M = PeakDistribution_Co_G3BP_M[row.names(PeakDistribution_Co_NLS_M), 'tags', drop = FALSE]
# 
# PeakDistribution_Co_G3BP_S = data.frame(data.frame(annotation = avgTagCounts_Not_mRNA$finalized_annotation, 
#                                                    AvgTagCounts = avgTagCounts_Not_mRNA$G3BP_E_S,
#                                                    BC = avgTagCounts_Not_mRNA$G3BP_E_S_BC)
#                                         %>% filter(BC >= BC_Threshold_E)
#                                         %>% group_by(annotation) 
#                                         %>% summarise(tags = sum(AvgTagCounts)), row.names = 1)
# PeakDistribution_Co_G3BP_S = rbind(PeakDistribution_Co_G3BP_S, c(0), c(0))
# row.names(PeakDistribution_Co_G3BP_S) = c(row.names(PeakDistribution_Co_G3BP_S)[1:9], 'scaRNA', 'miRNA')
# PeakDistribution_Co_G3BP_S = PeakDistribution_Co_G3BP_S[row.names(PeakDistribution_Co_NLS_M), 'tags', drop = FALSE]
# 
# colnames(PeakDistribution_Co_Input_M) = c('Co_M_Input')
# colnames(PeakDistribution_Co_Input_S) = c('Co_S_Input')
# colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
# colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
# colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')
# colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')
# colnames(PeakDistribution_Co_G3BP_M) = c('Co_M_G3BP')
# colnames(PeakDistribution_Co_G3BP_S) = c('Co_S_G3BP')
# 
# ## Tag Counts Distribution Stacked Bar Graph
# PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
#                                      PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
#                                      PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
#                                      PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
# PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)
# 
# PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
#   gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
#   select(Source, Freq, Annotation)
# 
# PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
# PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))
# 
# ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
#   geom_bar(position='stack', stat='identity') +
#   scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
#   ggtitle('non-mRNA Peak Density Distributions') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_brewer(palette = "Set3")
# 
# ## Fraction Distribution Stacked Bar Graph
# PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M/sum(PeakDistribution_Co_Input_M), PeakDistribution_Co_Input_S/sum(PeakDistribution_Co_Input_S), 
#                                      PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), 
#                                      PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S),
#                                      PeakDistribution_Co_G3BP_M/sum(PeakDistribution_Co_G3BP_M), PeakDistribution_Co_G3BP_S/sum(PeakDistribution_Co_G3BP_S))
# PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)
# 
# PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
#   gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
#   select(Source, Freq, Annotation)
# 
# PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
# PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))
# 
# ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
#   geom_bar(position='stack', stat='identity') +
#   scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
#   ggtitle('non-mRNA Peak Density Distributions by Fraction') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_brewer(palette = "Set3")
####################################################################################################################

## FIGURE2 Exploratory Stacked Bar Plots for Fraction vs CoCLIP Side-by-Side Peaks:
####################################################################################################################
PeakDistribution_F_Nuc_M = data.frame(table(Peak_F_Nuc_M$grouped_annotation), row.names = 1)
PeakDistribution_F_Cyt_M = data.frame(table(Peak_F_Cyt_M$grouped_annotation), row.names = 1)
PeakDistribution_Co_NLS_M = data.frame(table(Peak_Co_NLS_M$grouped_annotation), row.names = 1)
PeakDistribution_Co_NES_M = data.frame(table(Peak_Co_NES_M$grouped_annotation), row.names = 1)

PeakDistribution_F_Nuc_S = data.frame(table(Peak_F_Nuc_S$grouped_annotation), row.names = 1)
PeakDistribution_F_Cyt_S = data.frame(table(Peak_F_Cyt_S$grouped_annotation), row.names = 1)
PeakDistribution_Co_NLS_S = data.frame(table(Peak_Co_NLS_S$grouped_annotation), row.names = 1)
PeakDistribution_Co_NES_S = data.frame(table(Peak_Co_NES_S$grouped_annotation), row.names = 1)

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')

colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')
colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')

PeakDistribution_F_Nuc_M = PeakDistribution_F_Nuc_M %>% slice(-which(rownames(PeakDistribution_F_Nuc_M) == 'UnAn'))
PeakDistribution_F_Cyt_M = PeakDistribution_F_Cyt_M %>% slice(-which(rownames(PeakDistribution_F_Cyt_M) == 'UnAn'))
PeakDistribution_Co_NLS_M = PeakDistribution_Co_NLS_M %>% slice(-which(rownames(PeakDistribution_Co_NLS_M) == 'UnAn'))
PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M %>% slice(-which(rownames(PeakDistribution_Co_NES_M) == 'UnAn'))
PeakDistribution_F_Nuc_S = PeakDistribution_F_Nuc_S %>% slice(-which(rownames(PeakDistribution_F_Nuc_S) == 'UnAn'))
PeakDistribution_F_Cyt_S = PeakDistribution_F_Cyt_S %>% slice(-which(rownames(PeakDistribution_F_Cyt_S) == 'UnAn'))
PeakDistribution_Co_NLS_S = PeakDistribution_Co_NLS_S %>% slice(-which(rownames(PeakDistribution_Co_NLS_S) == 'UnAn'))
PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S %>% slice(-which(rownames(PeakDistribution_Co_NES_S) == 'UnAn'))

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Cyt_M, 
                                  PeakDistribution_Co_NLS_M, PeakDistribution_Co_NES_M,
                                  PeakDistribution_F_Nuc_S, PeakDistribution_F_Cyt_S,
                                  PeakDistribution_Co_NLS_S, PeakDistribution_Co_NES_S)

PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES') %>%
  select(Source, Freq, Annotation)

PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = c( 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'))
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K", 'UnAn'))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_M_Nuc' | Source == 'F_M_Cyt' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP')) +
  ggtitle('Mock Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_S_Nuc' | Source == 'F_S_Cyt' | Source == 'Co_S_NLS' | Source == 'Co_S_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP')) +
  ggtitle('Arsenite Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Fraction Distribution Stacked Bar Graph
PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), 
                                  PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M),
                                  PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S),
                                  PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S))
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES') %>%
  select(Source, Freq, Annotation)

PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = c( 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'))
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K", 'UnAn'))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_M_Nuc' | Source == 'F_M_Cyt' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP')) +
  ggtitle('Mock Peak Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_S_Nuc' | Source == 'F_S_Cyt' | Source == 'Co_S_NLS' | Source == 'Co_S_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP')) +
  ggtitle('Arsenite Peak Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

####################################################################################################################

## FIGURE2 Exploratory Stacked Bar Plots for Fraction vs CoCLIP Side-by-Side Peaks: Specific RNA biotypes ONLY
####################################################################################################################
## mRNA
Peak_F_Nuc_M_mRNA = Peak_F_Nuc_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_F_Cyt_M_mRNA = Peak_F_Cyt_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_NLS_M_mRNA = Peak_Co_NLS_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_NES_M_mRNA = Peak_Co_NES_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

Peak_F_Nuc_S_mRNA = Peak_F_Nuc_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_F_Cyt_S_mRNA = Peak_F_Cyt_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_NLS_S_mRNA = Peak_Co_NLS_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_NES_S_mRNA = Peak_Co_NES_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

PeakDistribution_F_Nuc_M = data.frame(table(Peak_F_Nuc_M_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_Cyt_M = data.frame(table(Peak_F_Cyt_M_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NLS_M = data.frame(table(Peak_Co_NLS_M_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NES_M = data.frame(table(Peak_Co_NES_M_mRNA$finalized_annotation), row.names = 1)

PeakDistribution_F_Nuc_S = data.frame(table(Peak_F_Nuc_S_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_Cyt_S = data.frame(table(Peak_F_Cyt_S_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NLS_S = data.frame(table(Peak_Co_NLS_S_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NES_S = data.frame(table(Peak_Co_NES_S_mRNA$finalized_annotation), row.names = 1)

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')

colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')
colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Cyt_M, 
                                  PeakDistribution_Co_NLS_M, PeakDistribution_Co_NES_M,
                                  PeakDistribution_F_Nuc_S, PeakDistribution_F_Cyt_S,
                                  PeakDistribution_Co_NLS_S, PeakDistribution_Co_NES_S)

PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES') %>%
  select(Source, Freq, Annotation)

PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = c( 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'))
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_M_Nuc' | Source == 'F_M_Cyt' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP')) +
  ggtitle('Mock mRNA Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_S_Nuc' | Source == 'F_S_Cyt' | Source == 'Co_S_NLS' | Source == 'Co_S_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP')) +
  ggtitle('Arsenite mRNA Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Fraction Distribution Stacked Bar Graph
PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), 
                                  PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M),
                                  PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S),
                                  PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S))
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES') %>%
  select(Source, Freq, Annotation)

PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = c( 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'))
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_M_Nuc' | Source == 'F_M_Cyt' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP')) +
  ggtitle('Mock mRNA Peak Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_S_Nuc' | Source == 'F_S_Cyt' | Source == 'Co_S_NLS' | Source == 'Co_S_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP')) +
  ggtitle('Arsenite mRNA Peak Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## ncRNA
Peak_F_Nuc_M_Not_mRNA = Peak_F_Nuc_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_F_Cyt_M_Not_mRNA = Peak_F_Cyt_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_NLS_M_Not_mRNA = Peak_Co_NLS_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_NES_M_Not_mRNA = Peak_Co_NES_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

Peak_F_Nuc_S_Not_mRNA = Peak_F_Nuc_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_F_Cyt_S_Not_mRNA = Peak_F_Cyt_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_NLS_S_Not_mRNA = Peak_Co_NLS_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_NES_S_Not_mRNA = Peak_Co_NES_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

PeakDistribution_F_Nuc_M = data.frame(table(Peak_F_Nuc_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_Nuc_M = rbind(PeakDistribution_F_Nuc_M, c(0))
row.names(PeakDistribution_F_Nuc_M) = c(row.names(PeakDistribution_F_Nuc_M)[1:10], 'scaRNA')

PeakDistribution_F_Cyt_M = data.frame(table(Peak_F_Cyt_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_Cyt_M = rbind(PeakDistribution_F_Cyt_M, c(0))
row.names(PeakDistribution_F_Cyt_M) = c(row.names(PeakDistribution_F_Cyt_M)[1:10], 'scaRNA')
PeakDistribution_F_Cyt_M = PeakDistribution_F_Cyt_M[row.names(PeakDistribution_F_Nuc_M), 'Freq', drop = FALSE]

PeakDistribution_Co_NLS_M = data.frame(table(Peak_Co_NLS_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NLS_M = rbind(PeakDistribution_Co_NLS_M, c(0))
row.names(PeakDistribution_Co_NLS_M) = c(row.names(PeakDistribution_Co_NLS_M)[1:10], 'scaRNA')
PeakDistribution_Co_NLS_M = PeakDistribution_Co_NLS_M[row.names(PeakDistribution_F_Nuc_M), 'Freq', drop = FALSE]

PeakDistribution_Co_NES_M = data.frame(table(Peak_Co_NES_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NES_M = rbind(PeakDistribution_Co_NES_M, c(0), c(0), c(0))
row.names(PeakDistribution_Co_NES_M) = c(row.names(PeakDistribution_Co_NES_M)[1:8], 'miRNA', 'nC_RI', 'scaRNA')
PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M[row.names(PeakDistribution_F_Nuc_M), 'Freq', drop = FALSE]

PeakDistribution_F_Nuc_S = data.frame(table(Peak_F_Nuc_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_Nuc_S = rbind(PeakDistribution_F_Nuc_S, c(0))
row.names(PeakDistribution_F_Nuc_S) = c(row.names(PeakDistribution_F_Nuc_S)[1:10], 'scaRNA')
PeakDistribution_F_Nuc_S = PeakDistribution_F_Nuc_S[row.names(PeakDistribution_F_Nuc_M), 'Freq', drop = FALSE]

PeakDistribution_F_Cyt_S = data.frame(table(Peak_F_Cyt_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_Cyt_S = rbind(PeakDistribution_F_Cyt_S, c(0), c(0))
row.names(PeakDistribution_F_Cyt_S) = c(row.names(PeakDistribution_F_Cyt_S)[1:9], 'miRNA', 'scaRNA')
PeakDistribution_F_Cyt_S = PeakDistribution_F_Cyt_S[row.names(PeakDistribution_F_Nuc_M), 'Freq', drop = FALSE]

PeakDistribution_Co_NLS_S = data.frame(table(Peak_Co_NLS_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NLS_S = rbind(PeakDistribution_Co_NLS_S, c(0))
row.names(PeakDistribution_Co_NLS_S) = c(row.names(PeakDistribution_Co_NLS_S)[1:10], 'scaRNA')
PeakDistribution_Co_NLS_S = PeakDistribution_Co_NLS_S[row.names(PeakDistribution_F_Nuc_M), 'Freq', drop = FALSE]

PeakDistribution_Co_NES_S = data.frame(table(Peak_Co_NES_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NES_S = rbind(PeakDistribution_Co_NES_S, c(0), c(0))
row.names(PeakDistribution_Co_NES_S) = c(row.names(PeakDistribution_Co_NES_S)[1:9], 'scaRNA', 'nC_RI')
PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S[row.names(PeakDistribution_F_Nuc_M), 'Freq', drop = FALSE]

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')

colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')
colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')


PeakDistribution_F_Nuc_M = PeakDistribution_F_Nuc_M %>% slice(-which(rownames(PeakDistribution_F_Nuc_M) == 'UnAn'))
PeakDistribution_F_Cyt_M = PeakDistribution_F_Cyt_M %>% slice(-which(rownames(PeakDistribution_F_Cyt_M) == 'UnAn'))
PeakDistribution_Co_NLS_M = PeakDistribution_Co_NLS_M %>% slice(-which(rownames(PeakDistribution_Co_NLS_M) == 'UnAn'))
PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M %>% slice(-which(rownames(PeakDistribution_Co_NES_M) == 'UnAn'))
PeakDistribution_F_Nuc_S = PeakDistribution_F_Nuc_S %>% slice(-which(rownames(PeakDistribution_F_Nuc_S) == 'UnAn'))
PeakDistribution_F_Cyt_S = PeakDistribution_F_Cyt_S %>% slice(-which(rownames(PeakDistribution_F_Cyt_S) == 'UnAn'))
PeakDistribution_Co_NLS_S = PeakDistribution_Co_NLS_S %>% slice(-which(rownames(PeakDistribution_Co_NLS_S) == 'UnAn'))
PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S %>% slice(-which(rownames(PeakDistribution_Co_NES_S) == 'UnAn'))

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Cyt_M, 
                                  PeakDistribution_Co_NLS_M, PeakDistribution_Co_NES_M,
                                  PeakDistribution_F_Nuc_S, PeakDistribution_F_Cyt_S,
                                  PeakDistribution_Co_NLS_S, PeakDistribution_Co_NES_S)

PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES') %>%
  select(Source, Freq, Annotation)

PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = c( 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'))
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_M_Nuc' | Source == 'F_M_Cyt' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP')) +
  ggtitle('Mock ncRNA Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_S_Nuc' | Source == 'F_S_Cyt' | Source == 'Co_S_NLS' | Source == 'Co_S_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP')) +
  ggtitle('Arsenite ncRNA Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Fraction Distribution Stacked Bar Graph
PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), 
                                  PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M),
                                  PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S),
                                  PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S))
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES') %>%
  select(Source, Freq, Annotation)

PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = c( 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'))
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_M_Nuc' | Source == 'F_M_Cyt' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP')) +
  ggtitle('Mock ncRNA Peak Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_S_Nuc' | Source == 'F_S_Cyt' | Source == 'Co_S_NLS' | Source == 'Co_S_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP')) +
  ggtitle('Arsenite ncRNA Peak Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

####################################################################################################################

## FIGURE2 Exploratory Stacked Bar Plots for Fraction vs CoCLIP Side-by-Side Tags:
####################################################################################################################
PeakDistribution_F_Nuc_M = data.frame(data.frame(annotation = Peak_F_Nuc_M$grouped_annotation, 
                                                 tagSum = rowSums(Peak_F_Nuc_M[, Nuc_F_M])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_Cyt_M = data.frame(data.frame(annotation = Peak_F_Cyt_M$grouped_annotation, 
                                                 tagSum = rowSums(Peak_F_Cyt_M[, Cyto_F_M])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_NLS_M = data.frame(data.frame(annotation = Peak_Co_NLS_M$grouped_annotation, 
                                                  tagSum = rowSums(Peak_Co_NLS_M[, c(NLS_E_M)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_NES_M = data.frame(data.frame(annotation = Peak_Co_NES_M$grouped_annotation, 
                                                  tagSum = rowSums(Peak_Co_NES_M[, c(NES_E_M)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_Nuc_S = data.frame(data.frame(annotation = Peak_F_Nuc_S$grouped_annotation, 
                                                 tagSum = rowSums(Peak_F_Nuc_S[, Nuc_F_S])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_Cyt_S = data.frame(data.frame(annotation = Peak_F_Cyt_S$grouped_annotation, 
                                                 tagSum = rowSums(Peak_F_Cyt_S[, Cyto_F_S])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_NLS_S = data.frame(data.frame(annotation = Peak_Co_NLS_S$grouped_annotation, 
                                                  tagSum = rowSums(Peak_Co_NLS_S[, c(NLS_E_S)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)


PeakDistribution_Co_NES_S = data.frame(data.frame(annotation = Peak_Co_NES_S$grouped_annotation, 
                                                  tagSum = rowSums(Peak_Co_NES_S[, c(NES_E_S)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')

colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')
colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')


## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Cyt_M, 
                                  PeakDistribution_Co_NLS_M, PeakDistribution_Co_NES_M,
                                  PeakDistribution_F_Nuc_S, PeakDistribution_F_Cyt_S,
                                  PeakDistribution_Co_NLS_S, PeakDistribution_Co_NES_S)

PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES') %>%
  select(Source, Freq, Annotation)

PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = c( 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'))
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K", 'UnAn'))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_M_Nuc' | Source == 'F_M_Cyt' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP')) +
  ggtitle('Mock Normalized Tag Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_S_Nuc' | Source == 'F_S_Cyt' | Source == 'Co_S_NLS' | Source == 'Co_S_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP')) +
  ggtitle('Arsenite Normalized Tag Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Fraction Distribution Stacked Bar Graph
PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), 
                                  PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M),
                                  PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S),
                                  PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S))
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES') %>%
  select(Source, Freq, Annotation)

PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = c( 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'))
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K", 'UnAn'))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_M_Nuc' | Source == 'F_M_Cyt' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP')) +
  ggtitle('Mock Normalized Tag Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_S_Nuc' | Source == 'F_S_Cyt' | Source == 'Co_S_NLS' | Source == 'Co_S_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP')) +
  ggtitle('Arsenite Normalized Tag Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

####################################################################################################################

## FIGURE2 Exploratory Stacked Bar Plots for Fraction vs CoCLIP Side-by-Side Tags: Specific RNA biotypes ONLY
####################################################################################################################
## mRNA
Peak_F_Nuc_M_mRNA = Peak_F_Nuc_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_F_Cyt_M_mRNA = Peak_F_Cyt_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_NLS_M_mRNA = Peak_Co_NLS_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_NES_M_mRNA = Peak_Co_NES_M %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

Peak_F_Nuc_S_mRNA = Peak_F_Nuc_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_F_Cyt_S_mRNA = Peak_F_Cyt_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_NLS_S_mRNA = Peak_Co_NLS_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')
Peak_Co_NES_S_mRNA = Peak_Co_NES_S %>% filter(finalized_annotation == "3'UTR" | finalized_annotation == "CDS" | finalized_annotation == "5'UTR" | finalized_annotation == "intron" | finalized_annotation == 'CDS_RI' | finalized_annotation == 'DS10K')

PeakDistribution_F_Nuc_M = data.frame(data.frame(annotation = Peak_F_Nuc_M_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_Nuc_M_mRNA[, Nuc_F_M])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_Cyt_M = data.frame(data.frame(annotation = Peak_F_Cyt_M_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_Cyt_M_mRNA[, Cyto_F_M])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_NLS_M = data.frame(data.frame(annotation = Peak_Co_NLS_M_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NLS_M_mRNA[, c(NLS_E_M)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_NES_M = data.frame(data.frame(annotation = Peak_Co_NES_M_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NES_M_mRNA[, c(NES_E_M)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_Nuc_S = data.frame(data.frame(annotation = Peak_F_Nuc_S_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_Nuc_S_mRNA[, Nuc_F_S])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_Cyt_S = data.frame(data.frame(annotation = Peak_F_Cyt_S_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_Cyt_S_mRNA[, Cyto_F_S])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_NLS_S = data.frame(data.frame(annotation = Peak_Co_NLS_S_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NLS_S_mRNA[, c(NLS_E_S)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_Co_NES_S = data.frame(data.frame(annotation = Peak_Co_NES_S_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NES_S_mRNA[, c(NES_E_S)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')

colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')
colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')


## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Cyt_M, 
                                  PeakDistribution_Co_NLS_M, PeakDistribution_Co_NES_M,
                                  PeakDistribution_F_Nuc_S, PeakDistribution_F_Cyt_S,
                                  PeakDistribution_Co_NLS_S, PeakDistribution_Co_NES_S)

PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES') %>%
  select(Source, Freq, Annotation)

PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = c( 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'))
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_M_Nuc' | Source == 'F_M_Cyt' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP')) +
  ggtitle('Mock mRNA Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_S_Nuc' | Source == 'F_S_Cyt' | Source == 'Co_S_NLS' | Source == 'Co_S_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP')) +
  ggtitle('Arsenite mRNA Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Fraction Distribution Stacked Bar Graph
PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), 
                                  PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M),
                                  PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S),
                                  PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S))
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES') %>%
  select(Source, Freq, Annotation)

PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = c( 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'))
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", "DS10K"))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_M_Nuc' | Source == 'F_M_Cyt' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP')) +
  ggtitle('Mock mRNA Peak Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_S_Nuc' | Source == 'F_S_Cyt' | Source == 'Co_S_NLS' | Source == 'Co_S_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP')) +
  ggtitle('Arsenite mRNA Peak Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## ncRNA
Peak_F_Nuc_M_Not_mRNA = Peak_F_Nuc_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_F_Cyt_M_Not_mRNA = Peak_F_Cyt_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_NLS_M_Not_mRNA = Peak_Co_NLS_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_NES_M_Not_mRNA = Peak_Co_NES_M %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

Peak_F_Nuc_S_Not_mRNA = Peak_F_Nuc_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_F_Cyt_S_Not_mRNA = Peak_F_Cyt_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_NLS_S_Not_mRNA = Peak_Co_NLS_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')
Peak_Co_NES_S_Not_mRNA = Peak_Co_NES_S %>% filter(finalized_annotation != "3'UTR" & finalized_annotation != "CDS" & finalized_annotation != "5'UTR" & finalized_annotation != "intron" & finalized_annotation != "CDS_RI" & finalized_annotation != 'DS10K')

PeakDistribution_F_Nuc_M = data.frame(data.frame(annotation = Peak_F_Nuc_M_Not_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_Nuc_M_Not_mRNA[, Nuc_F_M])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_F_Nuc_M = rbind(PeakDistribution_F_Nuc_M, c(0))
row.names(PeakDistribution_F_Nuc_M) = c(row.names(PeakDistribution_F_Nuc_M)[1:10], 'scaRNA')

PeakDistribution_F_Cyt_M = data.frame(data.frame(annotation = Peak_F_Cyt_M_Not_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_Cyt_M_Not_mRNA[, Cyto_F_M])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_F_Cyt_M = rbind(PeakDistribution_F_Cyt_M, c(0))
row.names(PeakDistribution_F_Cyt_M) = c(row.names(PeakDistribution_F_Cyt_M)[1:10], 'scaRNA')
PeakDistribution_F_Cyt_M = PeakDistribution_F_Cyt_M[row.names(PeakDistribution_F_Nuc_M), 'tags', drop = FALSE]

PeakDistribution_Co_NLS_M = data.frame(data.frame(annotation = Peak_Co_NLS_M_Not_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NLS_M_Not_mRNA[, c(NLS_E_M)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_NLS_M = rbind(PeakDistribution_Co_NLS_M, c(0))
row.names(PeakDistribution_Co_NLS_M) = c(row.names(PeakDistribution_Co_NLS_M)[1:10], 'scaRNA')
PeakDistribution_Co_NLS_M = PeakDistribution_Co_NLS_M[row.names(PeakDistribution_F_Nuc_M), 'tags', drop = FALSE]

PeakDistribution_Co_NES_M = data.frame(data.frame(annotation = Peak_Co_NES_M_Not_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NES_M_Not_mRNA[, c(NES_E_M)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_NES_M = rbind(PeakDistribution_Co_NES_M, c(0), c(0), c(0))
row.names(PeakDistribution_Co_NES_M) = c(row.names(PeakDistribution_Co_NES_M)[1:8], 'miRNA', 'nC_RI', 'scaRNA')
PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M[row.names(PeakDistribution_F_Nuc_M), 'tags', drop = FALSE]

PeakDistribution_F_Nuc_S = data.frame(data.frame(annotation = Peak_F_Nuc_S_Not_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_Nuc_S_Not_mRNA[, Nuc_F_S])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_F_Nuc_S = rbind(PeakDistribution_F_Nuc_S, c(0))
row.names(PeakDistribution_F_Nuc_S) = c(row.names(PeakDistribution_F_Nuc_S)[1:10], 'scaRNA')
PeakDistribution_F_Nuc_S = PeakDistribution_F_Nuc_S[row.names(PeakDistribution_F_Nuc_M), 'tags', drop = FALSE]

PeakDistribution_F_Cyt_S = data.frame(data.frame(annotation = Peak_F_Cyt_S_Not_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_Cyt_S_Not_mRNA[, Cyto_F_S])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_F_Cyt_S = rbind(PeakDistribution_F_Cyt_S, c(0), c(0))
row.names(PeakDistribution_F_Cyt_S) = c(row.names(PeakDistribution_F_Cyt_S)[1:9], 'miRNA', 'scaRNA')
PeakDistribution_F_Cyt_S = PeakDistribution_F_Cyt_S[row.names(PeakDistribution_F_Nuc_M), 'tags', drop = FALSE]

PeakDistribution_Co_NLS_S = data.frame(data.frame(annotation = Peak_Co_NLS_S_Not_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NLS_S_Not_mRNA[, c(NLS_E_S)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_NLS_S = rbind(PeakDistribution_Co_NLS_S, c(0))
row.names(PeakDistribution_Co_NLS_S) = c(row.names(PeakDistribution_Co_NLS_S)[1:10], 'scaRNA')
PeakDistribution_Co_NLS_S = PeakDistribution_Co_NLS_S[row.names(PeakDistribution_F_Nuc_M), 'tags', drop = FALSE]

PeakDistribution_Co_NES_S = data.frame(data.frame(annotation = Peak_Co_NES_S_Not_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NES_S_Not_mRNA[, c(NES_E_S)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_NES_S = rbind(PeakDistribution_Co_NES_S, c(0), c(0))
row.names(PeakDistribution_Co_NES_S) = c(row.names(PeakDistribution_Co_NES_S)[1:9], 'scaRNA', 'nC_RI')
PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S[row.names(PeakDistribution_F_Nuc_M), 'tags', drop = FALSE]

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_Co_NLS_M) = c('Co_M_NLS')
colnames(PeakDistribution_Co_NES_M) = c('Co_M_NES')

colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')
colnames(PeakDistribution_Co_NLS_S) = c('Co_S_NLS')
colnames(PeakDistribution_Co_NES_S) = c('Co_S_NES')


## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Cyt_M, 
                                  PeakDistribution_Co_NLS_M, PeakDistribution_Co_NES_M,
                                  PeakDistribution_F_Nuc_S, PeakDistribution_F_Cyt_S,
                                  PeakDistribution_Co_NLS_S, PeakDistribution_Co_NES_S)

PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES') %>%
  select(Source, Freq, Annotation)

PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = c( 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'))
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_M_Nuc' | Source == 'F_M_Cyt' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP')) +
  ggtitle('Mock ncRNA Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_S_Nuc' | Source == 'F_S_Cyt' | Source == 'Co_S_NLS' | Source == 'Co_S_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP')) +
  ggtitle('Arsenite ncRNA Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Fraction Distribution Stacked Bar Graph
PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), 
                                  PeakDistribution_Co_NLS_M/sum(PeakDistribution_Co_NLS_M), PeakDistribution_Co_NES_M/sum(PeakDistribution_Co_NES_M),
                                  PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S),
                                  PeakDistribution_Co_NLS_S/sum(PeakDistribution_Co_NLS_S), PeakDistribution_Co_NES_S/sum(PeakDistribution_Co_NES_S))
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES') %>%
  select(Source, Freq, Annotation)

PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = c( 'F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'))
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI', 'UnAn'))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_M_Nuc' | Source == 'F_M_Cyt' | Source == 'Co_M_NLS' | Source == 'Co_M_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP')) +
  ggtitle('Mock ncRNA Peak Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(PeakDistribution_combined %>% filter(Source == 'F_S_Nuc' | Source == 'F_S_Cyt' | Source == 'Co_S_NLS' | Source == 'Co_S_NES'), aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP')) +
  ggtitle('Arsenite ncRNA Peak Counts Distributions by Fraction') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

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

## FIGURE2 Scatter Plot of Nuclear Fraction vs CoCLIP:
####################################################################################################################
# Nuclear Mock:
Nuclear_Mock_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation != 'UnAn') & 
                                                      ((Nuc_F_M_BC >= BC_Threshold_F &
                                                         Nuc_F_M > median(NLS_E_M) * rowSum_Multiplier_F) & 
                                                      (NLS_E_M_BC >= BC_Threshold_E &
                                                         NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_F)))

Nuclear_Mock_Peaks_Filtered$grouped_annotation = factor(Nuclear_Mock_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(Nuclear_Mock_Peaks_Filtered, aes(x = log10(NLS_E_M), y = log10(Nuc_F_M), color = grouped_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log10(NLS CoCLIP)', y = 'log10(Nuclear Fraction)') +
  xlim(c(0, 4)) +
  ylim(c(0, 4)) +
  ggtitle(paste0('Mock HuR Normalized Tag Counts Per Peaks: NLS CoCLIP vs Nuclear FracCLIP (',  nrow(Nuclear_Mock_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

# Nuclear Mock - mRNA specific:
mRNA_Nuclear_Mock_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation == "5'UTR" | 
                                                            finalized_annotation == "3'UTR" | 
                                                            finalized_annotation == "CDS" | 
                                                            finalized_annotation == "intron" | 
                                                            finalized_annotation == "CDS_RI" |
                                                            finalized_annotation == "DS10K") &
                                                           (grouped_annotation != 'UnAn') & 
                                                           ((Nuc_F_M_BC >= BC_Threshold_F &
                                                               Nuc_F_M > median(NLS_E_M) * rowSum_Multiplier_F) & 
                                                              (NLS_E_M_BC >= BC_Threshold_E &
                                                                 NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_F)))

mRNA_Nuclear_Mock_Peaks_Filtered$finalized_annotation = factor(mRNA_Nuclear_Mock_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(mRNA_Nuclear_Mock_Peaks_Filtered, aes(x = log10(NLS_E_M), y = log10(Nuc_F_M), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log10(NLS CoCLIP)', y = 'log10(Nuclear Fraction)') + 
  xlim(c(0, 4)) +
  ylim(c(0, 4)) +
  ggtitle(paste0('Mock mRNA HuR Normalized Tag Counts Per Peaks: NLS CoCLIP vs Nuclear FracCLIP (',  nrow(mRNA_Nuclear_Mock_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

# Nuclear Mock - non-mRNA specific:
ncRNA_Nuclear_Mock_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation != "5'UTR" & 
                                                             finalized_annotation != "3'UTR" & 
                                                             finalized_annotation != "CDS" & 
                                                             finalized_annotation != "intron" &  
                                                             finalized_annotation != "CDS_RI" &
                                                             finalized_annotation != "DS10K" &
                                                             finalized_annotation != "UnAn") & 
                                                            ((Nuc_F_M_BC >= BC_Threshold_F &
                                                                Nuc_F_M > median(NLS_E_M) * rowSum_Multiplier_F) & 
                                                               (NLS_E_M_BC >= BC_Threshold_E &
                                                                  NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_F)))

ncRNA_Nuclear_Mock_Peaks_Filtered$finalized_annotation = factor(ncRNA_Nuclear_Mock_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(ncRNA_Nuclear_Mock_Peaks_Filtered, aes(x = log10(NLS_E_M), y = log10(Nuc_F_M), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log10(NLS CoCLIP)', y = 'log10(Nuclear Fraction)') +
  xlim(c(0, 4)) +
  ylim(c(0, 4)) +
  ggtitle(paste0('Mock ncRNA HuR Normalized Tag Counts Per Peaks: NLS CoCLIP vs Nuclear FracCLIP (',  nrow(ncRNA_Nuclear_Mock_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

# Nuclear Stress:
Nuclear_Stress_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation != 'UnAn') & 
                                                        ((Nuc_F_S_BC >= BC_Threshold_F &
                                                            Nuc_F_S > median(NLS_E_S) * rowSum_Multiplier_F) & 
                                                           (NLS_E_S_BC >= BC_Threshold_E &
                                                              NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_F)))

Nuclear_Stress_Peaks_Filtered$grouped_annotation = factor(Nuclear_Stress_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(Nuclear_Stress_Peaks_Filtered, aes(x = log10(NLS_E_S), y = log10(Nuc_F_S), color = grouped_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log10(NLS CoCLIP)', y = 'log10(Nuclear Fraction)') +
  xlim(c(0, 4)) +
  ylim(c(0, 4)) +
  ggtitle(paste0('Stress HuR Normalized Tag Counts Per Peaks: NLS CoCLIP vs Nuclear FracCLIP (',  nrow(Nuclear_Stress_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

# Nuclear Stress - mRNA specific:
mRNA_Nuclear_Stress_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation == "5'UTR" | 
                                                              finalized_annotation == "3'UTR" | 
                                                              finalized_annotation == "CDS" | 
                                                              finalized_annotation == "intron" | 
                                                              finalized_annotation == "CDS_RI" |
                                                              finalized_annotation == "DS10K") &
                                                             (grouped_annotation != 'UnAn') & 
                                                             ((Nuc_F_S_BC >= BC_Threshold_F &
                                                                 Nuc_F_S > median(NLS_E_S) * rowSum_Multiplier_F) & 
                                                                (NLS_E_S_BC >= BC_Threshold_E &
                                                                   NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_F)))

mRNA_Nuclear_Stress_Peaks_Filtered$finalized_annotation = factor(mRNA_Nuclear_Stress_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(mRNA_Nuclear_Stress_Peaks_Filtered, aes(x = log10(NLS_E_S), y = log10(Nuc_F_S), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log10(NLS CoCLIP)', y = 'log10(Nuclear Fraction)') +
  xlim(c(0, 4)) +
  ylim(c(0, 4)) +
  ggtitle(paste0('Stress mRNA HuR Normalized Tag Counts Per Peaks: NLS CoCLIP vs Nuclear FracCLIP (',  nrow(mRNA_Nuclear_Stress_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

# Nuclear Stress - non-mRNA specific:
ncRNA_Nuclear_Stress_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation != "5'UTR" & 
                                                               finalized_annotation != "3'UTR" & 
                                                               finalized_annotation != "CDS" & 
                                                               finalized_annotation != "intron" &  
                                                               finalized_annotation != "CDS_RI" &
                                                               finalized_annotation != "DS10K" &
                                                               finalized_annotation != "UnAn") & 
                                                              ((Nuc_F_S_BC >= BC_Threshold_F &
                                                                  Nuc_F_S > median(NLS_E_S) * rowSum_Multiplier_F) & 
                                                                 (NLS_E_S_BC >= BC_Threshold_E &
                                                                    NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_F)))

ncRNA_Nuclear_Stress_Peaks_Filtered$finalized_annotation = factor(ncRNA_Nuclear_Stress_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(ncRNA_Nuclear_Stress_Peaks_Filtered, aes(x = log10(NLS_E_S), y = log10(Nuc_F_S), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log10(NLS CoCLIP)', y = 'log10(Nuclear Fraction)') +
  xlim(c(0, 4)) +
  ylim(c(0, 4)) +
  ggtitle(paste0('Stress ncRNA HuR Normalized Tag Counts Per Peaks: NLS CoCLIP vs Nuclear FracCLIP (',  nrow(ncRNA_Nuclear_Stress_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

####################################################################################################################

## FIGURE2 Scatter Plot of Cytoplasm Fraction vs CoCLIP:
####################################################################################################################
# Cytoplasm Mock:
Cytoplasm_Mock_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation != 'UnAn') & 
                                                        ((Cyto_F_M_BC >= BC_Threshold_F &
                                                            Cyto_F_M > median(NES_E_M) * rowSum_Multiplier_F) & 
                                                           (NES_E_M_BC >= BC_Threshold_E &
                                                              NES_E_M > median(NES_E_M) * rowSum_Multiplier_F)))

# Cytoplasm_Mock_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation != 'UnAn') & 
#                                                         ((Cyto_F_M_BC >= BC_Threshold_F) &
#                                                            Cyto_F_M > median(Cyto_F_M) * rowSum_Multiplier_F & 
#                                                            (NES_E_M_BC >= 1) &
#                                                            NES_E_M > median(NES_E_M) * rowSum_Multiplier_F) & 
#                                                         ((Cyto_F_S_BC >= BC_Threshold_F) &
#                                                            Cyto_F_S > median(Cyto_F_S) * rowSum_Multiplier_F & 
#                                                            (NES_E_S_BC >= 1) &
#                                                            NES_E_S > median(NES_E_S) * rowSum_Multiplier_F))

Cytoplasm_Mock_Peaks_Filtered$grouped_annotation = factor(Cytoplasm_Mock_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(Cytoplasm_Mock_Peaks_Filtered, aes(x = log10(NES_E_M), y = log10(Cyto_F_M), color = grouped_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log10(NES Enrich)', y = 'log10(Cytoplasm Fraction)') + 
  xlim(c(0, 4)) +
  ylim(c(0, 4)) +
  ggtitle(paste0('Mock HuR Normalized Tag Counts Per Peaks: NES CoCLIP vs Cytoplasm FracCLIP (',  nrow(Cytoplasm_Mock_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

# Cytoplasm Mock - mRNA specific:
mRNA_Cytoplasm_Mock_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation == "5'UTR" | 
                                                              finalized_annotation == "3'UTR" | 
                                                              finalized_annotation == "CDS" | 
                                                              finalized_annotation == "intron" | 
                                                              finalized_annotation == "CDS_RI" |
                                                              finalized_annotation == "DS10K" &
                                                             grouped_annotation != 'UnAn') & 
                                                             ((Cyto_F_M_BC >= BC_Threshold_F &
                                                                 Cyto_F_M > median(NES_E_M) * rowSum_Multiplier_F) & 
                                                                (NES_E_M_BC >= BC_Threshold_E &
                                                                   NES_E_M > median(NES_E_M) * rowSum_Multiplier_F)))

mRNA_Cytoplasm_Mock_Peaks_Filtered$finalized_annotation = factor(mRNA_Cytoplasm_Mock_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(mRNA_Cytoplasm_Mock_Peaks_Filtered, aes(x = log10(NES_E_M), y = log10(Cyto_F_M), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log10(NES Enrich)', y = 'log10(Cytoplasm Fraction)') + 
  xlim(c(0, 4)) +
  ylim(c(0, 4)) +
  ggtitle(paste0('Mock mRNA HuR Normalized Tag Counts Per Peaks: NES CoCLIP vs Cytoplasm FracCLIP (',  nrow(mRNA_Cytoplasm_Mock_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

# Cytoplasm Mock - non-mRNA specific:
ncRNA_Cytoplasm_Mock_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation != "5'UTR" & 
                                                               finalized_annotation != "3'UTR" & 
                                                               finalized_annotation != "CDS" & 
                                                               finalized_annotation != "intron" &  
                                                               finalized_annotation != "CDS_RI" &
                                                               finalized_annotation != "DS10K" &
                                                               finalized_annotation != "UnAn") & 
                                                              ((Cyto_F_M_BC >= BC_Threshold_F &
                                                                  Cyto_F_M > median(NES_E_M) * rowSum_Multiplier_F) & 
                                                                 (NES_E_M_BC >= BC_Threshold_E &
                                                                    NES_E_M > median(NES_E_M) * rowSum_Multiplier_F)))

ncRNA_Cytoplasm_Mock_Peaks_Filtered$finalized_annotation = factor(ncRNA_Cytoplasm_Mock_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(ncRNA_Cytoplasm_Mock_Peaks_Filtered, aes(x = log10(NES_E_M), y = log10(Cyto_F_M), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log10(NES Enrich)', y = 'log10(Cytoplasm Fraction)') + 
  xlim(c(0, 4)) +
  ylim(c(0, 4)) +
  ggtitle(paste0('Mock ncRNA HuR Normalized Tag Counts Per Peaks: NES CoCLIP vs Cytoplasm FracCLIP (',  nrow(ncRNA_Cytoplasm_Mock_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

# Cytoplasm Stress:
Cytoplasm_Stress_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation != 'UnAn') & 
                                                          ((Cyto_F_S_BC >= BC_Threshold_F &
                                                              Cyto_F_S > median(NES_E_S) * rowSum_Multiplier_F) & 
                                                             (NES_E_S_BC >= BC_Threshold_E &
                                                                NES_E_S > median(NES_E_S) * rowSum_Multiplier_F)))

Cytoplasm_Stress_Peaks_Filtered$grouped_annotation = factor(Cytoplasm_Stress_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(Cytoplasm_Stress_Peaks_Filtered, aes(x = log10(NES_E_S), y = log10(Cyto_F_S), color = grouped_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log10(NES Enrich)', y = 'log10(Cytoplasm Fraction)') + 
  xlim(c(0, 4)) +
  ylim(c(0, 4)) +
  ggtitle(paste0('Stress HuR Normalized Tag Counts Per Peaks: NES CoCLIP vs Cytoplasm FracCLIP (',  nrow(Cytoplasm_Stress_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

# Cytoplasm Stress - mRNA specific:
mRNA_Cytoplasm_Stress_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation == "5'UTR" | 
                                                                finalized_annotation == "3'UTR" | 
                                                                finalized_annotation == "CDS" | 
                                                                finalized_annotation == "intron" | 
                                                                finalized_annotation == "CDS_RI" |
                                                                finalized_annotation == "DS10K" &
                                                               grouped_annotation != 'UnAn') & 
                                                               ((Cyto_F_S_BC >= BC_Threshold_F &
                                                                   Cyto_F_S > median(NES_E_S) * rowSum_Multiplier_F) & 
                                                                  (NES_E_S_BC >= BC_Threshold_E &
                                                                     NES_E_S > median(NES_E_S) * rowSum_Multiplier_F)))

mRNA_Cytoplasm_Stress_Peaks_Filtered$finalized_annotation = factor(mRNA_Cytoplasm_Stress_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(mRNA_Cytoplasm_Stress_Peaks_Filtered, aes(x = log10(NES_E_S), y = log10(Cyto_F_S), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log10(NES Enrich)', y = 'log10(Cytoplasm Fraction)') + 
  xlim(c(0, 4)) +
  ylim(c(0, 4)) +
  ggtitle(paste0('Stress mRNA HuR Normalized Tag Counts Per Peaks: NES CoCLIP vs Cytoplasm FracCLIP (',  nrow(mRNA_Cytoplasm_Stress_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

# Cytoplasm Stress - non-mRNA specific:
ncRNA_Cytoplasm_Stress_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation != "5'UTR" & 
                                                                 finalized_annotation != "3'UTR" & 
                                                                 finalized_annotation != "CDS" & 
                                                                 finalized_annotation != "intron" &  
                                                                 finalized_annotation != "CDS_RI" &
                                                                 finalized_annotation != "DS10K" &
                                                                 finalized_annotation != "UnAn") & 
                                                                ((Cyto_F_S_BC >= BC_Threshold_F &
                                                                    Cyto_F_S > median(NES_E_S) * rowSum_Multiplier_F) & 
                                                                   (NES_E_S_BC >= BC_Threshold_E &
                                                                      NES_E_S > median(NES_E_S) * rowSum_Multiplier_F)))

ncRNA_Cytoplasm_Stress_Peaks_Filtered$finalized_annotation = factor(ncRNA_Cytoplasm_Stress_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(ncRNA_Cytoplasm_Stress_Peaks_Filtered, aes(x = log10(NES_E_S), y = log10(Cyto_F_S), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log10(NES Enrich)', y = 'log10(Cytoplasm Fraction)') + 
  xlim(c(0, 4)) +
  ylim(c(0, 4)) +
  ggtitle(paste0('Stress ncRNA HuR Normalized Tag Counts Per Peaks: NES CoCLIP vs Cytoplasm FracCLIP (',  nrow(ncRNA_Cytoplasm_Stress_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

####################################################################################################################

## FIGURE2 Scatter Plot of Nuclear Fold Change Over Input Fraction vs CoCLIP:
####################################################################################################################
Mock_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != "UnAn") & 
                                                  ((Nuc_F_M_BC >= BC_Threshold_F &
                                                      Nuc_F_M > median(NLS_E_M) * rowSum_Multiplier_F & 
                                                      NLS_E_M_BC >= BC_Threshold_E &
                                                      NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_F & 
                                                      (NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I*3) | 
                                                     (Nuc_F_S_BC >= BC_Threshold_F &
                                                        Nuc_F_S > median(NLS_E_S) * rowSum_Multiplier_F &  
                                                        NLS_E_S_BC >= BC_Threshold_E &
                                                        NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_F & 
                                                        (NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I*3)))

data = data.frame(F_Nuc = Mock_Peaks_Filtered$Nuc_EvI_M, E_Nuc = Mock_Peaks_Filtered$NLS_EvI_M)
data$annotation = factor(Mock_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(E_Nuc), y = log2(F_Nuc), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Input)', y = 'log2(FracCLIP Nuclear/Input)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Mock HuR Peak Enrichment over Input: Nuclear CoCLIP vs FracCLIP (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

Mock_mRNA_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                        finalized_annotation == "3'UTR" | 
                                                        finalized_annotation == "CDS" | 
                                                        finalized_annotation == "intron" | 
                                                        finalized_annotation == "CDS_RI" |
                                                        finalized_annotation == "DS10K" &
                                                       grouped_annotation != 'UnAn') & 
                                                       ((Nuc_F_M_BC >= BC_Threshold_F &
                                                           Nuc_F_M > median(NLS_E_M) * rowSum_Multiplier_F & 
                                                           NLS_E_M_BC >= BC_Threshold_E &
                                                           NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_F & 
                                                           (NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I*3) | 
                                                          (Nuc_F_S_BC >= BC_Threshold_F &
                                                             Nuc_F_S > median(NLS_E_S) * rowSum_Multiplier_F &  
                                                             NLS_E_S_BC >= BC_Threshold_E &
                                                             NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_F & 
                                                             (NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I*3)))

data = data.frame(F_Nuc = Mock_mRNA_Peaks_Filtered$Nuc_EvI_M, E_Nuc = Mock_mRNA_Peaks_Filtered$NLS_EvI_M)
data$annotation = factor(Mock_mRNA_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(E_Nuc), y = log2(F_Nuc), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Input)', y = 'log2(FracCLIP Nuclear/Input)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Mock mRNA HuR Peak Enrichment over Input: Nuclear CoCLIP vs FracCLIP (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

Mock_ncRNA_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                         finalized_annotation != "3'UTR" & 
                                                         finalized_annotation != "CDS" & 
                                                         finalized_annotation != "intron" &  
                                                         finalized_annotation != "CDS_RI" &
                                                         finalized_annotation != "DS10K" &
                                                         finalized_annotation != "UnAn") & 
                                                        ((Nuc_F_M_BC >= BC_Threshold_F &
                                                            Nuc_F_M > median(NLS_E_M) * rowSum_Multiplier_F & 
                                                            NLS_E_M_BC >= BC_Threshold_E &
                                                            NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_F & 
                                                            (NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I*3) | 
                                                           (Nuc_F_S_BC >= BC_Threshold_F &
                                                              Nuc_F_S > median(NLS_E_S) * rowSum_Multiplier_F &  
                                                              NLS_E_S_BC >= BC_Threshold_E &
                                                              NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_F & 
                                                              (NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I*3)))

data = data.frame(F_Nuc = Mock_ncRNA_Peaks_Filtered$Nuc_EvI_M, E_Nuc = Mock_ncRNA_Peaks_Filtered$NLS_EvI_M)
data$annotation = factor(Mock_ncRNA_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(E_Nuc), y = log2(F_Nuc), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Input)', y = 'log2(FracCLIP Nuclear/Input)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Mock ncRNA HuR Peak Enrichment over Input: Nuclear CoCLIP vs FracCLIP (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

Arsenite_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != 'UnAn') & 
                                                      ((Nuc_F_M_BC >= BC_Threshold_F &
                                                          Nuc_F_M > median(NLS_E_M) * rowSum_Multiplier_F & 
                                                          NLS_E_M_BC >= BC_Threshold_E &
                                                          NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_F & 
                                                          (NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I*3) | 
                                                         (Nuc_F_S_BC >= BC_Threshold_F &
                                                            Nuc_F_S > median(NLS_E_S) * rowSum_Multiplier_F &  
                                                            NLS_E_S_BC >= BC_Threshold_E &
                                                            NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_F & 
                                                            (NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I*3)))

data = data.frame(F_Nuc = Arsenite_Peaks_Filtered$Nuc_EvI_S, E_Nuc = Arsenite_Peaks_Filtered$NLS_EvI_S)
data$annotation = factor(Arsenite_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(E_Nuc), y = log2(F_Nuc), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Input)', y = 'log2(FracCLIP Nuclear/Input)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Arsenite HuR Peak Enrichment over Input: Nuclear CoCLIP vs FracCLIP (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

Arsenite_mRNA_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                            finalized_annotation == "3'UTR" | 
                                                            finalized_annotation == "CDS" | 
                                                            finalized_annotation == "intron" | 
                                                            finalized_annotation == "CDS_RI" |
                                                            finalized_annotation == "DS10K" &
                                                           grouped_annotation != 'UnAn') & 
                                                           ((Nuc_F_M_BC >= BC_Threshold_F &
                                                               Nuc_F_M > median(NLS_E_M) * rowSum_Multiplier_F & 
                                                               NLS_E_M_BC >= BC_Threshold_E &
                                                               NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_F & 
                                                               (NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I*3) | 
                                                              (Nuc_F_S_BC >= BC_Threshold_F &
                                                                 Nuc_F_S > median(NLS_E_S) * rowSum_Multiplier_F &  
                                                                 NLS_E_S_BC >= BC_Threshold_E &
                                                                 NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_F & 
                                                                 (NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I*3)))

data = data.frame(F_Nuc = Arsenite_mRNA_Peaks_Filtered$Nuc_EvI_S, E_Nuc = Arsenite_mRNA_Peaks_Filtered$NLS_EvI_S)
data$annotation = factor(Arsenite_mRNA_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(E_Nuc), y = log2(F_Nuc), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Input)', y = 'log2(FracCLIP Nuclear/Input)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Arsenite mRNA HuR Peak Enrichment over Input: Nuclear CoCLIP vs FracCLIP (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

Arsenite_ncRNA_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                             finalized_annotation != "3'UTR" & 
                                                             finalized_annotation != "CDS" & 
                                                             finalized_annotation != "intron" &  
                                                             finalized_annotation != "CDS_RI" &
                                                             finalized_annotation != "DS10K" &
                                                             finalized_annotation != "UnAn") & 
                                                            ((Nuc_F_M_BC >= BC_Threshold_F &
                                                                Nuc_F_M > median(NLS_E_M) * rowSum_Multiplier_F & 
                                                                NLS_E_M_BC >= BC_Threshold_E &
                                                                NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_F & 
                                                                (NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I*3) | 
                                                               (Nuc_F_S_BC >= BC_Threshold_F &
                                                                  Nuc_F_S > median(NLS_E_S) * rowSum_Multiplier_F &  
                                                                  NLS_E_S_BC >= BC_Threshold_E &
                                                                  NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_F & 
                                                                  (NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I*3)))

data = data.frame(F_Nuc = Arsenite_ncRNA_Peaks_Filtered$Nuc_EvI_S, E_Nuc = Arsenite_ncRNA_Peaks_Filtered$NLS_EvI_S)
data$annotation = factor(Arsenite_ncRNA_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(E_Nuc), y = log2(F_Nuc), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Input)', y = 'log2(FracCLIP Nuclear/Input)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Arsenite ncRNA HuR Peak Enrichment over Input: Nuclear CoCLIP vs FracCLIP (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')


####################################################################################################################

## FIGURE2 Scatter Plot of Cytoplasm Fold Change Over Input Fraction vs CoCLIP:
####################################################################################################################
Mock_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != 'UnAn') & 
                                                  ((Cyto_F_M_BC >= BC_Threshold_F &
                                                      Cyto_F_M > median(NES_E_M) * rowSum_Multiplier_F & 
                                                     NES_E_M_BC >= BC_Threshold_E &
                                                        NES_E_M > median(NES_E_M) * rowSum_Multiplier_F & 
                                                     (NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I*3) | 
                                                  (Cyto_F_S_BC >= BC_Threshold_F &
                                                      Cyto_F_S > median(NES_E_S) * rowSum_Multiplier_F & 
                                                     NES_E_S_BC >= BC_Threshold_E &
                                                        NES_E_S > median(NES_E_S) * rowSum_Multiplier_F & 
                                                     (NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I*3)))

data = data.frame(F_Nuc = Mock_Peaks_Filtered$Cyto_EvI_M, E_Nuc = Mock_Peaks_Filtered$NES_EvI_M)
data$annotation = factor(Mock_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(E_Nuc), y = log2(F_Nuc), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Cytoplasm/Input)', y = 'log2(FracCLIP Cytoplasm/Input)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Mock HuR Peak Enrichment over Input: Cytoplasm CoCLIP vs FracCLIP (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

Mock_mRNA_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                        finalized_annotation == "3'UTR" | 
                                                        finalized_annotation == "CDS" | 
                                                        finalized_annotation == "intron" | 
                                                        finalized_annotation == "CDS_RI" |
                                                        finalized_annotation == "DS10K" &
                                                       grouped_annotation != 'UnAn') & 
                                                       ((Cyto_F_M_BC >= BC_Threshold_F &
                                                           Cyto_F_M > median(NES_E_M) * rowSum_Multiplier_F & 
                                                           NES_E_M_BC >= BC_Threshold_E &
                                                           NES_E_M > median(NES_E_M) * rowSum_Multiplier_F & 
                                                           (NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I*3) | 
                                                          (Cyto_F_S_BC >= BC_Threshold_F &
                                                             Cyto_F_S > median(NES_E_S) * rowSum_Multiplier_F & 
                                                             NES_E_S_BC >= BC_Threshold_E &
                                                             NES_E_S > median(NES_E_S) * rowSum_Multiplier_F & 
                                                             (NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I*3)))

data = data.frame(F_Nuc = Mock_mRNA_Peaks_Filtered$Cyto_EvI_M, E_Nuc = Mock_mRNA_Peaks_Filtered$NES_EvI_M)
data$annotation = factor(Mock_mRNA_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(E_Nuc), y = log2(F_Nuc), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Cytoplasm/Input)', y = 'log2(FracCLIP Cytoplasm/Input)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Mock mRNA HuR Peak Enrichment over Input: Cytoplasm CoCLIP vs FracCLIP (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

Mock_ncRNA_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                         finalized_annotation != "3'UTR" & 
                                                         finalized_annotation != "CDS" & 
                                                         finalized_annotation != "intron" &  
                                                         finalized_annotation != "CDS_RI" &
                                                         finalized_annotation != "DS10K" &
                                                         finalized_annotation != "UnAn") & 
                                                        ((Cyto_F_M_BC >= BC_Threshold_F &
                                                            Cyto_F_M > median(NES_E_M) * rowSum_Multiplier_F & 
                                                            NES_E_M_BC >= BC_Threshold_E &
                                                            NES_E_M > median(NES_E_M) * rowSum_Multiplier_F & 
                                                            (NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I*3) | 
                                                           (Cyto_F_S_BC >= BC_Threshold_F &
                                                              Cyto_F_S > median(NES_E_S) * rowSum_Multiplier_F & 
                                                              NES_E_S_BC >= BC_Threshold_E &
                                                              NES_E_S > median(NES_E_S) * rowSum_Multiplier_F & 
                                                              (NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I*3)))

data = data.frame(F_Nuc = Mock_ncRNA_Peaks_Filtered$Cyto_EvI_M, E_Nuc = Mock_ncRNA_Peaks_Filtered$NES_EvI_M)
data$annotation = factor(Mock_ncRNA_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(E_Nuc), y = log2(F_Nuc), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Cytoplasm/Input)', y = 'log2(FracCLIP Cytoplasm/Input)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Mock ncRNA HuR Peak Enrichment over Input: Cytoplasm CoCLIP vs FracCLIP (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')





Arsenite_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != 'UnAn') & 
                                                      ((Cyto_F_M_BC >= BC_Threshold_F &
                                                          Cyto_F_M > median(NES_E_M) * rowSum_Multiplier_F & 
                                                          NES_E_M_BC >= BC_Threshold_E &
                                                          NES_E_M > median(NES_E_M) * rowSum_Multiplier_F & 
                                                          (NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I*3) | 
                                                         (Cyto_F_S_BC >= BC_Threshold_F &
                                                            Cyto_F_S > median(NES_E_S) * rowSum_Multiplier_F & 
                                                            NES_E_S_BC >= BC_Threshold_E &
                                                            NES_E_S > median(NES_E_S) * rowSum_Multiplier_F & 
                                                            (NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I*3)))

data = data.frame(F_Nuc = Arsenite_Peaks_Filtered$Cyto_EvI_S, E_Nuc = Arsenite_Peaks_Filtered$NES_EvI_S)
data$annotation = factor(Arsenite_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(E_Nuc), y = log2(F_Nuc), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Cytoplasm/Input)', y = 'log2(FracCLIP Cytoplasm/Input)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Arsenite HuR Peak Enrichment over Input: Cytoplasm CoCLIP vs FracCLIP (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

Arsenite_mRNA_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                            finalized_annotation == "3'UTR" | 
                                                            finalized_annotation == "CDS" | 
                                                            finalized_annotation == "intron" | 
                                                            finalized_annotation == "CDS_RI" |
                                                            finalized_annotation == "DS10K" &
                                                           grouped_annotation != 'UnAn') & 
                                                           ((Cyto_F_M_BC >= BC_Threshold_F &
                                                               Cyto_F_M > median(NES_E_M) * rowSum_Multiplier_F & 
                                                               NES_E_M_BC >= BC_Threshold_E &
                                                               NES_E_M > median(NES_E_M) * rowSum_Multiplier_F & 
                                                               (NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I*3) | 
                                                              (Cyto_F_S_BC >= BC_Threshold_F &
                                                                 Cyto_F_S > median(NES_E_S) * rowSum_Multiplier_F & 
                                                                 NES_E_S_BC >= BC_Threshold_E &
                                                                 NES_E_S > median(NES_E_S) * rowSum_Multiplier_F & 
                                                                 (NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I*3)))

data = data.frame(F_Nuc = Arsenite_mRNA_Peaks_Filtered$Cyto_EvI_S, E_Nuc = Arsenite_mRNA_Peaks_Filtered$NES_EvI_S)
data$annotation = factor(Arsenite_mRNA_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(E_Nuc), y = log2(F_Nuc), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Cytoplasm/Input)', y = 'log2(FracCLIP Cytoplasm/Input)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Arsenite mRNA HuR Peak Enrichment over Input: Cytoplasm CoCLIP vs FracCLIP (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')

Arsenite_ncRNA_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                             finalized_annotation != "3'UTR" & 
                                                             finalized_annotation != "CDS" & 
                                                             finalized_annotation != "intron" &  
                                                             finalized_annotation != "CDS_RI" &
                                                             finalized_annotation != "DS10K" &
                                                             finalized_annotation != "UnAn") & 
                                                            ((Cyto_F_M_BC >= BC_Threshold_F &
                                                                Cyto_F_M > median(NES_E_M) * rowSum_Multiplier_F & 
                                                                NES_E_M_BC >= BC_Threshold_E &
                                                                NES_E_M > median(NES_E_M) * rowSum_Multiplier_F & 
                                                                (NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I*3) | 
                                                               (Cyto_F_S_BC >= BC_Threshold_F &
                                                                  Cyto_F_S > median(NES_E_S) * rowSum_Multiplier_F & 
                                                                  NES_E_S_BC >= BC_Threshold_E &
                                                                  NES_E_S > median(NES_E_S) * rowSum_Multiplier_F & 
                                                                  (NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I*3)))

data = data.frame(F_Nuc = Arsenite_ncRNA_Peaks_Filtered$Cyto_EvI_S, E_Nuc = Arsenite_ncRNA_Peaks_Filtered$NES_EvI_S)
data$annotation = factor(Arsenite_ncRNA_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(E_Nuc), y = log2(F_Nuc), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Cytoplasm/Input)', y = 'log2(FracCLIP Cytoplasm/Input)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Arsenite ncRNA HuR Peak Enrichment over Input: Cytoplasm CoCLIP vs FracCLIP (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_abline(linetype = 'dotted')


####################################################################################################################

## snoRNA Comparison
####################################################################################################################
snoRNA_Enrichment = peakRowSum %>% filter(finalized_annotation == 'snoRNA')
snoRNA_Enrichment = snoRNA_Enrichment %>% filter(Nuc_F_M_BC >= 1, 
                                                 Cyto_F_M_BC >= 1, 
                                                 NLS_E_M_BC >= 1, 
                                                 NES_E_M_BC >= 1)

snoRNA_E_M_NvC = snoRNA_Enrichment$NLS_E_M / snoRNA_Enrichment$NES_E_M
snoRNA_F_M_NvC = snoRNA_Enrichment$Nuc_F_M / snoRNA_Enrichment$Cyto_F_M
data = data.frame(log_snoRNA_E_M_NvC = log2(snoRNA_E_M_NvC), log_snoRNA_F_M_NvC = log2(snoRNA_F_M_NvC))

ggplot(data, aes(x = log_snoRNA_E_M_NvC, y = log_snoRNA_F_M_NvC)) +
  geom_point(col = 'black', pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP NLS/NES Mock)', y = 'log2(FractionCLIP Nuc/Cyto Mock)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle('snoRNA HuR Peaks: CoCLIP vs FracCLIP of Nuclear/Cytoplasm in Mock')


snoRNA_E_S_NvC = snoRNA_Enrichment$NLS_E_S / snoRNA_Enrichment$NES_E_S
snoRNA_F_S_NvC = snoRNA_Enrichment$Nuc_F_S / snoRNA_Enrichment$Cyto_F_S
data = data.frame(log_snoRNA_E_S_NvC = log2(snoRNA_E_S_NvC), log_snoRNA_F_S_NvC = log2(snoRNA_F_S_NvC))

ggplot(data, aes(x = log_snoRNA_E_S_NvC, y = log_snoRNA_F_S_NvC)) +
  geom_point(col = 'black', pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP NLS/NES Stress)', y = 'log2(FractionCLIP Nuc/Cyto Stress)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle('snoRNA HuR Peaks: CoCLIP vs FracCLIP of Nuclear/Cytoplasm in Stress')

####################################################################################################################

## Mock Vs Stress
####################################################################################################################
## Nuclear: 
NLS_Peaks_Filtered = peakEnrichment %>% filter(grouped_annotation != 'UnAn' & ((NLS_I_M_BC >= BC_Threshold_I & 
                                                  NLS_E_M_BC >= BC_Threshold_E) | 
                                                 (NLS_I_S_BC >= BC_Threshold_I & 
                                                 NLS_E_S_BC >= BC_Threshold_E)) &
                                                 I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                 E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_Peaks_Filtered$NLS_EvI_M, NLS_EvI_S = NLS_Peaks_Filtered$NLS_EvI_S)
data$annotation = factor(NLS_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(NLS_EvI_M), y = log2(NLS_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Mock Enriched/Input)', y = 'log2(Stress Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Nuclear HuR Peaks: Mock vs Stress of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Nuclear - mRNA specific:
NLS_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                  finalized_annotation == "3'UTR" | 
                                                  finalized_annotation == "CDS" | 
                                                  finalized_annotation == "intron" | 
                                                  finalized_annotation == "CDS_RI" |
                                                  finalized_annotation == "DS10K") & 
                                                 ((NLS_I_M_BC >= BC_Threshold_I & 
                                                   NLS_E_M_BC >= BC_Threshold_E) | 
                                                  (NLS_I_S_BC >= BC_Threshold_I & 
                                                     NLS_E_S_BC >= BC_Threshold_E)) &
                                                 I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                 E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_Peaks_Filtered$NLS_EvI_M, NLS_EvI_S = NLS_Peaks_Filtered$NLS_EvI_S)
data$annotation = factor(NLS_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(NLS_EvI_M), y = log2(NLS_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Mock Enriched/Input)', y = 'log2(Stress Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Nuclear HuR Peaks: Mock vs Stress of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Nuclear - non-mRNA specific:
NLS_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                  finalized_annotation != "3'UTR" & 
                                                  finalized_annotation != "CDS" & 
                                                  finalized_annotation != "intron" &  
                                                  finalized_annotation != "CDS_RI" &
                                                  finalized_annotation != "DS10K" &
                                                  finalized_annotation != "UnAn") & 
                                                 ((NLS_I_M_BC >= BC_Threshold_I & 
                                                     NLS_E_M_BC >= BC_Threshold_E) | 
                                                    (NLS_I_S_BC >= BC_Threshold_I & 
                                                       NLS_E_S_BC >= BC_Threshold_E)) &
                                                 I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                 E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_Peaks_Filtered$NLS_EvI_M, NLS_EvI_S = NLS_Peaks_Filtered$NLS_EvI_S)
data$annotation = factor(NLS_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(NLS_EvI_M), y = log2(NLS_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Mock Enriched/Input)', y = 'log2(Stress Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Nuclear HuR Peaks: Mock vs Stress of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Cytoplasm:
NES_Peaks_Filtered = peakEnrichment %>% filter(grouped_annotation != 'UnAn' & ((NES_I_M_BC >= BC_Threshold_I & 
                                                  NES_E_M_BC >= BC_Threshold_E) | 
                                                 (NES_I_S_BC >= BC_Threshold_I & 
                                                 NES_E_S_BC >= BC_Threshold_E)) &
                                                 I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                 E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_M = NES_Peaks_Filtered$NES_EvI_M, NES_EvI_S = NES_Peaks_Filtered$NES_EvI_S)
data$annotation = factor(NES_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(NES_EvI_M), y = log2(NES_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Mock Enriched/Input)', y = 'log2(Stress Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Cytoplasm HuR Peaks: Mock vs Stress of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Cytoplasm - mRNA specific:
NES_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                  finalized_annotation == "3'UTR" | 
                                                  finalized_annotation == "CDS" | 
                                                  finalized_annotation == "intron" | 
                                                  finalized_annotation == "CDS_RI" |
                                                  finalized_annotation == "DS10K") &
                                                 ((NES_I_M_BC >= BC_Threshold_I & 
                                                   NES_E_M_BC >= BC_Threshold_E) | 
                                                  (NES_I_S_BC >= BC_Threshold_I & 
                                                     NES_E_S_BC >= BC_Threshold_E)) &
                                                 I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                 E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_M = NES_Peaks_Filtered$NES_EvI_M, NES_EvI_S = NES_Peaks_Filtered$NES_EvI_S)
data$annotation = factor(NES_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(NES_EvI_M), y = log2(NES_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Mock Enriched/Input)', y = 'log2(Stress Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Cytoplasm HuR Peaks: Mock vs Stress of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Cytoplasm - non-mRNA specific:
NES_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                  finalized_annotation != "3'UTR" & 
                                                  finalized_annotation != "CDS" & 
                                                  finalized_annotation != "intron" &  
                                                  finalized_annotation != "CDS_RI" &
                                                  finalized_annotation != "DS10K" &
                                                  finalized_annotation != "UnAn") & 
                                                 ((NES_I_M_BC >= BC_Threshold_I & 
                                                     NES_E_M_BC >= BC_Threshold_E) | 
                                                    (NES_I_S_BC >= BC_Threshold_I & 
                                                       NES_E_S_BC >= BC_Threshold_E)) &
                                                 I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                 E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_M = NES_Peaks_Filtered$NES_EvI_M, NES_EvI_S = NES_Peaks_Filtered$NES_EvI_S)
data$annotation = factor(NES_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(NES_EvI_M), y = log2(NES_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Mock Enriched/Input)', y = 'log2(Stress Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Cytoplasm HuR Peaks: Mock vs Stress of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Stress Granule:
G3BP_Peaks_Filtered = peakEnrichment %>% filter(grouped_annotation != 'UnAn' & ((G3BP_I_M_BC >= BC_Threshold_I_SG & 
                                                   G3BP_E_M_BC >= BC_Threshold_E_SG) | 
                                                  (G3BP_I_S_BC >= BC_Threshold_I_SG & 
                                                  G3BP_E_S_BC >= BC_Threshold_E_SG)) &
                                                  I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                  E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(G3BP_EvI_M = G3BP_Peaks_Filtered$G3BP_EvI_M, G3BP_EvI_S = G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(G3BP_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(G3BP_EvI_M), y = log2(G3BP_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Mock Enriched/Input)', y = 'log2(Stress Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Stress Granule HuR Peaks: Mock vs Stress of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Stress Granule - mRNA specific:
G3BP_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                   finalized_annotation == "3'UTR" | 
                                                   finalized_annotation == "CDS" | 
                                                   finalized_annotation == "intron" | 
                                                   finalized_annotation == "CDS_RI" |
                                                   finalized_annotation == "DS10K") &
                                                  ((G3BP_I_M_BC >= BC_Threshold_I_SG & 
                                                    G3BP_E_M_BC >= BC_Threshold_E_SG) | 
                                                   (G3BP_I_S_BC >= BC_Threshold_I_SG & 
                                                      G3BP_E_S_BC >= BC_Threshold_E_SG)) &
                                                  I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                  E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(G3BP_EvI_M = G3BP_Peaks_Filtered$G3BP_EvI_M, G3BP_EvI_S = G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(G3BP_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(G3BP_EvI_M), y = log2(G3BP_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Mock Enriched/Input)', y = 'log2(Stress Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Stress Granule HuR Peaks: Mock vs Stress of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## Stress Granule - non-mRNA specific:
G3BP_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                   finalized_annotation != "3'UTR" & 
                                                   finalized_annotation != "CDS" & 
                                                   finalized_annotation != "intron" &  
                                                   finalized_annotation != "CDS_RI" &
                                                   finalized_annotation != "DS10K" &
                                                   finalized_annotation != "UnAn") & 
                                                  ((G3BP_I_M_BC >= BC_Threshold_I_SG & 
                                                    G3BP_E_M_BC >= BC_Threshold_E_SG) | 
                                                   (G3BP_I_S_BC >= BC_Threshold_I_SG & 
                                                      G3BP_E_S_BC >= BC_Threshold_E_SG)) &
                                                  I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                  E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(G3BP_EvI_M = G3BP_Peaks_Filtered$G3BP_EvI_M, G3BP_EvI_S = G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(G3BP_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(G3BP_EvI_M), y = log2(G3BP_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Mock Enriched/Input)', y = 'log2(Stress Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Stress Granule HuR Peaks: Mock vs Stress of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

####################################################################################################################

## Compartment Comparison
####################################################################################################################
# Mock NLS vs NES:
NLS_NES_Peaks_Filtered = peakEnrichment %>% filter(grouped_annotation != 'UnAn' & ((NLS_I_M_BC >= BC_Threshold_I & 
                                                      NLS_E_M_BC >= BC_Threshold_E) | 
                                                     (NES_I_M_BC >= BC_Threshold_I & 
                                                     NES_E_M_BC >= BC_Threshold_E)) &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                     E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_NES_Peaks_Filtered$NLS_EvI_M, NES_EvI_M = NLS_NES_Peaks_Filtered$NES_EvI_M)
data$annotation = factor(NLS_NES_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(NLS_EvI_M), y = log2(NES_EvI_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Cytoplasm Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: Nuclear vs Cytoplasm of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Mock NLS vs NES - mRNA specific:
NLS_NES_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                      finalized_annotation == "3'UTR" | 
                                                      finalized_annotation == "CDS" | 
                                                      finalized_annotation == "intron" | 
                                                      finalized_annotation == "CDS_RI" |
                                                      finalized_annotation == "DS10K") &
                                                     ((NLS_I_M_BC >= BC_Threshold_I & 
                                                       NLS_E_M_BC >= BC_Threshold_E) | 
                                                      (NES_I_M_BC >= BC_Threshold_I & 
                                                         NES_E_M_BC >= BC_Threshold_E)) &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                     E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_NES_Peaks_Filtered$NLS_EvI_M, NES_EvI_M = NLS_NES_Peaks_Filtered$NES_EvI_M)
data$annotation = factor(NLS_NES_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(NLS_EvI_M), y = log2(NES_EvI_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Cytoplasm Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: Nuclear vs Cytoplasm of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Mock NLS vs NES - non-mRNA specific:
NLS_NES_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                      finalized_annotation != "3'UTR" & 
                                                      finalized_annotation != "CDS" & 
                                                      finalized_annotation != "intron" &  
                                                      finalized_annotation != "CDS_RI" &
                                                      finalized_annotation != "DS10K" &
                                                      finalized_annotation != "UnAn") & 
                                                     ((NLS_I_M_BC >= BC_Threshold_I & 
                                                         NLS_E_M_BC >= BC_Threshold_E) | 
                                                        (NES_I_M_BC >= BC_Threshold_I & 
                                                           NES_E_M_BC >= BC_Threshold_E)) &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                     E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_NES_Peaks_Filtered$NLS_EvI_M, NES_EvI_M = NLS_NES_Peaks_Filtered$NES_EvI_M)
data$annotation = factor(NLS_NES_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(NLS_EvI_M), y = log2(NES_EvI_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Nuclear Enriched/Input)', y = 'log2(Cytoplasm Enriched/Input)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: Nuclear vs Cytoplasm of Enriched/Input (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

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

## FIGURE2 Fractionation and CoCLIP comparison
####################################################################################################################
# Mock Comparison
Mock_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != 'UnAn') & 
                                                  ((NLS_E_M_BC >= BC_Threshold_E & 
                                                      NLS_E_M >= median(NLS_E_M) * rowSum_Multiplier_E &
                                                      NES_E_M_BC >= BC_Threshold_E &
                                                      NES_E_M >= median(NES_E_M) * rowSum_Multiplier_E) |
                                                     (Nuc_F_M_BC >= BC_Threshold_F & 
                                                        Nuc_F_M >= median(Nuc_F_M) * rowSum_Multiplier_F & 
                                                        Cyto_F_M_BC >= BC_Threshold_F &
                                                        Cyto_F_M >= median(Cyto_F_M) * rowSum_Multiplier_F)))

data = data.frame(E_NvC_M = Mock_Peaks_Filtered$E_NvC_M, F_NvC_M = Mock_Peaks_Filtered$F_NvC_M)
data$annotation = factor(Mock_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(E_NvC_M), y = log2(F_NvC_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Cytoplasm)', y = 'log2(FracCLIP Nuclear/Cytoplasm)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: CoCLIP vs FracCLIP of Nuclear/Cytoplasm (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Mock Comparison - mRNA specific:
Mock_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                   finalized_annotation == "3'UTR" | 
                                                   finalized_annotation == "CDS" | 
                                                   finalized_annotation == "intron" | 
                                                   finalized_annotation == "CDS_RI" |
                                                   finalized_annotation == "DS10K") & 
                                                  ((NLS_E_M_BC >= BC_Threshold_E & 
                                                      NLS_E_M >= median(NLS_E_M) * rowSum_Multiplier_E &
                                                      NES_E_M_BC >= BC_Threshold_E &
                                                      NES_E_M >= median(NES_E_M) * rowSum_Multiplier_E) |
                                                     (Nuc_F_M_BC >= BC_Threshold_F & 
                                                        Nuc_F_M >= median(Nuc_F_M) * rowSum_Multiplier_F & 
                                                        Cyto_F_M_BC >= BC_Threshold_F &
                                                        Cyto_F_M >= median(Cyto_F_M) * rowSum_Multiplier_F)))

data = data.frame(E_NvC_M = Mock_Peaks_Filtered$E_NvC_M, F_NvC_M = Mock_Peaks_Filtered$F_NvC_M)
data$annotation = factor(Mock_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(E_NvC_M), y = log2(F_NvC_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Cytoplasm)', y = 'log2(FracCLIP Nuclear/Cytoplasm)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: CoCLIP vs FracCLIP of Nuclear/Cytoplasm (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Mock Comparison - non-mRNA specific:
Mock_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                   finalized_annotation != "3'UTR" & 
                                                   finalized_annotation != "CDS" & 
                                                   finalized_annotation != "intron" &  
                                                   finalized_annotation != "CDS_RI" &
                                                   finalized_annotation != "DS10K" &
                                                   finalized_annotation != "UnAn") & 
                                                  ((NLS_E_M_BC >= BC_Threshold_E & 
                                                      NLS_E_M >= median(NLS_E_M) * rowSum_Multiplier_E &
                                                      NES_E_M_BC >= BC_Threshold_E &
                                                      NES_E_M >= median(NES_E_M) * rowSum_Multiplier_E) |
                                                     (Nuc_F_M_BC >= BC_Threshold_F & 
                                                        Nuc_F_M >= median(Nuc_F_M) * rowSum_Multiplier_F & 
                                                        Cyto_F_M_BC >= BC_Threshold_F &
                                                        Cyto_F_M >= median(Cyto_F_M) * rowSum_Multiplier_F)))

data = data.frame(E_NvC_M = Mock_Peaks_Filtered$E_NvC_M, F_NvC_M = Mock_Peaks_Filtered$F_NvC_M)
data$annotation = factor(Mock_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(E_NvC_M), y = log2(F_NvC_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Cytoplasm)', y = 'log2(FracCLIP Nuclear/Cytoplasm)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: CoCLIP vs FracCLIP of Nuclear/Cytoplasm (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress Comparison
Stress_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != 'UnAn') & 
                                                    ((NLS_E_S_BC >= BC_Threshold_E & 
                                                        NLS_E_S >= median(NLS_E_S) * rowSum_Multiplier_E &
                                                        NES_E_S_BC >= BC_Threshold_E &
                                                        NES_E_S >= median(NES_E_S) * rowSum_Multiplier_E) |
                                                       (Nuc_F_S_BC >= BC_Threshold_F & 
                                                          Nuc_F_S >= median(Nuc_F_S) * rowSum_Multiplier_F & 
                                                          Cyto_F_S_BC >= BC_Threshold_F &
                                                          Cyto_F_S >= median(Cyto_F_S) * rowSum_Multiplier_F)))

data = data.frame(E_NvC_S= Stress_Peaks_Filtered$E_NvC_S, F_NvC_S = Stress_Peaks_Filtered$F_NvC_S)
data$annotation = factor(Stress_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(data, aes(x = log2(E_NvC_S), y = log2(F_NvC_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Cytoplasm)', y = 'log2(FracCLIP Nuclear/Cytoplasm)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Stress HuR Peaks: CoCLIP vs FracCLIP of Nuclear/Cytoplasm (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress Comparison - mRNA specific
Stress_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation == "5'UTR" | 
                                                     finalized_annotation == "3'UTR" | 
                                                     finalized_annotation == "CDS" | 
                                                     finalized_annotation == "intron" | 
                                                     finalized_annotation == "CDS_RI" |
                                                     finalized_annotation == "DS10K") & 
                                                    ((NLS_E_S_BC >= BC_Threshold_E & 
                                                        NLS_E_S >= median(NLS_E_S) * rowSum_Multiplier_E &
                                                        NES_E_S_BC >= BC_Threshold_E &
                                                        NES_E_S >= median(NES_E_S) * rowSum_Multiplier_E) |
                                                       (Nuc_F_S_BC >= BC_Threshold_F & 
                                                          Nuc_F_S >= median(Nuc_F_S) * rowSum_Multiplier_F & 
                                                          Cyto_F_S_BC >= BC_Threshold_F &
                                                          Cyto_F_S >= median(Cyto_F_S) * rowSum_Multiplier_F)))

data = data.frame(E_NvC_S= Stress_Peaks_Filtered$E_NvC_S, F_NvC_S = Stress_Peaks_Filtered$F_NvC_S)
data$annotation = factor(Stress_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(data, aes(x = log2(E_NvC_S), y = log2(F_NvC_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Cytoplasm)', y = 'log2(FracCLIP Nuclear/Cytoplasm)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Stress HuR Peaks: CoCLIP vs FracCLIP of Nuclear/Cytoplasm (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress Comparison - non-mRNA specific
Stress_Peaks_Filtered = peakEnrichment %>% filter((finalized_annotation != "5'UTR" & 
                                                     finalized_annotation != "3'UTR" & 
                                                     finalized_annotation != "CDS" & 
                                                     finalized_annotation != "intron" &  
                                                     finalized_annotation != "CDS_RI" &
                                                     finalized_annotation != "DS10K" &
                                                     finalized_annotation != "UnAn") & 
                                                    ((NLS_E_S_BC >= BC_Threshold_E & 
                                                        NLS_E_S >= median(NLS_E_S) * rowSum_Multiplier_E &
                                                        NES_E_S_BC >= BC_Threshold_E &
                                                        NES_E_S >= median(NES_E_S) * rowSum_Multiplier_E) |
                                                       (Nuc_F_S_BC >= BC_Threshold_F & 
                                                          Nuc_F_S >= median(Nuc_F_S) * rowSum_Multiplier_F & 
                                                          Cyto_F_S_BC >= BC_Threshold_F &
                                                          Cyto_F_S >= median(Cyto_F_S) * rowSum_Multiplier_F)))

data = data.frame(E_NvC_S= Stress_Peaks_Filtered$E_NvC_S, F_NvC_S = Stress_Peaks_Filtered$F_NvC_S)
data$annotation = factor(Stress_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(data, aes(x = log2(E_NvC_S), y = log2(F_NvC_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Cytoplasm)', y = 'log2(FracCLIP Nuclear/Cytoplasm)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Stress HuR Peaks: CoCLIP vs FracCLIP of Nuclear/Cytoplasm (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


####################################################################################################################

## Input and CoCLIP comparison
####################################################################################################################
# Nuclear Mock:
NLS_Mock_Peaks_Filtered = peakRowSum %>% filter(grouped_annotation != 'UnAn' & NLS_I_M_BC >= BC_Threshold_I & 
                                                      NLS_E_M_BC >= BC_Threshold_E  &
                                                      I_rowSum >= median(peakRowSum$I_rowSum)*2)
NLS_Mock_Peaks_Filtered$grouped_annotation = factor(NLS_Mock_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(NLS_Mock_Peaks_Filtered, aes(x = log10(NLS_E_M), y = log10(I_M), color = grouped_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'NLS Enrich', y = 'Input') + 
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Mock HuR Peaks: Nuclear CoCLIP vs Input (',  nrow(NLS_Mock_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Nuclear Mock - mRNA specific:
NLS_Mock_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation == "5'UTR" | 
                                                   finalized_annotation == "3'UTR" | 
                                                   finalized_annotation == "CDS" | 
                                                   finalized_annotation == "intron" | 
                                                   finalized_annotation == "CDS_RI" |
                                                   finalized_annotation == "DS10K") &
                                                  NLS_I_M_BC >= BC_Threshold_I & 
                                                  NLS_E_M_BC >= BC_Threshold_E  &
                                                  I_rowSum >= median(peakRowSum$I_rowSum)*2)
NLS_Mock_Peaks_Filtered$finalized_annotation = factor(NLS_Mock_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(NLS_Mock_Peaks_Filtered, aes(x = log10(NLS_E_M), y = log10(I_M), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'NLS Enrich', y = 'Input') + 
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Mock HuR Peaks: Nuclear CoCLIP vs Input (',  nrow(NLS_Mock_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Nuclear Mock - non-mRNA specific:
NLS_Mock_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation != "5'UTR" & 
                                                   finalized_annotation != "3'UTR" & 
                                                   finalized_annotation != "CDS" & 
                                                   finalized_annotation != "intron" &  
                                                   finalized_annotation != "CDS_RI" &
                                                   finalized_annotation != "DS10K" &
                                                   finalized_annotation != "UnAn") & 
                                                  NLS_I_M_BC >= BC_Threshold_I & 
                                                  NLS_E_M_BC >= BC_Threshold_E  &
                                                  I_rowSum >= median(peakRowSum$I_rowSum)*2)
NLS_Mock_Peaks_Filtered$finalized_annotation = factor(NLS_Mock_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(NLS_Mock_Peaks_Filtered, aes(x = log10(NLS_E_M), y = log10(I_M), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'NLS Enrich', y = 'Input') + 
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Mock HuR Peaks: Nuclear CoCLIP vs Input (',  nrow(NLS_Mock_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Cytoplasm Mock:
NES_Mock_Peaks_Filtered = peakRowSum %>% filter(grouped_annotation != 'UnAn' & NES_I_M_BC >= BC_Threshold_I & 
                                                      NES_E_M_BC >= BC_Threshold_E &
                                                      I_rowSum >= median(peakRowSum$I_rowSum)*2)
NES_Mock_Peaks_Filtered$grouped_annotation = factor(NES_Mock_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(NES_Mock_Peaks_Filtered, aes(x = log10(NES_E_M), y = log10(I_M), color = grouped_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'NES Enrich', y = 'Input') +
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Mock HuR Peaks: Cytoplasm CoCLIP vs Input (',  nrow(NES_Mock_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Cytoplasm Mock - mRNA specific:
NES_Mock_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation == "5'UTR" | 
                                                   finalized_annotation == "3'UTR" | 
                                                   finalized_annotation == "CDS" | 
                                                   finalized_annotation == "intron" | 
                                                   finalized_annotation == "CDS_RI" |
                                                   finalized_annotation == "DS10K") &
                                                  NES_I_M_BC >= BC_Threshold_I & 
                                                  NES_E_M_BC >= BC_Threshold_E &
                                                  I_rowSum >= median(peakRowSum$I_rowSum)*2)
NES_Mock_Peaks_Filtered$finalized_annotation = factor(NES_Mock_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(NES_Mock_Peaks_Filtered, aes(x = log10(NES_E_M), y = log10(I_M), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'NES Enrich', y = 'Input') +
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Mock HuR Peaks: Cytoplasm CoCLIP vs Input (',  nrow(NES_Mock_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Cytoplasm Mock - non-mRNA specific:
NES_Mock_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation != "5'UTR" & 
                                                   finalized_annotation != "3'UTR" & 
                                                   finalized_annotation != "CDS" & 
                                                   finalized_annotation != "intron" &  
                                                   finalized_annotation != "CDS_RI" &
                                                   finalized_annotation != "DS10K" &
                                                   finalized_annotation != "UnAn") & 
                                                  NES_I_M_BC >= BC_Threshold_I & 
                                                  NES_E_M_BC >= BC_Threshold_E &
                                                  I_rowSum >= median(peakRowSum$I_rowSum)*2)
NES_Mock_Peaks_Filtered$finalized_annotation = factor(NES_Mock_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(NES_Mock_Peaks_Filtered, aes(x = log10(NES_E_M), y = log10(I_M), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'NES Enrich', y = 'Input') +
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Mock HuR Peaks: Cytoplasm CoCLIP vs Input (',  nrow(NES_Mock_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress Granule Mock:
G3BP_Mock_Peaks_Filtered = peakRowSum %>% filter(grouped_annotation != 'UnAn' & G3BP_I_M_BC >= BC_Threshold_I_SG & 
                                                       G3BP_E_M_BC >= BC_Threshold_E_SG &
                                                       I_rowSum >= median(peakRowSum$I_rowSum)*2)
G3BP_Mock_Peaks_Filtered$grouped_annotation = factor(G3BP_Mock_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(G3BP_Mock_Peaks_Filtered, aes(x = log10(G3BP_E_M), y = log10(I_M), color = grouped_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'G3BP Enrich', y = 'Input') +
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Mock HuR Peaks: Stress Granule CoCLIP vs Input (',  nrow(G3BP_Mock_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress Granule Mock - mRNA specific:
G3BP_Mock_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation == "5'UTR" | 
                                                    finalized_annotation == "3'UTR" | 
                                                    finalized_annotation == "CDS" | 
                                                    finalized_annotation == "intron" | 
                                                    finalized_annotation == "CDS_RI" |
                                                    finalized_annotation == "DS10K") &
                                                   G3BP_I_M_BC >= BC_Threshold_I_SG & 
                                                   G3BP_E_M_BC >= BC_Threshold_E_SG &
                                                   I_rowSum >= median(peakRowSum$I_rowSum)*2)
G3BP_Mock_Peaks_Filtered$finalized_annotation = factor(G3BP_Mock_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(G3BP_Mock_Peaks_Filtered, aes(x = log10(G3BP_E_M), y = log10(I_M), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'G3BP Enrich', y = 'Input') +
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Mock HuR Peaks: Stress Granule CoCLIP vs Input (',  nrow(G3BP_Mock_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress Granule Mock - non-mRNA specific:
G3BP_Mock_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation != "5'UTR" & 
                                                    finalized_annotation != "3'UTR" & 
                                                    finalized_annotation != "CDS" & 
                                                    finalized_annotation != "intron" &  
                                                    finalized_annotation != "CDS_RI" &
                                                    finalized_annotation != "DS10K" &
                                                    finalized_annotation != "UnAn") & 
                                                   G3BP_I_M_BC >= BC_Threshold_I_SG & 
                                                   G3BP_E_M_BC >= BC_Threshold_E_SG &
                                                   I_rowSum >= median(peakRowSum$I_rowSum)*2)
G3BP_Mock_Peaks_Filtered$finalized_annotation = factor(G3BP_Mock_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(G3BP_Mock_Peaks_Filtered, aes(x = log10(G3BP_E_M), y = log10(I_M), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'G3BP Enrich', y = 'Input') +
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Mock HuR Peaks: Stress Granule CoCLIP vs Input (',  nrow(G3BP_Mock_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Nuclear Stress:
NLS_Stress_Peaks_Filtered = peakRowSum %>% filter(grouped_annotation != 'UnAn' & NLS_I_S_BC >= BC_Threshold_I & 
                                                        NLS_E_S_BC >= BC_Threshold_E &
                                                        I_rowSum >= median(peakRowSum$I_rowSum)*2)
NLS_Stress_Peaks_Filtered$grouped_annotation = factor(NLS_Stress_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(NLS_Stress_Peaks_Filtered, aes(x = log10(NLS_E_S), y = log10(I_S), color = grouped_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'NLS Enrich', y = 'Input') +
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Stress HuR Peaks: Nuclear CoCLIP vs Input (',  nrow(NLS_Stress_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Nuclear Stress - mRNA specific:
NLS_Stress_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation == "5'UTR" | 
                                                     finalized_annotation == "3'UTR" | 
                                                     finalized_annotation == "CDS" | 
                                                     finalized_annotation == "intron" | 
                                                     finalized_annotation == "CDS_RI" |
                                                     finalized_annotation == "DS10K") &
                                                    NLS_I_S_BC >= BC_Threshold_I & 
                                                    NLS_E_S_BC >= BC_Threshold_E &
                                                    I_rowSum >= median(peakRowSum$I_rowSum)*2)
NLS_Stress_Peaks_Filtered$finalized_annotation = factor(NLS_Stress_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(NLS_Stress_Peaks_Filtered, aes(x = log10(NLS_E_S), y = log10(I_S), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'NLS Enrich', y = 'Input') +
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Stress HuR Peaks: Nuclear CoCLIP vs Input (',  nrow(NLS_Stress_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Nuclear Stress - non-mRNA specific:
NLS_Stress_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation != "5'UTR" & 
                                                     finalized_annotation != "3'UTR" & 
                                                     finalized_annotation != "CDS" & 
                                                     finalized_annotation != "intron" &  
                                                     finalized_annotation != "CDS_RI" &
                                                     finalized_annotation != "DS10K" &
                                                     finalized_annotation != "UnAn") &
                                                    NLS_I_S_BC >= BC_Threshold_I & 
                                                    NLS_E_S_BC >= BC_Threshold_E &
                                                    I_rowSum >= median(peakRowSum$I_rowSum)*2)
NLS_Stress_Peaks_Filtered$finalized_annotation = factor(NLS_Stress_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(NLS_Stress_Peaks_Filtered, aes(x = log10(NLS_E_S), y = log10(I_S), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'NLS Enrich', y = 'Input') +
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Stress HuR Peaks: Nuclear CoCLIP vs Input (',  nrow(NLS_Stress_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Cytoplasm Stress:
NES_Stress_Peaks_Filtered = peakRowSum %>% filter(grouped_annotation != 'UnAn' & NES_I_S_BC >= BC_Threshold_I & 
                                                        NES_E_S_BC >= BC_Threshold_E &
                                                        I_rowSum >= median(peakRowSum$I_rowSum)*2)
NES_Stress_Peaks_Filtered$grouped_annotation = factor(NES_Stress_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(NES_Stress_Peaks_Filtered, aes(x = log10(NES_E_S), y = log10(I_S), color = grouped_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'NES Enrich', y = 'Input') +
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Stress HuR Peaks: Cytoplasm CoCLIP vs Input (',  nrow(NES_Stress_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Cytoplasm Stress - mRNA specific:
NES_Stress_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation == "5'UTR" | 
                                                     finalized_annotation == "3'UTR" | 
                                                     finalized_annotation == "CDS" | 
                                                     finalized_annotation == "intron" | 
                                                     finalized_annotation == "CDS_RI" |
                                                     finalized_annotation == "DS10K") &
                                                    NES_I_S_BC >= BC_Threshold_I & 
                                                    NES_E_S_BC >= BC_Threshold_E &
                                                    I_rowSum >= median(peakRowSum$I_rowSum)*2)
NES_Stress_Peaks_Filtered$finalized_annotation = factor(NES_Stress_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(NES_Stress_Peaks_Filtered, aes(x = log10(NES_E_S), y = log10(I_S), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'NES Enrich', y = 'Input') +
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Stress HuR Peaks: Cytoplasm CoCLIP vs Input (',  nrow(NES_Stress_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Cytoplasm Stress - non-mRNA specific:
NES_Stress_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation != "5'UTR" & 
                                                     finalized_annotation != "3'UTR" & 
                                                     finalized_annotation != "CDS" & 
                                                     finalized_annotation != "intron" &  
                                                     finalized_annotation != "CDS_RI" &
                                                     finalized_annotation != "DS10K" &
                                                     finalized_annotation != "UnAn") & 
                                                    NES_I_S_BC >= BC_Threshold_I & 
                                                    NES_E_S_BC >= BC_Threshold_E &
                                                    I_rowSum >= median(peakRowSum$I_rowSum)*2)
NES_Stress_Peaks_Filtered$finalized_annotation = factor(NES_Stress_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(NES_Stress_Peaks_Filtered, aes(x = log10(NES_E_S), y = log10(I_S), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'NES Enrich', y = 'Input') +
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Stress HuR Peaks: Cytoplasm CoCLIP vs Input (',  nrow(NES_Stress_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress Granule Stress:
G3BP_Stress_Peaks_Filtered = peakRowSum %>% filter(grouped_annotation != 'UnAn' & G3BP_I_S_BC >= BC_Threshold_I_SG & 
                                                         G3BP_E_S_BC >= BC_Threshold_E_SG &
                                                         I_rowSum >= median(peakRowSum$I_rowSum)*2)
G3BP_Stress_Peaks_Filtered$grouped_annotation = factor(G3BP_Stress_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "DS10K"))

ggplot(G3BP_Stress_Peaks_Filtered, aes(x = log10(G3BP_E_S), y = log10(I_S), color = grouped_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'G3BP Enrich', y = 'Input') +
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Stress HuR Peaks: Stress Granule CoCLIP vs Input (',  nrow(G3BP_Stress_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress Granule Stress - mRNA specific:
G3BP_Stress_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation == "5'UTR" | 
                                                      finalized_annotation == "3'UTR" | 
                                                      finalized_annotation == "CDS" | 
                                                      finalized_annotation == "intron" | 
                                                      finalized_annotation == "CDS_RI" |
                                                      finalized_annotation == "DS10K") &
                                                     G3BP_I_S_BC >= BC_Threshold_I_SG & 
                                                     G3BP_E_S_BC >= BC_Threshold_E_SG &
                                                     I_rowSum >= median(peakRowSum$I_rowSum)*2)
G3BP_Stress_Peaks_Filtered$finalized_annotation = factor(G3BP_Stress_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_RI", 'DS10K'))

ggplot(G3BP_Stress_Peaks_Filtered, aes(x = log10(G3BP_E_S), y = log10(I_S), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'G3BP Enrich', y = 'Input') +
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Stress HuR Peaks: Stress Granule CoCLIP vs Input (',  nrow(G3BP_Stress_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress Granule Stress - non-mRNA specific:
G3BP_Stress_Peaks_Filtered = peakRowSum %>% filter((finalized_annotation != "5'UTR" & 
                                                      finalized_annotation != "3'UTR" & 
                                                      finalized_annotation != "CDS" & 
                                                      finalized_annotation != "intron" &  
                                                      finalized_annotation != "CDS_RI" &
                                                      finalized_annotation != "DS10K" &
                                                      finalized_annotation != "UnAn") & 
                                                     G3BP_I_S_BC >= BC_Threshold_I_SG & 
                                                     G3BP_E_S_BC >= BC_Threshold_E_SG &
                                                     I_rowSum >= median(peakRowSum$I_rowSum)*2)
G3BP_Stress_Peaks_Filtered$finalized_annotation = factor(G3BP_Stress_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'nC_RI', 'TE', 'Other'))

ggplot(G3BP_Stress_Peaks_Filtered, aes(x = log10(G3BP_E_S), y = log10(I_S), color = finalized_annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'G3BP Enrich', y = 'Input') +
  xlim(c(0, 3)) +
  ylim(c(0, 3)) +
  ggtitle(paste0('Stress HuR Peaks: Stress Granule CoCLIP vs Input (',  nrow(G3BP_Stress_Peaks_Filtered), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
####################################################################################################################

