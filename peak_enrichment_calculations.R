## CoCLIP Analysis: 
## Peak Enrichment Calculation
## Written by Soon Yi
## Last Edit: 2023-08-31

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
peaksMatrix_PATH = 'L:/My Drive/CWRU/PhD/Luna Lab/1. coCLIP/Analysis/peaks/'    ## Use this for windows machine
# peaksMatrix_PATH = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/'
peaksMatrix_FILE = 'Combined_peakCoverage_groomed_normalized_annotated.txt'

peakMatrix = read_delim(paste0(peaksMatrix_PATH, peaksMatrix_FILE), show_col_types = FALSE)
peakMatrix = peakMatrix %>% mutate_at('TOTAL_BC', as.numeric)
# peakMatrix = peakMatrix %>% filter(annotation_count == 1)


## Filter by median
# peakCount_median = median(peakMatrix[, colnames(peakMatrix)[6:63]][peakMatrix[, colnames(peakMatrix)[6:63]] != 0], na.rm = TRUE)
# peakMatrix = peakMatrix %>% filter(if_any(all_of(comparison_vector), ~ . > peakCount_median ))

## Instead of filtering, add pseudocounts.
## This might be better because we will be filtering by BC down the road.
pseudoCount = min(peakMatrix[, colnames(peakMatrix)[6:63]][peakMatrix[, colnames(peakMatrix)[6:63]] != 0], na.rm = TRUE)
peakMatrix[, colnames(peakMatrix)[6:63]] = peakMatrix[, colnames(peakMatrix)[6:63]] + pseudoCount

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
peakMatrix$F_rowSum = rowSums(peakMatrix[, c(Nuc_F_M, Nuc_F_S, Cyto_F_M, Cyto_F_S)])
peakMatrix$I_rowSum = rowSums(peakMatrix[, c(NLS_I_M, NLS_I_S, NES_I_M, NES_I_S, G3BP_I_M, G3BP_I_S)])
peakMatrix$E_rowSum = rowSums(peakMatrix[, c(NLS_E_M, NLS_E_S, NES_E_M, NES_E_S, G3BP_E_M, G3BP_E_S)])

rowSum_columns = c('F_rowSum', 'I_rowSum', 'E_rowSum')

## Add deep intergenic to the finalized annotations:
peakMatrix = peakMatrix %>% mutate(annotation = ifelse(is.na(annotation), 'deep intergenic', annotation))
peakMatrix = peakMatrix %>% mutate(finalized_annotation = ifelse(is.na(finalized_annotation), 'deep intergenic', finalized_annotation))

####################################################################################################################

## PCA of Normalized Peaks:
####################################################################################################################
PCA_data = peakMatrix[, c(Nuc_F_M, Nuc_F_S, Cyto_F_M, Cyto_F_S, NLS_I_M, NLS_I_S, NES_I_M, NES_I_S, G3BP_I_M, G3BP_I_S, NLS_E_M, NLS_E_S, NES_E_M, NES_E_S, G3BP_E_M, G3BP_E_S)]
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

## BC Filter Criteria:
BC_Threshold_Fraction = 2
BC_Threshold_Input = 1
BC_Threshold_Input_SG = 2
BC_Threshold_CoCLIP = 2
BC_Threshold_CoCLIP_SG = 3

## Subset of Peaks for downstream analysis:
####################################################################################################################
Peak_all_M = (peakMatrix[, c(inert_columns, Nuc_F_M, Cyto_F_M, NLS_I_M, NES_I_M, G3BP_I_M, NLS_E_M, NES_E_M, G3BP_E_M, BC_columns)] 
              %>% filter((Nuc_F_M_BC >= BC_Threshold_Fraction & Cyto_F_M_BC >= BC_Threshold_Fraction) | 
                           (NLS_I_M_BC >= BC_Threshold_Input & NES_I_M_BC >= BC_Threshold_Input & G3BP_I_M_BC >= BC_Threshold_Input_SG) |
                           (NLS_E_M_BC >= BC_Threshold_CoCLIP & NES_E_M_BC >= BC_Threshold_CoCLIP & G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG)))

Peak_all_S = (peakMatrix[, c(inert_columns, Nuc_F_S, Cyto_F_S, NLS_I_S, NES_I_S, G3BP_I_S, NLS_E_S, NES_E_S, G3BP_E_S, BC_columns)] 
              %>% filter((Nuc_F_S_BC >= BC_Threshold_Fraction & Cyto_F_S_BC >= BC_Threshold_Fraction) |
                           (NLS_I_S_BC >= BC_Threshold_Input & NES_I_S_BC >= BC_Threshold_Input & G3BP_I_S_BC >= BC_Threshold_CoCLIP_SG) | 
                           (NLS_E_S_BC >= BC_Threshold_CoCLIP & NES_E_S_BC >= BC_Threshold_CoCLIP & G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG)))

Peak_all = (peakMatrix[, c(inert_columns, Nuc_F_S, Cyto_F_S, NLS_I_S, NES_I_S, G3BP_I_S, NLS_E_S, NES_E_S, G3BP_E_S, BC_columns)] 
            %>% filter(((Nuc_F_M_BC >= BC_Threshold_Fraction & Cyto_F_M_BC >= BC_Threshold_Fraction) | 
                          (NLS_I_M_BC >= BC_Threshold_Input & NES_I_M_BC >= BC_Threshold_Input & G3BP_I_M_BC >= BC_Threshold_CoCLIP_SG) | 
                          (NLS_E_M_BC >= BC_Threshold_CoCLIP & NES_E_M_BC >= BC_Threshold_CoCLIP & G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG)) | 
                         ((Nuc_F_S_BC >= BC_Threshold_Fraction & Cyto_F_S_BC >= BC_Threshold_Fraction) | 
                            (NLS_I_S_BC >= BC_Threshold_Input & NES_I_S_BC >= BC_Threshold_Input & G3BP_I_S_BC >= BC_Threshold_CoCLIP_SG) | 
                            NLS_E_S_BC >= BC_Threshold_CoCLIP & NES_E_S_BC >= BC_Threshold_CoCLIP & G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG)))


Peak_F_Nuc_M = (peakMatrix[, c(inert_columns, Nuc_F_M, BC_columns, rowSum_columns)] 
                %>% filter(Nuc_F_M_BC >= BC_Threshold_Fraction))

Peak_F_Nuc_S = (peakMatrix[, c(inert_columns, Nuc_F_S, BC_columns, rowSum_columns)] 
                %>% filter(Nuc_F_S_BC >= BC_Threshold_Fraction))

Peak_F_Cyt_M = (peakMatrix[, c(inert_columns, Cyto_F_M, BC_columns, rowSum_columns)] 
                %>% filter(Cyto_F_M_BC >= BC_Threshold_Fraction))

Peak_F_Cyt_S = (peakMatrix[, c(inert_columns, Cyto_F_S, BC_columns, rowSum_columns)] 
                %>% filter(Cyto_F_S_BC >= BC_Threshold_Fraction))

Peak_F_all_M = (peakMatrix[, c(inert_columns, Nuc_F_M, Cyto_F_M, BC_columns, rowSum_columns)] 
                %>% filter(Nuc_F_M_BC >= BC_Threshold_Fraction &
                           Cyto_F_M_BC >= BC_Threshold_Fraction))

Peak_F_all_S = (peakMatrix[, c(inert_columns, Nuc_F_S, Cyto_F_S, BC_columns, rowSum_columns)] 
                %>% filter(Nuc_F_S_BC >= BC_Threshold_Fraction & 
                           Cyto_F_S_BC >= BC_Threshold_Fraction))

Peak_Co_Input_M = (peakMatrix[, c(inert_columns, NLS_I_M, NES_I_M, G3BP_I_M, BC_columns, rowSum_columns)] 
                   %>% filter(NLS_I_M_BC >= BC_Threshold_Input & 
                              NES_I_M_BC >= BC_Threshold_Input & 
                              G3BP_I_M_BC >= BC_Threshold_Input_SG))

Peak_Co_Input_S = (peakMatrix[, c(inert_columns, NLS_I_S, NES_I_S, G3BP_I_S, BC_columns, rowSum_columns)] 
                   %>% filter(NLS_I_S_BC >= BC_Threshold_Input & 
                              NES_I_S_BC >= BC_Threshold_Input & 
                              G3BP_I_S_BC >= BC_Threshold_Input_SG))

Peak_Co_NLS_M = (peakMatrix[, c(inert_columns, NLS_E_M, BC_columns, rowSum_columns)] 
                 %>% filter(NLS_E_M_BC >= BC_Threshold_CoCLIP))

Peak_Co_NLS_S = (peakMatrix[, c(inert_columns, NLS_E_S, BC_columns, rowSum_columns)] 
                 %>% filter(NLS_E_S_BC >= BC_Threshold_CoCLIP))

Peak_Co_NES_M = (peakMatrix[, c(inert_columns, NES_E_M, BC_columns, rowSum_columns)] 
                 %>% filter(NES_E_M_BC >= BC_Threshold_CoCLIP))

Peak_Co_NES_S = (peakMatrix[, c(inert_columns, NES_E_S, BC_columns, rowSum_columns)] 
                 %>% filter(NES_E_S_BC >= BC_Threshold_CoCLIP))

Peak_Co_G3BP_M = (peakMatrix[, c(inert_columns, G3BP_E_M, BC_columns, rowSum_columns)] 
                  %>% filter(G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG))

Peak_Co_G3BP_S = (peakMatrix[, c(inert_columns, G3BP_E_S, BC_columns, rowSum_columns)] 
                  %>% filter(G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG))

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
Counts_perGene = Peak_all %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(TOTAL_TagCount), peakCounts = n())
Counts_perGene = Counts_perGene %>% mutate(tagDensity = tagCounts/peakCounts)

## All Mock tags
Counts_perGene_A_M= Peak_all_M %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(Nuc_F_M_1 + Nuc_F_M_2 + Nuc_F_M_3 + 
                                                                                                      Cyto_F_M_1 + Cyto_F_M_2 + Cyto_F_M_3 + 
                                                                                                      NLS_I_M_1 + NLS_I_M_2 + NES_I_M_1 + NES_I_M_2 + G3BP_I_M_1 + G3BP_I_M_2 + G3BP_I_M_3 + G3BP_I_M_4 +
                                                                                                      NLS_E_M_1 + NLS_E_M_2 + NLS_E_M_3 + NLS_E_M_4 + 
                                                                                                      NES_E_M_1 + NES_E_M_2 + NES_E_M_3 + NES_E_M_4 +
                                                                                                      G3BP_E_M_1 + G3BP_E_M_2 + G3BP_E_M_3 + G3BP_E_M_4 + G3BP_E_M_5 + G3BP_E_M_6), peakCounts = n())
Counts_perGene_A_M = Counts_perGene_A_M %>% mutate(tagDensity = tagCounts/peakCounts)

## All Stress tags
Counts_perGene_A_S= Peak_all_S %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(Nuc_F_S_1 + Nuc_F_S_2 + Nuc_F_S_3 +
                                                                                                      Cyto_F_S_1 + Cyto_F_S_2 + Cyto_F_S_3 +
                                                                                                      NLS_I_S_1 + NLS_I_S_2 + NES_I_S_1 + NES_I_S_2 + G3BP_I_S_1 + G3BP_I_S_2 + G3BP_I_S_3 + G3BP_I_S_4 +
                                                                                                      NLS_E_S_1 + NLS_E_S_2 + NLS_E_S_3 + NLS_E_S_4 +
                                                                                                      NES_E_S_1 + NES_E_S_2 + NES_E_S_3 + NES_E_S_4 +
                                                                                                      G3BP_E_S_1 + G3BP_E_S_2 + G3BP_E_S_3 + G3BP_E_S_4 + G3BP_E_S_5 + G3BP_E_S_6 + G3BP_E_S_7), peakCounts = n())
Counts_perGene_A_S = Counts_perGene_A_S %>% mutate(tagDensity = tagCounts/peakCounts)

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
## Mock vs Stress:
PeakDistribution_all_M = data.frame(table(Peak_all_M$grouped_annotation), row.names = 1)
PeakDistribution_all_S = data.frame(table(Peak_all_S$grouped_annotation), row.names = 1)
PeakDistribution_all = data.frame(table(Peak_all$grouped_annotation), row.names = 1)

colnames(PeakDistribution_all_M) = c('All_M')
colnames(PeakDistribution_all_S) = c('All_S')
colnames(PeakDistribution_all) = c('All_M+S')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_all_combined = cbind(PeakDistribution_all_M, 
                                      PeakDistribution_all_S, 
                                      PeakDistribution_all)
PeakDistribution_all_combined$Annotation = rownames(PeakDistribution_all_combined)

PeakDistribution_all_combined = PeakDistribution_all_combined %>%
  gather(key = "Source", value = "Freq", 'All_M', 'All_S', 'All_M+S') %>%
  select(Source, Freq, Annotation)

PeakDistribution_all_combined$Source = factor(PeakDistribution_all_combined$Source, levels = c('All_M+S', 'All_M', 'All_S'))
PeakDistribution_all_combined$Annotation = factor(PeakDistribution_all_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

ggplot(PeakDistribution_all_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Mock + Stress', 'Mock Only', 'Stress Only')) +
  ggtitle('All Raw Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## Fraction Distribution Stacked Bar Graph
PeakDistribution_all_combined = cbind(PeakDistribution_all_M/sum(PeakDistribution_all_M), 
                                      PeakDistribution_all_S/sum(PeakDistribution_all_S), 
                                      PeakDistribution_all/sum(PeakDistribution_all))
PeakDistribution_all_combined$Annotation = rownames(PeakDistribution_all_combined)

PeakDistribution_all_combined = PeakDistribution_all_combined %>%
  gather(key = "Source", value = "Freq", 'All_M', 'All_S', 'All_M+S') %>%
  select(Source, Freq, Annotation)

PeakDistribution_all_combined$Source = factor(PeakDistribution_all_combined$Source, levels = c('All_M+S', 'All_M', 'All_S'))
PeakDistribution_all_combined$Annotation = factor(PeakDistribution_all_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

ggplot(PeakDistribution_all_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Mock + Stress', 'Mock Only', 'Stress Only')) +
  ggtitle('All Peak Fractions Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")
####################################################################################################################

## Exploratory Stacked Bar Plots For Fractionation CLIP Peaks:
####################################################################################################################
## Mock vs Stress for Each Fraction
PeakDistribution_F_Nuc_M = data.frame(table(Peak_F_Nuc_M$grouped_annotation), row.names = 1)
PeakDistribution_F_Nuc_S = data.frame(table(Peak_F_Nuc_S$grouped_annotation), row.names = 1)

PeakDistribution_F_Cyt_M = data.frame(table(Peak_F_Cyt_M$grouped_annotation), row.names = 1)
PeakDistribution_F_Cyt_S = data.frame(table(Peak_F_Cyt_S$grouped_annotation), row.names = 1)

PeakDistribution_F_all_M = data.frame(table(Peak_F_all_M$grouped_annotation), row.names = 1)
PeakDistribution_F_all_S = data.frame(table(Peak_F_all_S$grouped_annotation), row.names = 1)

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')
colnames(PeakDistribution_F_all_M) = c('F_M_Nuc+Cyt')
colnames(PeakDistribution_F_all_S) = c('F_S_Nuc+Cyt')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Nuc_S, 
                                    PeakDistribution_F_Cyt_M, PeakDistribution_F_Cyt_S, 
                                    PeakDistribution_F_all_M, PeakDistribution_F_all_S)
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress', 'Both Fractions Mock','Both Fractions Stress')) +
  ggtitle('All Raw Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## Fraction Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), 
                                    PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S), 
                                    PeakDistribution_F_all_M/sum(PeakDistribution_F_all_M), PeakDistribution_F_all_S/sum(PeakDistribution_F_all_S))
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress', 'Both Fractions Mock','Both Fractions Stress')) +
  ggtitle('All Peak Fractions Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")
####################################################################################################################

## Exploratory Stacked Bar Plots For CoCLIP Peaks:
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
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  ggtitle('All Raw Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

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
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  ggtitle('All Peak Fractions Distributions') +
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

PeakDistribution_F_all_M = data.frame(data.frame(annotation = Peak_F_all_M$grouped_annotation, 
                                                 tagSum = rowSums(Peak_F_all_M[, c(Nuc_F_M, Cyto_F_M)])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_all_S = data.frame(data.frame(annotation = Peak_F_all_S$grouped_annotation, 
                                                 tagSum = rowSums(Peak_F_all_S[, c(Nuc_F_S, Cyto_F_S)])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')
colnames(PeakDistribution_F_all_M) = c('F_M_Nuc+Cyt')
colnames(PeakDistribution_F_all_S) = c('F_S_Nuc+Cyt')

## Tag Counts Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Nuc_S, 
                                    PeakDistribution_F_Cyt_M, PeakDistribution_F_Cyt_S, 
                                    PeakDistribution_F_all_M, PeakDistribution_F_all_S)
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress', 'Both Fractions Mock','Both Fractions Stress')) +
  ggtitle('All Normalized Tag Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## Fraction Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), 
                                    PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S), 
                                    PeakDistribution_F_all_M/sum(PeakDistribution_F_all_M), PeakDistribution_F_all_S/sum(PeakDistribution_F_all_S))
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress', 'Both Fractions Mock','Both Fractions Stress')) +
  ggtitle('All Normalized Tag Fractions Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

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
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  ggtitle('All Normalized Tag Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

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
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  ggtitle('All Normalized Tag Fractions Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")
####################################################################################################################

## Exploratory Stacked Bar Plots For Fractionation CLIP Peaks: Specific RNA biotypes ONLY
####################################################################################################################
## mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_F_Nuc_M_mRNA = Peak_F_Nuc_M %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')
Peak_F_Nuc_S_mRNA = Peak_F_Nuc_S %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')

Peak_F_Cyt_M_mRNA = Peak_F_Cyt_M %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')
Peak_F_Cyt_S_mRNA = Peak_F_Cyt_S %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')

Peak_F_all_M_mRNA = Peak_F_all_M %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')
Peak_F_all_S_mRNA = Peak_F_all_S %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')

PeakDistribution_F_Nuc_M = data.frame(table(Peak_F_Nuc_M_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_Nuc_S = data.frame(table(Peak_F_Nuc_S_mRNA$finalized_annotation), row.names = 1)

PeakDistribution_F_Cyt_M = data.frame(table(Peak_F_Cyt_M_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_Cyt_S = data.frame(table(Peak_F_Cyt_S_mRNA$finalized_annotation), row.names = 1)

PeakDistribution_F_all_M = data.frame(table(Peak_F_all_M_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_all_S = data.frame(table(Peak_F_all_S_mRNA$finalized_annotation), row.names = 1)

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')
colnames(PeakDistribution_F_all_M) = c('F_M_Nuc+Cyt')
colnames(PeakDistribution_F_all_S) = c('F_S_Nuc+Cyt')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Nuc_S, 
                                    PeakDistribution_F_Cyt_M, PeakDistribution_F_Cyt_S, 
                                    PeakDistribution_F_all_M, PeakDistribution_F_all_S)
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", "ncRNA_Retained_intron", 'downstream 10K'))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress', 'Both Fractions Mock','Both Fractions Stress')) +
  ggtitle('All Raw Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## Fraction Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), 
                                    PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S), 
                                    PeakDistribution_F_all_M/sum(PeakDistribution_F_all_M), PeakDistribution_F_all_S/sum(PeakDistribution_F_all_S))
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", "ncRNA_Retained_intron", 'downstream 10K'))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress', 'Both Fractions Mock','Both Fractions Stress')) +
  ggtitle('All Peak Fractions Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## non-mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_F_Nuc_M_Not_mRNA = Peak_F_Nuc_M %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != "downstream 10K")
Peak_F_Nuc_S_Not_mRNA = Peak_F_Nuc_S %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != "downstream 10K")

Peak_F_Cyt_M_Not_mRNA = Peak_F_Cyt_M %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != "downstream 10K")
Peak_F_Cyt_S_Not_mRNA = Peak_F_Cyt_S %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != "downstream 10K")

Peak_F_all_M_Not_mRNA = Peak_F_all_M %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != "downstream 10K")
Peak_F_all_S_Not_mRNA = Peak_F_all_S %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != "downstream 10K")

PeakDistribution_F_Nuc_M = data.frame(table(Peak_F_Nuc_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_Nuc_S = data.frame(table(Peak_F_Nuc_S_Not_mRNA$finalized_annotation), row.names = 1)

PeakDistribution_F_Cyt_M = data.frame(table(Peak_F_Cyt_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_Cyt_S = data.frame(table(Peak_F_Cyt_S_Not_mRNA$finalized_annotation), row.names = 1)

PeakDistribution_F_all_M = data.frame(table(Peak_F_all_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_F_all_S = data.frame(table(Peak_F_all_S_Not_mRNA$finalized_annotation), row.names = 1)

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')
colnames(PeakDistribution_F_all_M) = c('F_M_Nuc+Cyt')
colnames(PeakDistribution_F_all_S) = c('F_S_Nuc+Cyt')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Nuc_S, 
                                    PeakDistribution_F_Cyt_M, PeakDistribution_F_Cyt_S, 
                                    PeakDistribution_F_all_M, PeakDistribution_F_all_S)
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress', 'Both Fractions Mock','Both Fractions Stress')) +
  ggtitle('All Raw Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## Fraction Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), 
                                    PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S), 
                                    PeakDistribution_F_all_M/sum(PeakDistribution_F_all_M), PeakDistribution_F_all_S/sum(PeakDistribution_F_all_S))
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress', 'Both Fractions Mock','Both Fractions Stress')) +
  ggtitle('All Peak Fractions Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")
####################################################################################################################

## Exploratory Stacked Bar Plots For CoCLIP Peaks: Specific RNA biotypes ONLY
####################################################################################################################
## mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_Co_Input_M_mRNA = Peak_Co_Input_M %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')
Peak_Co_Input_S_mRNA = Peak_Co_Input_S %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')

Peak_Co_NLS_M_mRNA = Peak_Co_NLS_M %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')
Peak_Co_NLS_S_mRNA = Peak_Co_NLS_S %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')

Peak_Co_NES_M_mRNA = Peak_Co_NES_M %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')
Peak_Co_NES_S_mRNA = Peak_Co_NES_S %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')

Peak_Co_G3BP_M_mRNA = Peak_Co_G3BP_M %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')
Peak_Co_G3BP_S_mRNA = Peak_Co_G3BP_S %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')

PeakDistribution_Co_Input_M = data.frame(table(Peak_Co_Input_M_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_Input_S = data.frame(table(Peak_Co_Input_S_mRNA$finalized_annotation), row.names = 1)

PeakDistribution_Co_NLS_M = data.frame(table(Peak_Co_NLS_M_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NLS_S = data.frame(table(Peak_Co_NLS_S_mRNA$finalized_annotation), row.names = 1)

PeakDistribution_Co_NES_M = data.frame(table(Peak_Co_NES_M_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NES_M = rbind(PeakDistribution_Co_NES_M, c(0))
row.names(PeakDistribution_Co_NES_M) = c(row.names(PeakDistribution_Co_NES_M)[1:6], 'ncRNA_Retained_intron')
PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

PeakDistribution_Co_NES_S = data.frame(table(Peak_Co_NES_S_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NES_S = rbind(PeakDistribution_Co_NES_S, c(0))
row.names(PeakDistribution_Co_NES_S) = c(row.names(PeakDistribution_Co_NES_S)[1:6], 'ncRNA_Retained_intron')
PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

PeakDistribution_Co_G3BP_M = data.frame(table(Peak_Co_G3BP_M_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_G3BP_M = rbind(PeakDistribution_Co_G3BP_M, c(0))
row.names(PeakDistribution_Co_G3BP_M) = c(row.names(PeakDistribution_Co_G3BP_M)[1:6], 'ncRNA_Retained_intron')
PeakDistribution_Co_G3BP_M = PeakDistribution_Co_G3BP_M[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

PeakDistribution_Co_G3BP_S = data.frame(table(Peak_Co_G3BP_S_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_G3BP_S = rbind(PeakDistribution_Co_G3BP_S, c(0))
row.names(PeakDistribution_Co_G3BP_S) = c(row.names(PeakDistribution_Co_G3BP_S)[1:6], 'ncRNA_Retained_intron')
PeakDistribution_Co_G3BP_S = PeakDistribution_Co_G3BP_S[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

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
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", "ncRNA_Retained_intron", 'downstream 10K'))

ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  ggtitle('All Raw Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

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
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", "ncRNA_Retained_intron", 'downstream 10K'))

ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  ggtitle('All Peak Fractions Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## non-mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_Co_Input_M_Not_mRNA = Peak_Co_Input_M %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != 'downstream 10K')
Peak_Co_Input_S_Not_mRNA = Peak_Co_Input_S %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != 'downstream 10K')

Peak_Co_NLS_M_Not_mRNA = Peak_Co_NLS_M %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != 'downstream 10K')
Peak_Co_NLS_S_Not_mRNA = Peak_Co_NLS_S %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != 'downstream 10K')

Peak_Co_NES_M_Not_mRNA = Peak_Co_NES_M %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != 'downstream 10K')
Peak_Co_NES_S_Not_mRNA = Peak_Co_NES_S %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != 'downstream 10K')

Peak_Co_G3BP_M_Not_mRNA = Peak_Co_G3BP_M %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != 'downstream 10K')
Peak_Co_G3BP_S_Not_mRNA = Peak_Co_G3BP_S %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != 'downstream 10K')

PeakDistribution_Co_NLS_M = data.frame(table(Peak_Co_NLS_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NLS_S = data.frame(table(Peak_Co_NLS_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NLS_S = rbind(PeakDistribution_Co_NLS_S, c(0))
row.names(PeakDistribution_Co_NLS_S) = c(row.names(PeakDistribution_Co_NLS_S)[1:10], 'scaRNA')
PeakDistribution_Co_NLS_S = PeakDistribution_Co_NLS_S[row.names(PeakDistribution_Co_NLS_M), 'Freq', drop = FALSE]

PeakDistribution_Co_Input_M = data.frame(table(Peak_Co_Input_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_Input_M = rbind(PeakDistribution_Co_Input_M, c(0))
row.names(PeakDistribution_Co_Input_M) = c(row.names(PeakDistribution_Co_Input_M)[1:10], 'scaRNA')
PeakDistribution_Co_Input_M = PeakDistribution_Co_Input_M[row.names(PeakDistribution_Co_NLS_M), 'Freq', drop = FALSE]

PeakDistribution_Co_Input_S = data.frame(table(Peak_Co_Input_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_Input_S = rbind(PeakDistribution_Co_Input_S, c(0))
row.names(PeakDistribution_Co_Input_S) = c(row.names(PeakDistribution_Co_Input_S)[1:10], 'scaRNA')
PeakDistribution_Co_Input_S = PeakDistribution_Co_Input_S[row.names(PeakDistribution_Co_NLS_M), 'Freq', drop = FALSE]

PeakDistribution_Co_NES_M = data.frame(table(Peak_Co_NES_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NES_M = rbind(PeakDistribution_Co_NES_M, c(0), c(0), c(0))
row.names(PeakDistribution_Co_NES_M) = c(row.names(PeakDistribution_Co_NES_M)[1:8], 'miRNA', 'miscRNA', 'scaRNA')
PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

PeakDistribution_Co_NES_S = data.frame(table(Peak_Co_NES_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_NES_S = rbind(PeakDistribution_Co_NES_S, c(0))
row.names(PeakDistribution_Co_NES_S) = c(row.names(PeakDistribution_Co_NES_S)[1:10], 'scaRNA')
PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S[row.names(PeakDistribution_Co_NLS_M), 'Freq', drop = FALSE]

PeakDistribution_Co_G3BP_M = data.frame(table(Peak_Co_G3BP_M_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_G3BP_M = rbind(PeakDistribution_Co_G3BP_M, c(0), c(0), c(0))
row.names(PeakDistribution_Co_G3BP_M) = c(row.names(PeakDistribution_Co_G3BP_M)[1:8], 'miRNA', 'miscRNA', 'scaRNA')
PeakDistribution_Co_G3BP_M = PeakDistribution_Co_G3BP_M[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

PeakDistribution_Co_G3BP_S = data.frame(table(Peak_Co_G3BP_S_Not_mRNA$finalized_annotation), row.names = 1)
PeakDistribution_Co_G3BP_S = rbind(PeakDistribution_Co_G3BP_S, c(0), c(0), c(0))
row.names(PeakDistribution_Co_G3BP_S) = c(row.names(PeakDistribution_Co_G3BP_S)[1:8], 'miRNA', 'miscRNA', 'scaRNA')
PeakDistribution_Co_G3BP_S = PeakDistribution_Co_G3BP_S[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

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
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  ggtitle('All Raw Peak Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

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
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  ggtitle('All Peak Fractions Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")
####################################################################################################################

## Exploratory Stacked Bar Plots For Fractionation CLIP Tags: Specific RNA biotypes ONLY
####################################################################################################################
## mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_F_Nuc_M_mRNA = Peak_F_Nuc_M %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')
Peak_F_Nuc_S_mRNA = Peak_F_Nuc_S %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')

Peak_F_Cyt_M_mRNA = Peak_F_Cyt_M %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')
Peak_F_Cyt_S_mRNA = Peak_F_Cyt_S %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')

Peak_F_all_M_mRNA = Peak_F_all_M %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')
Peak_F_all_S_mRNA = Peak_F_all_S %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')

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

PeakDistribution_F_all_M = data.frame(data.frame(annotation = Peak_F_all_M_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_all_M_mRNA[, c(Nuc_F_M, Cyto_F_M)])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_all_S = data.frame(data.frame(annotation = Peak_F_all_S_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_all_S_mRNA[, c(Nuc_F_S, Cyto_F_S)])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')
colnames(PeakDistribution_F_all_M) = c('F_M_Nuc+Cyt')
colnames(PeakDistribution_F_all_S) = c('F_S_Nuc+Cyt')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Nuc_S, 
                                    PeakDistribution_F_Cyt_M, PeakDistribution_F_Cyt_S, 
                                    PeakDistribution_F_all_M, PeakDistribution_F_all_S)
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", "ncRNA_Retained_intron", 'downstream 10K'))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress', 'Both Fractions Mock','Both Fractions Stress')) +
  ggtitle('All Normalized Tag Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## Fraction Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), 
                                    PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S), 
                                    PeakDistribution_F_all_M/sum(PeakDistribution_F_all_M), PeakDistribution_F_all_S/sum(PeakDistribution_F_all_S))
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", "ncRNA_Retained_intron", 'downstream 10K'))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress', 'Both Fractions Mock','Both Fractions Stress')) +
  ggtitle('All Normalized Tag Fractions Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## non-mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_F_Nuc_M_Not_mRNA = Peak_F_Nuc_M %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != "downstream 10K")
Peak_F_Nuc_S_Not_mRNA = Peak_F_Nuc_S %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != "downstream 10K")

Peak_F_Cyt_M_Not_mRNA = Peak_F_Cyt_M %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != "downstream 10K")
Peak_F_Cyt_S_Not_mRNA = Peak_F_Cyt_S %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != "downstream 10K")

Peak_F_all_M_Not_mRNA = Peak_F_all_M %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != "downstream 10K")
Peak_F_all_S_Not_mRNA = Peak_F_all_S %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != "downstream 10K")

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

PeakDistribution_F_all_M = data.frame(data.frame(annotation = Peak_F_all_M_Not_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_all_M_Not_mRNA[, c(Nuc_F_M, Cyto_F_M)])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

PeakDistribution_F_all_S = data.frame(data.frame(annotation = Peak_F_all_S_Not_mRNA$finalized_annotation, 
                                                 tagSum = rowSums(Peak_F_all_S_Not_mRNA[, c(Nuc_F_S, Cyto_F_S)])) 
                                      %>% group_by(annotation) 
                                      %>% summarise(tags = sum(tagSum)), row.names = 1)

colnames(PeakDistribution_F_Nuc_M) = c('F_M_Nuc')
colnames(PeakDistribution_F_Nuc_S) = c('F_S_Nuc')
colnames(PeakDistribution_F_Cyt_M) = c('F_M_Cyt')
colnames(PeakDistribution_F_Cyt_S) = c('F_S_Cyt')
colnames(PeakDistribution_F_all_M) = c('F_M_Nuc+Cyt')
colnames(PeakDistribution_F_all_S) = c('F_S_Nuc+Cyt')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Nuc_S, 
                                    PeakDistribution_F_Cyt_M, PeakDistribution_F_Cyt_S, 
                                    PeakDistribution_F_all_M, PeakDistribution_F_all_S)
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress', 'Both Fractions Mock','Both Fractions Stress')) +
  ggtitle('All Normalized Tag Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## Fraction Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M/sum(PeakDistribution_F_Nuc_M), PeakDistribution_F_Nuc_S/sum(PeakDistribution_F_Nuc_S), 
                                    PeakDistribution_F_Cyt_M/sum(PeakDistribution_F_Cyt_M), PeakDistribution_F_Cyt_S/sum(PeakDistribution_F_Cyt_S), 
                                    PeakDistribution_F_all_M/sum(PeakDistribution_F_all_M), PeakDistribution_F_all_S/sum(PeakDistribution_F_all_S))
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress', 'Both Fractions Mock','Both Fractions Stress')) +
  ggtitle('All Normalized Tag Fractions Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")
####################################################################################################################

## Exploratory Stacked Bar Plots For CoCLIP Tags: Specific RNA biotypes ONLY
####################################################################################################################
## mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_Co_Input_M_mRNA = Peak_Co_Input_M %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')
Peak_Co_Input_S_mRNA = Peak_Co_Input_S %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')

Peak_Co_NLS_M_mRNA = Peak_Co_NLS_M %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')
Peak_Co_NLS_S_mRNA = Peak_Co_NLS_S %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')

Peak_Co_NES_M_mRNA = Peak_Co_NES_M %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')
Peak_Co_NES_S_mRNA = Peak_Co_NES_S %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')

Peak_Co_G3BP_M_mRNA = Peak_Co_G3BP_M %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')
Peak_Co_G3BP_S_mRNA = Peak_Co_G3BP_S %>% filter(grouped_annotation == "3'UTR" | grouped_annotation == "CDS" | grouped_annotation == "5'UTR" | grouped_annotation == "intron" | grouped_annotation == 'downstream 10K')

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
PeakDistribution_Co_NES_M = rbind(PeakDistribution_Co_NES_M, c(0))
row.names(PeakDistribution_Co_NES_M) = c(row.names(PeakDistribution_Co_NES_M)[1:6], 'ncRNA_Retained_intron')
PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M[row.names(PeakDistribution_Co_NLS_S), 'tags', drop = FALSE]

PeakDistribution_Co_NES_S = data.frame(data.frame(annotation = Peak_Co_NES_S_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NES_S_mRNA[, c(NES_E_S)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_NES_S = rbind(PeakDistribution_Co_NES_S, c(0))
row.names(PeakDistribution_Co_NES_S) = c(row.names(PeakDistribution_Co_NES_S)[1:6], 'ncRNA_Retained_intron')
PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S[row.names(PeakDistribution_Co_NLS_S), 'tags', drop = FALSE]

PeakDistribution_Co_G3BP_M = data.frame(data.frame(annotation = Peak_Co_G3BP_M_mRNA$finalized_annotation, 
                                                   tagSum = rowSums(Peak_Co_G3BP_M_mRNA[, c(G3BP_E_M)])) 
                                        %>% group_by(annotation) 
                                        %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_G3BP_M = rbind(PeakDistribution_Co_G3BP_M, c(0))
row.names(PeakDistribution_Co_G3BP_M) = c(row.names(PeakDistribution_Co_G3BP_M)[1:6], 'ncRNA_Retained_intron')
PeakDistribution_Co_G3BP_M = PeakDistribution_Co_G3BP_M[row.names(PeakDistribution_Co_NLS_S), 'tags', drop = FALSE]

PeakDistribution_Co_G3BP_S = data.frame(data.frame(annotation = Peak_Co_G3BP_S_mRNA$finalized_annotation, 
                                                   tagSum = rowSums(Peak_Co_G3BP_S_mRNA[, c(G3BP_E_S)])) 
                                        %>% group_by(annotation) 
                                        %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_G3BP_S = rbind(PeakDistribution_Co_G3BP_S, c(0))
row.names(PeakDistribution_Co_G3BP_S) = c(row.names(PeakDistribution_Co_G3BP_S)[1:6], 'ncRNA_Retained_intron')
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
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", "ncRNA_Retained_intron", 'downstream 10K'))

ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  ggtitle('All Normalized Tag Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

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
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", "ncRNA_Retained_intron", 'downstream 10K'))

ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  ggtitle('All Normalized Tag Fractions Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

## non-mRNA Features Only
## Mock vs Stress for Each Fraction
Peak_Co_Input_M_Not_mRNA = Peak_Co_Input_M %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != 'downstream 10K')
Peak_Co_Input_S_Not_mRNA = Peak_Co_Input_S %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != 'downstream 10K')

Peak_Co_NLS_M_Not_mRNA = Peak_Co_NLS_M %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != 'downstream 10K')
Peak_Co_NLS_S_Not_mRNA = Peak_Co_NLS_S %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != 'downstream 10K')

Peak_Co_NES_M_Not_mRNA = Peak_Co_NES_M %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != 'downstream 10K')
Peak_Co_NES_S_Not_mRNA = Peak_Co_NES_S %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != 'downstream 10K')

Peak_Co_G3BP_M_Not_mRNA = Peak_Co_G3BP_M %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != 'downstream 10K')
Peak_Co_G3BP_S_Not_mRNA = Peak_Co_G3BP_S %>% filter(grouped_annotation != "3'UTR" & grouped_annotation != "CDS" & grouped_annotation != "5'UTR" & grouped_annotation != "intron" & grouped_annotation != 'downstream 10K')

PeakDistribution_Co_NLS_M = data.frame(data.frame(annotation = Peak_Co_NLS_M_Not_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NLS_M_Not_mRNA[, c(NLS_E_M)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)

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
PeakDistribution_Co_Input_M = rbind(PeakDistribution_Co_Input_M, c(0))
row.names(PeakDistribution_Co_Input_M) = c(row.names(PeakDistribution_Co_Input_M)[1:10], 'scaRNA')
PeakDistribution_Co_Input_M = PeakDistribution_Co_Input_M[row.names(PeakDistribution_Co_NLS_M), 'tags', drop = FALSE]

PeakDistribution_Co_Input_S = data.frame(data.frame(annotation = Peak_Co_Input_S_Not_mRNA$finalized_annotation, 
                                                    tagSum = rowSums(Peak_Co_Input_S_Not_mRNA[, c(NLS_I_S, NES_I_S, G3BP_I_S)])) 
                                         %>% group_by(annotation) 
                                         %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_Input_S = rbind(PeakDistribution_Co_Input_S, c(0))
row.names(PeakDistribution_Co_Input_S) = c(row.names(PeakDistribution_Co_Input_S)[1:10], 'scaRNA')
PeakDistribution_Co_Input_S = PeakDistribution_Co_Input_S[row.names(PeakDistribution_Co_NLS_M), 'tags', drop = FALSE]

PeakDistribution_Co_NES_M = data.frame(data.frame(annotation = Peak_Co_NES_M_Not_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NES_M_Not_mRNA[, c(NES_E_M)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_NES_M = rbind(PeakDistribution_Co_NES_M, c(0), c(0), c(0))
row.names(PeakDistribution_Co_NES_M) = c(row.names(PeakDistribution_Co_NES_M)[1:8], 'miRNA', 'miscRNA', 'scaRNA')
PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M[row.names(PeakDistribution_Co_NLS_S), 'tags', drop = FALSE]

PeakDistribution_Co_NES_S = data.frame(data.frame(annotation = Peak_Co_NES_S_Not_mRNA$finalized_annotation, 
                                                  tagSum = rowSums(Peak_Co_NES_S_Not_mRNA[, c(NES_E_S)])) 
                                       %>% group_by(annotation) 
                                       %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_NES_S = rbind(PeakDistribution_Co_NES_S, c(0))
row.names(PeakDistribution_Co_NES_S) = c(row.names(PeakDistribution_Co_NES_S)[1:10], 'scaRNA')
PeakDistribution_Co_NES_S = PeakDistribution_Co_NES_S[row.names(PeakDistribution_Co_NLS_M), 'tags', drop = FALSE]

PeakDistribution_Co_G3BP_M = data.frame(data.frame(annotation = Peak_Co_G3BP_M_Not_mRNA$finalized_annotation, 
                                                   tagSum = rowSums(Peak_Co_G3BP_M_Not_mRNA[, c(G3BP_E_M)])) 
                                        %>% group_by(annotation) 
                                        %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_G3BP_M = rbind(PeakDistribution_Co_G3BP_M, c(0), c(0), c(0))
row.names(PeakDistribution_Co_G3BP_M) = c(row.names(PeakDistribution_Co_G3BP_M)[1:8], 'miRNA', 'miscRNA', 'scaRNA')
PeakDistribution_Co_G3BP_M = PeakDistribution_Co_G3BP_M[row.names(PeakDistribution_Co_NLS_S), 'tags', drop = FALSE]

PeakDistribution_Co_G3BP_S = data.frame(data.frame(annotation = Peak_Co_G3BP_S_Not_mRNA$finalized_annotation, 
                                                   tagSum = rowSums(Peak_Co_G3BP_S_Not_mRNA[, c(G3BP_E_S)])) 
                                        %>% group_by(annotation) 
                                        %>% summarise(tags = sum(tagSum)), row.names = 1)
PeakDistribution_Co_G3BP_S = rbind(PeakDistribution_Co_G3BP_S, c(0), c(0), c(0))
row.names(PeakDistribution_Co_G3BP_S) = c(row.names(PeakDistribution_Co_G3BP_S)[1:8], 'miRNA', 'miscRNA', 'scaRNA')
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
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  ggtitle('All Normalized Tag Counts Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

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
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  ggtitle('All Normalized Tag Fractions Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")
####################################################################################################################

## Start Building Enrichment Table:
####################################################################################################################
peakRowSum = peakMatrix[, c(inert_columns, rowSum_columns)]

peakRowSum = peakRowSum %>% mutate(Nuc_F_M = rowSums(peakMatrix[, Nuc_F_M])/length(Nuc_F_M) * 1e6)
peakRowSum = peakRowSum %>% mutate(Nuc_F_S = rowSums(peakMatrix[, Nuc_F_S])/length(Nuc_F_S) * 1e6)
peakRowSum = peakRowSum %>% mutate(Cyto_F_M = rowSums(peakMatrix[, Cyto_F_M])/length(Cyto_F_M) * 1e6)
peakRowSum = peakRowSum %>% mutate(Cyto_F_S = rowSums(peakMatrix[, Cyto_F_S])/length(Cyto_F_S) * 1e6)

peakRowSum = peakRowSum %>% mutate(NLS_I_M = rowSums(peakMatrix[, NLS_I_M])/length(NLS_I_M) * 1e6)
peakRowSum = peakRowSum %>% mutate(NES_I_M = rowSums(peakMatrix[, NES_I_M])/length(NES_I_M) * 1e6)
peakRowSum = peakRowSum %>% mutate(G3BP_I_M = rowSums(peakMatrix[, G3BP_I_M])/length(G3BP_I_M) * 1e6)
temp = rowSums(peakMatrix[, c(NLS_I_M, NES_I_M, G3BP_I_M)])/length(c(NLS_I_M, NES_I_M, G3BP_I_M)) * 1e6
peakRowSum = peakRowSum %>% mutate(I_M = temp)

peakRowSum = peakRowSum %>% mutate(NLS_I_S = rowSums(peakMatrix[, NLS_I_S])/length(NLS_I_S) * 1e6)
peakRowSum = peakRowSum %>% mutate(NES_I_S = rowSums(peakMatrix[, NES_I_S])/length(NES_I_S) * 1e6)
peakRowSum = peakRowSum %>% mutate(G3BP_I_S = rowSums(peakMatrix[, G3BP_I_S])/length(G3BP_I_S) * 1e6)
temp = rowSums(peakMatrix[, c(NLS_I_S, NES_I_S, G3BP_I_S)])/length(c(NLS_I_S, NES_I_S, G3BP_I_S)) * 1e6
peakRowSum = peakRowSum %>% mutate(I_S = temp)

peakRowSum = peakRowSum %>% mutate(NLS_E_M = rowSums(peakMatrix[, NLS_E_M])/length(NLS_E_M) * 1e6)
peakRowSum = peakRowSum %>% mutate(NLS_E_S = rowSums(peakMatrix[, NLS_E_S])/length(NLS_E_S) * 1e6)
peakRowSum = peakRowSum %>% mutate(NES_E_M = rowSums(peakMatrix[, NES_E_M])/length(NES_E_M) * 1e6)
peakRowSum = peakRowSum %>% mutate(NES_E_S = rowSums(peakMatrix[, NES_E_S])/length(NES_E_S) * 1e6)
peakRowSum = peakRowSum %>% mutate(G3BP_E_M = rowSums(peakMatrix[, G3BP_E_M])/length(G3BP_E_M) * 1e6)
peakRowSum = peakRowSum %>% mutate(G3BP_E_S = rowSums(peakMatrix[, G3BP_E_S])/length(G3BP_E_S) * 1e6)

peakRowSum = cbind(peakRowSum, peakMatrix[, BC_columns])

# write.table(peakRowSum, paste0(peaksMatrix_PATH, str_replace(peaksMatrix_FILE, ".txt", "_rowSum.txt")), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t', na = "")

peakEnrichment = peakMatrix[, c(inert_columns, BC_columns, rowSum_columns)]

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

# write.table(peakEnrichment, paste0(peaksMatrix_PATH, str_replace(peaksMatrix_FILE, ".txt", "_Enrichment.txt")), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t', na = "")
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
  ggtitle('snoRNA HuR Peaks: CoCLIP vs Fractionation CLIP of Nuclear/Cytoplasm in Mock')


snoRNA_E_S_NvC = snoRNA_Enrichment$NLS_E_S / snoRNA_Enrichment$NES_E_S
snoRNA_F_S_NvC = snoRNA_Enrichment$Nuc_F_S / snoRNA_Enrichment$Cyto_F_S
data = data.frame(log_snoRNA_E_S_NvC = log2(snoRNA_E_S_NvC), log_snoRNA_F_S_NvC = log2(snoRNA_F_S_NvC))

ggplot(data, aes(x = log_snoRNA_E_S_NvC, y = log_snoRNA_F_S_NvC)) +
  geom_point(col = 'black', pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP NLS/NES Stress)', y = 'log2(FractionCLIP Nuc/Cyto Stress)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle('snoRNA HuR Peaks: CoCLIP vs Fractionation CLIP of Nuclear/Cytoplasm in Stress')

####################################################################################################################

## Mock Vs Stress
####################################################################################################################
## Nuclear: 
NLS_Peaks_Filtered = peakEnrichment %>% filter(((NLS_I_M_BC >= BC_Threshold_Input & 
                                                  NLS_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                 (NLS_I_S_BC >= BC_Threshold_Input & 
                                                 NLS_E_S_BC >= BC_Threshold_CoCLIP)) &
                                                 I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                 E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_Peaks_Filtered$NLS_EvI_M, NLS_EvI_S = NLS_Peaks_Filtered$NLS_EvI_S)
data$annotation = factor(NLS_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

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
NLS_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation == "5'UTR" | 
                                                  grouped_annotation == "3'UTR" | 
                                                  grouped_annotation == "CDS" | 
                                                  grouped_annotation == "intron" |
                                                  grouped_annotation == "downstream 10K") & 
                                                 ((NLS_I_M_BC >= BC_Threshold_Input & 
                                                   NLS_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                  (NLS_I_S_BC >= BC_Threshold_Input & 
                                                     NLS_E_S_BC >= BC_Threshold_CoCLIP)) &
                                                 I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                 E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_Peaks_Filtered$NLS_EvI_M, NLS_EvI_S = NLS_Peaks_Filtered$NLS_EvI_S)
data$annotation = factor(NLS_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

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
NLS_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != "5'UTR" & 
                                                  grouped_annotation != "3'UTR" & 
                                                  grouped_annotation != "CDS" & 
                                                  grouped_annotation != "intron" & 
                                                  grouped_annotation != "downstream 10K") & 
                                                 ((NLS_I_M_BC >= BC_Threshold_Input & 
                                                     NLS_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                    (NLS_I_S_BC >= BC_Threshold_Input & 
                                                       NLS_E_S_BC >= BC_Threshold_CoCLIP)) &
                                                 I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                 E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_Peaks_Filtered$NLS_EvI_M, NLS_EvI_S = NLS_Peaks_Filtered$NLS_EvI_S)
data$annotation = factor(NLS_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

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
NES_Peaks_Filtered = peakEnrichment %>% filter(((NES_I_M_BC >= BC_Threshold_Input & 
                                                  NES_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                 (NES_I_S_BC >= BC_Threshold_Input & 
                                                 NES_E_S_BC >= BC_Threshold_CoCLIP)) &
                                                 I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                 E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_M = NES_Peaks_Filtered$NES_EvI_M, NES_EvI_S = NES_Peaks_Filtered$NES_EvI_S)
data$annotation = factor(NES_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

ggplot(data, aes(x = log2(NES_EvI_M), y = log2(NES_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3) +
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
NES_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation == "5'UTR" | 
                                                  grouped_annotation == "3'UTR" | 
                                                  grouped_annotation == "CDS" | 
                                                  grouped_annotation == "intron" |
                                                  grouped_annotation == "downstream 10K") &
                                                 ((NES_I_M_BC >= BC_Threshold_Input & 
                                                   NES_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                  (NES_I_S_BC >= BC_Threshold_Input & 
                                                     NES_E_S_BC >= BC_Threshold_CoCLIP)) &
                                                 I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                 E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_M = NES_Peaks_Filtered$NES_EvI_M, NES_EvI_S = NES_Peaks_Filtered$NES_EvI_S)
data$annotation = factor(NES_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

ggplot(data, aes(x = log2(NES_EvI_M), y = log2(NES_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3) +
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
NES_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != "5'UTR" & 
                                                  grouped_annotation != "3'UTR" & 
                                                  grouped_annotation != "CDS" & 
                                                  grouped_annotation != "intron" & 
                                                  grouped_annotation != "downstream 10K") &
                                                 ((NES_I_M_BC >= BC_Threshold_Input & 
                                                     NES_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                    (NES_I_S_BC >= BC_Threshold_Input & 
                                                       NES_E_S_BC >= BC_Threshold_CoCLIP)) &
                                                 I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                 E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_M = NES_Peaks_Filtered$NES_EvI_M, NES_EvI_S = NES_Peaks_Filtered$NES_EvI_S)
data$annotation = factor(NES_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

ggplot(data, aes(x = log2(NES_EvI_M), y = log2(NES_EvI_S), color = annotation)) +
  geom_point(pch = 16, size = 3) +
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
G3BP_Peaks_Filtered = peakEnrichment %>% filter(((G3BP_I_M_BC >= BC_Threshold_Input_SG & 
                                                   G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG) | 
                                                  (G3BP_I_S_BC >= BC_Threshold_Input_SG & 
                                                  G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG)) &
                                                  I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                  E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(G3BP_EvI_M = G3BP_Peaks_Filtered$G3BP_EvI_M, G3BP_EvI_S = G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(G3BP_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

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
G3BP_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation == "5'UTR" | 
                                                   grouped_annotation == "3'UTR" | 
                                                   grouped_annotation == "CDS" | 
                                                   grouped_annotation == "intron" |
                                                   grouped_annotation == "downstream 10K") &
                                                  ((G3BP_I_M_BC >= BC_Threshold_Input_SG & 
                                                    G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG) | 
                                                   (G3BP_I_S_BC >= BC_Threshold_Input_SG & 
                                                      G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG)) &
                                                  I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                  E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(G3BP_EvI_M = G3BP_Peaks_Filtered$G3BP_EvI_M, G3BP_EvI_S = G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(G3BP_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

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
G3BP_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != "5'UTR" & 
                                                   grouped_annotation != "3'UTR" & 
                                                   grouped_annotation != "CDS" & 
                                                   grouped_annotation != "intron" & 
                                                   grouped_annotation != "downstream 10K") &
                                                  ((G3BP_I_M_BC >= BC_Threshold_Input_SG & 
                                                    G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG) | 
                                                   (G3BP_I_S_BC >= BC_Threshold_Input_SG & 
                                                      G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG)) &
                                                  I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                  E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(G3BP_EvI_M = G3BP_Peaks_Filtered$G3BP_EvI_M, G3BP_EvI_S = G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(G3BP_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

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
NLS_NES_Peaks_Filtered = peakEnrichment %>% filter(((NLS_I_M_BC >= BC_Threshold_Input & 
                                                      NLS_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                     (NES_I_M_BC >= BC_Threshold_Input & 
                                                     NES_E_M_BC >= BC_Threshold_CoCLIP)) &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                     E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_NES_Peaks_Filtered$NLS_EvI_M, NES_EvI_M = NLS_NES_Peaks_Filtered$NES_EvI_M)
data$annotation = factor(NLS_NES_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

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
NLS_NES_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation == "5'UTR" | 
                                                      grouped_annotation == "3'UTR" | 
                                                      grouped_annotation == "CDS" | 
                                                      grouped_annotation == "intron" |
                                                      grouped_annotation == "downstream 10K") &
                                                     ((NLS_I_M_BC >= BC_Threshold_Input & 
                                                       NLS_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                      (NES_I_M_BC >= BC_Threshold_Input & 
                                                         NES_E_M_BC >= BC_Threshold_CoCLIP)) &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                     E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_NES_Peaks_Filtered$NLS_EvI_M, NES_EvI_M = NLS_NES_Peaks_Filtered$NES_EvI_M)
data$annotation = factor(NLS_NES_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

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
NLS_NES_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != "5'UTR" & 
                                                      grouped_annotation != "3'UTR" & 
                                                      grouped_annotation != "CDS" & 
                                                      grouped_annotation != "intron" & 
                                                      grouped_annotation != "downstream 10K") &
                                                     ((NLS_I_M_BC >= BC_Threshold_Input & 
                                                         NLS_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                        (NES_I_M_BC >= BC_Threshold_Input & 
                                                           NES_E_M_BC >= BC_Threshold_CoCLIP)) &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                     E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_NES_Peaks_Filtered$NLS_EvI_M, NES_EvI_M = NLS_NES_Peaks_Filtered$NES_EvI_M)
data$annotation = factor(NLS_NES_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

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
NLS_G3BP_Peaks_Filtered = peakEnrichment %>% filter(((NLS_I_M_BC >= BC_Threshold_Input & 
                                                       NLS_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                      (G3BP_I_M_BC >= BC_Threshold_Input_SG & 
                                                      G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_G3BP_Peaks_Filtered$NLS_EvI_M, G3BP_EvI_M = NLS_G3BP_Peaks_Filtered$G3BP_EvI_M)
data$annotation = factor(NLS_G3BP_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

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
NLS_G3BP_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation == "5'UTR" | 
                                                       grouped_annotation == "3'UTR" | 
                                                       grouped_annotation == "CDS" | 
                                                       grouped_annotation == "intron" |
                                                       grouped_annotation == "downstream 10K") &
                                                      ((NLS_I_M_BC >= BC_Threshold_Input & 
                                                        NLS_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                       (G3BP_I_M_BC >= BC_Threshold_Input_SG & 
                                                          G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_G3BP_Peaks_Filtered$NLS_EvI_M, G3BP_EvI_M = NLS_G3BP_Peaks_Filtered$G3BP_EvI_M)
data$annotation = factor(NLS_G3BP_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

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
NLS_G3BP_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != "5'UTR" & 
                                                       grouped_annotation != "3'UTR" & 
                                                       grouped_annotation != "CDS" & 
                                                       grouped_annotation != "intron" & 
                                                       grouped_annotation != "downstream 10K") &
                                                      ((NLS_I_M_BC >= BC_Threshold_Input & 
                                                          NLS_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                         (G3BP_I_M_BC >= BC_Threshold_Input_SG & 
                                                            G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_M = NLS_G3BP_Peaks_Filtered$NLS_EvI_M, G3BP_EvI_M = NLS_G3BP_Peaks_Filtered$G3BP_EvI_M)
data$annotation = factor(NLS_G3BP_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

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
NES_G3BP_Peaks_Filtered = peakEnrichment %>% filter(((NES_I_M_BC >= BC_Threshold_Input & 
                                                       NES_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                      (G3BP_I_M_BC >= BC_Threshold_Input_SG & 
                                                      G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_M = NES_G3BP_Peaks_Filtered$NES_EvI_M, G3BP_EvI_M = NES_G3BP_Peaks_Filtered$G3BP_EvI_M)
data$annotation = factor(NES_G3BP_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

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
NES_G3BP_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation == "5'UTR" | 
                                                       grouped_annotation == "3'UTR" | 
                                                       grouped_annotation == "CDS" | 
                                                       grouped_annotation == "intron" |
                                                       grouped_annotation == "downstream 10K") &
                                                      ((NES_I_M_BC >= BC_Threshold_Input & 
                                                        NES_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                       (G3BP_I_M_BC >= BC_Threshold_Input_SG & 
                                                          G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_M = NES_G3BP_Peaks_Filtered$NES_EvI_M, G3BP_EvI_M = NES_G3BP_Peaks_Filtered$G3BP_EvI_M)
data$annotation = factor(NES_G3BP_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

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
NES_G3BP_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != "5'UTR" & 
                                                       grouped_annotation != "3'UTR" & 
                                                       grouped_annotation != "CDS" & 
                                                       grouped_annotation != "intron" & 
                                                       grouped_annotation != "downstream 10K") &
                                                      ((NES_I_M_BC >= BC_Threshold_Input & 
                                                          NES_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                         (G3BP_I_M_BC >= BC_Threshold_Input_SG & 
                                                            G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_M = NES_G3BP_Peaks_Filtered$NES_EvI_M, G3BP_EvI_M = NES_G3BP_Peaks_Filtered$G3BP_EvI_M)
data$annotation = factor(NES_G3BP_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

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
NLS_NES_Peaks_Filtered = peakEnrichment %>% filter(((NLS_I_S_BC >= BC_Threshold_Input & 
                                                      NLS_E_S_BC >= BC_Threshold_CoCLIP) | 
                                                     (NES_I_S_BC >= BC_Threshold_Input & 
                                                     NES_E_S_BC >= BC_Threshold_CoCLIP)) &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                     E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_S = NLS_NES_Peaks_Filtered$NLS_EvI_S, NES_EvI_S = NLS_NES_Peaks_Filtered$NES_EvI_S)
data$annotation = factor(NLS_NES_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

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
NLS_NES_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation == "5'UTR" | 
                                                      grouped_annotation == "3'UTR" | 
                                                      grouped_annotation == "CDS" | 
                                                      grouped_annotation == "intron" |
                                                      grouped_annotation == "downstream 10K") &
                                                     ((NLS_I_S_BC >= BC_Threshold_Input & 
                                                       NLS_E_S_BC >= BC_Threshold_CoCLIP) | 
                                                      (NES_I_S_BC >= BC_Threshold_Input & 
                                                         NES_E_S_BC >= BC_Threshold_CoCLIP)) &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                     E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_S = NLS_NES_Peaks_Filtered$NLS_EvI_S, NES_EvI_S = NLS_NES_Peaks_Filtered$NES_EvI_S)
data$annotation = factor(NLS_NES_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

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
NLS_NES_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != "5'UTR" & 
                                                      grouped_annotation != "3'UTR" & 
                                                      grouped_annotation != "CDS" & 
                                                      grouped_annotation != "intron" & 
                                                      grouped_annotation != "downstream 10K") &
                                                     ((NLS_I_S_BC >= BC_Threshold_Input & 
                                                         NLS_E_S_BC >= BC_Threshold_CoCLIP) | 
                                                        (NES_I_S_BC >= BC_Threshold_Input & 
                                                           NES_E_S_BC >= BC_Threshold_CoCLIP)) &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                     E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_S = NLS_NES_Peaks_Filtered$NLS_EvI_S, NES_EvI_S = NLS_NES_Peaks_Filtered$NES_EvI_S)
data$annotation = factor(NLS_NES_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

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
NLS_G3BP_Peaks_Filtered = peakEnrichment %>% filter(((NLS_I_S_BC >= BC_Threshold_Input & 
                                                       NLS_E_S_BC >= BC_Threshold_CoCLIP) | 
                                                      (G3BP_I_S_BC >= BC_Threshold_Input_SG & 
                                                      G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_S = NLS_G3BP_Peaks_Filtered$NLS_EvI_S, G3BP_EvI_S = NLS_G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(NLS_G3BP_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

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
NLS_G3BP_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation == "5'UTR" | 
                                                       grouped_annotation == "3'UTR" | 
                                                       grouped_annotation == "CDS" | 
                                                       grouped_annotation == "intron" |
                                                       grouped_annotation == "downstream 10K") &
                                                      ((NLS_I_S_BC >= BC_Threshold_Input & 
                                                        NLS_E_S_BC >= BC_Threshold_CoCLIP) | 
                                                       (G3BP_I_S_BC >= BC_Threshold_Input_SG & 
                                                          G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_S = NLS_G3BP_Peaks_Filtered$NLS_EvI_S, G3BP_EvI_S = NLS_G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(NLS_G3BP_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

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
NLS_G3BP_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != "5'UTR" & 
                                                       grouped_annotation != "3'UTR" & 
                                                       grouped_annotation != "CDS" & 
                                                       grouped_annotation != "intron" & 
                                                       grouped_annotation != "downstream 10K") &
                                                      ((NLS_I_S_BC >= BC_Threshold_Input & 
                                                          NLS_E_S_BC >= BC_Threshold_CoCLIP) | 
                                                         (G3BP_I_S_BC >= BC_Threshold_Input_SG & 
                                                            G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NLS_EvI_S = NLS_G3BP_Peaks_Filtered$NLS_EvI_S, G3BP_EvI_S = NLS_G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(NLS_G3BP_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

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
NES_G3BP_Peaks_Filtered = peakEnrichment %>% filter(((NES_I_S_BC >= BC_Threshold_Input & 
                                                       NES_E_S_BC >= BC_Threshold_CoCLIP) | 
                                                      (G3BP_I_S_BC >= BC_Threshold_Input_SG & 
                                                      G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_S = NES_G3BP_Peaks_Filtered$NES_EvI_S, G3BP_EvI_S = NES_G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(NES_G3BP_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

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
NES_G3BP_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation == "5'UTR" | 
                                                       grouped_annotation == "3'UTR" | 
                                                       grouped_annotation == "CDS" | 
                                                       grouped_annotation == "intron" |
                                                       grouped_annotation == "downstream 10K") &
                                                      ((NES_I_S_BC >= BC_Threshold_Input & 
                                                        NES_E_S_BC >= BC_Threshold_CoCLIP) | 
                                                       (G3BP_I_S_BC >= BC_Threshold_Input_SG & 
                                                          G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_S = NES_G3BP_Peaks_Filtered$NES_EvI_S, G3BP_EvI_S = NES_G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(NES_G3BP_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

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
NES_G3BP_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != "5'UTR" & 
                                                       grouped_annotation != "3'UTR" & 
                                                       grouped_annotation != "CDS" & 
                                                       grouped_annotation != "intron" & 
                                                       grouped_annotation != "downstream 10K") &
                                                      ((NES_I_S_BC >= BC_Threshold_Input & 
                                                          NES_E_S_BC >= BC_Threshold_CoCLIP) | 
                                                         (G3BP_I_S_BC >= BC_Threshold_Input_SG & 
                                                            G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG)) &
                                                      I_rowSum >= median(peakEnrichment$I_rowSum)*2 &
                                                      E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(NES_EvI_S = NES_G3BP_Peaks_Filtered$NES_EvI_S, G3BP_EvI_S = NES_G3BP_Peaks_Filtered$G3BP_EvI_S)
data$annotation = factor(NES_G3BP_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

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
Mock_Peaks_Filtered = peakEnrichment %>% filter(((Cyto_F_M_BC >= BC_Threshold_Fraction & 
                                                  NES_I_M_BC >= BC_Threshold_Input) | 
                                                 (Nuc_F_M_BC >= BC_Threshold_Fraction & 
                                                  NLS_I_M_BC >= BC_Threshold_Input)) &
                                                 (F_rowSum >= median(peakEnrichment$F_rowSum)*2 &
                                                  I_rowSum >= median(peakEnrichment$I_rowSum)*2))

data = data.frame(Nuc_EvI_M = Mock_Peaks_Filtered$Nuc_EvI_M, Cyto_EvI_M = Mock_Peaks_Filtered$Cyto_EvI_M)
data$annotation = factor(Mock_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

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
Mock_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation == "5'UTR" | 
                                                   grouped_annotation == "3'UTR" | 
                                                   grouped_annotation == "CDS" | 
                                                   grouped_annotation == "intron" |
                                                   grouped_annotation == "downstream 10K") &
                                                  ((Cyto_F_M_BC >= BC_Threshold_Fraction & 
                                                    NES_I_M_BC >= BC_Threshold_Input) | 
                                                   (Nuc_F_M_BC >= BC_Threshold_Fraction & 
                                                      NLS_I_M_BC >= BC_Threshold_Input)) &
                                                  (F_rowSum >= median(peakEnrichment$F_rowSum)*2 &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*2))

data = data.frame(Nuc_EvI_M = Mock_Peaks_Filtered$Nuc_EvI_M, Cyto_EvI_M = Mock_Peaks_Filtered$Cyto_EvI_M)
data$annotation = factor(Mock_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

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
Mock_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != "5'UTR" & 
                                                   grouped_annotation != "3'UTR" & 
                                                   grouped_annotation != "CDS" & 
                                                   grouped_annotation != "intron" & 
                                                   grouped_annotation != "downstream 10K") &
                                                  ((Cyto_F_M_BC >= BC_Threshold_Fraction & 
                                                      NES_I_M_BC >= BC_Threshold_Input) | 
                                                     (Nuc_F_M_BC >= BC_Threshold_Fraction & 
                                                        NLS_I_M_BC >= BC_Threshold_Input)) &
                                                  (F_rowSum >= median(peakEnrichment$F_rowSum)*2 &
                                                     I_rowSum >= median(peakEnrichment$I_rowSum)*2))

data = data.frame(Nuc_EvI_M = Mock_Peaks_Filtered$Nuc_EvI_M, Cyto_EvI_M = Mock_Peaks_Filtered$Cyto_EvI_M)
data$annotation = factor(Mock_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

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
Stress_Peaks_Filtered = peakEnrichment %>% filter(((Cyto_F_S_BC >= BC_Threshold_Fraction & 
                                                    NES_I_S_BC >= BC_Threshold_Input) | 
                                                   (Nuc_F_S_BC >= BC_Threshold_Fraction & 
                                                    NLS_I_S_BC >= BC_Threshold_Input)) &
                                                   (F_rowSum >= median(peakEnrichment$F_rowSum)*2 &
                                                    I_rowSum >= median(peakEnrichment$I_rowSum)*2))

data = data.frame(Nuc_EvI_S = Stress_Peaks_Filtered$Nuc_EvI_S, Cyto_EvI_S = Stress_Peaks_Filtered$Cyto_EvI_S)
data$annotation = factor(Stress_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K")))

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
Stress_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation == "5'UTR" | 
                                                     grouped_annotation == "3'UTR" | 
                                                     grouped_annotation == "CDS" | 
                                                     grouped_annotation == "intron" |
                                                     grouped_annotation == "downstream 10K") &
                                                    ((Cyto_F_S_BC >= BC_Threshold_Fraction & 
                                                      NES_I_S_BC >= BC_Threshold_Input) | 
                                                     (Nuc_F_S_BC >= BC_Threshold_Fraction & 
                                                        NLS_I_S_BC >= BC_Threshold_Input)) &
                                                    (F_rowSum >= median(peakEnrichment$F_rowSum)*2 &
                                                       I_rowSum >= median(peakEnrichment$I_rowSum)*2))

data = data.frame(Nuc_EvI_S = Stress_Peaks_Filtered$Nuc_EvI_S, Cyto_EvI_S = Stress_Peaks_Filtered$Cyto_EvI_S)
data$annotation = factor(Stress_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

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
Stress_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != "5'UTR" & 
                                                     grouped_annotation != "3'UTR" & 
                                                     grouped_annotation != "CDS" & 
                                                     grouped_annotation != "intron" & 
                                                     grouped_annotation != "downstream 10K") &
                                                    ((Cyto_F_S_BC >= BC_Threshold_Fraction & 
                                                        NES_I_S_BC >= BC_Threshold_Input) | 
                                                       (Nuc_F_S_BC >= BC_Threshold_Fraction & 
                                                          NLS_I_S_BC >= BC_Threshold_Input)) &
                                                    (F_rowSum >= median(peakEnrichment$F_rowSum)*2 &
                                                       I_rowSum >= median(peakEnrichment$I_rowSum)*2))

data = data.frame(Nuc_EvI_S = Stress_Peaks_Filtered$Nuc_EvI_S, Cyto_EvI_S = Stress_Peaks_Filtered$Cyto_EvI_S)
data$annotation = factor(Stress_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

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

## Fractionation and CoCLIP comparison
####################################################################################################################
# Mock Comparison
Mock_Peaks_Filtered = peakEnrichment %>% filter(((NLS_E_M_BC >= BC_Threshold_CoCLIP & 
                                                  NES_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                 (Nuc_F_M_BC >= BC_Threshold_Fraction & 
                                                  Cyto_F_M_BC >= BC_Threshold_Fraction)) &
                                                  F_rowSum >= median(peakEnrichment$F_rowSum)*2 &
                                                  E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(E_NvC_M = Mock_Peaks_Filtered$E_NvC_M, F_NvC_M = Mock_Peaks_Filtered$F_NvC_M)
data$annotation = factor(Mock_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

ggplot(data, aes(x = log2(E_NvC_M), y = log2(F_NvC_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Cytoplasm)', y = 'log2(Fractionation CLIP Nuclear/Cytoplasm)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: CoCLIP vs Fractionation CLIP of Nuclear/Cytoplasm (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Mock Comparison - mRNA specific:
Mock_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation == "5'UTR" | 
                                                   grouped_annotation == "3'UTR" | 
                                                   grouped_annotation == "CDS" | 
                                                   grouped_annotation == "intron" |
                                                   grouped_annotation == "downstream 10K") &
                                                  ((NLS_E_M_BC >= BC_Threshold_CoCLIP & 
                                                    NES_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                   (Nuc_F_M_BC >= BC_Threshold_Fraction & 
                                                      Cyto_F_M_BC >= BC_Threshold_Fraction)) &
                                                  F_rowSum >= median(peakEnrichment$F_rowSum)*2 &
                                                  E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(E_NvC_M = Mock_Peaks_Filtered$E_NvC_M, F_NvC_M = Mock_Peaks_Filtered$F_NvC_M)
data$annotation = factor(Mock_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

ggplot(data, aes(x = log2(E_NvC_M), y = log2(F_NvC_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Cytoplasm)', y = 'log2(Fractionation CLIP Nuclear/Cytoplasm)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: CoCLIP vs Fractionation CLIP of Nuclear/Cytoplasm (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Mock Comparison - non-mRNA specific:
Mock_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != "5'UTR" & 
                                                   grouped_annotation != "3'UTR" & 
                                                   grouped_annotation != "CDS" & 
                                                   grouped_annotation != "intron" & 
                                                   grouped_annotation != "downstream 10K") &
                                                  ((NLS_E_M_BC >= BC_Threshold_CoCLIP & 
                                                      NES_E_M_BC >= BC_Threshold_CoCLIP) | 
                                                     (Nuc_F_M_BC >= BC_Threshold_Fraction & 
                                                        Cyto_F_M_BC >= BC_Threshold_Fraction)) &
                                                  F_rowSum >= median(peakEnrichment$F_rowSum)*2 &
                                                  E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(E_NvC_M = Mock_Peaks_Filtered$E_NvC_M, F_NvC_M = Mock_Peaks_Filtered$F_NvC_M)
data$annotation = factor(Mock_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

ggplot(data, aes(x = log2(E_NvC_M), y = log2(F_NvC_M), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Cytoplasm)', y = 'log2(Fractionation CLIP Nuclear/Cytoplasm)') +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  ggtitle(paste0('Mock HuR Peaks: CoCLIP vs Fractionation CLIP of Nuclear/Cytoplasm (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress Comparison
Stress_Peaks_Filtered = peakEnrichment %>% filter(((NLS_E_S_BC >= BC_Threshold_CoCLIP & 
                                                    NES_E_S_BC >= BC_Threshold_CoCLIP) | 
                                                   (Nuc_F_S_BC >= BC_Threshold_Fraction & 
                                                    Cyto_F_S_BC >= BC_Threshold_Fraction)) &
                                                    F_rowSum >= median(peakEnrichment$F_rowSum)*2 &
                                                    E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(E_NvC_S= Stress_Peaks_Filtered$E_NvC_S, F_NvC_S = Stress_Peaks_Filtered$F_NvC_S)
data$annotation = factor(Stress_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", 'snoRNA', "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

ggplot(data, aes(x = log2(E_NvC_S), y = log2(F_NvC_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Cytoplasm)', y = 'log2(Fractionation CLIP Nuclear/Cytoplasm)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Stress HuR Peaks: CoCLIP vs Fractionation CLIP of Nuclear/Cytoplasm (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress Comparison - mRNA specific
Stress_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation == "5'UTR" | 
                                                     grouped_annotation == "3'UTR" | 
                                                     grouped_annotation == "CDS" | 
                                                     grouped_annotation == "intron" |
                                                     grouped_annotation == "downstream 10K") &
                                                    ((NLS_E_S_BC >= BC_Threshold_CoCLIP & 
                                                      NES_E_S_BC >= BC_Threshold_CoCLIP) | 
                                                     (Nuc_F_S_BC >= BC_Threshold_Fraction & 
                                                        Cyto_F_S_BC >= BC_Threshold_Fraction)) &
                                                    F_rowSum >= median(peakEnrichment$F_rowSum)*2 &
                                                    E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(E_NvC_S= Stress_Peaks_Filtered$E_NvC_S, F_NvC_S = Stress_Peaks_Filtered$F_NvC_S)
data$annotation = factor(Stress_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

ggplot(data, aes(x = log2(E_NvC_S), y = log2(F_NvC_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Cytoplasm)', y = 'log2(Fractionation CLIP Nuclear/Cytoplasm)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Stress HuR Peaks: CoCLIP vs Fractionation CLIP of Nuclear/Cytoplasm (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# Stress Comparison - non-mRNA specific
Stress_Peaks_Filtered = peakEnrichment %>% filter((grouped_annotation != "5'UTR" & 
                                                     grouped_annotation != "3'UTR" & 
                                                     grouped_annotation != "CDS" & 
                                                     grouped_annotation != "intron" & 
                                                     grouped_annotation != "downstream 10K") &
                                                    ((NLS_E_S_BC >= BC_Threshold_CoCLIP & 
                                                        NES_E_S_BC >= BC_Threshold_CoCLIP) | 
                                                       (Nuc_F_S_BC >= BC_Threshold_Fraction & 
                                                          Cyto_F_S_BC >= BC_Threshold_Fraction)) &
                                                    F_rowSum >= median(peakEnrichment$F_rowSum)*2 &
                                                    E_rowSum >= median(peakEnrichment$E_rowSum)*2)

data = data.frame(E_NvC_S= Stress_Peaks_Filtered$E_NvC_S, F_NvC_S = Stress_Peaks_Filtered$F_NvC_S)
data$annotation = factor(Stress_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

ggplot(data, aes(x = log2(E_NvC_S), y = log2(F_NvC_S), color = annotation)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Nuclear/Cytoplasm)', y = 'log2(Fractionation CLIP Nuclear/Cytoplasm)') +
  xlim(c(-6, 6)) +
  ylim(c(-6, 6)) +
  ggtitle(paste0('Stress HuR Peaks: CoCLIP vs Fractionation CLIP of Nuclear/Cytoplasm (',  nrow(data), ' peaks)')) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

####################################################################################################################

## Input and CoCLIP comparison
####################################################################################################################
# Nuclear Mock:
NLS_Mock_Peaks_Filtered = peakRowSum %>% filter(NLS_I_M_BC >= BC_Threshold_Input & 
                                                      NLS_E_M_BC >= BC_Threshold_CoCLIP  &
                                                      I_rowSum >= median(peakRowSum$I_rowSum)*2)
NLS_Mock_Peaks_Filtered$grouped_annotation = factor(NLS_Mock_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", 'snoRNA', "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

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
NLS_Mock_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation == "5'UTR" | 
                                                   grouped_annotation == "3'UTR" | 
                                                   grouped_annotation == "CDS" | 
                                                   grouped_annotation == "intron" |
                                                   grouped_annotation == "downstream 10K") &
                                                  NLS_I_M_BC >= BC_Threshold_Input & 
                                                  NLS_E_M_BC >= BC_Threshold_CoCLIP  &
                                                  I_rowSum >= median(peakRowSum$I_rowSum)*2)
NLS_Mock_Peaks_Filtered$finalized_annotation = factor(NLS_Mock_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

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
NLS_Mock_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation != "5'UTR" & 
                                                   grouped_annotation != "3'UTR" & 
                                                   grouped_annotation != "CDS" & 
                                                   grouped_annotation != "intron" & 
                                                   grouped_annotation != "downstream 10K") &
                                                  NLS_I_M_BC >= BC_Threshold_Input & 
                                                  NLS_E_M_BC >= BC_Threshold_CoCLIP  &
                                                  I_rowSum >= median(peakRowSum$I_rowSum)*2)
NLS_Mock_Peaks_Filtered$finalized_annotation = factor(NLS_Mock_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

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
NES_Mock_Peaks_Filtered = peakRowSum %>% filter(NES_I_M_BC >= BC_Threshold_Input & 
                                                      NES_E_M_BC >= BC_Threshold_CoCLIP &
                                                      I_rowSum >= median(peakRowSum$I_rowSum)*2)
NES_Mock_Peaks_Filtered$grouped_annotation = factor(NES_Mock_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", 'snoRNA', "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

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
NES_Mock_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation == "5'UTR" | 
                                                   grouped_annotation == "3'UTR" | 
                                                   grouped_annotation == "CDS" | 
                                                   grouped_annotation == "intron" |
                                                   grouped_annotation == "downstream 10K") &
                                                  NES_I_M_BC >= BC_Threshold_Input & 
                                                  NES_E_M_BC >= BC_Threshold_CoCLIP &
                                                  I_rowSum >= median(peakRowSum$I_rowSum)*2)
NES_Mock_Peaks_Filtered$finalized_annotation = factor(NES_Mock_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

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
NES_Mock_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation != "5'UTR" & 
                                                   grouped_annotation != "3'UTR" & 
                                                   grouped_annotation != "CDS" & 
                                                   grouped_annotation != "intron" & 
                                                   grouped_annotation != "downstream 10K") &
                                                  NES_I_M_BC >= BC_Threshold_Input & 
                                                  NES_E_M_BC >= BC_Threshold_CoCLIP &
                                                  I_rowSum >= median(peakRowSum$I_rowSum)*2)
NES_Mock_Peaks_Filtered$finalized_annotation = factor(NES_Mock_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

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
G3BP_Mock_Peaks_Filtered = peakRowSum %>% filter(G3BP_I_M_BC >= BC_Threshold_Input_SG & 
                                                       G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG &
                                                       I_rowSum >= median(peakRowSum$I_rowSum)*2)
G3BP_Mock_Peaks_Filtered$grouped_annotation = factor(G3BP_Mock_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", 'snoRNA', "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

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
G3BP_Mock_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation == "5'UTR" | 
                                                    grouped_annotation == "3'UTR" | 
                                                    grouped_annotation == "CDS" | 
                                                    grouped_annotation == "intron" |
                                                    grouped_annotation == "downstream 10K") &
                                                   G3BP_I_M_BC >= BC_Threshold_Input_SG & 
                                                   G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG &
                                                   I_rowSum >= median(peakRowSum$I_rowSum)*2)
G3BP_Mock_Peaks_Filtered$finalized_annotation = factor(G3BP_Mock_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

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
G3BP_Mock_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation != "5'UTR" & 
                                                    grouped_annotation != "3'UTR" & 
                                                    grouped_annotation != "CDS" & 
                                                    grouped_annotation != "intron" & 
                                                    grouped_annotation != "downstream 10K") &
                                                   G3BP_I_M_BC >= BC_Threshold_Input_SG & 
                                                   G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG &
                                                   I_rowSum >= median(peakRowSum$I_rowSum)*2)
G3BP_Mock_Peaks_Filtered$finalized_annotation = factor(G3BP_Mock_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

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
NLS_Stress_Peaks_Filtered = peakRowSum %>% filter(NLS_I_S_BC >= BC_Threshold_Input & 
                                                        NLS_E_S_BC >= BC_Threshold_CoCLIP &
                                                        I_rowSum >= median(peakRowSum$I_rowSum)*2)
NLS_Stress_Peaks_Filtered$grouped_annotation = factor(NLS_Stress_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", 'snoRNA', "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

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
NLS_Stress_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation == "5'UTR" | 
                                                     grouped_annotation == "3'UTR" | 
                                                     grouped_annotation == "CDS" | 
                                                     grouped_annotation == "intron" |
                                                     grouped_annotation == "downstream 10K") &
                                                    NLS_I_S_BC >= BC_Threshold_Input & 
                                                    NLS_E_S_BC >= BC_Threshold_CoCLIP &
                                                    I_rowSum >= median(peakRowSum$I_rowSum)*2)
NLS_Stress_Peaks_Filtered$finalized_annotation = factor(NLS_Stress_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

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
NLS_Stress_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation != "5'UTR" & 
                                                     grouped_annotation != "3'UTR" & 
                                                     grouped_annotation != "CDS" & 
                                                     grouped_annotation != "intron" & 
                                                     grouped_annotation != "downstream 10K") &
                                                    NLS_I_S_BC >= BC_Threshold_Input & 
                                                    NLS_E_S_BC >= BC_Threshold_CoCLIP &
                                                    I_rowSum >= median(peakRowSum$I_rowSum)*2)
NLS_Stress_Peaks_Filtered$finalized_annotation = factor(NLS_Stress_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

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
NES_Stress_Peaks_Filtered = peakRowSum %>% filter(NES_I_S_BC >= BC_Threshold_Input & 
                                                        NES_E_S_BC >= BC_Threshold_CoCLIP &
                                                        I_rowSum >= median(peakRowSum$I_rowSum)*2)
NES_Stress_Peaks_Filtered$grouped_annotation = factor(NES_Stress_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", 'snoRNA', "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

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
NES_Stress_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation == "5'UTR" | 
                                                     grouped_annotation == "3'UTR" | 
                                                     grouped_annotation == "CDS" | 
                                                     grouped_annotation == "intron" |
                                                     grouped_annotation == "downstream 10K") &
                                                    NES_I_S_BC >= BC_Threshold_Input & 
                                                    NES_E_S_BC >= BC_Threshold_CoCLIP &
                                                    I_rowSum >= median(peakRowSum$I_rowSum)*2)
NES_Stress_Peaks_Filtered$finalized_annotation = factor(NES_Stress_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

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
NES_Stress_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation != "5'UTR" & 
                                                     grouped_annotation != "3'UTR" & 
                                                     grouped_annotation != "CDS" & 
                                                     grouped_annotation != "intron" & 
                                                     grouped_annotation != "downstream 10K") &
                                                    NES_I_S_BC >= BC_Threshold_Input & 
                                                    NES_E_S_BC >= BC_Threshold_CoCLIP &
                                                    I_rowSum >= median(peakRowSum$I_rowSum)*2)
NES_Stress_Peaks_Filtered$finalized_annotation = factor(NES_Stress_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

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
G3BP_Stress_Peaks_Filtered = peakRowSum %>% filter(G3BP_I_S_BC >= BC_Threshold_Input_SG & 
                                                         G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG &
                                                         I_rowSum >= median(peakRowSum$I_rowSum)*2)
G3BP_Stress_Peaks_Filtered$grouped_annotation = factor(G3BP_Stress_Peaks_Filtered$grouped_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", 'snoRNA', "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

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
G3BP_Stress_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation == "5'UTR" | 
                                                      grouped_annotation == "3'UTR" | 
                                                      grouped_annotation == "CDS" | 
                                                      grouped_annotation == "intron" |
                                                      grouped_annotation == "downstream 10K") &
                                                     G3BP_I_S_BC >= BC_Threshold_Input_SG & 
                                                     G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG &
                                                     I_rowSum >= median(peakRowSum$I_rowSum)*2)
G3BP_Stress_Peaks_Filtered$finalized_annotation = factor(G3BP_Stress_Peaks_Filtered$finalized_annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "CDS_Retained_intron", 'ncRNA_Retained_intron', 'downstream 10K'))

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
G3BP_Stress_Peaks_Filtered = peakRowSum %>% filter((grouped_annotation != "5'UTR" & 
                                                      grouped_annotation != "3'UTR" & 
                                                      grouped_annotation != "CDS" & 
                                                      grouped_annotation != "intron" & 
                                                      grouped_annotation != "downstream 10K") &
                                                     G3BP_I_S_BC >= BC_Threshold_Input_SG & 
                                                     G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG &
                                                     I_rowSum >= median(peakRowSum$I_rowSum)*2)
G3BP_Stress_Peaks_Filtered$finalized_annotation = factor(G3BP_Stress_Peaks_Filtered$finalized_annotation, levels = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'miscRNA', 'TE', 'deep intergenic', 'Other'))

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

