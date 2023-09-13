## CoCLIP Analysis: 
## Peak Processing for Plots
## Written by Soon Yi
## Last Edit: 2023-09-11

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
NES_I_M = c('NES_I_M_1', 'NES_I_M_2')
G3BP_I_M = c('G3BP_I_M_1', 'G3BP_I_M_2', 'G3BP_I_M_3', 'G3BP_I_M_4')

NLS_I_S = c('NLS_I_S_1', 'NLS_I_S_2')
NES_I_S = c('NES_I_S_1', 'NES_I_S_2')
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

## Add pseudocount:
pseudoCount = min(peaksMatrix[, colnames(peaksMatrix)[6:63]][peaksMatrix[, colnames(peaksMatrix)[6:63]] != 0], na.rm = TRUE)
peaksMatrix[, colnames(peaksMatrix)[6:63]] = peaksMatrix[, colnames(peaksMatrix)[6:63]] + pseudoCount

####################################################################################################################

## Custom Functions 
####################################################################################################################
## Filter peaks based on the designated criteria:
peakFilter = function(peak_matrix, sample_list, info_columns, BC_criteria, rowSum_criteria = NULL) {
  temp = peak_matrix[, c(info_columns, sample_list)]
  
  if (!is.null(rowSum_criteria)) {
    temp$rowAvg = rowSums(temp[, sample_list]) / length(sample_list)
    temp = temp[(rowSums(temp[, paste0(unique(sub("_[0-9]+$", "", sample_list)), '_BC')]) >= BC_criteria) & (temp$rowAvg > (median(temp$rowAvg) * rowSum_criteria)), ]
    temp$rowAvg = NULL
  } else {
    temp = temp[(rowSums(temp[, paste0(unique(sub("_[0-9]+$", "", sample_list)), '_BC')]) >= BC_criteria), ]
  }
  
  return(temp)
}

## Get peak counts per gene: 
peaksPerGene = function(peak_matrix, sample_list) {
  temp = peak_matrix %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(across(all_of(sample_list))), peakCounts = n())
  temp = temp %>% mutate(tagDensity = tagCounts/peakCounts)
  
  return(temp)
}
####################################################################################################################

## Filter Criteria:
####################################################################################################################
BC_Threshold_F = 2
BC_Threshold_I = 4
BC_Threshold_E = 2
BC_Threshold_E_SG = 3


rowSum_Multiplier_F = 4
rowSum_Multiplier_I = 2
rowSum_Multiplier_E = 2
####################################################################################################################

## Subset of Peaks for downstream analysis:
####################################################################################################################
## Fractionation CLIP
Peak_F_Nuc_M = peakFilter(peaksMatrix, Nuc_F_M, c(inert_columns, BC_columns), BC_Threshold_F, rowSum_Multiplier_F)
Peak_F_Nuc_S = peakFilter(peaksMatrix, Nuc_F_S, c(inert_columns, BC_columns), BC_Threshold_F, rowSum_Multiplier_F)
Peak_F_Cyt_M = peakFilter(peaksMatrix, Cyto_F_M, c(inert_columns, BC_columns), BC_Threshold_F, rowSum_Multiplier_F)
Peak_F_Cyt_S = peakFilter(peaksMatrix, Cyto_F_S, c(inert_columns, BC_columns), BC_Threshold_F, rowSum_Multiplier_F)

## CoCLIP
Peak_Co_Input_M = peakFilter(peaksMatrix, c(NLS_I_M, NES_I_M, G3BP_I_M), c(inert_columns, BC_columns), BC_Threshold_I, rowSum_Multiplier_I)
Peak_Co_Input_S = peakFilter(peaksMatrix, c(NLS_I_S, NES_I_S, G3BP_I_S), c(inert_columns, BC_columns), BC_Threshold_I, rowSum_Multiplier_I)
Peak_Co_NLS_M = peakFilter(peaksMatrix, NLS_E_M, c(inert_columns, BC_columns), BC_Threshold_E, rowSum_Multiplier_E)
Peak_Co_NLS_S = peakFilter(peaksMatrix, NLS_E_S, c(inert_columns, BC_columns), BC_Threshold_E, rowSum_Multiplier_E)
Peak_Co_NES_M = peakFilter(peaksMatrix, NES_E_M, c(inert_columns, BC_columns), BC_Threshold_E, rowSum_Multiplier_E)
Peak_Co_NES_S = peakFilter(peaksMatrix, NES_E_S, c(inert_columns, BC_columns), BC_Threshold_E, rowSum_Multiplier_E)
Peak_Co_G3BP_M = peakFilter(peaksMatrix, G3BP_E_M, c(inert_columns, BC_columns), BC_Threshold_E_SG, rowSum_Multiplier_E)
Peak_Co_G3BP_S = peakFilter(peaksMatrix, G3BP_E_S, c(inert_columns, BC_columns), BC_Threshold_E_SG, rowSum_Multiplier_E)
####################################################################################################################

## Peak Counts per Gene:
####################################################################################################################
## Fractionation CLIP
PPG_F_Nuc_M = peaksPerGene(Peak_F_Nuc_M, Nuc_F_M)
PPG_F_Nuc_S = peaksPerGene(Peak_F_Nuc_S, Nuc_F_S)
PPG_F_Cyto_M = peaksPerGene(Peak_F_Cyt_M, Cyto_F_M)
PPG_F_Cyto_S = peaksPerGene(Peak_F_Cyt_S, Cyto_F_S)

## CoCLIP
PPG_Co_Input_M = peaksPerGene(Peak_Co_Input_M, c(NLS_I_M, NES_I_M, G3BP_I_M))
PPG_Co_Input_S = peaksPerGene(Peak_Co_Input_S, c(NLS_I_S, NES_I_S, G3BP_I_S))
PPG_Co_NLS_M = peaksPerGene(Peak_Co_NLS_M, NLS_E_M)
PPG_Co_NLS_S = peaksPerGene(Peak_Co_NLS_S, NLS_E_S)
PPG_Co_NES_M = peaksPerGene(Peak_Co_NES_M, NES_E_M)
PPG_Co_NES_S = peaksPerGene(Peak_Co_NES_S, NES_E_S)
PPG_Co_G3BP_M = peaksPerGene(Peak_Co_G3BP_M, G3BP_E_M)
PPG_Co_G3BP_S = peaksPerGene(Peak_Co_G3BP_S, G3BP_E_S)
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
