## CoCLIP Analysis: 
## Peak Enrichment Calculation
## Written by Soon Yi
## Last Edit: 2023-08-25

library(stringr)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)

## Load peak matrix and clean up:
####################################################################################################################
peaksMatrix_PATH = 'L:/My Drive/CWRU/PhD/Luna Lab/1. coCLIP/Analysis/peaks/'  ## Use this for windows machine
# peaksMatrix_PATH = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/'
peaksMatrix_FILE = 'Combined_peakCoverage_groomed_normalized_annotated.txt'

peakMatrix = read_delim(paste0(peaksMatrix_PATH, peaksMatrix_FILE), show_col_types = FALSE)
peakMatrix = peakMatrix %>% mutate_at('TOTAL_BC', as.numeric)

## Filter by median
# peakCount_median = median(peakMatrix[, colnames(peakMatrix)[6:63]][peakMatrix[, colnames(peakMatrix)[6:63]] != 0], na.rm = TRUE)
# peakMatrix = peakMatrix %>% filter(if_any(all_of(comparison_vector), ~ . > peakCount_median ))

## Instead of filtering, add pseudocounts.
## This might be better because we will be filtering by BC down the road.
pseudoCount = min(peakMatrix[, colnames(peakMatrix)[6:63]][peakMatrix[, colnames(peakMatrix)[6:63]] != 0], na.rm = TRUE)
peakMatrix[, colnames(peakMatrix)[6:63]] = peakMatrix[, colnames(peakMatrix)[6:63]] + pseudoCount

## Column organization:
inert_columns = c('chrom', 'start', 'end', 'peak_names', 'score', 'strand', 
                  "gene", "annotation", "finalized_annotation", "grouped_annotation", "annotation_count")
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
G3BP_I_S = c('G3BP_I_S_1', 'G3BP_I_S_2', 'G3BP_I_S_3', 'G3BP_I_S_4')

NLS_E_M = c('NLS_E_M_1', 'NLS_E_M_2', 'NLS_E_M_3', 'NLS_E_M_4')
NLS_E_S = c('NLS_E_S_1', 'NLS_E_S_2', 'NLS_E_S_3', 'NLS_E_S_4')
NES_E_M = c('NES_E_M_1', 'NES_E_M_2', 'NES_E_M_3', 'NES_E_M_4')
NES_E_S = c('NES_E_S_1', 'NES_E_S_2', 'NES_E_S_3', 'NES_E_S_4')
G3BP_E_M = c('G3BP_E_M_1', 'G3BP_E_M_2', 'G3BP_E_M_3', 'G3BP_E_M_4', 'G3BP_E_M_5', 'G3BP_E_M_6')
G3BP_E_S = c('G3BP_E_S_1', 'G3BP_E_S_2', 'G3BP_E_S_3', 'G3BP_E_S_4', 'G3BP_E_S_5', 'G3BP_E_S_6', 'G3BP_E_S_7')
####################################################################################################################

## Initial Exploratory Plots:
####################################################################################################################
## All Mock vs Stress:
PeakDistribution_all_M = data.frame(table((peakMatrix[, c(inert_columns, Nuc_F_M, Cyto_F_M, NLS_I_M, NES_I_M, G3BP_I_M, NLS_E_M, NES_E_M, G3BP_E_M, BC_columns)] 
                                           %>% filter((Nuc_F_M_BC >= 2 & Cyto_F_M_BC >= 2) |  
                                                      (NLS_I_M_BC >= 2 & NES_I_M_BC >= 2 & G3BP_I_M_BC >= 3) |
                                                      NLS_E_M_BC >= 2 & NES_E_M_BC >= 2 & G3BP_E_M_BC >= 3))
                                          $grouped_annotation), row.names = 1)


PeakDistribution_all_S = data.frame(table((peakMatrix[, c(inert_columns, Nuc_F_S, Cyto_F_S, NLS_I_S, NES_I_S, G3BP_I_S, NLS_E_S, NES_E_S, G3BP_E_S, BC_columns)] 
                                           %>% filter((Nuc_F_S_BC >= 2 & Cyto_F_S_BC >= 2) |
                                                      (NLS_I_S_BC >= 2 & NES_I_S_BC >= 2 & G3BP_I_S_BC >= 3) | 
                                                      (NLS_E_S_BC >= 2 & NES_E_S_BC >= 2 & G3BP_E_S_BC >= 3)))
                                          $grouped_annotation), row.names = 1)

PeakDistribution_all = data.frame(table((peakMatrix[, c(inert_columns, Nuc_F_S, Cyto_F_S, NLS_I_S, NES_I_S, G3BP_I_S, NLS_E_S, NES_E_S, G3BP_E_S, BC_columns)] 
                                         %>% filter(((Nuc_F_M_BC >= 2 & Cyto_F_M_BC >= 2) | 
                                                       (NLS_I_M_BC >= 2 & NES_I_M_BC >= 2 & G3BP_I_M_BC >= 3) | 
                                                       (NLS_E_M_BC >= 2 & NES_E_M_BC >= 2 & G3BP_E_M_BC >= 3)) | 
                                                      ((Nuc_F_S_BC >= 2 & Cyto_F_S_BC >= 2) | 
                                                         (NLS_I_S_BC >= 2 & NES_I_S_BC >= 2 & G3BP_I_S_BC >= 3) | 
                                                         NLS_E_S_BC >= 2 & NES_E_S_BC >= 2 & G3BP_E_S_BC >= 3)))
                                        $grouped_annotation), row.names = 1)

colnames(PeakDistribution_all_M) = c('All_M')
colnames(PeakDistribution_all_S) = c('All_S')
colnames(PeakDistribution_all) = c('All_M+S')

## Raw Counts Distribution Stacked Bar Graph
PeakDistribution_all_combined = cbind(PeakDistribution_all_M, PeakDistribution_all_S, PeakDistribution_all)
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
  scale_fill_brewer(palette = "Set1")

## Fractionation CLIP Specific:
## Mock vs Stress
PeakDistribution_F_Nuc_M = data.frame(table((peakMatrix[, c(inert_columns, Nuc_F_M, BC_columns)] 
                                             %>% filter(Nuc_F_M_BC >= 2))
                                            $grouped_annotation), row.names = 1)

PeakDistribution_F_Cyt_M = data.frame(table((peakMatrix[, c(inert_columns, Cyto_F_M, BC_columns)] 
                                             %>% filter(Cyto_F_M_BC >= 2))
                                            $grouped_annotation), row.names = 1)
PeakDistribution_F_all_M = data.frame(table((peakMatrix[, c(inert_columns, Cyto_F_M, BC_columns)] 
                                             %>% filter(Nuc_F_M_BC >= 2 & Cyto_F_M_BC >= 2))
                                            $grouped_annotation), row.names = 1)


PeakDistribution_F_Nuc_S = table((peakMatrix[, c(inert_columns, Nuc_F_S, BC_columns)] %>% filter(Nuc_F_S_BC >= 2))$grouped_annotation)
PeakDistribution_F_Cyt_S = table((peakMatrix[, c(inert_columns, Cyto_F_S, BC_columns)] %>% filter(Cyto_F_S_BC >= 2))$grouped_annotation)
PeakDistribution_F_all_S

####################################################################################################################


## Build Enrichment Table:
####################################################################################################################
peakEnrichment = peakMatrix[, inert_columns]

peakEnrichment = peakEnrichment %>% mutate(Nuc_F_M = rowSums(peakMatrix[, Nuc_F_M])/length(Nuc_F_M))
peakEnrichment = peakEnrichment %>% mutate(Nuc_F_S = rowSums(peakMatrix[, Nuc_F_S])/length(Nuc_F_S))
peakEnrichment = peakEnrichment %>% mutate(Cyto_F_M = rowSums(peakMatrix[, Cyto_F_M])/length(Cyto_F_M))
peakEnrichment = peakEnrichment %>% mutate(Cyto_F_S = rowSums(peakMatrix[, Cyto_F_S])/length(Cyto_F_S))

peakEnrichment = peakEnrichment %>% mutate(NLS_I_M = rowSums(peakMatrix[, NLS_I_M])/length(NLS_I_M))
peakEnrichment = peakEnrichment %>% mutate(NLS_I_S = rowSums(peakMatrix[, NLS_I_S])/length(NLS_I_S))
peakEnrichment = peakEnrichment %>% mutate(NES_I_M = rowSums(peakMatrix[, NES_I_M])/length(NES_I_M))
peakEnrichment = peakEnrichment %>% mutate(NES_I_S = rowSums(peakMatrix[, NES_I_S])/length(NES_I_S))
peakEnrichment = peakEnrichment %>% mutate(G3BP_I_M = rowSums(peakMatrix[, G3BP_I_M])/length(G3BP_I_M))
peakEnrichment = peakEnrichment %>% mutate(G3BP_I_S = rowSums(peakMatrix[, G3BP_I_S])/length(G3BP_I_S))

peakEnrichment = peakEnrichment %>% mutate(NLS_E_M = rowSums(peakMatrix[, NLS_E_M])/length(NLS_E_M))
peakEnrichment = peakEnrichment %>% mutate(NLS_E_S = rowSums(peakMatrix[, NLS_E_S])/length(NLS_E_S))
peakEnrichment = peakEnrichment %>% mutate(NES_E_M = rowSums(peakMatrix[, NES_E_M])/length(NES_E_M))
peakEnrichment = peakEnrichment %>% mutate(NES_E_S = rowSums(peakMatrix[, NES_E_S])/length(NES_E_S))
peakEnrichment = peakEnrichment %>% mutate(G3BP_E_M = rowSums(peakMatrix[, G3BP_E_M])/length(G3BP_E_M))
peakEnrichment = peakEnrichment %>% mutate(G3BP_E_S = rowSums(peakMatrix[, G3BP_E_S])/length(G3BP_E_S))

peakEnrichment = cbind(peakEnrichment, peakMatrix[, BC_columns])
# write.table(peakEnrichment, paste0(peaksMatrix_PATH, str_replace(peaksMatrix_FILE, ".txt", "_enrichment.txt")), 
# quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t', na = "")
####################################################################################################################

## snoRNA Comparison
####################################################################################################################
snoRNA_Enrichment = peakEnrichment %>% filter(finalized_annotation == 'snoRNA')
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


Input_BC_Threshold = 1 # out of 2
CoCLIP_BC_Threshold = 2 # out of 4

SG_Input_BC_Threshold = 2 # out of 4
SG_CoCLIP_BC_Threshold = 3 # out of 6/7

## Nuclear Comparison 
####################################################################################################################
NuclearPeaks_Filtered = peakEnrichment %>% filter(NLS_I_M_BC >= Input_BC_Threshold, 
                                                  NLS_I_S_BC >= Input_BC_Threshold, 
                                                  NLS_E_M_BC >= CoCLIP_BC_Threshold, 
                                                  NLS_E_S_BC >= CoCLIP_BC_Threshold)

# Mock vs Stress of Enriched/Input
Nu_M_IvE = NuclearPeaks_Filtered$NLS_E_M / NuclearPeaks_Filtered$NLS_I_M
Nu_S_IvE = NuclearPeaks_Filtered$NLS_E_S / NuclearPeaks_Filtered$NLS_I_S
data = data.frame(log_Nu_M_IvE = log2(Nu_M_IvE), log_Nu_S_IvE = log2(Nu_S_IvE))

ggplot(data, aes(x = log_Nu_M_IvE, y = log_Nu_S_IvE)) +
  geom_point(col = 'black', pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Mock Enriched/Input)', y = 'log2(Stress Enriched/Input)') +
  xlim(c(-5, 5)) +
  ylim(c(-5, 5)) +
  ggtitle('Nuclear HuR Peaks: Mock vs Stress of Enriched/Input')


# CoCLIP vs Input of Mock/Stress
Nu_E_MvS = NuclearPeaks_Filtered$NLS_E_M / NuclearPeaks_Filtered$NLS_E_S
Nu_I_MvS = NuclearPeaks_Filtered$NLS_I_M / NuclearPeaks_Filtered$NLS_I_S
data = data.frame(log_Nu_E_MvS = log2(Nu_E_MvS), log_Nu_I_MvS = log2(Nu_I_MvS))

ggplot(data, aes(x = log_Nu_E_MvS, y = log_Nu_I_MvS)) +
  geom_point(col = 'black', pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Mock/Stress)', y = 'log2(Input Mock/Stress)') +
  xlim(c(-5, 5)) +
  ylim(c(-5, 5)) +
  ggtitle('Nuclear HuR Peaks: CoCLIP vs Input of Mock/Stress')

####################################################################################################################

## Cytoplasm Comparison 
####################################################################################################################
CytoPeaks_Filtered= peakEnrichment %>% filter(NES_I_M_BC >= Input_BC_Threshold, 
                                              NES_I_S_BC >= Input_BC_Threshold, 
                                              NES_E_M_BC >= CoCLIP_BC_Threshold, 
                                              NES_E_S_BC >= CoCLIP_BC_Threshold)

# Mock vs Stress of Enriched/Input
Cy_M_IvE = CytoPeaks_Filtered$NES_E_M / CytoPeaks_Filtered$NES_I_M
Cy_S_IvE = CytoPeaks_Filtered$NES_E_S / CytoPeaks_Filtered$NES_I_S
data = data.frame(log_Cy_M_IvE = log2(Cy_M_IvE), log_Cy_S_IvE = log2(Cy_S_IvE))

ggplot(data, aes(x = log_Cy_M_IvE, y = log_Cy_S_IvE)) +
  geom_point(col = 'black', pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Mock Enriched/Input)', y = 'log2(Stress Enriched/Input)') +
  xlim(c(-5, 5)) +
  ylim(c(-5, 5)) +
  ggtitle('Cytoplasm HuR Peaks: Mock vs Stress of Enriched/Input')


# CoCLIP vs Input of Mock/Stress
Cy_E_MvS = CytoPeaks_Filtered$NES_E_M / CytoPeaks_Filtered$NES_E_S
Cy_I_MvS = CytoPeaks_Filtered$NES_I_M / CytoPeaks_Filtered$NES_I_S
data = data.frame(log_Cy_E_MvS = log2(Cy_E_MvS), log_Cy_I_MvS = log2(Cy_I_MvS))

ggplot(data, aes(x = log_Cy_E_MvS, y = log_Cy_I_MvS)) +
  geom_point(col = 'black', pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Mock/Stress)', y = 'log2(Input Mock/Stress)') +
  xlim(c(-5, 5)) +
  ylim(c(-5, 5)) +
  ggtitle('Cytoplasm HuR Peaks: CoCLIP vs Input of Mock/Stress')

####################################################################################################################

## Stress Granule Comparison 
####################################################################################################################
SGPeaks_Filtered= peakEnrichment %>% filter(G3BP_I_M_BC >= SG_Input_BC_Threshold, 
                                            G3BP_I_S_BC >= SG_Input_BC_Threshold, 
                                            G3BP_E_M_BC >= SG_CoCLIP_BC_Threshold, 
                                            G3BP_E_S_BC >= SG_CoCLIP_BC_Threshold)

# Mock vs Stress of Enriched/Input
SG_M_IvE = SGPeaks_Filtered$G3BP_E_M / SGPeaks_Filtered$G3BP_I_M
SG_S_IvE = SGPeaks_Filtered$G3BP_E_S / SGPeaks_Filtered$G3BP_I_S
data = data.frame(log_SG_M_IvE = log2(SG_M_IvE), log_SG_S_IvE = log2(SG_S_IvE))

ggplot(data, aes(x = log_SG_M_IvE, y = log_SG_S_IvE)) +
  geom_point(col = 'black', pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(Mock Enriched/Input)', y = 'log2(Stress Enriched/Input)') +
  xlim(c(-10, 10)) +
  ylim(c(-10, 10)) +
  ggtitle('Stress Granule HuR Peaks: Mock vs Stress of Enriched/Input')


# CoCLIP vs Input of Mock/Stress
SG_E_MvS = SGPeaks_Filtered$G3BP_E_M / SGPeaks_Filtered$G3BP_E_S
SG_I_MvS = SGPeaks_Filtered$G3BP_I_M / SGPeaks_Filtered$G3BP_I_S
data = data.frame(log_SG_E_MvS = log2(SG_E_MvS), log_SG_I_MvS = log2(SG_I_MvS))

ggplot(data, aes(x = log_SG_E_MvS, y = log_SG_I_MvS)) +
  geom_point(col = 'black', pch = 16, size = 3, alpha = 0.5) +
  labs(x = 'log2(CoCLIP Mock/Stress)', y = 'log2(Input Mock/Stress)') +
  xlim(c(-10, 10)) +
  ylim(c(-10, 10)) +
  ggtitle('Stress Granule HuR Peaks: CoCLIP vs Input of Mock/Stress')

####################################################################################################################


