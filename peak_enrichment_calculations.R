## CoCLIP Analysis: 
## Peak Enrichment Calculation
## Written by Soon Yi
## Last Edit: 2023-08-28

library(stringr)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)



## Load peak matrix and clean up:
####################################################################################################################
# peaksMatrix_PATH = 'L:/My Drive/CWRU/PhD/Luna Lab/1. coCLIP/Analysis/peaks/'    ## Use this for windows machine
peaksMatrix_PATH = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/'
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

## Subset of Peaks for downstream analysis:
####################################################################################################################
BC_Threshold_Fraction = 2
BC_Threshold_Input = 1
BC_Threshold_Input_SG = 2
BC_Threshold_CoCLIP = 2
BC_Threshold_CoCLIP_SG = 3

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


Peak_F_Nuc_M = (peakMatrix[, c(inert_columns, Nuc_F_M, BC_columns)] 
                %>% filter(Nuc_F_M_BC >= BC_Threshold_Fraction))

Peak_F_Nuc_S = (peakMatrix[, c(inert_columns, Nuc_F_S, BC_columns)] 
                %>% filter(Nuc_F_S_BC >= BC_Threshold_Fraction))

Peak_F_Cyt_M = (peakMatrix[, c(inert_columns, Cyto_F_M, BC_columns)] 
                %>% filter(Cyto_F_M_BC >= BC_Threshold_Fraction))

Peak_F_Cyt_S = (peakMatrix[, c(inert_columns, Cyto_F_S, BC_columns)] 
                %>% filter(Cyto_F_S_BC >= BC_Threshold_Fraction))

Peak_F_all_M = (peakMatrix[, c(inert_columns, Nuc_F_M, Cyto_F_M, BC_columns)] 
                %>% filter(Nuc_F_M_BC >= BC_Threshold_Fraction & Cyto_F_M_BC >= BC_Threshold_Fraction))

Peak_F_all_S = (peakMatrix[, c(inert_columns, Nuc_F_S, Cyto_F_S, BC_columns)] 
                %>% filter(Nuc_F_S_BC >= BC_Threshold_Fraction & Cyto_F_S_BC >= BC_Threshold_Fraction))


Peak_Co_Input_M = (peakMatrix[, c(inert_columns, NLS_I_M, NES_I_M, G3BP_I_M, BC_columns)] 
                   %>% filter(NLS_I_M_BC >= BC_Threshold_Input & NES_I_M_BC >= BC_Threshold_Input & G3BP_I_M_BC >= BC_Threshold_Input_SG))

Peak_Co_Input_S = (peakMatrix[, c(inert_columns, NLS_I_S, NES_I_S, G3BP_I_S, BC_columns)] 
                   %>% filter(NLS_I_S_BC >= BC_Threshold_Input & NES_I_S_BC >= BC_Threshold_Input & G3BP_I_S_BC >= BC_Threshold_Input_SG))

Peak_Co_NLS_M = (peakMatrix[, c(inert_columns, NLS_E_M, BC_columns)] 
                 %>% filter(NLS_E_M_BC >= BC_Threshold_CoCLIP))

Peak_Co_NLS_S = (peakMatrix[, c(inert_columns, NLS_E_S, BC_columns)] 
                 %>% filter(NLS_E_S_BC >= BC_Threshold_CoCLIP))

Peak_Co_NES_M = (peakMatrix[, c(inert_columns, NES_E_M, BC_columns)] 
                 %>% filter(NES_E_M_BC >= BC_Threshold_CoCLIP))

Peak_Co_NES_S = (peakMatrix[, c(inert_columns, NES_E_S, BC_columns)] 
                %>% filter(NES_E_S_BC >= BC_Threshold_CoCLIP))

Peak_Co_G3BP_M = (peakMatrix[, c(inert_columns, G3BP_E_M, BC_columns)] 
                  %>% filter(G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG))

Peak_Co_G3BP_S = (peakMatrix[, c(inert_columns, G3BP_E_S, BC_columns)] 
                  %>% filter(G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG))

# These outputs are in peaks format suitable for Homer motif analysis:
write.table(Peak_all_M[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_all_M.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
write.table(Peak_all_S[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_all_S.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
write.table(Peak_all[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_all.txt'), quote = F, sep = '\t', row.names = F, col.names = F)

write.table(Peak_F_Nuc_M[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_F_Nuc_M.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
write.table(Peak_F_Nuc_S[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_F_Nuc_S.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
write.table(Peak_F_Cyt_M[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_F_Cyt_M.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
write.table(Peak_F_Cyt_S[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_F_Cyt_S.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
write.table(Peak_F_all_M[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_F_all_M.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
write.table(Peak_F_all_S[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_F_all_S.txt'), quote = F, sep = '\t', row.names = F, col.names = F)

write.table(Peak_Co_Input_M[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_Co_Input_M.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
write.table(Peak_Co_Input_S[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_Co_Input_S.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
write.table(Peak_Co_NLS_M[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_Co_NLS_M.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
write.table(Peak_Co_NLS_S[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_Co_NLS_S.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
write.table(Peak_Co_NES_M[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_Co_NES_M.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
write.table(Peak_Co_NES_S[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_Co_NES_S.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
write.table(Peak_Co_G3BP_M[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_Co_G3BP_M.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
write.table(Peak_Co_G3BP_S[, c('peak_names', 'chrom', 'start', 'end', 'strand')], paste0(peaksMatrix_PATH, '/Peak_Subsets/Peak_Co_G3BP_S.txt'), quote = F, sep = '\t', row.names = F, col.names = F)

####################################################################################################################


## Exploratory Stacked Bar Plots Across All Samples:
####################################################################################################################
## Mock vs Stress:
PeakDistribution_all_M = data.frame(table(Peak_all_M$grouped_annotation), row.names = 1)
PeakDistribution_all_S = data.frame(table(Peak_all_S$grouped_annotation), row.names = 1)
PeakDistribution_all = data.frame(table(Peak_all$grouped_annotation), row.names = 1)

colnames(PeakDistribution_all_M) = c('All_M')
colnames(PeakDistribution_all_S) = c('All_S')
colnames(PeakDistribution_all) = c('All_M+S')

## Raw Counts Distribution Stacked Bar Graph
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

## Exploratory Stacked Bar Plots For Fractionation CLIP Samples:
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

## Raw Counts Distribution Stacked Bar Graph
PeakDistribution_F_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Nuc_S, 
                                    PeakDistribution_F_Cyt_M, PeakDistribution_F_Cyt_S, 
                                    PeakDistribution_F_all_M, PeakDistribution_F_all_S)
PeakDistribution_F_combined$Annotation = rownames(PeakDistribution_F_combined)

PeakDistribution_F_combined = PeakDistribution_F_combined %>%
  gather(key = "Source", value = "Freq", 'F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt') %>%
  select(Source, Freq, Annotation)

PeakDistribution_F_combined$Source = factor(PeakDistribution_F_combined$Source, levels = c('F_M_Nuc', 'F_S_Nuc', 'F_M_Cyt', 'F_S_Cyt', 'F_M_Nuc+Cyt', 'F_S_Nuc+Cyt'))
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress', 'Both Fraction Mock','Both Fraction Stress')) +
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
PeakDistribution_F_combined$Annotation = factor(PeakDistribution_F_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

ggplot(PeakDistribution_F_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Nuclear Mock', 'Nuclear Stress', 'Cyto Mock', 'Cyto Stress', 'Both Fraction Mock','Both Fraction Stress')) +
  ggtitle('All Peak Fractions Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")
####################################################################################################################

## Exploratory Stacked Bar Plots For CoCLIP Samples:
####################################################################################################################
## Mock vs Stress for Each Fraction
PeakDistribution_Co_Input_M = data.frame(table(Peak_Co_Input_M$grouped_annotation), row.names = 1)
PeakDistribution_Co_Input_S = data.frame(table(Peak_Co_Input_S$grouped_annotation), row.names = 1)

PeakDistribution_Co_NLS_M = data.frame(table(Peak_Co_NLS_M$grouped_annotation), row.names = 1)
PeakDistribution_Co_NLS_S = data.frame(table(Peak_Co_NLS_S$grouped_annotation), row.names = 1)

PeakDistribution_Co_NES_M = data.frame(table(Peak_Co_NES_M$grouped_annotation), row.names = 1)
# Use this when we are filtering to peaks with only a single annotation:
# PeakDistribution_Co_NES_M = rbind(PeakDistribution_Co_NES_M, c(0), c(0))
# row.names(PeakDistribution_Co_NES_M) = c(row.names(PeakDistribution_Co_NES_M)[1:7], 'CDS', 'TE')
# PeakDistribution_Co_NES_M = PeakDistribution_Co_NES_M[row.names(PeakDistribution_Co_NLS_S), 'Freq', drop = FALSE]

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

## Raw Counts Distribution Stacked Bar Graph
PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
                                     PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
                                     PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
                                     PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>%
  gather(key = "Source", value = "Freq", 'Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP') %>%
  select(Source, Freq, Annotation)

PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP'))
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

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
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = c("5'UTR", "CDS", "3'UTR", "intron", "ncRNA", "TE", "Other", "deep intergenic", "downstream 10K"))

ggplot(PeakDistribution_Co_combined, aes(fill = Annotation, y=Freq, x=Source)) + 
  geom_bar(position='stack', stat='identity') +
  scale_x_discrete(labels= c('Input Mock', 'Input Stress', 'NLS Mock', 'NLS Stress', 'NES Mock', 'NES Stress', 'G3BP Mock', 'G3BP Stress')) +
  ggtitle('All Peak Fractions Distributions') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")
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

## Nuclear Comparison 
####################################################################################################################
NuclearPeaks_Filtered = peakEnrichment %>% filter(NLS_I_M_BC >= BC_Threshold_Input, 
                                                  NLS_I_S_BC >= BC_Threshold_Input, 
                                                  NLS_E_M_BC >= BC_Threshold_CoCLIP, 
                                                  NLS_E_S_BC >= BC_Threshold_CoCLIP)

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
CytoPeaks_Filtered= peakEnrichment %>% filter(NES_I_M_BC >= BC_Threshold_Input, 
                                              NES_I_S_BC >= BC_Threshold_Input, 
                                              NES_E_M_BC >= BC_Threshold_CoCLIP, 
                                              NES_E_S_BC >= BC_Threshold_CoCLIP)

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
SGPeaks_Filtered= peakEnrichment %>% filter(G3BP_I_M_BC >= BC_Threshold_Input_SG, 
                                            G3BP_I_S_BC >= BC_Threshold_Input_SG, 
                                            G3BP_E_M_BC >= BC_Threshold_CoCLIP_SG, 
                                            G3BP_E_S_BC >= BC_Threshold_CoCLIP_SG)

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


