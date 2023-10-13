## CoCLIP Analysis: 
## Peak Processing for Stackd Bar Graphs
## Written by Soon Yi
## Last Edit: 2023-10-12

library(stringr)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)

## Load peak matrix and clean up:
####################################################################################################################
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

# ## Add row sum columns for further filtering:
# peaksMatrix$F_rowSum = rowSums(peaksMatrix[, c(Nuc_F_M, Nuc_F_S, Cyto_F_M, Cyto_F_S)])
# peaksMatrix$I_rowSum = rowSums(peaksMatrix[, c(NLS_I_M, NLS_I_S, NES_I_M, NES_I_S, G3BP_I_M, G3BP_I_S)])
# peaksMatrix$E_rowSum = rowSums(peaksMatrix[, c(NLS_E_M, NLS_E_S, NES_E_M, NES_E_S, G3BP_E_M, G3BP_E_S)])
# rowSum_columns = c('F_rowSum', 'I_rowSum', 'E_rowSum')

## Add pseudocount:
pseudoCount = min(peaksMatrix[, colnames(peaksMatrix)[6:63]][peaksMatrix[, colnames(peaksMatrix)[6:63]] != 0], na.rm = TRUE)
peaksMatrix[, colnames(peaksMatrix)[6:63]] = peaksMatrix[, colnames(peaksMatrix)[6:63]] + pseudoCount

####################################################################################################################

## Custom Functions 
####################################################################################################################
## Filter peaks based on the designated criteria:
filterPeakMatrix = function(peak_matrix, sample_list, info_columns, BC_criteria, rowSum_criteria = NULL) {
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

## Filter peaks by annotation:
filterPeaksByAnnotation = function(peak_matrix, annotation_column, list_of_annotation) {
  temp = peak_matrix %>% filter({{ annotation_column }} %in% list_of_annotation)
  return(temp)
}

## Get Annotation counts:
countAnnotation = function(peak_matrix, annotation_column, new_column_name = NULL, annotation_to_skip = NULL, fraction = NULL) {
  temp = data.frame(table(peak_matrix[, annotation_column]), row.names = 1)
  if(!is.null(new_column_name)) {
    colnames(temp) = new_column_name
  }
  
  if(!is.null(annotation_to_skip)) {
    temp = temp[rownames(temp) != annotation_to_skip, , drop = FALSE]
  }
  
  if(!is.null(fraction)) {
    temp = temp/sum(temp)
  }
  
  return(temp)
}

## Fill Annotation counts if anything is mixing:
fillAnnotation = function(annotation_counts, annotation_list) {
  colnames((annotation_counts))
  
  temp = data.frame(Sample = numeric(length(annotation_list)))
  rownames(temp) = annotation_list
  temp2 = merge(temp, annotation_counts, by = "row.names", all = TRUE)
  temp2[is.na(temp2)] = 0 
  temp2$Sample = NULL
  
  rownames(temp2) = temp2$Row.names
  temp2 = temp2[ncRNA_List, -1, drop = FALSE]
  
  return(temp2)
}

## Plot Stacked bar:
plotStackedBar = function(annotation_counts, sample_list, sample_label, title, y_lim = NULL, y_tick = NULL) {
  plot = ggplot(annotation_counts %>% filter(Source %in% sample_list), aes(fill = Annotation, y=Freq, x=Source)) + 
    geom_bar(position='stack', stat='identity') +
    scale_x_discrete(labels = sample_label) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_brewer(palette = "Set3") +
    theme_bw() + 
    theme(axis.text = element_text(size=14), 
          axis.title = element_text(size=14, face = 'bold'), 
          legend.text = element_text(size=14))
  
  if(!is.null(y_lim)) {
    plot = plot + ylim(y_lim)
  }
  
  if (!is.null(y_tick)) {
    plot = plot + scale_y_continuous(breaks = seq(0, y_lim[2], by=y_tick), limits=c(0, y_lim[2]))
  }
  
  return(plot)
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
Peak_F_Nuc_M = filterPeakMatrix(peaksMatrix, Nuc_F_M, c(inert_columns, BC_columns), BC_Threshold_F, rowSum_Multiplier_F)
Peak_F_Nuc_S = filterPeakMatrix(peaksMatrix, Nuc_F_S, c(inert_columns, BC_columns), BC_Threshold_F, rowSum_Multiplier_F)
Peak_F_Cyt_M = filterPeakMatrix(peaksMatrix, Cyto_F_M, c(inert_columns, BC_columns), BC_Threshold_F, rowSum_Multiplier_F)
Peak_F_Cyt_S = filterPeakMatrix(peaksMatrix, Cyto_F_S, c(inert_columns, BC_columns), BC_Threshold_F, rowSum_Multiplier_F)

## CoCLIP
Peak_Co_Input_M = filterPeakMatrix(peaksMatrix, c(NLS_I_M, NES_I_M, G3BP_I_M), c(inert_columns, BC_columns), BC_Threshold_I, rowSum_Multiplier_I)
Peak_Co_Input_S = filterPeakMatrix(peaksMatrix, c(NLS_I_S, NES_I_S, G3BP_I_S), c(inert_columns, BC_columns), BC_Threshold_I, rowSum_Multiplier_I)
Peak_Co_NLS_M = filterPeakMatrix(peaksMatrix, NLS_E_M, c(inert_columns, BC_columns), BC_Threshold_E, rowSum_Multiplier_E)
Peak_Co_NLS_S = filterPeakMatrix(peaksMatrix, NLS_E_S, c(inert_columns, BC_columns), BC_Threshold_E, rowSum_Multiplier_E)
Peak_Co_NES_M = filterPeakMatrix(peaksMatrix, NES_E_M, c(inert_columns, BC_columns), BC_Threshold_E, rowSum_Multiplier_E)
Peak_Co_NES_S = filterPeakMatrix(peaksMatrix, NES_E_S, c(inert_columns, BC_columns), BC_Threshold_E, rowSum_Multiplier_E)
Peak_Co_G3BP_M = filterPeakMatrix(peaksMatrix, G3BP_E_M, c(inert_columns, BC_columns), BC_Threshold_E_SG, rowSum_Multiplier_E)
Peak_Co_G3BP_S = filterPeakMatrix(peaksMatrix, G3BP_E_S, c(inert_columns, BC_columns), BC_Threshold_E_SG, rowSum_Multiplier_E)
####################################################################################################################

## Global Variables:
####################################################################################################################
CoCLIP_List = c('Co_M_Input', 'Co_S_Input', 'Co_M_NLS', 'Co_S_NLS', 'Co_M_NES', 'Co_S_NES', 'Co_M_G3BP', 'Co_S_G3BP')
F_CoCLIP_List = c('F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES', 'F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES')

All_Annotation_List = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", 'ncRNA', "TE", "Other", "DS10K", 'UnAn')
mRNA_List = c("5'UTR", "CDS", "3'UTR", "intron", 'CDS_RI', 'DS10K')
ncRNA_List = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI')
####################################################################################################################

## FIGURE1/5 Mock vs Stress Stacked Bar Plots For CoCLIP Peaks:
####################################################################################################################
## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_Co_Input_M = countAnnotation(Peak_Co_Input_M, 'grouped_annotation', 'Co_M_Input', 'UnAn')
PeakDistribution_Co_Input_S = countAnnotation(Peak_Co_Input_S, 'grouped_annotation', 'Co_S_Input', 'UnAn')
PeakDistribution_Co_NLS_M = countAnnotation(Peak_Co_NLS_M, 'grouped_annotation', 'Co_M_NLS', 'UnAn')
PeakDistribution_Co_NLS_S = countAnnotation(Peak_Co_NLS_S, 'grouped_annotation', 'Co_S_NLS', 'UnAn')
PeakDistribution_Co_NES_M = countAnnotation(Peak_Co_NES_M, 'grouped_annotation', 'Co_M_NES', 'UnAn')
PeakDistribution_Co_NES_S = countAnnotation(Peak_Co_NES_S, 'grouped_annotation', 'Co_S_NES', 'UnAn')
PeakDistribution_Co_G3BP_M = countAnnotation(Peak_Co_G3BP_M, 'grouped_annotation', 'Co_M_G3BP', 'UnAn')
PeakDistribution_Co_G3BP_S = countAnnotation(Peak_Co_G3BP_S, 'grouped_annotation', 'Co_S_G3BP', 'UnAn')

PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
                                     PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
                                     PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
                                     PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>% gather(key = "Source", value = "Freq", CoCLIP_List) %>% select(Source, Freq, Annotation)
PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = CoCLIP_List)
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = All_Annotation_List)

plotStackedBar(PeakDistribution_Co_combined, c('Co_M_Input', 'Co_M_NLS', 'Co_M_NES'), c('Input', 'NLS', 'NES'), 'CoCLIP Mock Peak Distribution', y_lim = c(0, 10000), y_tick = 2000) 
plotStackedBar(PeakDistribution_Co_combined, c('Co_M_Input', 'Co_M_G3BP'), c('Input', 'G3BP'), 'CoCLIP Mock Peak Distribution', y_lim = c(0, 30000), y_tick = 5000)
plotStackedBar(PeakDistribution_Co_combined, c('Co_S_Input', 'Co_S_G3BP'), c('Input', 'G3BP'), 'CoCLIP Arsenite Peak Distribution', y_lim = c(0, 30000), y_tick = 5000) 

## Fraction Distribution Stacked Bar Graph
PeakDistribution_Co_Input_M = countAnnotation(Peak_Co_Input_M, 'grouped_annotation', 'Co_M_Input', 'UnAn', fraction = TRUE)
PeakDistribution_Co_Input_S = countAnnotation(Peak_Co_Input_S, 'grouped_annotation', 'Co_S_Input', 'UnAn', fraction = TRUE)
PeakDistribution_Co_NLS_M = countAnnotation(Peak_Co_NLS_M, 'grouped_annotation', 'Co_M_NLS', 'UnAn', fraction = TRUE)
PeakDistribution_Co_NLS_S = countAnnotation(Peak_Co_NLS_S, 'grouped_annotation', 'Co_S_NLS', 'UnAn', fraction = TRUE)
PeakDistribution_Co_NES_M = countAnnotation(Peak_Co_NES_M, 'grouped_annotation', 'Co_M_NES', 'UnAn', fraction = TRUE)
PeakDistribution_Co_NES_S = countAnnotation(Peak_Co_NES_S, 'grouped_annotation', 'Co_S_NES', 'UnAn', fraction = TRUE)
PeakDistribution_Co_G3BP_M = countAnnotation(Peak_Co_G3BP_M, 'grouped_annotation', 'Co_M_G3BP', 'UnAn', fraction = TRUE)
PeakDistribution_Co_G3BP_S = countAnnotation(Peak_Co_G3BP_S, 'grouped_annotation', 'Co_S_G3BP', 'UnAn', fraction = TRUE)

PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
                                     PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
                                     PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
                                     PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>% gather(key = "Source", value = "Freq", CoCLIP_List) %>% select(Source, Freq, Annotation)
PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = CoCLIP_List)
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = All_Annotation_List)

plotStackedBar(PeakDistribution_Co_combined, c('Co_M_Input', 'Co_M_NLS', 'Co_M_NES'), c('Input', 'NLS', 'NES'), 'CoCLIP Mock Peak Distribution By Fraction') 
plotStackedBar(PeakDistribution_Co_combined, c('Co_M_Input', 'Co_M_G3BP'), c('Input', 'G3BP'), 'CoCLIP Mock Peak Distribution By Fraction') 
plotStackedBar(PeakDistribution_Co_combined, c('Co_S_Input', 'Co_S_G3BP'), c('Input', 'G3BP'), 'CoCLIP Arsenite Peak Distribution By Fraction') 

####################################################################################################################

## FIGURE1/5 Mock vs Stress Stacked Bar Plots For CoCLIP mRNA Peaks:
####################################################################################################################
## Filter to mRNA Peaks
Peak_Co_Input_M_mRNA = filterPeaksByAnnotation(Peak_Co_Input_M, finalized_annotation, mRNA_List)
Peak_Co_Input_S_mRNA = filterPeaksByAnnotation(Peak_Co_Input_S, finalized_annotation, mRNA_List)
Peak_Co_NLS_M_mRNA = filterPeaksByAnnotation(Peak_Co_NLS_M, finalized_annotation, mRNA_List)
Peak_Co_NLS_S_mRNA = filterPeaksByAnnotation(Peak_Co_NLS_S, finalized_annotation, mRNA_List)
Peak_Co_NES_M_mRNA = filterPeaksByAnnotation(Peak_Co_NES_M, finalized_annotation, mRNA_List)
Peak_Co_NES_S_mRNA = filterPeaksByAnnotation(Peak_Co_NES_S, finalized_annotation, mRNA_List)
Peak_Co_G3BP_M_mRNA = filterPeaksByAnnotation(Peak_Co_G3BP_M, finalized_annotation, mRNA_List)
Peak_Co_G3BP_S_mRNA = filterPeaksByAnnotation(Peak_Co_G3BP_S, finalized_annotation, mRNA_List)

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_Co_Input_M = countAnnotation(Peak_Co_Input_M_mRNA, 'finalized_annotation', 'Co_M_Input', 'UnAn')
PeakDistribution_Co_Input_S = countAnnotation(Peak_Co_Input_S_mRNA, 'finalized_annotation', 'Co_S_Input', 'UnAn')
PeakDistribution_Co_NLS_M = countAnnotation(Peak_Co_NLS_M_mRNA, 'finalized_annotation', 'Co_M_NLS', 'UnAn')
PeakDistribution_Co_NLS_S = countAnnotation(Peak_Co_NLS_S_mRNA, 'finalized_annotation', 'Co_S_NLS', 'UnAn')
PeakDistribution_Co_NES_M = countAnnotation(Peak_Co_NES_M_mRNA, 'finalized_annotation', 'Co_M_NES', 'UnAn')
PeakDistribution_Co_NES_S = countAnnotation(Peak_Co_NES_S_mRNA, 'finalized_annotation', 'Co_S_NES', 'UnAn')
PeakDistribution_Co_G3BP_M = countAnnotation(Peak_Co_G3BP_M_mRNA, 'finalized_annotation', 'Co_M_G3BP', 'UnAn')
PeakDistribution_Co_G3BP_S = countAnnotation(Peak_Co_G3BP_S_mRNA, 'finalized_annotation', 'Co_S_G3BP', 'UnAn')

PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
                                     PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
                                     PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
                                     PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>% gather(key = "Source", value = "Freq", CoCLIP_List) %>% select(Source, Freq, Annotation)
PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = CoCLIP_List)
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = mRNA_List)

## break 2000 up to 10K
plotStackedBar(PeakDistribution_Co_combined, c('Co_M_Input', 'Co_M_NLS', 'Co_M_NES'), c('Input', 'NLS', 'NES'), 'CoCLIP Mock mRNA Peak Distribution', y_lim = c(0, 10000), y_tick = 2000) 
plotStackedBar(PeakDistribution_Co_combined, c('Co_M_Input', 'Co_M_G3BP'), c('Input', 'G3BP'), 'CoCLIP Mock mRNA Peak Distribution', y_lim = c(0, 25000), y_tick = 5000) 
plotStackedBar(PeakDistribution_Co_combined, c('Co_S_Input', 'Co_S_G3BP'), c('Input', 'G3BP'), 'CoCLIP Arsenite mRNA Peak Distribution', y_lim = c(0, 25000), y_tick = 5000) 

## Fraction Distribution Stacked Bar Graph
PeakDistribution_Co_Input_M = countAnnotation(Peak_Co_Input_M_mRNA, 'finalized_annotation', 'Co_M_Input', 'UnAn', fraction = T)
PeakDistribution_Co_Input_S = countAnnotation(Peak_Co_Input_S_mRNA, 'finalized_annotation', 'Co_S_Input', 'UnAn', fraction = T)
PeakDistribution_Co_NLS_M = countAnnotation(Peak_Co_NLS_M_mRNA, 'finalized_annotation', 'Co_M_NLS', 'UnAn', fraction = T)
PeakDistribution_Co_NLS_S = countAnnotation(Peak_Co_NLS_S_mRNA, 'finalized_annotation', 'Co_S_NLS', 'UnAn', fraction = T)
PeakDistribution_Co_NES_M = countAnnotation(Peak_Co_NES_M_mRNA, 'finalized_annotation', 'Co_M_NES', 'UnAn', fraction = T)
PeakDistribution_Co_NES_S = countAnnotation(Peak_Co_NES_S_mRNA, 'finalized_annotation', 'Co_S_NES', 'UnAn', fraction = T)
PeakDistribution_Co_G3BP_M = countAnnotation(Peak_Co_G3BP_M_mRNA, 'finalized_annotation', 'Co_M_G3BP', 'UnAn', fraction = T)
PeakDistribution_Co_G3BP_S = countAnnotation(Peak_Co_G3BP_S_mRNA, 'finalized_annotation', 'Co_S_G3BP', 'UnAn', fraction = T)

PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
                                     PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
                                     PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
                                     PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>% gather(key = "Source", value = "Freq", CoCLIP_List) %>% select(Source, Freq, Annotation)
PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = CoCLIP_List)
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = mRNA_List)

plotStackedBar(PeakDistribution_Co_combined, c('Co_M_Input', 'Co_M_NLS', 'Co_M_NES'), c('Input', 'NLS', 'NES'), 'CoCLIP Mock mRNA Peak Distribution by Fraction') 
plotStackedBar(PeakDistribution_Co_combined, c('Co_M_Input', 'Co_M_G3BP'), c('Input', 'G3BP'), 'CoCLIP Mock mRNA Peak Distribution by Fraction') 
plotStackedBar(PeakDistribution_Co_combined, c('Co_S_Input', 'Co_S_G3BP'), c('Input', 'G3BP'), 'CoCLIP Arsenite mRNA Peak Distribution by Fraction') 
####################################################################################################################

## FIGURE1/5 Mock vs Stress Stacked Bar Plots For CoCLIP ncRNA Peaks:
####################################################################################################################
## Filter to ncRNA Peaks
Peak_Co_Input_M_ncRNA = filterPeaksByAnnotation(Peak_Co_Input_M, finalized_annotation, ncRNA_List)
Peak_Co_Input_S_ncRNA = filterPeaksByAnnotation(Peak_Co_Input_S, finalized_annotation, ncRNA_List)
Peak_Co_NLS_M_ncRNA = filterPeaksByAnnotation(Peak_Co_NLS_M, finalized_annotation, ncRNA_List)
Peak_Co_NLS_S_ncRNA = filterPeaksByAnnotation(Peak_Co_NLS_S, finalized_annotation, ncRNA_List)
Peak_Co_NES_M_ncRNA = filterPeaksByAnnotation(Peak_Co_NES_M, finalized_annotation, ncRNA_List)
Peak_Co_NES_S_ncRNA = filterPeaksByAnnotation(Peak_Co_NES_S, finalized_annotation, ncRNA_List)
Peak_Co_G3BP_M_ncRNA = filterPeaksByAnnotation(Peak_Co_G3BP_M, finalized_annotation, ncRNA_List)
Peak_Co_G3BP_S_ncRNA = filterPeaksByAnnotation(Peak_Co_G3BP_S, finalized_annotation, ncRNA_List)

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_Co_Input_M = countAnnotation(Peak_Co_Input_M_ncRNA, 'finalized_annotation', 'Co_M_Input', 'UnAn')
PeakDistribution_Co_Input_S = countAnnotation(Peak_Co_Input_S_ncRNA, 'finalized_annotation', 'Co_S_Input', 'UnAn')
PeakDistribution_Co_NLS_M = countAnnotation(Peak_Co_NLS_M_ncRNA, 'finalized_annotation', 'Co_M_NLS', 'UnAn')
PeakDistribution_Co_NLS_S = countAnnotation(Peak_Co_NLS_S_ncRNA, 'finalized_annotation', 'Co_S_NLS', 'UnAn')
PeakDistribution_Co_NES_M = countAnnotation(Peak_Co_NES_M_ncRNA, 'finalized_annotation', 'Co_M_NES', 'UnAn')
PeakDistribution_Co_NES_S = countAnnotation(Peak_Co_NES_S_ncRNA, 'finalized_annotation', 'Co_S_NES', 'UnAn')
PeakDistribution_Co_G3BP_M = countAnnotation(Peak_Co_G3BP_M_ncRNA, 'finalized_annotation', 'Co_M_G3BP', 'UnAn')
PeakDistribution_Co_G3BP_S = countAnnotation(Peak_Co_G3BP_S_ncRNA, 'finalized_annotation', 'Co_S_G3BP', 'UnAn')

PeakDistribution_Co_Input_M = fillAnnotation(PeakDistribution_Co_Input_M, ncRNA_List)
PeakDistribution_Co_Input_S = fillAnnotation(PeakDistribution_Co_Input_S, ncRNA_List)
PeakDistribution_Co_NLS_M = fillAnnotation(PeakDistribution_Co_NLS_M, ncRNA_List)
PeakDistribution_Co_NLS_S = fillAnnotation(PeakDistribution_Co_NLS_S, ncRNA_List)
PeakDistribution_Co_NES_M = fillAnnotation(PeakDistribution_Co_NES_M, ncRNA_List)
PeakDistribution_Co_NES_S = fillAnnotation(PeakDistribution_Co_NES_S, ncRNA_List)
PeakDistribution_Co_G3BP_M = fillAnnotation(PeakDistribution_Co_G3BP_M, ncRNA_List)
PeakDistribution_Co_G3BP_S = fillAnnotation(PeakDistribution_Co_G3BP_S, ncRNA_List)

PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
                                     PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
                                     PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
                                     PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>% gather(key = "Source", value = "Freq", CoCLIP_List) %>% select(Source, Freq, Annotation)
PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = CoCLIP_List)
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = ncRNA_List)

plotStackedBar(PeakDistribution_Co_combined, c('Co_M_Input', 'Co_M_NLS', 'Co_M_NES'), c('Input', 'NLS', 'NES'), 'CoCLIP Mock ncRNA Peak Distribution')
plotStackedBar(PeakDistribution_Co_combined, c('Co_M_Input', 'Co_M_G3BP'), c('Input', 'G3BP'), 'CoCLIP Mock ncRNA Peak Distribution', y_lim = c(0, 6000), y_tick = 2000) 
plotStackedBar(PeakDistribution_Co_combined, c('Co_S_Input', 'Co_S_G3BP'), c('Input', 'G3BP'), 'CoCLIP Arsenite ncRNA Peak Distribution', y_lim = c(0, 6000), y_tick = 2000) 

## Fraction Distribution Stacked Bar Graph
PeakDistribution_Co_Input_M = countAnnotation(Peak_Co_Input_M_ncRNA, 'finalized_annotation', 'Co_M_Input', 'UnAn', fraction = T)
PeakDistribution_Co_Input_S = countAnnotation(Peak_Co_Input_S_ncRNA, 'finalized_annotation', 'Co_S_Input', 'UnAn', fraction = T)
PeakDistribution_Co_NLS_M = countAnnotation(Peak_Co_NLS_M_ncRNA, 'finalized_annotation', 'Co_M_NLS', 'UnAn', fraction = T)
PeakDistribution_Co_NLS_S = countAnnotation(Peak_Co_NLS_S_ncRNA, 'finalized_annotation', 'Co_S_NLS', 'UnAn', fraction = T)
PeakDistribution_Co_NES_M = countAnnotation(Peak_Co_NES_M_ncRNA, 'finalized_annotation', 'Co_M_NES', 'UnAn', fraction = T)
PeakDistribution_Co_NES_S = countAnnotation(Peak_Co_NES_S_ncRNA, 'finalized_annotation', 'Co_S_NES', 'UnAn', fraction = T)
PeakDistribution_Co_G3BP_M = countAnnotation(Peak_Co_G3BP_M_ncRNA, 'finalized_annotation', 'Co_M_G3BP', 'UnAn', fraction = T)
PeakDistribution_Co_G3BP_S = countAnnotation(Peak_Co_G3BP_S_ncRNA, 'finalized_annotation', 'Co_S_G3BP', 'UnAn', fraction = T)

PeakDistribution_Co_Input_M = fillAnnotation(PeakDistribution_Co_Input_M, ncRNA_List)
PeakDistribution_Co_Input_S = fillAnnotation(PeakDistribution_Co_Input_S, ncRNA_List)
PeakDistribution_Co_NLS_M = fillAnnotation(PeakDistribution_Co_NLS_M, ncRNA_List)
PeakDistribution_Co_NLS_S = fillAnnotation(PeakDistribution_Co_NLS_S, ncRNA_List)
PeakDistribution_Co_NES_M = fillAnnotation(PeakDistribution_Co_NES_M, ncRNA_List)
PeakDistribution_Co_NES_S = fillAnnotation(PeakDistribution_Co_NES_S, ncRNA_List)
PeakDistribution_Co_G3BP_M = fillAnnotation(PeakDistribution_Co_G3BP_M, ncRNA_List)
PeakDistribution_Co_G3BP_S = fillAnnotation(PeakDistribution_Co_G3BP_S, ncRNA_List)

PeakDistribution_Co_combined = cbind(PeakDistribution_Co_Input_M, PeakDistribution_Co_Input_S, 
                                     PeakDistribution_Co_NLS_M, PeakDistribution_Co_NLS_S, 
                                     PeakDistribution_Co_NES_M, PeakDistribution_Co_NES_S,
                                     PeakDistribution_Co_G3BP_M, PeakDistribution_Co_G3BP_S)
PeakDistribution_Co_combined$Annotation = rownames(PeakDistribution_Co_combined)

PeakDistribution_Co_combined = PeakDistribution_Co_combined %>% gather(key = "Source", value = "Freq", CoCLIP_List) %>% select(Source, Freq, Annotation)
PeakDistribution_Co_combined$Source = factor(PeakDistribution_Co_combined$Source, levels = CoCLIP_List)
PeakDistribution_Co_combined$Annotation = factor(PeakDistribution_Co_combined$Annotation, levels = ncRNA_List)

plotStackedBar(PeakDistribution_Co_combined, c('Co_M_Input', 'Co_M_NLS', 'Co_M_NES'), c('Input', 'NLS', 'NES'), 'CoCLIP Mock ncRNA Peak Distribution by Fraction')
plotStackedBar(PeakDistribution_Co_combined, c('Co_M_Input', 'Co_M_G3BP'), c('Input', 'G3BP'), 'CoCLIP Mock ncRNA Peak Distribution by Fraction') 
plotStackedBar(PeakDistribution_Co_combined, c('Co_S_Input', 'Co_S_G3BP'), c('Input', 'G3BP'), 'CoCLIP Arsenite ncRNA Peak Distribution by Fraction') 
####################################################################################################################

## FIGURE2 Stacked Bar Plots for Fraction vs CoCLIP Side-by-Side Peaks:
####################################################################################################################
## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_F_Nuc_M = countAnnotation(Peak_F_Nuc_M, 'grouped_annotation', 'F_M_Nuc', 'UnAn')
PeakDistribution_F_Cyt_M = countAnnotation(Peak_F_Cyt_M, 'grouped_annotation', 'F_M_Cyt', 'UnAn')
PeakDistribution_Co_NLS_M = countAnnotation(Peak_Co_NLS_M, 'grouped_annotation', 'Co_M_NLS', 'UnAn')
PeakDistribution_Co_NES_M = countAnnotation(Peak_Co_NES_M, 'grouped_annotation', 'Co_M_NES', 'UnAn')
PeakDistribution_F_Nuc_S = countAnnotation(Peak_F_Nuc_S, 'grouped_annotation', 'F_S_Nuc', 'UnAn')
PeakDistribution_F_Cyt_S = countAnnotation(Peak_F_Cyt_S, 'grouped_annotation', 'F_S_Cyt', 'UnAn')
PeakDistribution_Co_NLS_S = countAnnotation(Peak_Co_NLS_S, 'grouped_annotation', 'Co_S_NLS', 'UnAn')
PeakDistribution_Co_NES_S = countAnnotation(Peak_Co_NES_S, 'grouped_annotation', 'Co_S_NES', 'UnAn')

PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Cyt_M, 
                                  PeakDistribution_Co_NLS_M, PeakDistribution_Co_NES_M,
                                  PeakDistribution_F_Nuc_S, PeakDistribution_F_Cyt_S,
                                  PeakDistribution_Co_NLS_S, PeakDistribution_Co_NES_S)
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>% gather(key = "Source", value = "Freq", F_CoCLIP_List) %>% select(Source, Freq, Annotation)
PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = F_CoCLIP_List)
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = All_Annotation_List)

plotStackedBar(PeakDistribution_combined, 
               c('F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES'), 
               c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP'), 
               'Mock Peak Counts Distributions') 

plotStackedBar(PeakDistribution_combined, 
               c('F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'), 
               c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP'), 
               'Arsenite Peak Counts Distributions') 

## Fraction Distribution Stacked Bar Graph
PeakDistribution_F_Nuc_M = countAnnotation(Peak_F_Nuc_M, 'grouped_annotation', 'F_M_Nuc', 'UnAn', fraction = T)
PeakDistribution_F_Cyt_M = countAnnotation(Peak_F_Cyt_M, 'grouped_annotation', 'F_M_Cyt', 'UnAn', fraction = T)
PeakDistribution_Co_NLS_M = countAnnotation(Peak_Co_NLS_M, 'grouped_annotation', 'Co_M_NLS', 'UnAn', fraction = T)
PeakDistribution_Co_NES_M = countAnnotation(Peak_Co_NES_M, 'grouped_annotation', 'Co_M_NES', 'UnAn', fraction = T)
PeakDistribution_F_Nuc_S = countAnnotation(Peak_F_Nuc_S, 'grouped_annotation', 'F_S_Nuc', 'UnAn', fraction = T)
PeakDistribution_F_Cyt_S = countAnnotation(Peak_F_Cyt_S, 'grouped_annotation', 'F_S_Cyt', 'UnAn', fraction = T)
PeakDistribution_Co_NLS_S = countAnnotation(Peak_Co_NLS_S, 'grouped_annotation', 'Co_S_NLS', 'UnAn', fraction = T)
PeakDistribution_Co_NES_S = countAnnotation(Peak_Co_NES_S, 'grouped_annotation', 'Co_S_NES', 'UnAn', fraction = T)

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Cyt_M, 
                                  PeakDistribution_Co_NLS_M, PeakDistribution_Co_NES_M,
                                  PeakDistribution_F_Nuc_S, PeakDistribution_F_Cyt_S,
                                  PeakDistribution_Co_NLS_S, PeakDistribution_Co_NES_S)
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>% gather(key = "Source", value = "Freq", F_CoCLIP_List) %>% select(Source, Freq, Annotation)
PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = F_CoCLIP_List)
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = All_Annotation_List)

plotStackedBar(PeakDistribution_combined, 
               c('F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES'), 
               c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP'), 
               'Mock Peak Counts Distributions by Fraction') 

plotStackedBar(PeakDistribution_combined, 
               c('F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'), 
               c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP'), 
               'Arsenite Peak Counts Distributions by Fraction')

####################################################################################################################

## FIGURE2 Stacked Bar Plots for Fraction vs CoCLIP Side-by-Side mRNA Peaks
####################################################################################################################
## Filter to mRNA Peaks
Peak_F_Nuc_M_mRNA = filterPeaksByAnnotation(Peak_F_Nuc_M, finalized_annotation, mRNA_List)
Peak_F_Cyt_M_mRNA = filterPeaksByAnnotation(Peak_F_Cyt_M, finalized_annotation, mRNA_List)
Peak_Co_NLS_M_mRNA = filterPeaksByAnnotation(Peak_Co_NLS_M, finalized_annotation, mRNA_List)
Peak_Co_NES_M_mRNA = filterPeaksByAnnotation(Peak_Co_NES_M, finalized_annotation, mRNA_List)
Peak_F_Nuc_S_mRNA = filterPeaksByAnnotation(Peak_F_Nuc_S, finalized_annotation, mRNA_List)
Peak_F_Cyt_S_mRNA = filterPeaksByAnnotation(Peak_F_Cyt_S, finalized_annotation, mRNA_List)
Peak_Co_NLS_S_mRNA = filterPeaksByAnnotation(Peak_Co_NLS_S, finalized_annotation, mRNA_List)
Peak_Co_NES_S_mRNA = filterPeaksByAnnotation(Peak_Co_NES_S, finalized_annotation, mRNA_List)

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_F_Nuc_M = countAnnotation(Peak_F_Nuc_M_mRNA, 'finalized_annotation', 'F_M_Nuc', 'UnAn')
PeakDistribution_F_Cyt_M = countAnnotation(Peak_F_Cyt_M_mRNA, 'finalized_annotation', 'F_M_Cyt', 'UnAn')
PeakDistribution_Co_NLS_M = countAnnotation(Peak_Co_NLS_M_mRNA, 'finalized_annotation', 'Co_M_NLS', 'UnAn')
PeakDistribution_Co_NES_M = countAnnotation(Peak_Co_NES_M_mRNA, 'finalized_annotation', 'Co_M_NES', 'UnAn')
PeakDistribution_F_Nuc_S = countAnnotation(Peak_F_Nuc_S_mRNA, 'finalized_annotation', 'F_S_Nuc', 'UnAn')
PeakDistribution_F_Cyt_S = countAnnotation(Peak_F_Cyt_S_mRNA, 'finalized_annotation', 'F_S_Cyt', 'UnAn')
PeakDistribution_Co_NLS_S = countAnnotation(Peak_Co_NLS_S_mRNA, 'finalized_annotation', 'Co_S_NLS', 'UnAn')
PeakDistribution_Co_NES_S = countAnnotation(Peak_Co_NES_S_mRNA, 'finalized_annotation', 'Co_S_NES', 'UnAn')

PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Cyt_M, 
                                  PeakDistribution_Co_NLS_M, PeakDistribution_Co_NES_M,
                                  PeakDistribution_F_Nuc_S, PeakDistribution_F_Cyt_S,
                                  PeakDistribution_Co_NLS_S, PeakDistribution_Co_NES_S)
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>% gather(key = "Source", value = "Freq", F_CoCLIP_List) %>% select(Source, Freq, Annotation)
PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = F_CoCLIP_List)
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = mRNA_List)

plotStackedBar(PeakDistribution_combined, 
               c('F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES'), 
               c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP'), 
               'Mock mRNA Peak Counts Distributions') 

plotStackedBar(PeakDistribution_combined, 
               c('F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'), 
               c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP'), 
               'Arsenite mRNA Peak Counts Distributions') 

## Fraction Distribution Stacked Bar Graph
PeakDistribution_F_Nuc_M = countAnnotation(Peak_F_Nuc_M_mRNA, 'finalized_annotation', 'F_M_Nuc', 'UnAn', fraction = T)
PeakDistribution_F_Cyt_M = countAnnotation(Peak_F_Cyt_M_mRNA, 'finalized_annotation', 'F_M_Cyt', 'UnAn', fraction = T)
PeakDistribution_Co_NLS_M = countAnnotation(Peak_Co_NLS_M_mRNA, 'finalized_annotation', 'Co_M_NLS', 'UnAn', fraction = T)
PeakDistribution_Co_NES_M = countAnnotation(Peak_Co_NES_M_mRNA, 'finalized_annotation', 'Co_M_NES', 'UnAn', fraction = T)
PeakDistribution_F_Nuc_S = countAnnotation(Peak_F_Nuc_S_mRNA, 'finalized_annotation', 'F_S_Nuc', 'UnAn', fraction = T)
PeakDistribution_F_Cyt_S = countAnnotation(Peak_F_Cyt_S_mRNA, 'finalized_annotation', 'F_S_Cyt', 'UnAn', fraction = T)
PeakDistribution_Co_NLS_S = countAnnotation(Peak_Co_NLS_S_mRNA, 'finalized_annotation', 'Co_S_NLS', 'UnAn', fraction = T)
PeakDistribution_Co_NES_S = countAnnotation(Peak_Co_NES_S_mRNA, 'finalized_annotation', 'Co_S_NES', 'UnAn', fraction = T)

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Cyt_M, 
                                  PeakDistribution_Co_NLS_M, PeakDistribution_Co_NES_M,
                                  PeakDistribution_F_Nuc_S, PeakDistribution_F_Cyt_S,
                                  PeakDistribution_Co_NLS_S, PeakDistribution_Co_NES_S)
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>% gather(key = "Source", value = "Freq", F_CoCLIP_List) %>% select(Source, Freq, Annotation)
PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = F_CoCLIP_List)
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = mRNA_List)

plotStackedBar(PeakDistribution_combined, 
               c('F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES'), 
               c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP'), 
               'Mock mRNA Peak Counts Distributions by Fraction') 

plotStackedBar(PeakDistribution_combined, 
               c('F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'), 
               c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP'), 
               'Arsenite mRNA Peak Counts Distributions by Fraction')
####################################################################################################################

## FIGURE2 Stacked Bar Plots for Fraction vs CoCLIP Side-by-Side ncRNA Peaks
####################################################################################################################
## Filter to ncRNA Peaks
Peak_F_Nuc_M_ncRNA = filterPeaksByAnnotation(Peak_F_Nuc_M, finalized_annotation, ncRNA_List)
Peak_F_Cyt_M_ncRNA = filterPeaksByAnnotation(Peak_F_Cyt_M, finalized_annotation, ncRNA_List)
Peak_Co_NLS_M_ncRNA = filterPeaksByAnnotation(Peak_Co_NLS_M, finalized_annotation, ncRNA_List)
Peak_Co_NES_M_ncRNA = filterPeaksByAnnotation(Peak_Co_NES_M, finalized_annotation, ncRNA_List)
Peak_F_Nuc_S_ncRNA = filterPeaksByAnnotation(Peak_F_Nuc_S, finalized_annotation, ncRNA_List)
Peak_F_Cyt_S_ncRNA = filterPeaksByAnnotation(Peak_F_Cyt_S, finalized_annotation, ncRNA_List)
Peak_Co_NLS_S_ncRNA = filterPeaksByAnnotation(Peak_Co_NLS_S, finalized_annotation, ncRNA_List)
Peak_Co_NES_S_ncRNA = filterPeaksByAnnotation(Peak_Co_NES_S, finalized_annotation, ncRNA_List)

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_F_Nuc_M = countAnnotation(Peak_F_Nuc_M_ncRNA, 'finalized_annotation', 'F_M_Nuc', 'UnAn')
PeakDistribution_F_Cyt_M = countAnnotation(Peak_F_Cyt_M_ncRNA, 'finalized_annotation', 'F_M_Cyt', 'UnAn')
PeakDistribution_Co_NLS_M = countAnnotation(Peak_Co_NLS_M_ncRNA, 'finalized_annotation', 'Co_M_NLS', 'UnAn')
PeakDistribution_Co_NES_M = countAnnotation(Peak_Co_NES_M_ncRNA, 'finalized_annotation', 'Co_M_NES', 'UnAn')
PeakDistribution_F_Nuc_S = countAnnotation(Peak_F_Nuc_S_ncRNA, 'finalized_annotation', 'F_S_Nuc', 'UnAn')
PeakDistribution_F_Cyt_S = countAnnotation(Peak_F_Cyt_S_ncRNA, 'finalized_annotation', 'F_S_Cyt', 'UnAn')
PeakDistribution_Co_NLS_S = countAnnotation(Peak_Co_NLS_S_ncRNA, 'finalized_annotation', 'Co_S_NLS', 'UnAn')
PeakDistribution_Co_NES_S = countAnnotation(Peak_Co_NES_S_ncRNA, 'finalized_annotation', 'Co_S_NES', 'UnAn')

PeakDistribution_F_Nuc_M = fillAnnotation(PeakDistribution_F_Nuc_M, ncRNA_List)
PeakDistribution_F_Cyt_M = fillAnnotation(PeakDistribution_F_Cyt_M, ncRNA_List)
PeakDistribution_Co_NLS_M = fillAnnotation(PeakDistribution_Co_NLS_M, ncRNA_List)
PeakDistribution_Co_NES_M = fillAnnotation(PeakDistribution_Co_NES_M, ncRNA_List)
PeakDistribution_F_Nuc_S = fillAnnotation(PeakDistribution_F_Nuc_S, ncRNA_List)
PeakDistribution_F_Cyt_S = fillAnnotation(PeakDistribution_F_Cyt_S, ncRNA_List)
PeakDistribution_Co_NLS_S = fillAnnotation(PeakDistribution_Co_NLS_S, ncRNA_List)
PeakDistribution_Co_NES_S = fillAnnotation(PeakDistribution_Co_NES_S, ncRNA_List)

PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Cyt_M, 
                                  PeakDistribution_Co_NLS_M, PeakDistribution_Co_NES_M,
                                  PeakDistribution_F_Nuc_S, PeakDistribution_F_Cyt_S,
                                  PeakDistribution_Co_NLS_S, PeakDistribution_Co_NES_S)
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>% gather(key = "Source", value = "Freq", F_CoCLIP_List) %>% select(Source, Freq, Annotation)
PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = F_CoCLIP_List)
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = ncRNA_List)

plotStackedBar(PeakDistribution_combined, 
               c('F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES'), 
               c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP'), 
               'Mock ncRNA Peak Counts Distributions') 

plotStackedBar(PeakDistribution_combined, 
               c('F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'), 
               c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP'), 
               'Arsenite mRNA Peak Counts Distributions') 

## Fraction Distribution Stacked Bar Graph
PeakDistribution_F_Nuc_M = countAnnotation(Peak_F_Nuc_M_ncRNA, 'finalized_annotation', 'F_M_Nuc', 'UnAn', fraction = T)
PeakDistribution_F_Cyt_M = countAnnotation(Peak_F_Cyt_M_ncRNA, 'finalized_annotation', 'F_M_Cyt', 'UnAn', fraction = T)
PeakDistribution_Co_NLS_M = countAnnotation(Peak_Co_NLS_M_ncRNA, 'finalized_annotation', 'Co_M_NLS', 'UnAn', fraction = T)
PeakDistribution_Co_NES_M = countAnnotation(Peak_Co_NES_M_ncRNA, 'finalized_annotation', 'Co_M_NES', 'UnAn', fraction = T)
PeakDistribution_F_Nuc_S = countAnnotation(Peak_F_Nuc_S_ncRNA, 'finalized_annotation', 'F_S_Nuc', 'UnAn', fraction = T)
PeakDistribution_F_Cyt_S = countAnnotation(Peak_F_Cyt_S_ncRNA, 'finalized_annotation', 'F_S_Cyt', 'UnAn', fraction = T)
PeakDistribution_Co_NLS_S = countAnnotation(Peak_Co_NLS_S_ncRNA, 'finalized_annotation', 'Co_S_NLS', 'UnAn', fraction = T)
PeakDistribution_Co_NES_S = countAnnotation(Peak_Co_NES_S_ncRNA, 'finalized_annotation', 'Co_S_NES', 'UnAn', fraction = T)

PeakDistribution_F_Nuc_M = fillAnnotation(PeakDistribution_F_Nuc_M, ncRNA_List)
PeakDistribution_F_Cyt_M = fillAnnotation(PeakDistribution_F_Cyt_M, ncRNA_List)
PeakDistribution_Co_NLS_M = fillAnnotation(PeakDistribution_Co_NLS_M, ncRNA_List)
PeakDistribution_Co_NES_M = fillAnnotation(PeakDistribution_Co_NES_M, ncRNA_List)
PeakDistribution_F_Nuc_S = fillAnnotation(PeakDistribution_F_Nuc_S, ncRNA_List)
PeakDistribution_F_Cyt_S = fillAnnotation(PeakDistribution_F_Cyt_S, ncRNA_List)
PeakDistribution_Co_NLS_S = fillAnnotation(PeakDistribution_Co_NLS_S, ncRNA_List)
PeakDistribution_Co_NES_S = fillAnnotation(PeakDistribution_Co_NES_S, ncRNA_List)

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_combined = cbind(PeakDistribution_F_Nuc_M, PeakDistribution_F_Cyt_M, 
                                  PeakDistribution_Co_NLS_M, PeakDistribution_Co_NES_M,
                                  PeakDistribution_F_Nuc_S, PeakDistribution_F_Cyt_S,
                                  PeakDistribution_Co_NLS_S, PeakDistribution_Co_NES_S)
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>% gather(key = "Source", value = "Freq", F_CoCLIP_List) %>% select(Source, Freq, Annotation)
PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = F_CoCLIP_List)
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = ncRNA_List)

plotStackedBar(PeakDistribution_combined, 
               c('F_M_Nuc', 'F_M_Cyt', 'Co_M_NLS', 'Co_M_NES'), 
               c('Mock\nNuclear\nFracCLIP', 'Mock\nCytoplasm\nFracCLIP', 'Mock\nNLS\nCoCLIP', 'Mock\nNES\nCoCLIP'), 
               'Mock ncRNA Peak Counts Distributions by Fraction') 

plotStackedBar(PeakDistribution_combined, 
               c('F_S_Nuc', 'F_S_Cyt', 'Co_S_NLS', 'Co_S_NES'), 
               c('Arsenite\nNuclear\nFracCLIP', 'Arsenite\nCytoplasm\nFracCLIP', 'Arsenite\nNLS\nCoCLIP', 'Arsenite\nNES\nCoCLIP'), 
               'Arsenite mRNA Peak Counts Distributions by Fraction')

####################################################################################################################