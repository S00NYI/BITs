## CoCLIP Analysis: 
## Peak Processing for Enrichment Analysis
## Written by Soon Yi
## Last Edit: 2023-09-13

library(stringr)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)

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

## Get peak counts per gene: 
peaksPerGene = function(peak_matrix, sample_list) {
  temp = peak_matrix %>% group_by(gene, external_gene_name) %>% summarise(tagCounts = sum(across(all_of(sample_list))), peakCounts = n())
  temp = temp %>% mutate(tagDensity = tagCounts/peakCounts)
  
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
plotStackedBar = function(annotation_counts, sample_list, sample_label, title) {
  ggplot(annotation_counts %>% filter(Source %in% sample_list), aes(fill = Annotation, y=Freq, x=Source)) + 
    geom_bar(position='stack', stat='identity') +
    scale_x_discrete(labels = sample_label) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_brewer(palette = "Set3") +
    theme_bw() + 
    theme(axis.text = element_text(size=14), 
          axis.title = element_text(size=14, face = 'bold'), 
          legend.text = element_text(size=14))
}

## Plot Scatter 1:
plotScatter1 = function(peak_matrix, annotation_level, x_axis, y_axis, x_label, y_label, x_lim, y_lim, title) {
  ggplot(peak_matrix, aes(x = ({{x_axis}}), y = ({{y_axis}}), color = {{annotation_level}})) +
    geom_point(pch = 16, size = 3, alpha = 0.5) +
    labs(x = x_label, y = y_label) +
    xlim(x_lim) +
    ylim(y_lim) +
    ggtitle(title) +
    scale_fill_brewer(palette = "Set3") +
    theme_bw() + 
    theme(axis.text = element_text(size=14), 
          axis.title = element_text(size=14, face = 'bold'), 
          legend.text = element_text(size=14)) + 
    geom_abline(linetype = 'dotted')
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

All_Annotation_List = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", 'ncRNA', "TE", "Other", "DS10K")
mRNA_List = c("5'UTR", "CDS", "3'UTR", "intron", 'CDS_RI', 'DS10K')
ncRNA_List = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI')
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
Peak_Nuc_M = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                         (Nuc_F_M_BC >= BC_Threshold_F & Nuc_F_M > median(Nuc_F_M) * rowSum_Multiplier_F) &
                                         (NLS_E_M_BC >= BC_Threshold_E & NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_E))
Peak_Nuc_M$grouped_annotation = factor(Peak_Nuc_M$grouped_annotation, levels = All_Annotation_List)

plotScatter1(Peak_Nuc_M, grouped_annotation,
             log10(NLS_E_M), log10(Nuc_F_M), 
             'log10(NLS CoCLIP)', 'log10(Nuclear Fraction)',
             c(0, 4), c(0, 4), 
             paste0('Mock HuR Normalized Tag Counts Per Peaks: NLS CoCLIP vs Nuclear FracCLIP (',  nrow(Peak_Nuc_M), ' peaks)'))

# Nuclear Mock - mRNA specific:
Peak_Nuc_M_mRNA = Peak_Nuc_M %>% filter((finalized_annotation == "5'UTR" | 
                                           finalized_annotation == "3'UTR" | 
                                           finalized_annotation == "CDS" | 
                                           finalized_annotation == "intron" | 
                                           finalized_annotation == "CDS_RI" |
                                           finalized_annotation == "DS10K") &
                                          (grouped_annotation != 'UnAn'))
Peak_Nuc_M_mRNA$finalized_annotation = factor(Peak_Nuc_M_mRNA$finalized_annotation, levels = mRNA_List)

plotScatter1(Peak_Nuc_M_mRNA, finalized_annotation,
             log10(NLS_E_M), log10(Nuc_F_M), 
             'log10(NLS CoCLIP)', 'log10(Nuclear Fraction)',
             c(0, 4), c(0, 4), 
             paste0('Mock mRNA HuR Normalized Tag Counts Per Peaks: NLS CoCLIP vs Nuclear FracCLIP (',  nrow(Peak_Nuc_M_mRNA), ' peaks)'))

# Nuclear Mock - ncRNA specific:
Peak_Nuc_M_ncRNA = Peak_Nuc_M %>% filter((finalized_annotation != "5'UTR" & 
                                            finalized_annotation != "3'UTR" & 
                                            finalized_annotation != "CDS" & 
                                            finalized_annotation != "intron" &  
                                            finalized_annotation != "CDS_RI" &
                                            finalized_annotation != "DS10K" &
                                            finalized_annotation != "UnAn"))
Peak_Nuc_M_ncRNA$finalized_annotation = factor(Peak_Nuc_M_ncRNA$finalized_annotation, levels = ncRNA_List)

plotScatter1(Peak_Nuc_M_ncRNA, finalized_annotation,
             log10(NLS_E_M), log10(Nuc_F_M), 
             'log10(NLS CoCLIP)', 'log10(Nuclear Fraction)',
             c(0, 4), c(0, 4), 
             paste0('Mock ncRNA HuR Normalized Tag Counts Per Peaks: NLS CoCLIP vs Nuclear FracCLIP (',  nrow(Peak_Nuc_M_ncRNA), ' peaks)'))


# Nuclear Stress:
Peak_Nuc_S = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                         (Nuc_F_S_BC >= BC_Threshold_F & Nuc_F_S > median(Nuc_F_S) * rowSum_Multiplier_F) &
                                         (NLS_E_S_BC >= BC_Threshold_E & NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_E))
Peak_Nuc_S$grouped_annotation = factor(Peak_Nuc_S$grouped_annotation, levels = All_Annotation_List)

plotScatter1(Peak_Nuc_S, grouped_annotation,
             log10(NLS_E_S), log10(Nuc_F_S), 
             'log10(NLS CoCLIP)', 'log10(Nuclear Fraction)',
             c(0, 4), c(0, 4), 
             paste0('Stress HuR Normalized Tag Counts Per Peaks: NLS CoCLIP vs Nuclear FracCLIP (',  nrow(Peak_Nuc_M), ' peaks)'))

# Nuclear Stress - mRNA specific:
Peak_Nuc_S_mRNA = Peak_Nuc_S %>% filter((finalized_annotation == "5'UTR" | 
                                           finalized_annotation == "3'UTR" | 
                                           finalized_annotation == "CDS" | 
                                           finalized_annotation == "intron" | 
                                           finalized_annotation == "CDS_RI" |
                                           finalized_annotation == "DS10K") &
                                          (grouped_annotation != 'UnAn'))
Peak_Nuc_S_mRNA$finalized_annotation = factor(Peak_Nuc_S_mRNA$finalized_annotation, levels = mRNA_List)

plotScatter1(Peak_Nuc_S_mRNA, finalized_annotation,
             log10(NLS_E_S), log10(Nuc_F_S), 
             'log10(NLS CoCLIP)', 'log10(Nuclear Fraction)',
             c(0, 4), c(0, 4), 
             paste0('Stress mRNA HuR Normalized Tag Counts Per Peaks: NLS CoCLIP vs Nuclear FracCLIP (',  nrow(Peak_Nuc_M_mRNA), ' peaks)'))

# Nuclear Stress - ncRNA specific:
Peak_Nuc_S_ncRNA = Peak_Nuc_S %>% filter((finalized_annotation != "5'UTR" & 
                                            finalized_annotation != "3'UTR" & 
                                            finalized_annotation != "CDS" & 
                                            finalized_annotation != "intron" &  
                                            finalized_annotation != "CDS_RI" &
                                            finalized_annotation != "DS10K" &
                                            finalized_annotation != "UnAn"))
Peak_Nuc_S_ncRNA$finalized_annotation = factor(Peak_Nuc_S_ncRNA$finalized_annotation, levels = ncRNA_List)

plotScatter1(Peak_Nuc_S_ncRNA, finalized_annotation,
             log10(NLS_E_S), log10(Nuc_F_S), 
             'log10(NLS CoCLIP)', 'log10(Nuclear Fraction)',
             c(0, 4), c(0, 4), 
             paste0('Stress ncRNA HuR Normalized Tag Counts Per Peaks: NLS CoCLIP vs Nuclear FracCLIP (',  nrow(Peak_Nuc_M_ncRNA), ' peaks)'))

####################################################################################################################

## FIGURE2 Scatter Plot of Cytoplasm Fraction vs CoCLIP:
####################################################################################################################
# Cytoplasm Mock:
Peak_Cyt_M = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                         (Cyto_F_M_BC >= BC_Threshold_F & Cyto_F_M > median(Cyto_F_M) * rowSum_Multiplier_F) &
                                         (NES_E_M_BC >= BC_Threshold_E & NES_E_M > median(NES_E_M) * rowSum_Multiplier_E))
Peak_Cyt_M$grouped_annotation = factor(Peak_Cyt_M$grouped_annotation, levels = All_Annotation_List)

plotScatter1(Peak_Cyt_M, grouped_annotation,
             log10(NES_E_M), log10(Cyto_F_M), 
             'log10(NES CoCLIP)', 'log10(Cytoplasm Fraction)',
             c(0, 4), c(0, 4), 
             paste0('Mock HuR Normalized Tag Counts Per Peaks: NES CoCLIP vs Cytoplasm FracCLIP (',  nrow(Peak_Cyt_M), ' peaks)'))

# Cytoplasm Mock - mRNA specific:
Peak_Cyt_M_mRNA = Peak_Cyt_M %>% filter((finalized_annotation == "5'UTR" | 
                                           finalized_annotation == "3'UTR" | 
                                           finalized_annotation == "CDS" | 
                                           finalized_annotation == "intron" | 
                                           finalized_annotation == "CDS_RI" |
                                           finalized_annotation == "DS10K") &
                                          (grouped_annotation != 'UnAn'))
Peak_Cyt_M_mRNA$finalized_annotation = factor(Peak_Cyt_M_mRNA$finalized_annotation, levels = mRNA_List)

plotScatter1(Peak_Cyt_M_mRNA, finalized_annotation,
             log10(NES_E_M), log10(Cyto_F_M), 
             'log10(NES CoCLIP)', 'log10(Cytoplasm Fraction)',
             c(0, 4), c(0, 4), 
             paste0('Mock mRNA HuR Normalized Tag Counts Per Peaks: NES CoCLIP vs Cytoplasm FracCLIP (',  nrow(Peak_Cyt_M_mRNA), ' peaks)'))

# Cytoplasm Mock - ncRNA specific:
Peak_Cyt_M_ncRNA = Peak_Cyt_M %>% filter((finalized_annotation != "5'UTR" & 
                                            finalized_annotation != "3'UTR" & 
                                            finalized_annotation != "CDS" & 
                                            finalized_annotation != "intron" &  
                                            finalized_annotation != "CDS_RI" &
                                            finalized_annotation != "DS10K" &
                                            finalized_annotation != "UnAn"))
Peak_Cyt_M_ncRNA$finalized_annotation = factor(Peak_Cyt_M_ncRNA$finalized_annotation, levels = ncRNA_List)

plotScatter1(Peak_Cyt_M_ncRNA, finalized_annotation,
             log10(NES_E_M), log10(Cyto_F_M), 
             'log10(NES CoCLIP)', 'log10(Cytoplasm Fraction)',
             c(0, 4), c(0, 4), 
             paste0('Mock ncRNA HuR Normalized Tag Counts Per Peaks: NES CoCLIP vs Cytoplasm FracCLIP (',  nrow(Peak_Cyt_M_ncRNA), ' peaks)'))


# Cytoplasm Stress:
Peak_Cyt_S = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                         (Cyto_F_S_BC >= BC_Threshold_F & Cyto_F_S > median(Cyto_F_S) * rowSum_Multiplier_F) &
                                         (NES_E_S_BC >= BC_Threshold_E & NES_E_S > median(NES_E_S) * rowSum_Multiplier_E))
Peak_Cyt_S$grouped_annotation = factor(Peak_Cyt_S$grouped_annotation, levels = All_Annotation_List)

plotScatter1(Peak_Cyt_S, grouped_annotation, 
             log10(NES_E_S), log10(Cyto_F_S), 
             'log10(NES CoCLIP)', 'log10(Cytoplasm Fraction)',
             c(0, 4), c(0, 4), 
             paste0('Stress HuR Normalized Tag Counts Per Peaks: NES CoCLIP vs Cytoplasm FracCLIP (',  nrow(Peak_Cyt_M), ' peaks)'))

# Cytoplasm Stress - mRNA specific:
Peak_Cyt_S_mRNA = Peak_Cyt_S %>% filter((finalized_annotation == "5'UTR" | 
                                           finalized_annotation == "3'UTR" | 
                                           finalized_annotation == "CDS" | 
                                           finalized_annotation == "intron" | 
                                           finalized_annotation == "CDS_RI" |
                                           finalized_annotation == "DS10K") &
                                          (grouped_annotation != 'UnAn'))
Peak_Cyt_S_mRNA$finalized_annotation = factor(Peak_Cyt_S_mRNA$finalized_annotation, levels = mRNA_List)

plotScatter1(Peak_Cyt_S_mRNA, finalized_annotation,
             log10(NES_E_S), log10(Cyto_F_S), 
             'log10(NES CoCLIP)', 'log10(Cytoplasm Fraction)',
             c(0, 4), c(0, 4), 
             paste0('Stress mRNA HuR Normalized Tag Counts Per Peaks: NES CoCLIP vs Cytoplasm FracCLIP (',  nrow(Peak_Cyt_S_mRNA), ' peaks)'))

# Cytoplasm Stress - ncRNA specific:
Peak_Cyt_S_ncRNA = Peak_Cyt_S %>% filter((finalized_annotation != "5'UTR" & 
                                            finalized_annotation != "3'UTR" & 
                                            finalized_annotation != "CDS" & 
                                            finalized_annotation != "intron" &  
                                            finalized_annotation != "CDS_RI" &
                                            finalized_annotation != "DS10K" &
                                            finalized_annotation != "UnAn"))
Peak_Cyt_S_ncRNA$finalized_annotation = factor(Peak_Cyt_S_ncRNA$finalized_annotation, levels = ncRNA_List)

plotScatter1(Peak_Cyt_S_ncRNA, finalized_annotation,
             log10(NES_E_S), log10(Cyto_F_S), 
             'log10(NES CoCLIP)', 'log10(Cytoplasm Fraction)',
             c(0, 4), c(0, 4), 
             paste0('Stress ncRNA HuR Normalized Tag Counts Per Peaks: NES CoCLIP vs Cytoplasm FracCLIP (',  nrow(Peak_Cyt_S_ncRNA), ' peaks)'))
####################################################################################################################

## FIGURE2 Scatter Plot of Nuclear Fold Change Over Input Fraction vs CoCLIP:
####################################################################################################################
x_lim = c(-6, 6)
y_lim = x_lim

# Nuclear:
Peak_Nuc_FC = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                            ((Nuc_F_M_BC >= BC_Threshold_F & Nuc_F_M > median(Nuc_F_M) * rowSum_Multiplier_F &
                                            NLS_E_M_BC >= BC_Threshold_E & NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_E &
                                            (NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I & I_M > median(I_M) * rowSum_Multiplier_I) &
                                            (Nuc_F_S_BC >= BC_Threshold_F & Nuc_F_S > median(Nuc_F_S) * rowSum_Multiplier_F &
                                               NLS_E_S_BC >= BC_Threshold_E & NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_E &
                                               (NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I & I_S > median(I_S) * rowSum_Multiplier_I)))

Peak_Nuc_FC$grouped_annotation = factor(Peak_Nuc_FC$grouped_annotation, levels = All_Annotation_List)

plotScatter1(Peak_Nuc_FC, grouped_annotation,
             log2(NLS_EvI_M), log2(Nuc_EvI_M),
             'log2(CoCLIP Nuclear/Input)', 'log2(FracCLIP Nuclear/Input)',
             x_lim, y_lim, 
             paste0('Mock HuR Peak Enrichment over Input: Nuclear CoCLIP vs FracCLIP (',  nrow(Peak_Nuc_FC), ' peaks)'))

plotScatter1(Peak_Nuc_FC, grouped_annotation,
             log2(NLS_EvI_S), log2(Nuc_EvI_S),
             'log2(CoCLIP Nuclear/Input)', 'log2(FracCLIP Nuclear/Input)',
             x_lim, y_lim, 
             paste0('Stress HuR Peak Enrichment over Input: Nuclear CoCLIP vs FracCLIP (',  nrow(Peak_Nuc_FC), ' peaks)'))

# Nuclear - mRNA specific:
Peak_Nuc_FC_mRNA = Peak_Nuc_FC %>% filter((finalized_annotation == "5'UTR" | 
                                                  finalized_annotation == "3'UTR" | 
                                                  finalized_annotation == "CDS" | 
                                                  finalized_annotation == "intron" | 
                                                  finalized_annotation == "CDS_RI" |
                                                  finalized_annotation == "DS10K"))

Peak_Nuc_FC_mRNA$finalized_annotation = factor(Peak_Nuc_FC_mRNA$finalized_annotation, levels = mRNA_List)

plotScatter1(Peak_Nuc_FC_mRNA, finalized_annotation,
             log2(NLS_EvI_M), log2(Nuc_EvI_M),
             'log2(CoCLIP Nuclear/Input)', 'log2(FracCLIP Nuclear/Input)',
             x_lim, y_lim, 
             paste0('Mock HuR mRNA Peak Enrichment over Input: Nuclear CoCLIP vs FracCLIP (',  nrow(Peak_Nuc_FC_mRNA), ' peaks)'))

plotScatter1(Peak_Nuc_FC_mRNA, finalized_annotation,
             log2(NLS_EvI_S), log2(Nuc_EvI_S),
             'log2(CoCLIP Nuclear/Input)', 'log2(FracCLIP Nuclear/Input)',
             x_lim, y_lim, 
             paste0('Stress HuR mRNA Peak Enrichment over Input: Nuclear CoCLIP vs FracCLIP (',  nrow(Peak_Nuc_FC_mRNA), ' peaks)'))

# Nuclear - ncRNA specific:
Peak_Nuc_FC_ncRNA = Peak_Nuc_FC %>% filter((finalized_annotation != "5'UTR" & 
                                                  finalized_annotation != "3'UTR" & 
                                                  finalized_annotation != "CDS" & 
                                                  finalized_annotation != "intron" &  
                                                  finalized_annotation != "CDS_RI" &
                                                  finalized_annotation != "DS10K" &
                                                  finalized_annotation != "UnAn"))

Peak_Nuc_FC_ncRNA$finalized_annotation = factor(Peak_Nuc_FC_ncRNA$finalized_annotation, levels = ncRNA_List)

plotScatter1(Peak_Nuc_FC_ncRNA, finalized_annotation,
             log2(NLS_EvI_M), log2(Nuc_EvI_M),
             'log2(CoCLIP Nuclear/Input)', 'log2(FracCLIP Nuclear/Input)',
             x_lim, y_lim, 
             paste0('Mock HuR ncRNA Peak Enrichment over Input: Nuclear CoCLIP vs FracCLIP (',  nrow(Peak_Nuc_FC_ncRNA), ' peaks)'))

plotScatter1(Peak_Nuc_FC_ncRNA, finalized_annotation,
             log2(NLS_EvI_S), log2(Nuc_EvI_S),
             'log2(CoCLIP Nuclear/Input)', 'log2(FracCLIP Nuclear/Input)',
             x_lim, y_lim, 
             paste0('Stress HuR ncRNA Peak Enrichment over Input: Nuclear CoCLIP vs FracCLIP (',  nrow(Peak_Nuc_FC_ncRNA), ' peaks)'))
####################################################################################################################

## FIGURE2 Scatter Plot of Cytoplasm Fold Change Over Input Fraction vs CoCLIP:
####################################################################################################################
x_lim = c(-8, 8)
y_lim = x_lim

# Cytoplasm:
Peak_Cyt_FC = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                          ((Cyto_F_M_BC >= BC_Threshold_F & Cyto_F_M > median(Cyto_F_M) * rowSum_Multiplier_F &
                                              NES_E_M_BC >= BC_Threshold_E & NES_E_M > median(NES_E_M) * rowSum_Multiplier_E &
                                              (NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I & I_M > median(I_M) * rowSum_Multiplier_I) &
                                             (Cyto_F_S_BC >= BC_Threshold_F & Cyto_F_S > median(Cyto_F_S) * rowSum_Multiplier_F &
                                                NES_E_S_BC >= BC_Threshold_E & NES_E_S > median(NES_E_S) * rowSum_Multiplier_E &
                                                (NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I & I_S > median(I_S) * rowSum_Multiplier_I)))

Peak_Cyt_FC$grouped_annotation = factor(Peak_Cyt_FC$grouped_annotation, levels = All_Annotation_List)

plotScatter1(Peak_Cyt_FC, grouped_annotation,
             log2(NES_EvI_M), log2(Cyto_EvI_M),
             'log2(CoCLIP Cytoplasm/Input)', 'log2(FracCLIP Cytoplasm/Input)',
             x_lim, y_lim, 
             paste0('Mock HuR Peak Enrichment over Input: Cytoplasm CoCLIP vs FracCLIP (',  nrow(Peak_Cyt_FC), ' peaks)'))

plotScatter1(Peak_Cyt_FC, grouped_annotation,
             log2(NES_EvI_S), log2(Cyto_EvI_S),
             'log2(CoCLIP Nuclear/Input)', 'log2(FracCLIP Nuclear/Input)',
             x_lim, y_lim, 
             paste0('Stress HuR Peak Enrichment over Input: Cytoplasm CoCLIP vs FracCLIP (',  nrow(Peak_Cyt_FC), ' peaks)'))

# Cytoplasm - mRNA specific:
Peak_Cyt_FC_mRNA = Peak_Cyt_FC %>% filter((finalized_annotation == "5'UTR" | 
                                             finalized_annotation == "3'UTR" | 
                                             finalized_annotation == "CDS" | 
                                             finalized_annotation == "intron" | 
                                             finalized_annotation == "CDS_RI" |
                                             finalized_annotation == "DS10K"))

Peak_Cyt_FC_mRNA$finalized_annotation = factor(Peak_Cyt_FC_mRNA$finalized_annotation, levels = mRNA_List)

plotScatter1(Peak_Cyt_FC_mRNA, finalized_annotation,
             log2(NES_EvI_M), log2(Cyto_EvI_M),
             'log2(CoCLIP Cytoplasm/Input)', 'log2(FracCLIP Cytoplasm/Input)',
             x_lim, y_lim, 
             paste0('Mock HuR mRNA Peak Enrichment over Input: Cytoplasm CoCLIP vs FracCLIP (',  nrow(Peak_Cyt_FC_mRNA), ' peaks)'))

plotScatter1(Peak_Cyt_FC_mRNA, finalized_annotation,
             log2(NES_EvI_S), log2(Cyto_EvI_S),
             'log2(CoCLIP Cytoplasm/Input)', 'log2(FracCLIP Cytoplasm/Input)',
             x_lim, y_lim, 
             paste0('Stress HuR mRNA Peak Enrichment over Input: Cytoplasm CoCLIP vs FracCLIP (',  nrow(Peak_Cyt_FC_mRNA), ' peaks)'))

# Cytoplasm Mock - ncRNA specific:
Peak_Cyt_FC_ncRNA = Peak_Cyt_FC %>% filter((finalized_annotation != "5'UTR" & 
                                              finalized_annotation != "3'UTR" & 
                                              finalized_annotation != "CDS" & 
                                              finalized_annotation != "intron" &  
                                              finalized_annotation != "CDS_RI" &
                                              finalized_annotation != "DS10K" &
                                              finalized_annotation != "UnAn"))

Peak_Cyt_FC_ncRNA$finalized_annotation = factor(Peak_Cyt_FC_ncRNA$finalized_annotation, levels = ncRNA_List)

plotScatter1(Peak_Cyt_FC_ncRNA, finalized_annotation,
             log2(NES_EvI_M), log2(Cyto_EvI_M),
             'log2(CoCLIP Nuclear/Input)', 'log2(FracCLIP Nuclear/Input)',
             x_lim, y_lim, 
             paste0('Mock HuR ncRNA Peak Enrichment over Input: Nuclear CoCLIP vs FracCLIP (',  nrow(Peak_Cyt_FC_ncRNA), ' peaks)'))

plotScatter1(Peak_Cyt_FC_ncRNA, finalized_annotation,
             log2(NES_EvI_S), log2(Cyto_EvI_S),
             'log2(CoCLIP Nuclear/Input)', 'log2(FracCLIP Nuclear/Input)',
             x_lim, y_lim, 
             paste0('Stress HuR ncRNA Peak Enrichment over Input: Nuclear CoCLIP vs FracCLIP (',  nrow(Peak_Cyt_FC_ncRNA), ' peaks)'))

####################################################################################################################

## FIGURE2 Fractionation and CoCLIP comparison
####################################################################################################################
Peak_Nuc_M = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                         (Nuc_F_M_BC >= BC_Threshold_F & Nuc_F_M > median(Nuc_F_M) * rowSum_Multiplier_F) &
                                         (NLS_E_M_BC >= BC_Threshold_E & NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_E))

Peak_Cyt_M = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                         (Cyto_F_M_BC >= BC_Threshold_F & Cyto_F_M > median(Cyto_F_M) * rowSum_Multiplier_F) &
                                         (NES_E_M_BC >= BC_Threshold_E & NES_E_M > median(NES_E_M) * rowSum_Multiplier_E))

Peak_Nuc_FC = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                          ((Nuc_F_M_BC >= BC_Threshold_F & Nuc_F_M > median(Nuc_F_M) * rowSum_Multiplier_F &
                                              NLS_E_M_BC >= BC_Threshold_E & NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_E &
                                              (NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I & I_M > median(I_M) * rowSum_Multiplier_I) |
                                             (Nuc_F_S_BC >= BC_Threshold_F & Nuc_F_S > median(Nuc_F_S) * rowSum_Multiplier_F &
                                                NLS_E_S_BC >= BC_Threshold_E & NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_E &
                                                (NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I & I_S > median(I_S) * rowSum_Multiplier_I)))

Peak_Cyt_FC = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                          ((Cyto_F_M_BC >= BC_Threshold_F & Cyto_F_M > median(Cyto_F_M) * rowSum_Multiplier_F &
                                              NES_E_M_BC >= BC_Threshold_E & NES_E_M > median(NES_E_M) * rowSum_Multiplier_E &
                                              (NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I & I_M > median(I_M) * rowSum_Multiplier_I) |
                                             (Cyto_F_S_BC >= BC_Threshold_F & Cyto_F_S > median(Cyto_F_S) * rowSum_Multiplier_F &
                                                NES_E_S_BC >= BC_Threshold_E & NES_E_S > median(NES_E_S) * rowSum_Multiplier_E &
                                                (NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I & I_S > median(I_S) * rowSum_Multiplier_I)))


Peak_NvC_M

Peak_NvC_S

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