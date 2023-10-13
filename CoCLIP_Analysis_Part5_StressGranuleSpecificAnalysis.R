## CoCLIP Analysis: 
## Peak Processing for Enrichment Analysis
## Written by Soon Yi
## Last Edit: 2023-09-18

library(stringr)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(scales)

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

## Plot Scatter:
plotScatter = function(peak_matrix, annotation_level, x_axis, y_axis, x_label = NULL, y_label = NULL, x_lim = NULL, y_lim = NULL, title = NULL, diag = FALSE) {
  plot = ggplot(peak_matrix, aes(x = ({{x_axis}}), y = ({{y_axis}}), color = {{annotation_level}})) +
    geom_point(pch = 16, size = 3, alpha = 0.5) +
    scale_fill_brewer(palette = "Set3") +
    theme_bw() + 
    theme(axis.text = element_text(size=14), 
          axis.title = element_text(size=14, face = 'bold'), 
          legend.text = element_text(size=14))
  
  if (!is.null(title)) {
    plot = plot + ggtitle(title)
  }
  if (!is.null(x_label)) {
    plot = plot + labs(x = x_label)
  }
  if (!is.null(y_label)) {
    plot = plot + labs(y = y_label)
  }
  if (!is.null(x_lim)) {
    plot = plot + xlim(x_lim)
  }
  if (!is.null(y_lim)) {
    plot = plot + ylim(y_lim)
  }
  if (diag == TRUE) {
    plot = plot + geom_abline(linetype = 'dotted') 
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

peakEnrichment = peakEnrichment %>% mutate(NLS_SvM = peakRowSum$NLS_E_S / peakRowSum$NLS_E_M)
peakEnrichment = peakEnrichment %>% mutate(NES_SvM = peakRowSum$NES_E_S / peakRowSum$NES_E_M)
peakEnrichment = peakEnrichment %>% mutate(G3BP_SvM = peakRowSum$G3BP_E_S / peakRowSum$G3BP_E_M)

peakEnrichment = cbind(peakEnrichment, peakRowSum[, colnames(peakRowSum)[17:34]])

# write.table(peakEnrichment, paste0(peaksMatrix_PATH, str_replace(peaksMatrix_FILE, ".txt", "_Enrichment.txt")), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t', na = "")
####################################################################################################################

x_lim = c(-7, 7)
y_lim = x_lim

## FIGURE3 Mock G3BP E/I vs Mock Nucleus E/I 
####################################################################################################################
Peak_G3BP_Nuc_M_EvI = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                                   (((G3BP_E_M_BC >= BC_Threshold_E_SG & G3BP_E_M > median(G3BP_E_M) * rowSum_Multiplier_E) & 
                                                       (G3BP_I_M_BC >= BC_Threshold_E & G3BP_I_M > median(G3BP_I_M) * rowSum_Multiplier_I)) | 
                                                      ((NLS_E_M_BC >= BC_Threshold_E & NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_E) & 
                                                         (NLS_I_M_BC >= 1 & NLS_I_M > median(NLS_I_M) * rowSum_Multiplier_I))))

Peak_G3BP_Nuc_M_EvI$grouped_annotation = factor(Peak_G3BP_Nuc_M_EvI$grouped_annotation, levels = All_Annotation_List)

plotScatter(Peak_G3BP_Nuc_M_EvI, grouped_annotation,
            log2(NLS_EvI_M), log2(G3BP_EvI_M),
            x_label = 'log2(CoCLIP Mock Nuclear/Input)', 'log2(CoCLIP Mock SG/Input)',
            x_lim, y_lim,
            title = paste0('HuR Peaks: Mock SG vs Nuclear Enriched over Input (',  nrow(Peak_G3BP_Nuc_M_EvI), ' peaks)'))

Peak_G3BP_Nuc_M_EvI_mRNA = Peak_G3BP_Nuc_M_EvI %>% filter((finalized_annotation == "5'UTR" | 
                                                               finalized_annotation == "3'UTR" | 
                                                               finalized_annotation == "CDS" | 
                                                               finalized_annotation == "intron" | 
                                                               finalized_annotation == "CDS_RI" |
                                                               finalized_annotation == "DS10K"))
Peak_G3BP_Nuc_M_EvI_mRNA$finalized_annotation = factor(Peak_G3BP_Nuc_M_EvI_mRNA$finalized_annotation, levels = mRNA_List)

plotScatter(Peak_G3BP_Nuc_M_EvI_mRNA, finalized_annotation,
            log2(NLS_EvI_M), log2(G3BP_EvI_M),
            x_label = 'log2(CoCLIP Mock Nuclear/Input)', 'log2(CoCLIP Mock SG/Input)',
            x_lim, y_lim,
            paste0('HuR mRNA Peaks: Mock SG vs Nuclear Enriched over Input (',  nrow(Peak_G3BP_Nuc_M_EvI_mRNA), ' peaks)'))

Peak_G3BP_Nuc_M_EvI_ncRNA = Peak_G3BP_Nuc_M_EvI %>% filter((finalized_annotation != "5'UTR" & 
                                                                finalized_annotation != "3'UTR" & 
                                                                finalized_annotation != "CDS" & 
                                                                finalized_annotation != "intron" &  
                                                                finalized_annotation != "CDS_RI" &
                                                                finalized_annotation != "DS10K" &
                                                                finalized_annotation != "UnAn"))
Peak_G3BP_Nuc_M_EvI_ncRNA$finalized_annotation = factor(Peak_G3BP_Nuc_M_EvI_ncRNA$finalized_annotation, levels = ncRNA_List)

plotScatter(Peak_G3BP_Nuc_M_EvI_ncRNA, finalized_annotation,
            log2(NLS_EvI_M), log2(G3BP_EvI_M),
            x_label = 'log2(CoCLIP Mock Nuclear/Input)', 'log2(CoCLIP Mock SG/Input)',
            x_lim, y_lim,
            paste0('HuR ncRNA Peaks: Mock SG vs Nuclear Enriched over Input (',  nrow(Peak_G3BP_Nuc_M_EvI_ncRNA), ' peaks)'))
#################################################################################################################### 

## FIGURE3 Stress G3BP E/I vs Stress Nucleus E/I 
####################################################################################################################
Peak_G3BP_Nuc_S_EvI = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                                   (((G3BP_E_S_BC >= BC_Threshold_E_SG & G3BP_E_S > median(G3BP_E_S) * rowSum_Multiplier_E) & 
                                                       (G3BP_I_S_BC >= BC_Threshold_E & G3BP_I_S > median(G3BP_I_S) * rowSum_Multiplier_I)) | 
                                                      ((NLS_E_S_BC >= BC_Threshold_E & NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_E) & 
                                                         (NLS_I_S_BC >= 1 & NLS_I_S > median(NLS_I_S) * rowSum_Multiplier_I))))

Peak_G3BP_Nuc_S_EvI$grouped_annotation = factor(Peak_G3BP_Nuc_S_EvI$grouped_annotation, levels = All_Annotation_List)

plotScatter(Peak_G3BP_Nuc_S_EvI, grouped_annotation,
            log2(NLS_EvI_S), log2(G3BP_EvI_S),
            x_label = 'log2(CoCLIP Stress Nuclear/Input)', 'log2(CoCLIP Stress SG/Input)',
            x_lim, y_lim,
            title = paste0('HuR Peaks: Stress SG vs Nuclear Enriched over Input (',  nrow(Peak_G3BP_Nuc_S_EvI), ' peaks)'))

Peak_G3BP_Nuc_S_EvI_mRNA = Peak_G3BP_Nuc_S_EvI %>% filter((finalized_annotation == "5'UTR" | 
                                                               finalized_annotation == "3'UTR" | 
                                                               finalized_annotation == "CDS" | 
                                                               finalized_annotation == "intron" | 
                                                               finalized_annotation == "CDS_RI" |
                                                               finalized_annotation == "DS10K"))
Peak_G3BP_Nuc_S_EvI_mRNA$finalized_annotation = factor(Peak_G3BP_Nuc_S_EvI_mRNA$finalized_annotation, levels = mRNA_List)

plotScatter(Peak_G3BP_Nuc_S_EvI_mRNA, finalized_annotation,
            log2(NLS_EvI_S), log2(G3BP_EvI_S),
            x_label = 'log2(CoCLIP Stress Nuclear/Input)', 'log2(CoCLIP Stress SG/Input)',
            x_lim, y_lim,
            paste0('HuR mRNA Peaks: Stress SG vs Nuclear Enriched over Input (',  nrow(Peak_G3BP_Nuc_S_EvI_mRNA), ' peaks)'))

Peak_G3BP_Nuc_S_EvI_ncRNA = Peak_G3BP_Nuc_S_EvI %>% filter((finalized_annotation != "5'UTR" & 
                                                                finalized_annotation != "3'UTR" & 
                                                                finalized_annotation != "CDS" & 
                                                                finalized_annotation != "intron" &  
                                                                finalized_annotation != "CDS_RI" &
                                                                finalized_annotation != "DS10K" &
                                                                finalized_annotation != "UnAn"))
Peak_G3BP_Nuc_S_EvI_ncRNA$finalized_annotation = factor(Peak_G3BP_Nuc_S_EvI_ncRNA$finalized_annotation, levels = ncRNA_List)

plotScatter(Peak_G3BP_Nuc_S_EvI_ncRNA, finalized_annotation,
            log2(NLS_EvI_S), log2(G3BP_EvI_S),
            x_label = 'log2(CoCLIP Stress Nuclear/Input)', 'log2(CoCLIP Stress SG/Input)',
            x_lim, y_lim,
            paste0('HuR ncRNA Peaks: Stress SG vs Nuclear Enriched over Input (',  nrow(Peak_G3BP_Nuc_S_EvI_ncRNA), ' peaks)'))
####################################################################################################################

## FIGURE3 Mock G3BP E/I vs Mock Cyto E/I 
####################################################################################################################
Peak_G3BP_Cyto_M_EvI = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                                 (((G3BP_E_M_BC >= BC_Threshold_E_SG & G3BP_E_M > median(G3BP_E_M) * rowSum_Multiplier_E) & 
                                                     (G3BP_I_M_BC >= BC_Threshold_E & G3BP_I_M > median(G3BP_I_M) * rowSum_Multiplier_I)) | 
                                                    ((NES_E_M_BC >= BC_Threshold_E & NES_E_M > median(NES_E_M) * rowSum_Multiplier_E) & 
                                                       (NES_I_M_BC >= 1 & NES_I_M > median(NES_I_M) * rowSum_Multiplier_I))))

Peak_G3BP_Cyto_M_EvI$grouped_annotation = factor(Peak_G3BP_Cyto_M_EvI$grouped_annotation, levels = All_Annotation_List)

plotScatter(Peak_G3BP_Cyto_M_EvI, grouped_annotation,
            log2(NES_EvI_M), log2(G3BP_EvI_M),
            x_label = 'log2(CoCLIP Mock Cytoplasm/Input)', 'log2(CoCLIP Mock SG/Input)',
            x_lim, y_lim,
            title = paste0('HuR Peaks: Mock SG vs Cytoplasm Enriched over Input (',  nrow(Peak_G3BP_Cyto_M_EvI), ' peaks)'))

Peak_G3BP_Cyto_M_EvI_mRNA = Peak_G3BP_Cyto_M_EvI %>% filter((finalized_annotation == "5'UTR" | 
                                                           finalized_annotation == "3'UTR" | 
                                                           finalized_annotation == "CDS" | 
                                                           finalized_annotation == "intron" | 
                                                           finalized_annotation == "CDS_RI" |
                                                           finalized_annotation == "DS10K"))
Peak_G3BP_Cyto_M_EvI_mRNA$finalized_annotation = factor(Peak_G3BP_Cyto_M_EvI_mRNA$finalized_annotation, levels = mRNA_List)

plotScatter(Peak_G3BP_Cyto_M_EvI_mRNA, finalized_annotation,
            log2(NES_EvI_M), log2(G3BP_EvI_M),
            x_label = 'log2(CoCLIP Mock Cytoplasm/Input)', 'log2(CoCLIP Mock SG/Input)',
            x_lim, y_lim,
            paste0('HuR mRNA Peaks: Mock SG vs Cytoplasm Enriched over Input (',  nrow(Peak_G3BP_Cyto_M_EvI_mRNA), ' peaks)'))

Peak_G3BP_Cyto_M_EvI_ncRNA = Peak_G3BP_Cyto_M_EvI %>% filter((finalized_annotation != "5'UTR" & 
                                                            finalized_annotation != "3'UTR" & 
                                                            finalized_annotation != "CDS" & 
                                                            finalized_annotation != "intron" &  
                                                            finalized_annotation != "CDS_RI" &
                                                            finalized_annotation != "DS10K" &
                                                            finalized_annotation != "UnAn"))
Peak_G3BP_Cyto_M_EvI_ncRNA$finalized_annotation = factor(Peak_G3BP_Cyto_M_EvI_ncRNA$finalized_annotation, levels = ncRNA_List)

plotScatter(Peak_G3BP_Cyto_M_EvI_ncRNA, finalized_annotation,
            log2(NES_EvI_M), log2(G3BP_EvI_M),
            x_label = 'log2(CoCLIP Mock Cytoplasm/Input)', 'log2(CoCLIP Mock SG/Input)',
            x_lim, y_lim,
            paste0('HuR ncRNA Peaks: Mock SG vs Cytoplasm Enriched over Input (',  nrow(Peak_G3BP_Cyto_M_EvI_ncRNA), ' peaks)'))
#################################################################################################################### 

## FIGURE3 Stress G3BP E/I vs Stress Cyto E/I 
####################################################################################################################
Peak_G3BP_Cyto_S_EvI = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                             (((G3BP_E_S_BC >= BC_Threshold_E_SG & G3BP_E_S > median(G3BP_E_S) * rowSum_Multiplier_E) & 
                                                 (G3BP_I_S_BC >= BC_Threshold_E & G3BP_I_S > median(G3BP_I_S) * rowSum_Multiplier_I)) | 
                                                ((NES_E_S_BC >= BC_Threshold_E & NES_E_S > median(NES_E_S) * rowSum_Multiplier_E) & 
                                                   (NES_I_S_BC >= 1 & NES_I_S > median(NES_I_S) * rowSum_Multiplier_I))))

Peak_G3BP_Cyto_S_EvI$grouped_annotation = factor(Peak_G3BP_Cyto_S_EvI$grouped_annotation, levels = All_Annotation_List)

plotScatter(Peak_G3BP_Cyto_S_EvI, grouped_annotation,
            log2(NES_EvI_S), log2(G3BP_EvI_S),
            x_label = 'log2(CoCLIP Stress Cytoplasm/Input)', 'log2(CoCLIP Stress SG/Input)',
            x_lim, y_lim,
            title = paste0('HuR Peaks: Stress SG vs Cytoplasm Enriched over Input (',  nrow(Peak_G3BP_Cyto_S_EvI), ' peaks)'))

Peak_G3BP_Cyto_S_EvI_mRNA = Peak_G3BP_Cyto_S_EvI %>% filter((finalized_annotation == "5'UTR" | 
                                                   finalized_annotation == "3'UTR" | 
                                                   finalized_annotation == "CDS" | 
                                                   finalized_annotation == "intron" | 
                                                   finalized_annotation == "CDS_RI" |
                                                   finalized_annotation == "DS10K"))
Peak_G3BP_Cyto_S_EvI_mRNA$finalized_annotation = factor(Peak_G3BP_Cyto_S_EvI_mRNA$finalized_annotation, levels = mRNA_List)

plotScatter(Peak_G3BP_Cyto_S_EvI_mRNA, finalized_annotation,
            log2(NES_EvI_S), log2(G3BP_EvI_S),
            x_label = 'log2(CoCLIP Stress Cytoplasm/Input)', 'log2(CoCLIP Stress SG/Input)',
            x_lim, y_lim,
            paste0('HuR mRNA Peaks: Stress SG vs Cytoplasm Enriched over Input (',  nrow(Peak_G3BP_Cyto_S_EvI_mRNA), ' peaks)'))

Peak_G3BP_Cyto_S_EvI_ncRNA = Peak_G3BP_Cyto_S_EvI %>% filter((finalized_annotation != "5'UTR" & 
                                                    finalized_annotation != "3'UTR" & 
                                                    finalized_annotation != "CDS" & 
                                                    finalized_annotation != "intron" &  
                                                    finalized_annotation != "CDS_RI" &
                                                    finalized_annotation != "DS10K" &
                                                    finalized_annotation != "UnAn"))
Peak_G3BP_Cyto_S_EvI_ncRNA$finalized_annotation = factor(Peak_G3BP_Cyto_S_EvI_ncRNA$finalized_annotation, levels = ncRNA_List)

plotScatter(Peak_G3BP_Cyto_S_EvI_ncRNA, finalized_annotation,
            log2(NES_EvI_S), log2(G3BP_EvI_S),
            x_label = 'log2(CoCLIP Stress Cytoplasm/Input)', 'log2(CoCLIP Stress SG/Input)',
            x_lim, y_lim,
            paste0('HuR ncRNA Peaks: Stress SG vs Cytoplasm Enriched over Input (',  nrow(Peak_G3BP_Cyto_S_EvI_ncRNA), ' peaks)'))
####################################################################################################################
