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
library(scales)

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

## Plot cumulative distribution:
plot_CD = function(count_table, y_data, colormaps, linetypes) {
  plot_data = count_table %>% select(counts, {{y_data}})
  plot_data_long = plot_data %>% pivot_longer(cols = {{y_data}}, names_to = "Sample", values_to = "nCounts")
  plot_data_long$Sample = factor(plot_data_long$Sample, levels = unique(plot_data_long$Sample)) # Use unique() to ensure correct order of levels
  
  plot = ggplot(plot_data_long, aes(x = counts, y = nCounts, color = Sample, linetype = Sample)) +
    geom_line() +
    theme_minimal() +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = 'bold')) +
    scale_x_continuous(trans = 'log2') +
    scale_color_manual(values = colormaps) +
    scale_linetype_manual(values = linetypes)
  
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

## FIGURE SG to Nuclear and/or Cytoplasm comparison
####################################################################################################################
x_lim = c(-8, 8)
y_lim = x_lim

Peak_SG_M = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                        (((G3BP_E_M_BC >= BC_Threshold_E_SG & G3BP_E_M > median(G3BP_E_M) * rowSum_Multiplier_E) &
                                            ((NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I & I_M > median(I_M) * rowSum_Multiplier_I)) |
                                           ((NLS_E_M_BC >= BC_Threshold_E & NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_E) &
                                              (NES_E_M_BC >= BC_Threshold_E & NES_E_M > median(NES_E_M) * rowSum_Multiplier_E))))

Peak_SG_M$grouped_annotation = factor(Peak_SG_M$grouped_annotation, levels = All_Annotation_List)

plotScatter(Peak_SG_M, grouped_annotation,
            log2(E_NvC_M), log2(G3BP_EvI_M),
            x_label = 'log2(CoCLIP Nuclear/Cytoplasm)', 'log2(CoCLIP Stress Granule/Input)',
            x_lim, y_lim,
            paste0('Mock HuR Peaks: Stress Granule CoCLIP Enrichment (',  nrow(Peak_SG_M), ' peaks)'))


Peak_SG_M_mRNA = Peak_SG_M %>% filter((finalized_annotation == "5'UTR" | 
                                         finalized_annotation == "3'UTR" | 
                                         finalized_annotation == "CDS" | 
                                         finalized_annotation == "intron" | 
                                         finalized_annotation == "CDS_RI" |
                                         finalized_annotation == "DS10K"))
Peak_SG_M_mRNA$finalized_annotation = factor(Peak_SG_M_mRNA$finalized_annotation, levels = mRNA_List)

plotScatter(Peak_SG_M_mRNA, finalized_annotation,
            log2(E_NvC_M), log2(G3BP_EvI_M),
            x_label = 'log2(CoCLIP Nuclear/Cytoplasm)', 'log2(CoCLIP Stress Granule/Input)',
            x_lim, y_lim,
            paste0('Mock HuR mRNA Peaks: Stress Granule CoCLIP Enrichment (',  nrow(Peak_SG_M_mRNA), ' peaks)'))

Peak_SG_M_ncRNA = Peak_SG_M %>% filter((finalized_annotation != "5'UTR" & 
                                          finalized_annotation != "3'UTR" & 
                                          finalized_annotation != "CDS" & 
                                          finalized_annotation != "intron" &  
                                          finalized_annotation != "CDS_RI" &
                                          finalized_annotation != "DS10K" &
                                          finalized_annotation != "UnAn"))
Peak_SG_M_ncRNA$finalized_annotation = factor(Peak_SG_M_ncRNA$finalized_annotation, levels = ncRNA_List)

plotScatter(Peak_SG_M_ncRNA, finalized_annotation,
            log2(E_NvC_M), log2(G3BP_EvI_M),
            x_label = 'log2(CoCLIP Nuclear/Cytoplasm)', 'log2(CoCLIP Stress Granule/Input)',
            x_lim, y_lim,
            paste0('Mock HuR ncRNA Peaks: Stress Granule CoCLIP Enrichment (',  nrow(Peak_SG_M_ncRNA), ' peaks)'))

Peak_SG_S = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                        (((G3BP_E_S_BC >= BC_Threshold_E_SG & G3BP_E_S > median(G3BP_E_S) * rowSum_Multiplier_E) &
                                            ((NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I & I_S > median(I_S) * rowSum_Multiplier_I)) |
                                           ((NLS_E_S_BC >= BC_Threshold_E & NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_E) &
                                              (NES_E_S_BC >= BC_Threshold_E & NES_E_S > median(NES_E_S) * rowSum_Multiplier_E))))

# Peak_SG_S = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
#                                         (((G3BP_E_M_BC >= BC_Threshold_E_SG & G3BP_E_M > median(G3BP_E_M) * rowSum_Multiplier_E) &
#                                             ((NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I & I_M > median(I_M) * rowSum_Multiplier_I)) |
#                                            ((NLS_E_M_BC >= BC_Threshold_E & NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_E) &
#                                               (NES_E_M_BC >= BC_Threshold_E & NES_E_M > median(NES_E_M) * rowSum_Multiplier_E))) &
#                                         (((G3BP_E_S_BC >= BC_Threshold_E_SG & G3BP_E_S > median(G3BP_E_S) * rowSum_Multiplier_E) &
#                                             ((NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I & I_S > median(I_S) * rowSum_Multiplier_I)) |
#                                            ((NLS_E_S_BC >= BC_Threshold_E & NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_E) &
#                                               (NES_E_S_BC >= BC_Threshold_E & NES_E_S > median(NES_E_S) * rowSum_Multiplier_E))))

Peak_SG_S$grouped_annotation = factor(Peak_SG_S$grouped_annotation, levels = All_Annotation_List)

plotScatter(Peak_SG_S, grouped_annotation,
            log2(E_NvC_S), log2(G3BP_EvI_S),
            x_label = 'log2(CoCLIP Nuclear/Cytoplasm)', 'log2(CoCLIP Stress Granule/Input)',
            x_lim, y_lim,
            paste0('Arsenite HuR Peaks: Stress Granule CoCLIP Enrichment (',  nrow(Peak_SG_S), ' peaks)'))

Peak_SG_S_mRNA = Peak_SG_S %>% filter((finalized_annotation == "5'UTR" | 
                                         finalized_annotation == "3'UTR" | 
                                         finalized_annotation == "CDS" | 
                                         finalized_annotation == "intron" | 
                                         finalized_annotation == "CDS_RI" |
                                         finalized_annotation == "DS10K"))
Peak_SG_S_mRNA$finalized_annotation = factor(Peak_SG_S_mRNA$finalized_annotation, levels = mRNA_List)

plotScatter(Peak_SG_S_mRNA, finalized_annotation,
            log2(E_NvC_S), log2(G3BP_EvI_S),
            x_label = 'log2(CoCLIP Nuclear/Cytoplasm)', 'log2(CoCLIP Stress Granule/Input)',
            x_lim, y_lim,
            paste0('Arsenite HuR mRNA Peaks: Stress Granule CoCLIP Enrichment (',  nrow(Peak_SG_S_mRNA), ' peaks)'))

Peak_SG_S_ncRNA = Peak_SG_S %>% filter((finalized_annotation != "5'UTR" & 
                                          finalized_annotation != "3'UTR" & 
                                          finalized_annotation != "CDS" & 
                                          finalized_annotation != "intron" &  
                                          finalized_annotation != "CDS_RI" &
                                          finalized_annotation != "DS10K" &
                                          finalized_annotation != "UnAn"))
Peak_SG_S_ncRNA$finalized_annotation = factor(Peak_SG_S_ncRNA$finalized_annotation, levels = ncRNA_List)

plotScatter(Peak_SG_S_ncRNA, finalized_annotation,
            log2(E_NvC_S), log2(G3BP_EvI_S),
            x_label = 'log2(CoCLIP Nuclear/Cytoplasm)', 'log2(CoCLIP Stress Granule/Input)',
            x_lim, y_lim,
            paste0('Arsenite HuR ncRNA Peaks: Stress Granule CoCLIP Enrichment (',  nrow(Peak_SG_S_ncRNA), ' peaks)'))

####################################################################################################################

## FIGURE SG to Nuclear and/or Cytoplasm comparison (SG Mock v Stress) vs (Mock Nuclear v Cytoplasm)
####################################################################################################################
Peak_SG_SvM_M = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                            (((G3BP_E_M_BC >= BC_Threshold_E_SG & G3BP_E_M > median(G3BP_E_M) * rowSum_Multiplier_E) & 
                                                (G3BP_E_S_BC >= BC_Threshold_E_SG & G3BP_E_S > median(G3BP_E_S) * rowSum_Multiplier_E)) |
                                               ((NLS_E_M_BC >= BC_Threshold_E & NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_E) & 
                                                  (NES_E_M_BC >= BC_Threshold_E & NES_E_M > median(NES_E_M) * rowSum_Multiplier_E))))

Peak_SG_SvM_M$grouped_annotation = factor(Peak_SG_SvM_M$grouped_annotation, levels = All_Annotation_List)

plotScatter(Peak_SG_SvM_M, grouped_annotation,
            log2(E_NvC_M), log2(G3BP_SvM),
            x_label = 'log2(CoCLIP Nuclear/Cytoplasm)', 'log2(CoCLIP Stress Granule Stress/Mock)',
            x_lim, y_lim,
            paste0('Mock HuR Peaks: Stress Granule CoCLIP Enrichment (',  nrow(Peak_SG_SvM_M), ' peaks)'))


Peak_SG_SvM_M_mRNA = Peak_SG_SvM_M %>% filter((finalized_annotation == "5'UTR" | 
                                                 finalized_annotation == "3'UTR" | 
                                                 finalized_annotation == "CDS" | 
                                                 finalized_annotation == "intron" | 
                                                 finalized_annotation == "CDS_RI" |
                                                 finalized_annotation == "DS10K"))
Peak_SG_SvM_M_mRNA$finalized_annotation = factor(Peak_SG_SvM_M_mRNA$finalized_annotation, levels = mRNA_List)

plotScatter(Peak_SG_SvM_M_mRNA, finalized_annotation,
            log2(E_NvC_M), log2(G3BP_SvM),
            x_label = 'log2(CoCLIP Nuclear/Cytoplasm)', 'log2(CoCLIP Stress Granule Stress/Mock)',
            x_lim, y_lim,
            paste0('Mock HuR mRNA Peaks: Stress Granule CoCLIP Enrichment (',  nrow(Peak_SG_SvM_M_mRNA), ' peaks)'))

Peak_SG_SvM_M_ncRNA = Peak_SG_SvM_M %>% filter((finalized_annotation != "5'UTR" & 
                                                  finalized_annotation != "3'UTR" & 
                                                  finalized_annotation != "CDS" & 
                                                  finalized_annotation != "intron" &  
                                                  finalized_annotation != "CDS_RI" &
                                                  finalized_annotation != "DS10K" &
                                                  finalized_annotation != "UnAn"))
Peak_SG_SvM_M_ncRNA$finalized_annotation = factor(Peak_SG_SvM_M_ncRNA$finalized_annotation, levels = ncRNA_List)

plotScatter(Peak_SG_SvM_M_ncRNA, finalized_annotation,
            log2(E_NvC_M), log2(G3BP_SvM),
            x_label = 'log2(CoCLIP Nuclear/Cytoplasm)', 'log2(CoCLIP Stress Granule Stress/Mock)',
            x_lim, y_lim,
            paste0('Mock HuR ncRNA Peaks: Stress Granule CoCLIP Enrichment (',  nrow(Peak_SG_SvM_M_ncRNA), ' peaks)'))

####################################################################################################################

## FIGURE SG to Nuclear and/or Cytoplasm comparison (SG Mock v Stress) vs (Stress Nuclear v Cytoplasm)
####################################################################################################################
Peak_SG_SvM_S = peakEnrichment %>% filter((grouped_annotation != 'UnAn') &
                                            (((G3BP_E_M_BC >= BC_Threshold_E_SG & G3BP_E_M > median(G3BP_E_M) * rowSum_Multiplier_E) & 
                                                (G3BP_E_S_BC >= BC_Threshold_E_SG & G3BP_E_S > median(G3BP_E_S) * rowSum_Multiplier_E)) |
                                               ((NLS_E_S_BC >= BC_Threshold_E & NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_E) & 
                                                  (NES_E_S_BC >= BC_Threshold_E & NES_E_S > median(NES_E_S) * rowSum_Multiplier_E))))

Peak_SG_SvM_S$grouped_annotation = factor(Peak_SG_SvM_S$grouped_annotation, levels = All_Annotation_List)

plotScatter(Peak_SG_SvM_S, grouped_annotation,
            log2(E_NvC_S), log2(G3BP_SvM),
            x_label = 'log2(CoCLIP Nuclear/Cytoplasm)', 'log2(CoCLIP Stress Granule Stress/Mock)',
            x_lim, y_lim,
            paste0('Arsenite HuR Peaks: Stress Granule CoCLIP Enrichment (',  nrow(Peak_SG_SvM), ' peaks)'))


Peak_SG_SvM_S_mRNA = Peak_SG_SvM_S %>% filter((finalized_annotation == "5'UTR" | 
                                                 finalized_annotation == "3'UTR" | 
                                                 finalized_annotation == "CDS" | 
                                                 finalized_annotation == "intron" | 
                                                 finalized_annotation == "CDS_RI" |
                                                 finalized_annotation == "DS10K"))
Peak_SG_SvM_S_mRNA$finalized_annotation = factor(Peak_SG_SvM_S_mRNA$finalized_annotation, levels = mRNA_List)

plotScatter(Peak_SG_SvM_S_mRNA, finalized_annotation,
            log2(E_NvC_S), log2(G3BP_SvM),
            x_label = 'log2(CoCLIP Nuclear/Cytoplasm)', 'log2(CoCLIP Stress Granule Stress/Mock)',
            x_lim, y_lim,
            paste0('Arsenite HuR mRNA Peaks: Stress Granule CoCLIP Enrichment (',  nrow(Peak_SG_SvM_S_mRNA), ' peaks)'))

Peak_SG_SvM_S_ncRNA = Peak_SG_SvM_S %>% filter((finalized_annotation != "5'UTR" & 
                                                  finalized_annotation != "3'UTR" & 
                                                  finalized_annotation != "CDS" & 
                                                  finalized_annotation != "intron" &  
                                                  finalized_annotation != "CDS_RI" &
                                                  finalized_annotation != "DS10K" &
                                                  finalized_annotation != "UnAn"))
Peak_SG_SvM_S_ncRNA$finalized_annotation = factor(Peak_SG_SvM_S_ncRNA$finalized_annotation, levels = ncRNA_List)

plotScatter(Peak_SG_SvM_S_ncRNA, finalized_annotation,
            log2(E_NvC_S), log2(G3BP_SvM),
            x_label = 'log2(CoCLIP Nuclear/Cytoplasm)', 'log2(CoCLIP Stress Granule Stress/Mock)',
            x_lim, y_lim,
            paste0('Arsenite HuR ncRNA Peaks: Stress Granule CoCLIP Enrichment (',  nrow(Peak_SG_SvM_S_ncRNA), ' peaks)'))

####################################################################################################################





## SG Enrichment and Motif High Order Analysis
####################################################################################################################
baseDir = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/motifEnrichment/motifs/counts'
# baseDir = 'L:/.shortcut-targets-by-id/13hY9t_p6eUdnvP-c2OClHbzyeWMOruD_/CoCLIP_HuR_Paper/Data/Homer_Outputs/motifs/counts'
setwd(baseDir)
countFiles = list.files(baseDir)
countFiles = countFiles[grep('MotifCounts.txt', countFiles)]

## G3BP Mock
countFile = countFiles[4]
motifCounts = read.delim(countFile, header = T)
motifCounts_M_AAAAA = motifCounts %>% filter(Motif.Name == '0-AAAAA')
motifCounts_M_UUUUU = motifCounts %>% filter(Motif.Name == '0-UUUUU')

motifCounts_perPeak_M_AAAAA = data.frame(table(motifCounts_M_AAAAA$PositionID))
motifCounts_perPeak_M_UUUUU = data.frame(table(motifCounts_M_UUUUU$PositionID))


## G3BP Arsenite
countFile = countFiles[3]
motifCounts = read.delim(countFile, header = T)
motifCounts_S_AAAAA = motifCounts %>% filter(Motif.Name == '0-AAAAA')
motifCounts_S_UUUUU = motifCounts %>% filter(Motif.Name == '0-UUUUU')

motifCounts_perPeak_S_AAAAA = data.frame(table(motifCounts_S_AAAAA$PositionID))
motifCounts_perPeak_S_UUUUU = data.frame(table(motifCounts_S_UUUUU$PositionID))

## Plot Cumulative Distribution

## freq_table will look like:
## counts sample1 sample2
## counts: individual motif counts per peaks
## sample1: freq of the motif counts per peaks in sample1

freq_table = data.frame(
  counts = 0:50,
  M_AAAAA = rep(NA, 51),
  M_UUUUU = rep(NA, 51),
  S_AAAAA = rep(NA, 51),
  S_UUUUU = rep(NA, 51))

freq_table$M_AAAAA[freq_table$counts %in% rownames(table(motifCounts_perPeak_M_AAAAA$Freq))] = table(motifCounts_perPeak_M_AAAAA$Freq)
freq_table$M_UUUUU[freq_table$counts %in% rownames(table(motifCounts_perPeak_M_UUUUU$Freq))] = table(motifCounts_perPeak_M_UUUUU$Freq)
freq_table$S_AAAAA[freq_table$counts %in% rownames(table(motifCounts_perPeak_S_AAAAA$Freq))] = table(motifCounts_perPeak_S_AAAAA$Freq)
freq_table$S_UUUUU[freq_table$counts %in% rownames(table(motifCounts_perPeak_S_UUUUU$Freq))] = table(motifCounts_perPeak_S_UUUUU$Freq)
freq_table[is.na(freq_table)] = 0

cumulative_table = cbind(freq_table$counts, cumsum(freq_table[, c('M_AAAAA', 'M_UUUUU', 'S_AAAAA', 'S_UUUUU')]))
colnames(cumulative_table)[1] = 'counts'

normed_table = cbind(freq_table$counts, data.frame(sapply(cumsum(freq_table[, c('M_AAAAA', 'M_UUUUU', 'S_AAAAA', 'S_UUUUU')]), rescale)))
colnames(normed_table)[1] = 'counts'

plot_CD(normed_table, y_data = c('M_AAAAA', 'M_UUUUU', 'S_AAAAA', 'S_UUUUU'), colormaps = c('blue', 'blue', 'red', 'red'), linetypes = c(1, 2, 1, 2))

####################################################################################################################

freq_table = data.frame(counts = 0:50)
freq_table = freq_table %>% mutate(temp5A = rep(NA, 51))

for (countFile in countFiles[c(6, 10, 8, 4, 5, 9, 7, 3)]) {
  freq_table = freq_table %>% mutate(temp5A = rep(NA, 51))
  freq_table = freq_table %>% mutate(temp5U = rep(NA, 51))
  
  motifCounts = read.delim(countFile, header = T)
  temp5A_Counts = motifCounts %>% filter(Motif.Name == '0-AAAAA')
  temp5U_Counts = motifCounts %>% filter(Motif.Name == '0-UUUUU')
  temp5A_CountsPerPeak = data.frame(table(temp5A_Counts$PositionID))
  temp5U_CountsPerPeak = data.frame(table(temp5U_Counts$PositionID))
  
  freq_table$temp5A[freq_table$counts %in% rownames(table(temp5A_CountsPerPeak$Freq))] = table(temp5A_CountsPerPeak$Freq)
  rownames(table(temp5A_CountsPerPeak$Freq))[length(rownames(table(temp5A_CountsPerPeak$Freq)))]
  freq_table$temp5U[freq_table$counts %in% rownames(table(temp5U_CountsPerPeak$Freq))] = table(temp5U_CountsPerPeak$Freq)
  rownames(table(temp5U_CountsPerPeak$Freq))[length(rownames(table(temp5U_CountsPerPeak$Freq)))]
  
  freq_table = freq_table %>% rename(!!paste0(str_split(countFile, '\\.')[[1]][1], '_AAAAA') := temp5A)
  freq_table = freq_table %>% rename(!!paste0(str_split(countFile, '\\.')[[1]][1], '_UUUUU') := temp5U)
}

freq_table[is.na(freq_table)] = 0

cumulative_table = cbind(freq_table$counts, cumsum(freq_table[, colnames(freq_table)[2:ncol(freq_table)]]))
colnames(cumulative_table)[1] = 'counts'

normed_table = cbind(freq_table$counts, data.frame(sapply(cumsum(freq_table[, colnames(freq_table)[2:ncol(freq_table)]]), rescale)))
colnames(normed_table)[1] = 'counts'

plot_CD(normed_table, y_data = c('CoCLIP_NLS_Mock_AAAAA', 'CoCLIP_NLS_Mock_UUUUU', 'CoCLIP_NLS_Arsenite_AAAAA', 'CoCLIP_NLS_Arsenite_UUUUU'), colormaps = c('blue', 'blue', 'red', 'red'), linetypes = c(1, 2, 1, 2))

plot_CD(normed_table, y_data = c('CoCLIP_NES_Mock_AAAAA', 'CoCLIP_NES_Mock_UUUUU', 'CoCLIP_NES_Arsenite_AAAAA', 'CoCLIP_NES_Arsenite_UUUUU'), colormaps = c('blue', 'blue', 'red', 'red'), linetypes = c(1, 2, 1, 2))
plot_CD(normed_table, y_data = c('CoCLIP_NES_Mock_AAAAA', 'CoCLIP_NES_Mock_UUUUU'), colormaps = c('blue', 'blue'), linetypes = c(1, 2))
plot_CD(normed_table, y_data = c('CoCLIP_NES_Arsenite_AAAAA', 'CoCLIP_NES_Arsenite_UUUUU'), colormaps = c('red', 'red'), linetypes = c(1, 2))

plot_CD(normed_table, y_data = c('CoCLIP_G3BP_Mock_AAAAA', 'CoCLIP_G3BP_Mock_UUUUU', 'CoCLIP_G3BP_Arsenite_AAAAA', 'CoCLIP_G3BP_Arsenite_UUUUU'), colormaps = c('blue', 'blue', 'red', 'red'), linetypes = c(1, 2, 1, 2))
plot_CD(normed_table, y_data = c('CoCLIP_G3BP_Mock_AAAAA', 'CoCLIP_G3BP_Mock_UUUUU'), colormaps = c('blue', 'blue'), linetypes = c(1, 2))
plot_CD(normed_table, y_data = c('CoCLIP_G3BP_Arsenite_AAAAA', 'CoCLIP_G3BP_Arsenite_UUUUU'), colormaps = c('red', 'red'), linetypes = c(1, 2))




## Plot peaks containing the motifs:
# Peak_M_AAAAA = motifCounts_perPeak_M_AAAAA$Var1
# Peak_M_UUUUU = motifCounts_perPeak_M_UUUUU$Var1
# Peak_SG_M_AAAAA = Peak_SG_M %>% filter(peak_names %in% Peak_M_AAAAA)
# Peak_SG_M_UUUUU = Peak_SG_M %>% filter(peak_names %in% Peak_M_UUUUU)
# 
# plotScatter(Peak_SG_M_AAAAA, grouped_annotation,
#             log2(E_NvC_M), log2(G3BP_EvI_M),
#             x_label = 'log2(CoCLIP Nuclear/Cytoplasm)', 'log2(CoCLIP Stress Granule/Input)',
#             x_lim, y_lim,
#             paste0('Mock HuR Peaks: Stress Granule CoCLIP Enrichment (',  nrow(Peak_SG_M), ' peaks)'))

## Plot peaks containing the motifs:
# Peak_S_AAAAA = motifCounts_perPeak_S_AAAAA$Var1
# Peak_S_UUUUU = motifCounts_perPeak_S_UUUUU$Var1
# Peak_SG_S_AAAAA = Peak_SG_S %>% filter(peak_names %in% Peak_S_AAAAA)
# Peak_SG_S_UUUUU = Peak_SG_S %>% filter(peak_names %in% Peak_S_UUUUU)
# 
# plotScatter(Peak_SG_S_AAAAA, grouped_annotation,
#             log2(E_NvC_S), log2(G3BP_EvI_S),
#             x_label = 'log2(CoCLIP Nuclear/Cytoplasm)', 'log2(CoCLIP Stress Granule/Input)',
#             x_lim, y_lim,
#             paste0('Arsenite HuR Peaks: Stress Granule CoCLIP Enrichment (',  nrow(Peak_SG_S), ' peaks)'))


