## CoCLIP Analysis: 
## Peak Processing for Gene Level Analysis:
## Written by Soon Yi
## Last Edit: 2023-09-23

library(stringr)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(scales)
library(biomaRt)
library(ggsignif)
library(fgsea)
library(pheatmap)
library(RColorBrewer)


## Load peak matrix and clean up:
####################################################################################################################
# peaksMatrix_PATH = 'L:/My Drive/CWRU/PhD/Luna Lab/1. coCLIP/Analysis/peaks/'    ## Use this for windows machine
peaksMatrix_PATH = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/'
peaksMatrix_FILE = 'Combined_peakCoverage_groomed_normalized_annotated.txt'
# peaksMatrix_FILE = 'Combined_peakCoverage_groomed_annotated.txt'

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
# peaksMatrix = peaksMatrix %>% mutate_at(c(Nuc_F_M, Nuc_F_S, Cyto_F_M, Cyto_F_S, NLS_I_M, NLS_I_S, NES_I_M, NES_I_S, G3BP_I_M, G3BP_I_S, NLS_E_M, NLS_E_S, NES_E_M, NES_E_S, G3BP_E_M, G3BP_E_S), as.double)
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
plot_CD = function(count_table, y_data, colormaps, linetypes, linesize) {
  plot_data = count_table %>% select(counts, {{y_data}})
  plot_data_long = plot_data %>% pivot_longer(cols = {{y_data}}, names_to = "Sample", values_to = "nCounts")
  plot_data_long$Sample = factor(plot_data_long$Sample, levels = unique(plot_data_long$Sample)) # Use unique() to ensure correct order of levels
  
  plot = ggplot(plot_data_long, aes(x = counts, y = nCounts, color = Sample, linetype = Sample)) +
    geom_line(linewidth = linesize) +
    theme_minimal() +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = 'bold')) +
    scale_x_continuous(trans = 'log2') +
    scale_color_manual(values = colormaps) +
    scale_linetype_manual(values = linetypes)
  
  return(plot)
}

motifCounts = function(peak_matrix, motifs) {
  result_container = list()
  for (motif in motifs) {
    org_motif = motif
    motif = gsub("U", "T", motif)
    motif_positions = sapply(motif, function(motif) gregexpr(paste0("(?=", motif, ")"), peak_matrix$sequence, perl=TRUE))
    motif_positions = lapply(motif_positions, function(x) { attributes(x) <- NULL; x })
    counts = numeric(length(motif_positions))
    for (i in 1:length(motif_positions)) {
      sublist = motif_positions[[i]]
      if (length(sublist) == 1 && sublist[[1]] == -1) {
        counts[i] = 0
      } else {
        counts[i] = length(sublist)
      }
    }
    result_container[[org_motif]] = counts
  }
  result_df = data.frame(result_container)
  colnames(result_df) = motifs
  return(result_df)
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



## Gene Level Fold Change: Nuclear Mock vs Arsenite
####################################################################################################################
mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

Gene_Co_NLS_M = peakRowSum[, c('gene', 'external_gene_name', 'peak_names', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
Gene_Co_NLS_M = Gene_Co_NLS_M %>% filter(peak_names %in% Peak_Co_NLS_M$peak_names) %>% filter(!is.na(gene))
Gene_Co_NLS_M$peak_names = NULL
Gene_Co_NLS_M = Gene_Co_NLS_M %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_Co_NLS_M = Gene_Co_NLS_M[, c('gene', 'external_gene_name', 'NLS_E_M')]

geneLength = getBM(attributes = c("ensembl_gene_id", "transcript_length"),
                   filters = "ensembl_gene_id",
                   values = Gene_Co_NLS_M$gene,
                   mart = mart.hs)
Gene_Co_NLS_M = cbind(Gene_Co_NLS_M, Length = geneLength$transcript_length[match(Gene_Co_NLS_M$gene, geneLength$ensembl_gene_id)])
Gene_Co_NLS_M = Gene_Co_NLS_M %>% mutate(across(all_of('NLS_E_M'), ~ . / Length*1e6))
Gene_Co_NLS_M$Length = NULL

Gene_Co_NLS_S = peakRowSum[, c('gene', 'external_gene_name', 'peak_names', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
Gene_Co_NLS_S = Gene_Co_NLS_S %>% filter(peak_names %in% Peak_Co_NLS_S$peak_names) %>% filter(!is.na(gene))
Gene_Co_NLS_S$peak_names = NULL
Gene_Co_NLS_S = Gene_Co_NLS_S %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_Co_NLS_S = Gene_Co_NLS_S[, c('gene', 'external_gene_name', 'NLS_E_S')]

geneLength = getBM(attributes = c("ensembl_gene_id", "transcript_length"),
                   filters = "ensembl_gene_id",
                   values = Gene_Co_NLS_S$gene,
                   mart = mart.hs)
Gene_Co_NLS_S = cbind(Gene_Co_NLS_S, Length = geneLength$transcript_length[match(Gene_Co_NLS_S$gene, geneLength$ensembl_gene_id)])
Gene_Co_NLS_S = Gene_Co_NLS_S %>% mutate(across(all_of('NLS_E_S'), ~ . / Length*1e6))
Gene_Co_NLS_S$Length = NULL

geneOverlap = intersect(Gene_Co_NLS_M$gene, Gene_Co_NLS_S$gene)
Gene_Co_NLS_M_OL = Gene_Co_NLS_M %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Co_NLS_S_OL = Gene_Co_NLS_S %>% filter(gene %in% geneOverlap) %>% arrange(gene)

Gene_Co_NLS_SvM = Gene_Co_NLS_M_OL[, c('gene', 'external_gene_name')]
Gene_Co_NLS_SvM = Gene_Co_NLS_SvM %>% mutate(log2FoldChange = log2(Gene_Co_NLS_S_OL$NLS_E_S/Gene_Co_NLS_M_OL$NLS_E_M))

####################################################################################################################

## Gene Level Fold Change: Cytoplasm Mock vs Arsenite
####################################################################################################################
mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

Gene_Co_NES_M = peakRowSum[, c('gene', 'external_gene_name', 'peak_names', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
Gene_Co_NES_M = Gene_Co_NES_M %>% filter(peak_names %in% Peak_Co_NES_M$peak_names) %>% filter(!is.na(gene))
Gene_Co_NES_M$peak_names = NULL
Gene_Co_NES_M = Gene_Co_NES_M %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_Co_NES_M = Gene_Co_NES_M[, c('gene', 'external_gene_name', 'NES_E_M')]

geneLength = getBM(attributes = c("ensembl_gene_id", "transcript_length"),
                   filters = "ensembl_gene_id",
                   values = Gene_Co_NES_M$gene,
                   mart = mart.hs)
Gene_Co_NES_M = cbind(Gene_Co_NES_M, Length = geneLength$transcript_length[match(Gene_Co_NES_M$gene, geneLength$ensembl_gene_id)])
Gene_Co_NES_M = Gene_Co_NES_M %>% mutate(across(all_of('NES_E_M'), ~ . / Length*1e6))
Gene_Co_NES_M$Length = NULL

Gene_Co_NES_S = peakRowSum[, c('gene', 'external_gene_name', 'peak_names', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
Gene_Co_NES_S = Gene_Co_NES_S %>% filter(peak_names %in% Peak_Co_NES_S$peak_names) %>% filter(!is.na(gene))
Gene_Co_NES_S$peak_names = NULL
Gene_Co_NES_S = Gene_Co_NES_S %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_Co_NES_S = Gene_Co_NES_S[, c('gene', 'external_gene_name', 'NES_E_S')]

geneLength = getBM(attributes = c("ensembl_gene_id", "transcript_length"),
                   filters = "ensembl_gene_id",
                   values = Gene_Co_NES_S$gene,
                   mart = mart.hs)
Gene_Co_NES_S = cbind(Gene_Co_NES_S, Length = geneLength$transcript_length[match(Gene_Co_NES_S$gene, geneLength$ensembl_gene_id)])
Gene_Co_NES_S = Gene_Co_NES_S %>% mutate(across(all_of('NES_E_S'), ~ . / Length*1e6))
Gene_Co_NES_S$Length = NULL

geneOverlap = intersect(Gene_Co_NES_M$gene, Gene_Co_NES_S$gene)
Gene_Co_NES_M_OL = Gene_Co_NES_M %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Co_NES_S_OL = Gene_Co_NES_S %>% filter(gene %in% geneOverlap) %>% arrange(gene)

Gene_Co_NES_SvM = Gene_Co_NES_M_OL[, c('gene', 'external_gene_name')]
Gene_Co_NES_SvM = Gene_Co_NES_SvM %>% mutate(log2FoldChange = log2(Gene_Co_NES_S_OL$NES_E_S/Gene_Co_NES_M_OL$NES_E_M))

####################################################################################################################

## Gene Level Fold Change: SGranule Mock vs Arsenite
####################################################################################################################
mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

Gene_Co_G3BP_M = peakRowSum[, c('gene', 'external_gene_name', 'peak_names', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
Gene_Co_G3BP_M = Gene_Co_G3BP_M %>% filter(peak_names %in% Peak_Co_G3BP_M$peak_names) %>% filter(!is.na(gene))
Gene_Co_G3BP_M$peak_names = NULL
Gene_Co_G3BP_M = Gene_Co_G3BP_M %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_Co_G3BP_M = Gene_Co_G3BP_M[, c('gene', 'external_gene_name', 'G3BP_E_M')]

geneLength = getBM(attributes = c("ensembl_gene_id", "transcript_length"),
                   filters = "ensembl_gene_id",
                   values = Gene_Co_G3BP_M$gene,
                   mart = mart.hs)
Gene_Co_G3BP_M = cbind(Gene_Co_G3BP_M, Length = geneLength$transcript_length[match(Gene_Co_G3BP_M$gene, geneLength$ensembl_gene_id)])
Gene_Co_G3BP_M = Gene_Co_G3BP_M %>% mutate(across(all_of('G3BP_E_M'), ~ . / Length*1e6))
Gene_Co_G3BP_M$Length = NULL

Gene_Co_G3BP_S = peakRowSum[, c('gene', 'external_gene_name', 'peak_names', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
Gene_Co_G3BP_S = Gene_Co_G3BP_S %>% filter(peak_names %in% Peak_Co_G3BP_S$peak_names) %>% filter(!is.na(gene))
Gene_Co_G3BP_S$peak_names = NULL
Gene_Co_G3BP_S = Gene_Co_G3BP_S %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_Co_G3BP_S = Gene_Co_G3BP_S[, c('gene', 'external_gene_name', 'G3BP_E_S')]

geneLength = getBM(attributes = c("ensembl_gene_id", "transcript_length"),
                   filters = "ensembl_gene_id",
                   values = Gene_Co_G3BP_S$gene,
                   mart = mart.hs)
Gene_Co_G3BP_S = cbind(Gene_Co_G3BP_S, Length = geneLength$transcript_length[match(Gene_Co_G3BP_S$gene, geneLength$ensembl_gene_id)])
Gene_Co_G3BP_S = Gene_Co_G3BP_S %>% mutate(across(all_of('G3BP_E_S'), ~ . / Length*1e6))
Gene_Co_G3BP_S$Length = NULL

geneOverlap = intersect(Gene_Co_G3BP_M$gene, Gene_Co_G3BP_S$gene)
Gene_Co_G3BP_M_OL = Gene_Co_G3BP_M %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Co_G3BP_S_OL = Gene_Co_G3BP_S %>% filter(gene %in% geneOverlap) %>% arrange(gene)

Gene_Co_G3BP_SvM = Gene_Co_G3BP_M_OL[, c('gene', 'external_gene_name')]
Gene_Co_G3BP_SvM = Gene_Co_G3BP_SvM %>% mutate(log2FoldChange = log2(Gene_Co_G3BP_S_OL$G3BP_E_S/Gene_Co_G3BP_M_OL$G3BP_E_M))

####################################################################################################################

## Gene Level Fold Change: Other Comparisons
####################################################################################################################
## Nuclear Mock --> Cytoplasm Arsenite
geneOverlap = intersect(Gene_Co_NLS_M$gene, Gene_Co_NES_S$gene)
Gene_Co_NLS_M_OL = Gene_Co_NLS_M %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Co_NES_S_OL = Gene_Co_NES_S %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Co_NES_SvNLS_M = Gene_Co_NLS_M_OL[, c('gene', 'external_gene_name')]
Gene_Co_NES_SvNLS_M = Gene_Co_NES_SvNLS_M %>% mutate(log2FoldChange = log2(Gene_Co_NES_S_OL$NES_E_S/Gene_Co_NLS_M_OL$NLS_E_M))

## Nuclear Mock --> SGranule Arsenite
geneOverlap = intersect(Gene_Co_NLS_M$gene, Gene_Co_G3BP_S$gene)
Gene_Co_NLS_M_OL = Gene_Co_NLS_M %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Co_G3BP_S_OL = Gene_Co_G3BP_S %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Co_G3BP_SvNLS_M = Gene_Co_NLS_M_OL[, c('gene', 'external_gene_name')]
Gene_Co_G3BP_SvNLS_M = Gene_Co_G3BP_SvNLS_M %>% mutate(log2FoldChange = log2(Gene_Co_G3BP_S_OL$G3BP_E_S/Gene_Co_NLS_M_OL$NLS_E_M))

## Cytoplasm Mock --> Nuclear Arsenite
geneOverlap = intersect(Gene_Co_NES_M$gene, Gene_Co_NLS_S$gene)
Gene_Co_NES_M_OL = Gene_Co_NES_M %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Co_NLS_S_OL = Gene_Co_NLS_S %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Co_NLS_SvNES_M = Gene_Co_NES_M_OL[, c('gene', 'external_gene_name')]
Gene_Co_NLS_SvNES_M = Gene_Co_NLS_SvNES_M %>% mutate(log2FoldChange = log2(Gene_Co_NLS_S_OL$NLS_E_S/Gene_Co_NES_M_OL$NES_E_M))

## Cytoplasm Mock --> SGranule Arsenite
geneOverlap = intersect(Gene_Co_NES_M$gene, Gene_Co_G3BP_S$gene)
Gene_Co_NES_M_OL = Gene_Co_NES_M %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Co_G3BP_S_OL = Gene_Co_G3BP_S %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Co_G3BP_SvNES_M = Gene_Co_NES_M_OL[, c('gene', 'external_gene_name')]
Gene_Co_G3BP_SvNES_M = Gene_Co_G3BP_SvNES_M %>% mutate(log2FoldChange = log2(Gene_Co_G3BP_S_OL$G3BP_E_S/Gene_Co_NES_M_OL$NES_E_M))

####################################################################################################################

## RNASeq Fold Changes:
####################################################################################################################
## Read in DESeq results:
RNASeq_PATH = '/Users/soonyi/Desktop/Genomics/RNASeq/RNASeq_APEXCelllines/'
RNASeq_FILE_All = 'All_DESeq.tsv'
RNASeq_FILE_Nuclear = 'Nuclear_DESeq.tsv'
RNASeq_FILE_Cytoplasm = 'Cytoplasm_DESeq.tsv'
RNASeq_FILE_SGranule = 'SGranule_DESeq.tsv'

## All APEX Arsenite/Mock Comparison
All_DESeq = read_delim(paste0(RNASeq_PATH, RNASeq_FILE_All), show_col_types = FALSE)
geneIDs = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                filters = "external_gene_name",
                values = All_DESeq$gene[!grepl("^ENSG", All_DESeq$gene)],
                mart = mart.hs)

All_DESeq = cbind(All_DESeq, ENSEMBL = geneIDs$ensembl_gene_id[match(All_DESeq$gene, geneIDs$external_gene_name)])
All_DESeq$ENSEMBL[grep("^ENSG", All_DESeq$gene)] = All_DESeq$gene[grep("^ENSG", All_DESeq$gene)]
All_DESeq = All_DESeq[, c('ENSEMBL', 'gene', colnames(All_DESeq)[2:7])]
All_DESeq$gene[grepl("^ENSG", All_DESeq$gene)] = NA

## NLS-APEX Arsenite/Mock Comparison
Nuclear_DESeq = read_delim(paste0(RNASeq_PATH, RNASeq_FILE_Nuclear), show_col_types = FALSE)
geneIDs = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                filters = "external_gene_name",
                values = Nuclear_DESeq$gene[!grepl("^ENSG", Nuclear_DESeq$gene)],
                mart = mart.hs)

Nuclear_DESeq = cbind(Nuclear_DESeq, ENSEMBL = geneIDs$ensembl_gene_id[match(Nuclear_DESeq$gene, geneIDs$external_gene_name)])
Nuclear_DESeq$ENSEMBL[grep("^ENSG", Nuclear_DESeq$gene)] = Nuclear_DESeq$gene[grep("^ENSG", Nuclear_DESeq$gene)]
Nuclear_DESeq = Nuclear_DESeq[, c('ENSEMBL', 'gene', colnames(Nuclear_DESeq)[2:7])]
Nuclear_DESeq$gene[grepl("^ENSG", Nuclear_DESeq$gene)] = NA

## NES-APEX Arsenite/Mock Comparison
Cytoplasm_DESeq = read_delim(paste0(RNASeq_PATH, RNASeq_FILE_Cytoplasm), show_col_types = FALSE)
geneIDs = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                filters = "external_gene_name",
                values = Cytoplasm_DESeq$gene[!grepl("^ENSG", Cytoplasm_DESeq$gene)],
                mart = mart.hs)

Cytoplasm_DESeq = cbind(Cytoplasm_DESeq, ENSEMBL = geneIDs$ensembl_gene_id[match(Cytoplasm_DESeq$gene, geneIDs$external_gene_name)])
Cytoplasm_DESeq$ENSEMBL[grep("^ENSG", Cytoplasm_DESeq$gene)] = Cytoplasm_DESeq$gene[grep("^ENSG", Cytoplasm_DESeq$gene)]
Cytoplasm_DESeq = Cytoplasm_DESeq[, c('ENSEMBL', 'gene', colnames(Cytoplasm_DESeq)[2:7])]
Cytoplasm_DESeq$gene[grepl("^ENSG", Cytoplasm_DESeq$gene)] = NA

## G3BP-APEX Arsenite/Mock Comparison
SGranule_DESeq = read_delim(paste0(RNASeq_PATH, RNASeq_FILE_SGranule), show_col_types = FALSE)
geneIDs = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                filters = "external_gene_name",
                values = SGranule_DESeq$gene[!grepl("^ENSG", SGranule_DESeq$gene)],
                mart = mart.hs)

SGranule_DESeq = cbind(SGranule_DESeq, ENSEMBL = geneIDs$ensembl_gene_id[match(SGranule_DESeq$gene, geneIDs$external_gene_name)])
SGranule_DESeq$ENSEMBL[grep("^ENSG", SGranule_DESeq$gene)] = SGranule_DESeq$gene[grep("^ENSG", SGranule_DESeq$gene)]
SGranule_DESeq = SGranule_DESeq[, c('ENSEMBL', 'gene', colnames(SGranule_DESeq)[2:7])]
SGranule_DESeq$gene[grepl("^ENSG", SGranule_DESeq$gene)] = NA
####################################################################################################################

## CoCLIP vs RNASeq 
####################################################################################################################
## SGranule S v Nucleus M:
All_DESeqSub = All_DESeq %>% filter(!is.na(padj))

geneOverlap = intersect(All_DESeqSub$ENSEMBL, Gene_Co_G3BP_SvNLS_M$gene)
All_DESeq_OL = All_DESeqSub %>% filter(ENSEMBL %in% geneOverlap) %>% arrange(ENSEMBL)
Gene_Co_G3BP_SvNLS_M_OL = Gene_Co_G3BP_SvNLS_M %>% filter(gene %in% geneOverlap) %>% arrange(gene)

plotData = data.frame(Gene = All_DESeq_OL$gene,
                      RNASeq = All_DESeq_OL$log2FoldChange, 
                      CoCLIP = Gene_Co_G3BP_SvNLS_M_OL$log2FoldChange,
                      pAdj = All_DESeq_OL$padj)

ggplot(plotData, aes(x = RNASeq, y = CoCLIP,)) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  geom_density_2d(linewidth = 0.5, colour = 'skyblue2') + 
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) +
  ggtitle(paste0('SG_ars vs Nuc_mock (', nrow(plotData), ' genes)')) +
  xlim(c(-4, 4)) +
  ylim(c(-8, 8))

## SGranule S v Cytoplasm M:
geneOverlap = intersect(All_DESeqSub$ENSEMBL, Gene_Co_G3BP_SvNES_M$gene)
All_DESeq_OL = All_DESeqSub %>% filter(ENSEMBL %in% geneOverlap) %>% arrange(ENSEMBL)
Gene_Co_G3BP_SvNES_M_OL = Gene_Co_G3BP_SvNES_M %>% filter(gene %in% geneOverlap) %>% arrange(gene)

plotData = data.frame(Gene = All_DESeq_OL$gene,
                      RNASeq = All_DESeq_OL$log2FoldChange, 
                      CoCLIP = Gene_Co_G3BP_SvNES_M_OL$log2FoldChange,
                      pAdj = All_DESeq_OL$padj)

ggplot(plotData, aes(x = RNASeq, y = CoCLIP)) +
  geom_point(pch = 16, size = 3, alpha = 0.6) +
  geom_density_2d(linewidth = 0.5, colour = 'darkseagreen2') + 
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) +
  ggtitle(paste0('SG_ars vs Cyto_mock (', nrow(plotData), ' genes)')) +
  xlim(c(-4, 4)) +
  ylim(c(-8, 8))

## SGranule S v SGranule M:
geneOverlap = intersect(All_DESeqSub$ENSEMBL, Gene_Co_G3BP_SvM$gene)
All_DESeq_OL = All_DESeqSub %>% filter(ENSEMBL %in% geneOverlap) %>% arrange(ENSEMBL)
Gene_Co_G3BP_SvM_OL = Gene_Co_G3BP_SvM %>% filter(gene %in% geneOverlap) %>% arrange(gene)

plotData = data.frame(Gene = All_DESeq_OL$gene,
                      RNASeq = All_DESeq_OL$log2FoldChange, 
                      CoCLIP = Gene_Co_G3BP_SvM_OL$log2FoldChange,
                      pAdj = All_DESeq_OL$padj)

ggplot(plotData, aes(x = RNASeq, y = CoCLIP)) +
  geom_point(pch = 16, size = 3, alpha = 0.6) +
  geom_density_2d(linewidth = 0.5, colour = 'salmon') + 
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) +
  ggtitle(paste0('SG_ars vs SG_mock (', nrow(plotData), ' genes)')) +
  xlim(c(-4, 4)) +
  ylim(c(-8, 8))


## Nuclear S v M:
# Nuclear_DESeqSub = Nuclear_DESeq %>% filter(!is.na(pvalue))
# geneOverlap = intersect(Nuclear_DESeqSub$ENSEMBL, Gene_Co_NLS_SvM$gene)
# Nuclear_DESeq_OL = Nuclear_DESeqSub %>% filter(ENSEMBL %in% geneOverlap) %>% arrange(ENSEMBL)
# Gene_Co_NLS_SvM_OL = Gene_Co_NLS_SvM %>% filter(gene %in% geneOverlap) %>% arrange(gene)
# 
# plotData = data.frame(Gene = Nuclear_DESeq_OL$gene,
#                       RNASeq = Nuclear_DESeq_OL$log2FoldChange, 
#                       CoCLIP = Gene_Co_NLS_SvM_OL$log2FoldChange,
#                       pval = Nuclear_DESeq_OL$pvalue)
# 
# ggplot(plotData, aes(x = RNASeq, y = CoCLIP)) +
#   geom_point(pch = 16, size = 3, alpha = 0.2) +
#   # scale_color_distiller(palette = "Greens") +
#   theme_bw() + 
#   theme(axis.text = element_text(size=14), 
#         axis.title = element_text(size=14, face = 'bold'), 
#         legend.text = element_text(size=14)) +
#   ggtitle(c('Nuclear S/M Gene Level Enrichment RNASeq vs CoCLIP ')) +
#   xlim(c(-5, 5)) +
#   ylim(c(-5, 5))

## Cytoplasm S v M:
# geneOverlap = intersect(Cytoplasm_DESeq$ENSEMBL, Gene_Co_NES_SvM$gene)
# Cytoplasm_DESeq_OL = Cytoplasm_DESeq %>% filter(ENSEMBL %in% geneOverlap) %>% arrange(ENSEMBL)
# Gene_Co_NES_SvM_OL = Gene_Co_NES_SvM %>% filter(gene %in% geneOverlap) %>% arrange(gene)
# 
# plotData = data.frame(Gene = Cytoplasm_DESeq_OL$gene,
#                       RNASeq = Cytoplasm_DESeq_OL$log2FoldChange, 
#                       CoCLIP = Gene_Co_NES_SvM_OL$log2FoldChange)
# 
# ggplot(plotData, aes(x = RNASeq, y = CoCLIP)) +
#   geom_point(pch = 16, size = 3, alpha = 0.5) +
#   scale_fill_brewer(palette = "Set3") +
#   theme_bw() + 
#   theme(axis.text = element_text(size=14), 
#         axis.title = element_text(size=14, face = 'bold'), 
#         legend.text = element_text(size=14)) +
#   ggtitle(c('Cytoplasm S/M Gene Level Enrichment RNASeq vs CoCLIP ')) +
#   xlim(c(-5, 5)) +
#   ylim(c(-5, 5))

## SGranule S v M:
# geneOverlap = intersect(SGranule_DESeq$ENSEMBL, Gene_Co_G3BP_SvM$gene)
# SGranule_DESeq_OL = SGranule_DESeq %>% filter(ENSEMBL %in% geneOverlap) %>% arrange(ENSEMBL)
# Gene_Co_G3BP_SvM_OL = Gene_Co_G3BP_SvM %>% filter(gene %in% geneOverlap) %>% arrange(gene)
# 
# plotData = data.frame(Gene = SGranule_DESeq_OL$gene,
#                       RNASeq = SGranule_DESeq_OL$log2FoldChange, 
#                       CoCLIP = Gene_Co_G3BP_SvM_OL$log2FoldChange)
# 
# ggplot(plotData, aes(x = RNASeq, y = CoCLIP)) +
#   geom_density_2d() + 
#   geom_point(pch = 16, size = 3, alpha = 0.5) +
#   # scale_fill_brewer(palette = "Set3") +
#   theme_bw() + 
#   theme(axis.text = element_text(size=14), 
#         axis.title = element_text(size=14, face = 'bold'), 
#         legend.text = element_text(size=14)) +
#   ggtitle(c('SGranule S/M Gene Level Enrichment RNASeq vs CoCLIP ')) +
#   xlim(c(-5, 5)) +
#   ylim(c(-10, 10))
####################################################################################################################

## RNASeq Internal Comparison
####################################################################################################################
RNASeq_Pairwise = RNASeq[, c('ENSEMBL', 'gene')]
RNASeq_Pairwise = RNASeq_Pairwise %>% mutate(Nuclear_Mock = (RNASeq$S9 + RNASeq$S10)/2)
RNASeq_Pairwise = RNASeq_Pairwise %>% mutate(Cytoplasm_Mock = (RNASeq$S1 + RNASeq$S2)/2)
RNASeq_Pairwise = RNASeq_Pairwise %>% mutate(SG_Mock = (RNASeq$S13 + RNASeq$S14)/2)
RNASeq_Pairwise = RNASeq_Pairwise %>% mutate(Nuclear_Arse = (RNASeq$S11 + RNASeq$S12)/2)
RNASeq_Pairwise = RNASeq_Pairwise %>% mutate(Cytoplasm_Arse = (RNASeq$S3 + RNASeq$S4)/2)
RNASeq_Pairwise = RNASeq_Pairwise %>% mutate(SG_Arse = (RNASeq$S15 + RNASeq$S16)/2)

ggplot(RNASeq_Pairwise, aes(x = log10(Nuclear_Mock), y = log10(Cytoplasm_Mock))) +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) +
  xlim(c(0, 10)) +
  ylim(c(0, 10)) +
  # scale_y_continuous(trans = 'log10') +
  # scale_x_continuous(trans = 'log10') +
  geom_abline(color = 'red')

####################################################################################################################



## Localization Specific Genes
####################################################################################################################
Gene_Input_M = distinct(data.frame(ENSEMBL = Peak_Co_Input_M$gene, SYMBOL = Peak_Co_Input_M$external_gene_name))
Gene_Input_S = distinct(data.frame(ENSEMBL = Peak_Co_Input_S$gene, SYMBOL = Peak_Co_Input_S$external_gene_name))

Gene_Nuclear_M = distinct(data.frame(ENSEMBL = Peak_Co_NLS_M$gene, SYMBOL = Peak_Co_NLS_M$external_gene_name))
Gene_Cytoplasm_M = distinct(data.frame(ENSEMBL = Peak_Co_NES_M$gene, SYMBOL = Peak_Co_NES_M$external_gene_name))
Gene_SGranule_M = distinct(data.frame(ENSEMBL = Peak_Co_G3BP_M$gene, SYMBOL = Peak_Co_G3BP_M$external_gene_name))

Gene_Nuclear_S = distinct(data.frame(ENSEMBL = Peak_Co_NLS_S$gene, SYMBOL = Peak_Co_NLS_S$external_gene_name))
Gene_Cytoplasm_S = distinct(data.frame(ENSEMBL = Peak_Co_NES_S$gene, SYMBOL = Peak_Co_NES_S$external_gene_name))
Gene_SGranule_S = distinct(data.frame(ENSEMBL = Peak_Co_G3BP_S$gene, SYMBOL = Peak_Co_G3BP_S$external_gene_name))

Gene_Nuclear_Only_M = Gene_Nuclear_M %>% anti_join(Gene_Cytoplasm_M) %>% anti_join(Gene_SGranule_M)               # 1
Gene_Cytoplasm_Only_M = Gene_Cytoplasm_M %>% anti_join(Gene_Nuclear_M) %>% anti_join(Gene_SGranule_M)             # 2
Gene_SGranule_Only_M = Gene_SGranule_M %>% anti_join(Gene_Nuclear_M) %>% anti_join(Gene_Cytoplasm_M)              # 3

Gene_Nuclear_Only_S = Gene_Nuclear_S %>% anti_join(Gene_Cytoplasm_S) %>% anti_join(Gene_SGranule_S)               # 4
Gene_Cytoplasm_Only_S = Gene_Cytoplasm_S %>% anti_join(Gene_Nuclear_S) %>% anti_join(Gene_SGranule_S)             # 5
Gene_SGranule_Only_S = Gene_SGranule_S %>% anti_join(Gene_Nuclear_S) %>% anti_join(Gene_Cytoplasm_S)              # 6

Gene_Input_Both_MS = Gene_Input_M %>% intersect(Gene_Input_S)
Gene_Nuclear_Both_MS = Gene_Nuclear_Only_M %>% intersect(Gene_Nuclear_Only_S)                                     # 7
Gene_Cytoplasm_Both_MS = Gene_Cytoplasm_Only_M %>% intersect(Gene_Cytoplasm_Only_S)                               # 8
Gene_SGranule_Both_MS = Gene_SGranule_Only_M %>% intersect(Gene_SGranule_Only_S)                                  # 9

Gene_Input_Only_M_Only = Gene_Input_M %>% anti_join(Gene_Input_Both_MS)
Gene_Input_Only_S_Only = Gene_Input_S %>% anti_join(Gene_Input_Both_MS)
Gene_Nuclear_Only_M_Only = Gene_Nuclear_Only_M %>% anti_join(Gene_Nuclear_Both_MS)                                # 10                
Gene_Nuclear_Only_S_Only = Gene_Nuclear_Only_S %>% anti_join(Gene_Nuclear_Both_MS)                                # 11
Gene_Cytoplasm_Only_M_Only = Gene_Cytoplasm_Only_M %>% anti_join(Gene_Cytoplasm_Both_MS)                          # 12   
Gene_Cytoplasm_Only_S_Only = Gene_Cytoplasm_Only_S %>% anti_join(Gene_Cytoplasm_Both_MS)                          # 13
Gene_SGranule_Only_M_Only = Gene_SGranule_Only_M %>% anti_join(Gene_SGranule_Both_MS)                             # 14
Gene_SGranule_Only_S_Only = Gene_SGranule_Only_S %>% anti_join(Gene_SGranule_Both_MS)                             # 15
####################################################################################################################

## Transcript length comparison 
####################################################################################################################
# mart.hs = useMart("ensembl", host = "https://useast.ensembl.org", dataset = "hsapiens_gene_ensembl")
# mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

groups = c('Input', 'Nuclear', 'Cytoplasm', 'SGranule')

gene_Lengths = list()
gene_Lengths_longest = list()
gene_Lengths_tsl1 = list()
transcript_Counts = list()
transcript_Counts_avg = list()
transcript_Counts_med = list()

for (group in groups) {
  mock_DataName = paste0('Gene_', group, '_Only_M_Only')
  arse_DataName = paste0('Gene_', group, '_Only_S_Only')
  
  mock_ListName = paste0('M_', group)
  arse_ListName = paste0('S_', group)
  
  mock_Data = get(mock_DataName)
  arse_Data = get(arse_DataName)
  
  mock_Length = getBM(attributes = c("ensembl_gene_id", "transcript_length", "transcript_tsl", "gene_biotype"),
                     filters = "ensembl_gene_id",
                     values = mock_Data$ENSEMBL,
                     mart = mart.hs)
  
  arse_Length = getBM(attributes = c("ensembl_gene_id", "transcript_length", "transcript_tsl", "gene_biotype"),
                     filters = "ensembl_gene_id",
                     values = arse_Data$ENSEMBL,
                     mart = mart.hs)
  
  ## filter by biotype
  # mock_Length = mock_Length %>% filter(gene_biotype == 'protein_coding')
  # arse_Length = arse_Length %>% filter(gene_biotype == 'protein_coding')
  
  mock_Length$transcript_tsl = as.integer(str_extract(mock_Length$transcript_tsl, "(?<=tsl)\\d+"))
  arse_Length$transcript_tsl = as.integer(str_extract(arse_Length$transcript_tsl, "(?<=tsl)\\d+"))
  
  # mock_Length = mock_Length %>% filter(transcript_tsl != "" & transcript_tsl != 'tslNA')
  # arse_Length = arse_Length %>% filter(transcript_tsl != "" & transcript_tsl != 'tslNA')
  
  mock_Length = mock_Length %>% filter(transcript_tsl != "")
  arse_Length = arse_Length %>% filter(transcript_tsl != "")
  
  gene_Lengths[mock_ListName] = list(mock_Length)
  gene_Lengths[arse_ListName] = list(arse_Length)
  
  gene_Lengths_longest[mock_ListName] = list(mock_Length %>% group_by(ensembl_gene_id) %>% slice(which.max(transcript_length)))
  gene_Lengths_longest[arse_ListName] = list(arse_Length %>% group_by(ensembl_gene_id) %>% slice(which.max(transcript_length)))
  
  gene_Lengths_tsl1[mock_ListName] = list(mock_Length %>% group_by(ensembl_gene_id) %>% slice(which(transcript_tsl == 1 | is.na(transcript_tsl))))
  gene_Lengths_tsl1[arse_ListName] = list(arse_Length %>% group_by(ensembl_gene_id) %>% slice(which(transcript_tsl == 1 | is.na(transcript_tsl))))
  
  transcript_Counts[mock_ListName] = list(mock_Length %>% group_by(ensembl_gene_id) %>% summarize(transcript_counts = n_distinct(transcript_length)))
  transcript_Counts[arse_ListName] = list(arse_Length %>% group_by(ensembl_gene_id) %>% summarize(transcript_counts = n_distinct(transcript_length)))
  
  transcript_Counts_avg[mock_ListName] = mean((mock_Length %>% group_by(ensembl_gene_id) %>% summarize(transcript_counts = n_distinct(transcript_length)))$transcript_counts)
  transcript_Counts_avg[arse_ListName] = mean((arse_Length %>% group_by(ensembl_gene_id) %>% summarize(transcript_counts = n_distinct(transcript_length)))$transcript_counts)
  
  transcript_Counts_med[mock_ListName] = median((mock_Length %>% group_by(ensembl_gene_id) %>% summarize(transcript_counts = n_distinct(transcript_length)))$transcript_counts)
  transcript_Counts_med[arse_ListName] = median((arse_Length %>% group_by(ensembl_gene_id) %>% summarize(transcript_counts = n_distinct(transcript_length)))$transcript_counts)
  
}

# Create a new data frame to store the selected data
plotSelections = c('M_Input', 'M_Nuclear', 'M_Cytoplasm', 'M_SGranule')
plotData = data.frame()

# Loop through the selected lists and extract the data
for (plotSelection in plotSelections) {
  DataSubset = data.frame(matrix(NA, nrow = length(gene_Lengths_tsl1[[plotSelection]]$transcript_length), ncol = 2))
  colnames(DataSubset) = c('transcript_lengths', 'sample_source')
  DataSubset$transcript_lengths = gene_Lengths_tsl1[[plotSelection]]$transcript_length
  DataSubset$sample_source = plotSelection
    
  plotData = rbind(plotData, DataSubset)
}

ggplot(plotData, aes(x = factor(sample_source, levels = plotSelections), y = transcript_lengths, fill = sample_source)) +
  geom_boxplot(notch = T, outlier.shape = NA) +
  labs(title = 'Mock Protein Coding Transcripts Length Comparison', x = "Samples", y = "Transcript Length") + 
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold')) +
  scale_fill_manual(values = c('ivory4', 'skyblue', 'darkseagreen2', 'salmon')) +
  scale_x_discrete(labels = c('Input', 'Nuclear', 'Cytoplasm', 'SGranule')) +
  ylim(c(0, 15000)) +
  geom_signif(comparisons = list(c("M_SGranule", "M_Input")),
              test = c('wilcox.test'),
              map_signif_level = TRUE,
              textsize = 5, y_position = 11000, tip_length = 0) +
  geom_signif(comparisons = list(c("M_SGranule", "M_Nuclear")),
              test = c('wilcox.test'),
              map_signif_level = TRUE,
              textsize = 5, y_position = 10000, tip_length = 0) +
  geom_signif(comparisons = list(c("M_SGranule", "M_Cytoplasm")),
              test = c('wilcox.test'),
              map_signif_level = TRUE,
              textsize = 5, y_position = 9000, tip_length = 0)
  
ks.test(plotData[plotData$sample_source == 'M_SGranule', ]$transcript_lengths, plotData[plotData$sample_source == 'M_Input', ]$transcript_lengths)$p.value
ks.test(plotData[plotData$sample_source == 'M_SGranule', ]$transcript_lengths, plotData[plotData$sample_source == 'M_Nuclear', ]$transcript_lengths)$p.value
ks.test(plotData[plotData$sample_source == 'M_SGranule', ]$transcript_lengths, plotData[plotData$sample_source == 'M_Cytoplasm', ]$transcript_lengths)$p.value

wilcox.test(plotData[plotData$sample_source == 'M_SGranule', ]$transcript_lengths, plotData[plotData$sample_source == 'M_Input', ]$transcript_lengths)$p.value
wilcox.test(plotData[plotData$sample_source == 'M_SGranule', ]$transcript_lengths, plotData[plotData$sample_source == 'M_Nuclear', ]$transcript_lengths)$p.value
wilcox.test(plotData[plotData$sample_source == 'M_SGranule', ]$transcript_lengths, plotData[plotData$sample_source == 'M_Cytoplasm', ]$transcript_lengths)$p.value


# Create a new data frame to store the selected data
plotData = data.frame()
plotSelections = c('S_Input', 'S_Nuclear', 'S_Cytoplasm', 'S_SGranule')

# Loop through the selected lists and extract the data
for (plotSelection in plotSelections) {
  DataSubset = data.frame(matrix(NA, nrow = length(gene_Lengths_tsl1[[plotSelection]]$transcript_length), ncol = 2))
  colnames(DataSubset) = c('transcript_lengths', 'sample_source')
  DataSubset$transcript_lengths = gene_Lengths_tsl1[[plotSelection]]$transcript_length
  DataSubset$sample_source = plotSelection
  
  plotData = rbind(plotData, DataSubset)
}

ggplot(plotData, aes(x = factor(sample_source, levels = plotSelections), y = transcript_lengths, fill = sample_source)) +
  geom_boxplot(notch = T, outlier.shape = NA) +
  labs(title = 'Stress Protein Coding Transcripts Length Comparison', x = "Samples", y = "Transcript Length") + 
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold')) +
  scale_fill_manual(values = c('ivory4', 'skyblue', 'darkseagreen2', 'salmon')) +
  scale_x_discrete(labels = c('Input', 'Nuclear', 'Cytoplasm', 'SGranule')) +
  ylim(c(0, 15000)) +
  geom_signif(comparisons = list(c("S_SGranule", "S_Input")),
              test = c('wilcox.test'),
              map_signif_level = TRUE,
              textsize = 5, y_position = 11000, tip_length = 0) +
  geom_signif(comparisons = list(c("S_SGranule", "S_Nuclear")),
              test = c('wilcox.test'),
              map_signif_level = TRUE,
              textsize = 5, y_position = 10000, tip_length = 0) +
  geom_signif(comparisons = list(c("S_SGranule", "S_Cytoplasm")),
              test = c('wilcox.test'),
              map_signif_level = TRUE,
              textsize = 5, y_position = 9000, tip_length = 0)


ks.test(plotData[plotData$sample_source == 'S_SGranule', ]$transcript_lengths, plotData[plotData$sample_source == 'S_Input', ]$transcript_lengths)$p.value
ks.test(plotData[plotData$sample_source == 'S_SGranule', ]$transcript_lengths, plotData[plotData$sample_source == 'S_Nuclear', ]$transcript_lengths)$p.value
ks.test(plotData[plotData$sample_source == 'S_SGranule', ]$transcript_lengths, plotData[plotData$sample_source == 'S_Cytoplasm', ]$transcript_lengths)$p.value

wilcox.test(plotData[plotData$sample_source == 'S_SGranule', ]$transcript_lengths, plotData[plotData$sample_source == 'S_Input', ]$transcript_lengths)$p.value
wilcox.test(plotData[plotData$sample_source == 'S_SGranule', ]$transcript_lengths, plotData[plotData$sample_source == 'S_Nuclear', ]$transcript_lengths)$p.value
wilcox.test(plotData[plotData$sample_source == 'S_SGranule', ]$transcript_lengths, plotData[plotData$sample_source == 'S_Cytoplasm', ]$transcript_lengths)$p.value

####################################################################################################################





## Load gene sets for GSEA:
####################################################################################################################
H_SET_GMT = gmtPathways('/Users/soonyi/Desktop/Genomics/Annotations/GSEA/h.all.v2023.1.Hs.symbols.gmt')
C2_SET_GMT = gmtPathways('/Users/soonyi/Desktop/Genomics/Annotations/GSEA/c2.all.v2023.1.Hs.symbols.gmt')
C5_SET_GMT = gmtPathways('/Users/soonyi/Desktop/Genomics/Annotations/GSEA/c5.all.v2023.1.Hs.symbols.gmt')
####################################################################################################################

## Make Gene Enrichment Table for GSEA:
####################################################################################################################
# geneEnrichment = peakRowSum %>% filter(grouped_annotation == "intron")
# geneEnrichment = geneEnrichment[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
geneEnrichment = peakRowSum[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
geneEnrichment = geneEnrichment %>% filter(!is.na(gene))
geneEnrichment = geneEnrichment %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()

geneEnrichment = geneEnrichment %>% mutate(NLS_EvI_M = geneEnrichment$NLS_E_M / geneEnrichment$NLS_I_M)
geneEnrichment = geneEnrichment %>% mutate(NES_EvI_M = geneEnrichment$NES_E_M / geneEnrichment$NES_I_M)
geneEnrichment = geneEnrichment %>% mutate(G3BP_EvI_M = geneEnrichment$G3BP_E_M / geneEnrichment$G3BP_I_M)
geneEnrichment = geneEnrichment %>% mutate(NLS_EvI_S = geneEnrichment$NLS_E_S / geneEnrichment$NLS_I_S)
geneEnrichment = geneEnrichment %>% mutate(NES_EvI_S = geneEnrichment$NES_E_S / geneEnrichment$NES_I_S)
geneEnrichment = geneEnrichment %>% mutate(G3BP_EvI_S = geneEnrichment$G3BP_E_S / geneEnrichment$G3BP_I_S)

geneEnrichment = geneEnrichment %>% mutate(Nuc_EvI_M = geneEnrichment$Nuc_F_M / geneEnrichment$I_M)
geneEnrichment = geneEnrichment %>% mutate(Cyto_EvI_M = geneEnrichment$Cyto_F_M / geneEnrichment$I_M)
geneEnrichment = geneEnrichment %>% mutate(Nuc_EvI_S = geneEnrichment$Nuc_F_S / geneEnrichment$I_S)
geneEnrichment = geneEnrichment %>% mutate(Cyto_EvI_S = geneEnrichment$Cyto_F_S / geneEnrichment$I_S)

geneEnrichment = geneEnrichment %>% mutate(E_NvC_M = geneEnrichment$NLS_E_M / geneEnrichment$NES_E_M)
geneEnrichment = geneEnrichment %>% mutate(F_NvC_M = geneEnrichment$Nuc_F_M / geneEnrichment$Cyto_F_M)
geneEnrichment = geneEnrichment %>% mutate(E_NvC_S = geneEnrichment$NLS_E_S / geneEnrichment$NES_E_S)
geneEnrichment = geneEnrichment %>% mutate(F_NvC_S = geneEnrichment$Nuc_F_S / geneEnrichment$Cyto_F_S)

geneEnrichment = geneEnrichment %>% mutate(Input_SvM = geneEnrichment$I_S / geneEnrichment$I_M)
geneEnrichment = geneEnrichment %>% mutate(NLS_SvM = geneEnrichment$NLS_E_S / geneEnrichment$NLS_E_M)
geneEnrichment = geneEnrichment %>% mutate(NES_SvM = geneEnrichment$NES_E_S / geneEnrichment$NES_E_M)
geneEnrichment = geneEnrichment %>% mutate(G3BP_SvM = geneEnrichment$G3BP_E_S / geneEnrichment$G3BP_E_M)

####################################################################################################################

## Location Specific GSEA:
####################################################################################################################
Gene_INP_M = geneEnrichment %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_INP_M = Gene_INP_M %>% mutate(Gene_INP_M = ecdf(I_M)(I_M)*100)
Gene_INP_M = (Gene_INP_M %>% filter(!is.na(external_gene_name)))[, c('external_gene_name', 'Gene_INP_M')] 

Gene_INP_S = geneEnrichment %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_INP_S = Gene_INP_S %>% mutate(Gene_INP_S = ecdf(I_S)(I_S)*100)
Gene_INP_S = (Gene_INP_S %>% filter(!is.na(external_gene_name)))[, c('external_gene_name', 'Gene_INP_S')] 

Gene_NLS_M = geneEnrichment %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_NLS_M = Gene_NLS_M %>% mutate(Gene_NLS_M = ecdf(NLS_E_M)(NLS_E_M)*100)
Gene_NLS_M = (Gene_NLS_M %>% filter(!is.na(external_gene_name)))[, c('external_gene_name', 'Gene_NLS_M')] 

Gene_NLS_S = geneEnrichment %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_NLS_S = Gene_NLS_S %>% mutate(Gene_NLS_S = ecdf(NLS_E_S)(NLS_E_S)*100)
Gene_NLS_S = (Gene_NLS_S %>% filter(!is.na(external_gene_name)))[, c('external_gene_name', 'Gene_NLS_S')] 

Gene_NES_M = geneEnrichment %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_NES_M = Gene_NES_M %>% mutate(Gene_NES_M = ecdf(NES_E_M)(NES_E_M)*100)
Gene_NES_M = (Gene_NES_M %>% filter(!is.na(external_gene_name)))[, c('external_gene_name', 'Gene_NES_M')] 

Gene_NES_S = geneEnrichment %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_NES_S = Gene_NES_S %>% mutate(Gene_NES_S = ecdf(NES_E_S)(NES_E_S)*100)
Gene_NES_S = (Gene_NES_S %>% filter(!is.na(external_gene_name)))[, c('external_gene_name', 'Gene_NES_S')] 

Gene_G3BP_M = geneEnrichment %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_G3BP_M = Gene_G3BP_M %>% mutate(Gene_G3BP_M = ecdf(G3BP_E_M)(G3BP_E_M)*100)
Gene_G3BP_M = (Gene_G3BP_M %>% filter(!is.na(external_gene_name)))[, c('external_gene_name', 'Gene_G3BP_M')] 

Gene_G3BP_S = geneEnrichment %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_G3BP_S = Gene_G3BP_S %>% mutate(Gene_G3BP_S = ecdf(G3BP_E_S)(G3BP_E_S)*100)
Gene_G3BP_S = (Gene_G3BP_S %>% filter(!is.na(external_gene_name)))[, c('external_gene_name', 'Gene_G3BP_S')] 

## Make GSEA object:
GS_Pathways = c('H_SET_GMT', 'C2_SET_GMT', 'C5_SET_GMT')
GS_Pathway = get(GS_Pathways[3])

fgsea_Gene_INP_M = fgsea(pathways = GS_Pathway,
                      stats = setNames(as.vector(Gene_INP_M$Gene_INP_M), Gene_INP_M$external_gene_name),
                      scoreType = 'pos')

fgsea_Gene_NLS_M = fgsea(pathways = GS_Pathway,
                      stats = setNames(as.vector(Gene_NLS_M$Gene_NLS_M), Gene_NLS_M$external_gene_name),
                      scoreType = 'pos')

fgsea_Gene_NES_M = fgsea(pathways = GS_Pathway,
                         stats = setNames(as.vector(Gene_NES_M$Gene_NES_M), Gene_NES_M$external_gene_name),
                         scoreType = 'pos')

fgsea_Gene_G3BP_M = fgsea(pathways = GS_Pathway,
                         stats = setNames(as.vector(Gene_G3BP_M$Gene_G3BP_M), Gene_G3BP_M$external_gene_name),
                         scoreType = 'pos')

fgsea_Gene_INP_S = fgsea(pathways = GS_Pathway,
                      stats = setNames(as.vector(Gene_INP_S$Gene_INP_S), Gene_INP_S$external_gene_name),
                      scoreType = 'pos')

fgsea_Gene_NLS_S = fgsea(pathways = GS_Pathway,
                         stats = setNames(as.vector(Gene_NLS_S$Gene_NLS_S), Gene_NLS_S$external_gene_name),
                         scoreType = 'pos')

fgsea_Gene_NES_S = fgsea(pathways = GS_Pathway,
                         stats = setNames(as.vector(Gene_NES_S$Gene_NES_S), Gene_NES_S$external_gene_name),
                         scoreType = 'pos')

fgsea_Gene_G3BP_S = fgsea(pathways = GS_Pathway,
                          stats = setNames(as.vector(Gene_G3BP_S$Gene_G3BP_S), Gene_G3BP_S$external_gene_name),
                          scoreType = 'pos')

## Calculate Score that combines both enrichment score and adjusted p_values
fgsea_Gene_INP_M = fgsea_Gene_INP_M %>% mutate(newScore = NES * -log10(padj))
fgsea_Gene_INP_S = fgsea_Gene_INP_S %>% mutate(newScore = NES * -log10(padj))
fgsea_Gene_NLS_M = fgsea_Gene_NLS_M %>% mutate(newScore = NES * -log10(padj))
fgsea_Gene_NLS_S = fgsea_Gene_NLS_S %>% mutate(newScore = NES * -log10(padj))
fgsea_Gene_NES_M = fgsea_Gene_NES_M %>% mutate(newScore = NES * -log10(padj))
fgsea_Gene_NES_S = fgsea_Gene_NES_S %>% mutate(newScore = NES * -log10(padj))
fgsea_Gene_G3BP_M = fgsea_Gene_G3BP_M %>% mutate(newScore = NES * -log10(padj))
fgsea_Gene_G3BP_S = fgsea_Gene_G3BP_S %>% mutate(newScore = NES * -log10(padj))

## Filter by GO Terms:
# fgsea_Gene_INP_M = fgsea_Gene_INP_M %>% mutate(ontology = unlist(lapply(str_split(fgsea_Gene_INP_M$pathway, '_'), function(x) (str_sub(x[1], 1, 2)))))
# fgsea_Gene_INP_M = fgsea_Gene_INP_M %>% filter(ontology == 'GO')
# fgsea_Gene_INP_M = fgsea_Gene_INP_M %>% mutate(GO_Subset = unlist(lapply(str_split(fgsea_Gene_INP_M$pathway, '_'), function(x) (str_sub(x[1], 3, 4)))))
# 
# fgsea_Gene_INP_M_GO_BP = fgsea_Gene_INP_M %>% filter(GO_Subset == 'BP')
# fgsea_Gene_INP_M_GO_CC = fgsea_Gene_INP_M %>% filter(GO_Subset == 'CC')
# fgsea_Gene_INP_M_GO_MF = fgsea_Gene_INP_M %>% filter(GO_Subset == 'MF')


## Plot GSEA Results:
topPathwaysUp = fgsea_Gene_INP_M[ES > 0][head(order(-newScore), n = 20), pathway]
topPathwaysDown = fgsea_Gene_INP_M[ES < 0][head(order(-newScore), n = 20), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], 
              setNames(as.vector(Gene_INP_S$Gene_INP_S), Gene_INP_S$external_gene_name), 
              fgsea_Gene_INP_M,
              gseaParam = 0.5)

topPathwaysUp = fgsea_Gene_NLS_M[ES > 0][head(order(-newScore), n = 20), pathway]
topPathwaysDown = fgsea_Gene_NLS_M[ES < 0][head(order(newScore), n = 20), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], 
              setNames(as.vector(Gene_NLS_M$Gene_NLS_M), Gene_NLS_M$external_gene_name), 
              fgsea_Gene_NLS_M,
              gseaParam = 0.5)

topPathwaysUp = fgsea_Gene_NES_M[ES > 0][head(order(-newScore), n = 20), pathway]
topPathwaysDown = fgsea_Gene_NES_M[ES < 0][head(order(-newScore), n = 20), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], 
              setNames(as.vector(Gene_NES_M$Gene_NES_M), Gene_NES_M$external_gene_name), 
              fgsea_Gene_NES_M,
              gseaParam = 0.5)

topPathwaysUp = fgsea_Gene_G3BP_M[ES > 0][head(order(-newScore), n = 20), pathway]
topPathwaysDown = fgsea_Gene_G3BP_M[ES < 0][head(order(-newScore), n = 20), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], 
              setNames(as.vector(Gene_G3BP_M$Gene_G3BP_M), Gene_G3BP_M$external_gene_name), 
              fgsea_Gene_G3BP_M,
              gseaParam = 0.5)

topPathwaysUp = fgsea_Gene_INP_S[ES > 0][head(order(-newScore), n = 20), pathway]
topPathwaysDown = fgsea_Gene_INP_S[ES < 0][head(order(-newScore), n = 20), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], 
              setNames(as.vector(Gene_INP_S$Gene_INP_S), Gene_INP_S$external_gene_name), 
              fgsea_Gene_INP_S,
              gseaParam = 0.5)

topPathwaysUp = fgsea_Gene_NLS_S[ES > 0][head(order(-newScore), n = 20), pathway]
topPathwaysDown = fgsea_Gene_NLS_S[ES < 0][head(order(-newScore), n = 20), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], 
              setNames(as.vector(Gene_NLS_S$Gene_NLS_S), Gene_NLS_S$external_gene_name), 
              fgsea_Gene_NLS_S,
              gseaParam = 0.5)

topPathwaysUp = fgsea_Gene_NES_S[ES > 0][head(order(-newScore), n = 20), pathway]
topPathwaysDown = fgsea_Gene_NES_S[ES < 0][head(order(-newScore), n = 20), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], 
              setNames(as.vector(Gene_NES_S$Gene_NES_S), Gene_NES_S$external_gene_name), 
              fgsea_Gene_NES_S,
              gseaParam = 0.5)

topPathwaysUp = fgsea_Gene_G3BP_S[ES > 0][head(order(-newScore), n = 20), pathway]
topPathwaysDown = fgsea_Gene_G3BP_S[ES < 0][head(order(-newScore), n = 20), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], 
              setNames(as.vector(Gene_G3BP_S$Gene_G3BP_S), Gene_G3BP_S$external_gene_name), 
              fgsea_Gene_G3BP_S,
              gseaParam = 0.5)




samples = c('INP_M', 'NLS_M', 'NES_M', 'G3BP_M', 'INP_S', 'NLS_S', 'NES_S', 'G3BP_S')

for (sample in samples) {
  if (sample == 'INP_M') {
    fgsea_All = fgsea_Gene_INP_M[, c('pathway')]
    fgsea_All$INP_M = ecdf(fgsea_Gene_INP_M$newScore[match(fgsea_All$pathway, fgsea_Gene_INP_M$pathway)])(fgsea_Gene_INP_M$newScore[match(fgsea_All$pathway, fgsea_Gene_INP_M$pathway)])*100
    # fgsea_All$INP_M = fgsea_Gene_INP_M$newScore
  } else {
    fgsea_sample = get(paste0('fgsea_Gene_', sample))
    fgsea_All = cbind(fgsea_All, sample = ecdf(fgsea_sample$newScore[match(fgsea_All$pathway, fgsea_sample$pathway)])(fgsea_sample$newScore[match(fgsea_All$pathway, fgsea_sample$pathway)])*100)
    # fgsea_All = cbind(fgsea_All, sample = fgsea_sample$newScore)
  }
}
colnames(fgsea_All) = c('pathway', samples)
fgsea_All = data.frame(fgsea_All)


fgsea_All = fgsea_All %>% mutate(ontology = unlist(lapply(str_split(fgsea_All$pathway, '_'), function(x) (str_sub(x[1], 1, 2)))))
fgsea_All_GO = fgsea_All %>% filter(ontology == 'GO')
fgsea_All_GO = fgsea_All_GO %>% mutate(GO_Subset = unlist(lapply(str_split(fgsea_All_GO$pathway, '_'), function(x) (str_sub(x[1], 3, 4)))))

fgsea_All_GO_BP = fgsea_All_GO %>% filter(GO_Subset == 'BP')
fgsea_All_GO_CC = fgsea_All_GO %>% filter(GO_Subset == 'CC')
fgsea_All_GO_MF = fgsea_All_GO %>% filter(GO_Subset == 'MF')




samplesOrder = c('INP_M', 'NLS_M', 'NES_M', 'G3BP_M', 'INP_S', 'NLS_S', 'NES_S', 'G3BP_S')
samplesOrder = c('NLS_M', 'NES_M', 'G3BP_M', 'NLS_S', 'NES_S', 'G3BP_S')
samplesOrder = c('INP_M', 'NLS_M', 'NES_M', 'G3BP_M')
samplesOrder = c('INP_S', 'NLS_S', 'NES_S', 'G3BP_S')

row_cluster = T
col_cluster = T

CorrMatrix = cor(fgsea_All[, samplesOrder])
CorrMatrix = matrix(round(CorrMatrix, 2), nrow = length(samplesOrder))
colnames(CorrMatrix) = samplesOrder
rownames(CorrMatrix) = samplesOrder
sampleCorr = pheatmap(CorrMatrix, cluster_rows = row_cluster, cluster_cols = col_cluster, color = colorRampPalette(brewer.pal(9, "RdYlBu"))(100))



breaks = seq(0.85, 1, length.out = 101)

CorrMatrix = cor(fgsea_All_GO[, samplesOrder])
CorrMatrix = matrix(round(CorrMatrix, 2), nrow = length(samplesOrder))
colnames(CorrMatrix) = samplesOrder
rownames(CorrMatrix) = samplesOrder
sampleCorr = pheatmap(CorrMatrix, cluster_rows = row_cluster, cluster_cols = col_cluster, color = rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(100)))

CorrMatrix = cor(fgsea_All_GO_BP[, samplesOrder])
CorrMatrix = matrix(round(CorrMatrix, 2), nrow = length(samplesOrder))
colnames(CorrMatrix) = samplesOrder
rownames(CorrMatrix) = samplesOrder
sampleCorr = pheatmap(CorrMatrix, cluster_rows = row_cluster, cluster_cols = col_cluster, color = rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(100)))
plot(sampleCorr$tree_row)

CorrMatrix = cor(fgsea_All_GO_CC[, samplesOrder])
CorrMatrix = matrix(round(CorrMatrix, 2), nrow = length(samplesOrder))
colnames(CorrMatrix) = samplesOrder
rownames(CorrMatrix) = samplesOrder
sampleCorr = pheatmap(CorrMatrix, cluster_rows = row_cluster, cluster_cols = col_cluster, color = rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(100)))
plot(sampleCorr$tree_row)

CorrMatrix = cor(fgsea_All_GO_MF[, samplesOrder])
CorrMatrix = matrix(round(CorrMatrix, 2), nrow = length(samplesOrder))
colnames(CorrMatrix) = samplesOrder
rownames(CorrMatrix) = samplesOrder
sampleCorr = pheatmap(CorrMatrix, cluster_rows = row_cluster, cluster_cols = col_cluster, color = rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(100)))




ALL_heatmap = pheatmap(fgsea_All[, samplesOrder], 
                       cluster_rows = row_cluster, 
                       cluster_cols = col_cluster,
                       cutree_rows = 3,
                       color = colorRampPalette(brewer.pal(9, "RdYlBu"))(100))

GO_heatmap = pheatmap(fgsea_All_GO[, samplesOrder], 
                      cluster_rows = row_cluster, 
                      cluster_cols = F,
                      color = colorRampPalette(brewer.pal(9, "RdYlBu"))(100))

BP_heatmap = pheatmap(fgsea_All_GO_BP[, samplesOrder], cluster_rows = row_cluster, cluster_cols = col_cluster, color = colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
CC_heatmap = pheatmap(fgsea_All_GO_CC[, samplesOrder], cluster_rows = row_cluster, cluster_cols = col_cluster, color = colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
MF_heatmap = pheatmap(fgsea_All_GO_MF[, samplesOrder], cluster_rows = row_cluster, cluster_cols = col_cluster, color = colorRampPalette(brewer.pal(9, "RdYlBu"))(100))

GO_cluster = cbind(fgsea_All_GO, cluster = cutree(GO_heatmap$tree_row, k = 3))
BP_cluster = cbind(fgsea_All_GO_BP, cluster = cutree(BP_heatmap$tree_row, k = nrow(BP_heatmap$tree_row$merge)))
CC_cluster = cbind(fgsea_All_GO_CC, cluster = cutree(CC_heatmap$tree_row, k = nrow(CC_heatmap$tree_row$merge)))
MF_cluster = cbind(fgsea_All_GO_MF, cluster = cutree(MF_heatmap$tree_row, k = nrow(MF_heatmap$tree_row$merge)))

####################################################################################################################


get the list of significant terms
get the unique list across all and set it as column 1
make columns Input NLS NES SG for Mock and Stress
Fill in the rows, if not 0


















#########################################################################################################################
#######################################################DEPRECATED########################################################
#########################################################################################################################



## Stress Over Mock GSEA:
####################################################################################################################
Input_SvM = (geneEnrichment %>% filter(!is.na(external_gene_name)))[, c('gene', 'external_gene_name', 'Input_SvM')] 
NLS_SvM = (geneEnrichment %>% filter(!is.na(external_gene_name)))[, c('gene', 'external_gene_name', 'NLS_SvM')] 
NES_SvM = (geneEnrichment %>% filter(!is.na(external_gene_name)))[, c('gene', 'external_gene_name', 'NES_SvM')]
G3BP_SvM = (geneEnrichment %>% filter(!is.na(external_gene_name)))[, c('gene', 'external_gene_name', 'G3BP_SvM')]

Input_SvM = Input_SvM %>%
  filter(!(external_gene_name %in% external_gene_name[duplicated(external_gene_name)])) %>%
  filter(!(Input_SvM %in% Input_SvM[duplicated(Input_SvM)]))

NLS_SvM = NLS_SvM %>% 
  filter(!(external_gene_name %in% external_gene_name[duplicated(external_gene_name)])) %>%
  filter(!(NLS_SvM %in% NLS_SvM[duplicated(NLS_SvM)]))

NES_SvM = NES_SvM %>% 
  filter(!(external_gene_name %in% external_gene_name[duplicated(external_gene_name)])) %>%
  filter(!(NES_SvM %in% NES_SvM[duplicated(NES_SvM)]))

G3BP_SvM = G3BP_SvM %>% 
  filter(!(external_gene_name %in% external_gene_name[duplicated(external_gene_name)])) %>%
  filter(!(G3BP_SvM %in% G3BP_SvM[duplicated(G3BP_SvM)]))

## Make GSEA object:
GS_Pathways = c('H_SET_GMT', 'C2_SET_GMT', 'C5_SET_GMT')
GS_Pathway = get(GS_Pathways[3])

fgsea_Input_SvM = fgsea(pathways = GS_Pathway,
                      stats = setNames(as.vector(Input_SvM$Input_SvM), Input_SvM$external_gene_name),
                      nPermSimple = 10000)

fgsea_NLS_SvM = fgsea(pathways = GS_Pathway,
                      stats = setNames(as.vector(NLS_SvM$NLS_SvM), NLS_SvM$external_gene_name),
                      nPermSimple = 10000)

fgsea_NES_SvM = fgsea(pathways = GS_Pathway,
                      stats = setNames(as.vector(NES_SvM$NES_SvM), NES_SvM$external_gene_name),
                      nPermSimple = 10000)

fgsea_G3BP_SvM = fgsea(pathways = GS_Pathway,
                      stats = setNames(as.vector(G3BP_SvM$G3BP_SvM), G3BP_SvM$external_gene_name),
                      nPermSimple = 10000)

## Plot GSEA Results:
topPathwaysUp = fgsea_Input_SvM[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown = fgsea_Input_SvM[ES < 0][head(order(pval), n = 10), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], 
              setNames(as.vector(Input_SvM$Input_SvM), Input_SvM$external_gene_name), 
              fgsea_Input_SvM,
              gseaParam = 0.5)

topPathwaysUp = fgsea_NLS_SvM[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown = fgsea_NLS_SvM[ES < 0][head(order(pval), n = 10), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], 
              setNames(as.vector(NLS_SvM$NLS_SvM), NLS_SvM$external_gene_name), 
              fgsea_NLS_SvM,
              gseaParam = 0.5)

topPathwaysUp = fgsea_NES_SvM[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown = fgsea_NES_SvM[ES < 0][head(order(pval), n = 10), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], 
              setNames(as.vector(NES_SvM$NES_SvM), NES_SvM$external_gene_name), 
              fgsea_NES_SvM, 
              gseaParam = 0.5)

topPathwaysUp = fgsea_G3BP_SvM[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown = fgsea_G3BP_SvM[ES < 0][head(order(pval), n = 10), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], setNames(as.vector(G3BP_SvM$G3BP_SvM), G3BP_SvM$external_gene_name), 
              fgsea_G3BP_SvM,
              gseaParam=0.5)

####################################################################################################################

## Enriched Over Input GSEA:
####################################################################################################################
# NLS_EvI_M = peakRowSum %>% filter((NLS_I_M_BC >= 1 & NLS_I_M > median(NLS_I_M) * rowSum_Multiplier_I) & 
#                                     (NLS_E_M_BC >= BC_Threshold_E & NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_E))
NLS_EvI_M = peakRowSum[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
# NLS_EvI_M = NLS_EvI_M %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
NLS_EvI_M = NLS_EvI_M %>% mutate(NLS_EvI_M = log(NLS_E_M / NLS_I_M))
NLS_EvI_M = NLS_EvI_M %>% group_by(gene, external_gene_name) %>% slice_max(order_by = NLS_EvI_M, n = 1) %>% ungroup()
NLS_EvI_M = (NLS_EvI_M %>% filter(!is.na(external_gene_name)))[, c('gene', 'external_gene_name', 'NLS_EvI_M')] 

# NES_EvI_M = peakRowSum %>% filter((NES_I_M_BC >= 1 & NES_I_M > median(NES_I_M) * rowSum_Multiplier_I) & 
#                                     (NES_E_M_BC >= BC_Threshold_E & NES_E_M > median(NES_E_M) * rowSum_Multiplier_E))
NES_EvI_M = peakRowSum[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
# NES_EvI_M = NES_EvI_M %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
NES_EvI_M = NES_EvI_M %>% mutate(NES_EvI_M = log(NES_E_M / NES_I_M))
NES_EvI_M = NES_EvI_M %>% group_by(gene, external_gene_name) %>% slice_max(order_by = NES_EvI_M, n = 1) %>% ungroup()
NES_EvI_M = (NES_EvI_M %>% filter(!is.na(external_gene_name)))[, c('gene', 'external_gene_name', 'NES_EvI_M')] 

# G3BP_EvI_M = peakRowSum %>% filter((G3BP_I_M_BC >= 2 & G3BP_I_M > median(G3BP_I_M) * rowSum_Multiplier_I) & 
#                                     (G3BP_E_M_BC >= BC_Threshold_E_SG & G3BP_E_M > median(G3BP_E_M) * rowSum_Multiplier_E))
G3BP_EvI_M = peakRowSum[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
# G3BP_EvI_M = G3BP_EvI_M %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
G3BP_EvI_M = G3BP_EvI_M %>% mutate(G3BP_EvI_M = log(G3BP_E_M / G3BP_I_M))
G3BP_EvI_M = G3BP_EvI_M %>% group_by(gene, external_gene_name) %>% slice_max(order_by = G3BP_EvI_M, n = 1) %>% ungroup()
G3BP_EvI_M = (G3BP_EvI_M %>% filter(!is.na(external_gene_name)))[, c('gene', 'external_gene_name', 'G3BP_EvI_M')] 

# NLS_EvI_S = peakRowSum %>% filter((NLS_I_S_BC >= 1 & NLS_I_S > median(NLS_I_S) * rowSum_Multiplier_I) & 
#                                     (NLS_E_S_BC >= BC_Threshold_E & NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_E))
NLS_EvI_S = peakRowSum[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
# NLS_EvI_S = NLS_EvI_S %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
NLS_EvI_S = NLS_EvI_S %>% mutate(NLS_EvI_S = log(NLS_E_S / NLS_I_S))
NLS_EvI_S = NLS_EvI_S %>% group_by(gene, external_gene_name) %>% slice_max(order_by = NLS_EvI_S, n = 1) %>% ungroup()
NLS_EvI_S = (NLS_EvI_S %>% filter(!is.na(external_gene_name)))[, c('gene', 'external_gene_name', 'NLS_EvI_S')] 

# NES_EvI_S = peakRowSum %>% filter((NES_I_S_BC >= 1 & NES_I_S > median(NES_I_S) * rowSum_Multiplier_I) & 
#                                     (NES_E_S_BC >= BC_Threshold_E & NES_E_S > median(NES_E_S) * rowSum_Multiplier_E))
NES_EvI_S = peakRowSum[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
# NES_EvI_S = NES_EvI_S %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
NES_EvI_S = NES_EvI_S %>% mutate(NES_EvI_S = log(NES_E_S / NES_I_S))
NES_EvI_S = NES_EvI_S %>% group_by(gene, external_gene_name) %>% slice_max(order_by = NES_EvI_S, n = 1) %>% ungroup()
NES_EvI_S = (NES_EvI_S %>% filter(!is.na(external_gene_name)))[, c('gene', 'external_gene_name', 'NES_EvI_S')] 

# G3BP_EvI_S = peakRowSum %>% filter((G3BP_I_S_BC >= 2 & G3BP_I_S > median(G3BP_I_S) * rowSum_Multiplier_I) & 
#                                     (G3BP_E_S_BC >= BC_Threshold_E_SG & G3BP_E_S > median(G3BP_E_S) * rowSum_Multiplier_E))
G3BP_EvI_S = peakRowSum[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
# G3BP_EvI_S = G3BP_EvI_S %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
G3BP_EvI_S = G3BP_EvI_S %>% mutate(G3BP_EvI_S = log(G3BP_E_S / G3BP_I_S))
G3BP_EvI_S = G3BP_EvI_S %>% group_by(gene, external_gene_name) %>% slice_max(order_by = G3BP_EvI_S, n = 1) %>% ungroup()
G3BP_EvI_S = (G3BP_EvI_S %>% filter(!is.na(external_gene_name)))[, c('gene', 'external_gene_name', 'G3BP_EvI_S')] 


# NLS_EvI_M = NLS_EvI_M %>% 
#   filter(!(external_gene_name %in% external_gene_name[duplicated(external_gene_name)])) %>%
#   filter(!(NLS_EvI_M %in% NLS_EvI_M[duplicated(NLS_EvI_M)]))
# 
# NES_EvI_M = NES_EvI_M %>% 
#   filter(!(external_gene_name %in% external_gene_name[duplicated(external_gene_name)])) %>%
#   filter(!(NES_EvI_M %in% NES_EvI_M[duplicated(NES_EvI_M)]))
# 
# G3BP_EvI_M = G3BP_EvI_M %>% 
#   filter(!(external_gene_name %in% external_gene_name[duplicated(external_gene_name)])) %>%
#   filter(!(G3BP_EvI_M %in% G3BP_EvI_M[duplicated(G3BP_EvI_M)]))
# 
# NLS_EvI_S = NLS_EvI_S %>% 
#   filter(!(external_gene_name %in% external_gene_name[duplicated(external_gene_name)])) %>%
#   filter(!(NLS_EvI_S %in% NLS_EvI_S[duplicated(NLS_EvI_S)]))
# 
# NES_EvI_S = NES_EvI_S %>% 
#   filter(!(external_gene_name %in% external_gene_name[duplicated(external_gene_name)])) %>%
#   filter(!(NES_EvI_S %in% NES_EvI_S[duplicated(NES_EvI_S)]))
# 
# G3BP_EvI_S = G3BP_EvI_S %>% 
#   filter(!(external_gene_name %in% external_gene_name[duplicated(external_gene_name)])) %>%
#   filter(!(G3BP_EvI_S %in% G3BP_EvI_S[duplicated(G3BP_EvI_S)]))


## Make GSEA object:
GS_Pathways = c('H_SET_GMT', 'C2_SET_GMT', 'C5_SET_GMT')
GS_Pathway = get(GS_Pathways[3])

fgsea_NLS_EvI_M = fgsea(pathways = GS_Pathway,
                      stats = setNames(as.vector(NLS_EvI_M$NLS_EvI_M), NLS_EvI_M$external_gene_name), eps = 0, nPermSimple = 10000)

fgsea_NES_EvI_M = fgsea(pathways = GS_Pathway,
                      stats = setNames(as.vector(NES_EvI_M$NES_EvI_M), NES_EvI_M$external_gene_name), eps = 0, nPermSimple = 10000)

fgsea_G3BP_EvI_M = fgsea(pathways = GS_Pathway,
                       stats = setNames(as.vector(G3BP_EvI_M$G3BP_EvI_M), G3BP_EvI_M$external_gene_name), eps = 0, nPermSimple = 10000)

fgsea_NLS_EvI_S = fgsea(pathways = GS_Pathway,
                        stats = setNames(as.vector(NLS_EvI_S$NLS_EvI_S), NLS_EvI_S$external_gene_name), eps = 0, nPermSimple = 10000)

fgsea_NES_EvI_S = fgsea(pathways = GS_Pathway,
                        stats = setNames(as.vector(NES_EvI_S$NES_EvI_S), NES_EvI_S$external_gene_name), eps = 0, nPermSimple = 10000)

fgsea_G3BP_EvI_S = fgsea(pathways = GS_Pathway,
                         stats = setNames(as.vector(G3BP_EvI_S$G3BP_EvI_S), G3BP_EvI_S$external_gene_name), eps = 0, nPermSimple = 10000)

## Plot GSEA Results:
topPathwaysUp = fgsea_NLS_EvI_M[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown = fgsea_NLS_EvI_M[ES < 0][head(order(pval), n = 10), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], 
              setNames(as.vector(NLS_EvI_M$NLS_EvI_M), NLS_EvI_M$external_gene_name), 
              fgsea_NLS_EvI_M)

topPathwaysUp = fgsea_NES_EvI_M[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown = fgsea_NES_EvI_M[ES < 0][head(order(pval), n = 10), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], 
              setNames(as.vector(NES_EvI_M$NES_EvI_M), NES_EvI_M$external_gene_name), 
              fgsea_NES_EvI_M)

topPathwaysUp = fgsea_G3BP_EvI_M[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown = fgsea_G3BP_EvI_M[ES < 0][head(order(pval), n = 10), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], setNames(as.vector(G3BP_EvI_M$G3BP_EvI_M), G3BP_EvI_M$external_gene_name), 
              fgsea_G3BP_EvI_M)

topPathwaysUp = fgsea_NLS_EvI_S[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown = fgsea_NLS_EvI_S[ES < 0][head(order(pval), n = 10), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], 
              setNames(as.vector(NLS_EvI_S$NLS_EvI_S), NLS_EvI_S$external_gene_name), 
              fgsea_NLS_EvI_S)

topPathwaysUp = fgsea_NES_EvI_S[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown = fgsea_NES_EvI_S[ES < 0][head(order(pval), n = 10), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], 
              setNames(as.vector(NES_EvI_S$NES_EvI_S), NES_EvI_S$external_gene_name), 
              fgsea_NES_EvI_S)

topPathwaysUp = fgsea_G3BP_EvI_S[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown = fgsea_G3BP_EvI_S[ES < 0][head(order(pval), n = 10), pathway]
topPathways = c(topPathwaysUp, 
                rev(topPathwaysDown))
plotGseaTable(GS_Pathway[topPathways], setNames(as.vector(G3BP_EvI_S$G3BP_EvI_S), G3BP_EvI_S$external_gene_name), 
              fgsea_G3BP_EvI_S)

####################################################################################################################

## Write File for PantherDB Analysis:
####################################################################################################################
NLS_EvI_M = peakRowSum %>% filter((NLS_I_M_BC >= 1 & NLS_I_M > median(NLS_I_M) * rowSum_Multiplier_I) &
                                    (NLS_E_M_BC >= BC_Threshold_E & NLS_E_M > median(NLS_E_M) * rowSum_Multiplier_E))
NLS_EvI_M = NLS_EvI_M[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
NLS_EvI_M = NLS_EvI_M %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
NLS_EvI_M = NLS_EvI_M %>% mutate(NLS_EvI_M = ecdf(NLS_E_M / NLS_I_M)(NLS_E_M / NLS_I_M)*100)
NLS_EvI_M = (NLS_EvI_M %>% filter(!is.na(external_gene_name)))[, c('gene', 'NLS_EvI_M')] 

NES_EvI_M = peakRowSum %>% filter((NES_I_M_BC >= 1 & NES_I_M > median(NES_I_M) * rowSum_Multiplier_I) &
                                    (NES_E_M_BC >= BC_Threshold_E & NES_E_M > median(NES_E_M) * rowSum_Multiplier_E))
NES_EvI_M = NES_EvI_M[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
NES_EvI_M = NES_EvI_M %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
NES_EvI_M = NES_EvI_M %>% mutate(NES_EvI_M = ecdf(NES_E_M / NES_I_M)(NES_E_M / NES_I_M)*100)
NES_EvI_M = (NES_EvI_M %>% filter(!is.na(external_gene_name)))[, c('gene', 'NES_EvI_M')] 

G3BP_EvI_M = peakRowSum %>% filter((G3BP_I_M_BC >= 2 & G3BP_I_M > median(G3BP_I_M) * rowSum_Multiplier_I) &
                                    (G3BP_E_M_BC >= BC_Threshold_E_SG & G3BP_E_M > median(G3BP_E_M) * rowSum_Multiplier_E))
G3BP_EvI_M = G3BP_EvI_M[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
G3BP_EvI_M = G3BP_EvI_M %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
G3BP_EvI_M = G3BP_EvI_M %>% mutate(G3BP_EvI_M = ecdf(G3BP_E_M / G3BP_I_M)(G3BP_E_M / G3BP_I_M)*100)
G3BP_EvI_M = (G3BP_EvI_M %>% filter(!is.na(external_gene_name)))[, c('gene', 'G3BP_EvI_M')] 

NLS_EvI_S = peakRowSum %>% filter((NLS_I_S_BC >= 1 & NLS_I_S > median(NLS_I_S) * rowSum_Multiplier_I) &
                                    (NLS_E_S_BC >= BC_Threshold_E & NLS_E_S > median(NLS_E_S) * rowSum_Multiplier_E))
NLS_EvI_S = NLS_EvI_S[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
NLS_EvI_S = NLS_EvI_S %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
NLS_EvI_S = NLS_EvI_S %>% mutate(NLS_EvI_S = ecdf(NLS_E_S / NLS_I_S)(NLS_E_S / NLS_I_S)*100)
NLS_EvI_S = (NLS_EvI_S %>% filter(!is.na(external_gene_name)))[, c('gene', 'NLS_EvI_S')] 

NES_EvI_S = peakRowSum %>% filter((NES_I_S_BC >= 1 & NES_I_S > median(NES_I_S) * rowSum_Multiplier_I) &
                                    (NES_E_S_BC >= BC_Threshold_E & NES_E_S > median(NES_E_S) * rowSum_Multiplier_E))
NES_EvI_S = NES_EvI_S[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
NES_EvI_S = NES_EvI_S %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
NES_EvI_S = NES_EvI_S %>% mutate(NES_EvI_S = ecdf(NES_E_S / NES_I_S)(NES_E_S / NES_I_S)*100)
NES_EvI_S = (NES_EvI_S %>% filter(!is.na(external_gene_name)))[, c('gene', 'NES_EvI_S')] 

G3BP_EvI_S = peakRowSum %>% filter((G3BP_I_S_BC >= 2 & G3BP_I_S > median(G3BP_I_S) * rowSum_Multiplier_I) &
                                    (G3BP_E_S_BC >= BC_Threshold_E_SG & G3BP_E_S > median(G3BP_E_S) * rowSum_Multiplier_E))
G3BP_EvI_S = G3BP_EvI_S[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
G3BP_EvI_S = G3BP_EvI_S %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
G3BP_EvI_S = G3BP_EvI_S %>% mutate(G3BP_EvI_S = ecdf(G3BP_E_S / G3BP_I_S)(G3BP_E_S / G3BP_I_S)*100)
G3BP_EvI_S = (G3BP_EvI_S %>% filter(!is.na(external_gene_name)))[, c('gene', 'G3BP_EvI_S')] 


write.table(NLS_EvI_M, './NLS_EvI_M.txt', sep = '\t', quote = F, row.names = F, col.names = F)
write.table(NES_EvI_M, './NES_EvI_M.txt', sep = '\t', quote = F, row.names = F, col.names = F)
write.table(G3BP_EvI_M, './G3BP_EvI_M.txt', sep = '\t', quote = F, row.names = F, col.names = F)
write.table(NLS_EvI_S, './NLS_EvI_S.txt', sep = '\t', quote = F, row.names = F, col.names = F)
write.table(NES_EvI_S, './NES_EvI_S.txt', sep = '\t', quote = F, row.names = F, col.names = F)
write.table(G3BP_EvI_S, './G3BP_EvI_S.txt', sep = '\t', quote = F, row.names = F, col.names = F)
####################################################################################################################

## Make Gene Enrichment Table:
####################################################################################################################
# geneEnrichment = peakRowSum %>% filter(grouped_annotation == "intron")
# geneEnrichment = geneEnrichment[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
geneEnrichment = peakRowSum[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
geneEnrichment = geneEnrichment %>% filter(!is.na(gene))
geneEnrichment = geneEnrichment %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()

## For each gene, pick annotation with highest frequency as the representative annotation:
AnnotationPerGene = peakRowSum %>%
  group_by(gene, grouped_annotation) %>%
  summarise(count = n()) %>%
  arrange(gene, desc(count)) %>%
  filter(!any(duplicated(count))) %>%
  group_by(gene) %>%
  filter(rank(desc(count)) == 1) 

## Use BiomaRt to gt legnth information of genes.
## For each gene, use longest transcript with TSL 1 or NA as representative length:
mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
geneLength = getBM(attributes = c("ensembl_gene_id", "transcript_length", "transcript_tsl", "gene_biotype"),
                   filters = "ensembl_gene_id",
                   values = geneEnrichment$gene,
                   mart = mart.hs)

geneLength = geneLength %>% group_by(ensembl_gene_id) %>% slice(which.max(transcript_length))
geneLength$transcript_tsl = as.integer(str_extract(geneLength$transcript_tsl, "(?<=tsl)\\d+"))

## Append length information to the enrichment table:
geneEnrichment = cbind(geneEnrichment, Length = geneLength$transcript_length[match(geneEnrichment$gene, geneLength$ensembl_gene_id)])
geneEnrichment = geneEnrichment %>% mutate(across(all_of(colnames(peakRowSum)[17:34]), ~ . / Length*1e6))

## Append representative peak annotation:
geneEnrichment = geneEnrichment %>% filter(gene %in% AnnotationPerGene$gene)
geneEnrichment = cbind(geneEnrichment, MFA = AnnotationPerGene$grouped_annotation[match(geneEnrichment$gene, AnnotationPerGene$gene)])
geneEnrichment = geneEnrichment %>% filter(MFA != 'UnAn')

geneEnrichment = geneEnrichment %>% mutate(NLS_EvI_M = geneEnrichment$NLS_E_M / geneEnrichment$NLS_I_M)
geneEnrichment = geneEnrichment %>% mutate(NES_EvI_M = geneEnrichment$NES_E_M / geneEnrichment$NES_I_M)
geneEnrichment = geneEnrichment %>% mutate(G3BP_EvI_M = geneEnrichment$G3BP_E_M / geneEnrichment$G3BP_I_M)
geneEnrichment = geneEnrichment %>% mutate(NLS_EvI_S = geneEnrichment$NLS_E_S / geneEnrichment$NLS_I_S)
geneEnrichment = geneEnrichment %>% mutate(NES_EvI_S = geneEnrichment$NES_E_S / geneEnrichment$NES_I_S)
geneEnrichment = geneEnrichment %>% mutate(G3BP_EvI_S = geneEnrichment$G3BP_E_S / geneEnrichment$G3BP_I_S)

geneEnrichment = geneEnrichment %>% mutate(Nuc_EvI_M = geneEnrichment$Nuc_F_M / geneEnrichment$I_M)
geneEnrichment = geneEnrichment %>% mutate(Cyto_EvI_M = geneEnrichment$Cyto_F_M / geneEnrichment$I_M)
geneEnrichment = geneEnrichment %>% mutate(Nuc_EvI_S = geneEnrichment$Nuc_F_S / geneEnrichment$I_S)
geneEnrichment = geneEnrichment %>% mutate(Cyto_EvI_S = geneEnrichment$Cyto_F_S / geneEnrichment$I_S)

geneEnrichment = geneEnrichment %>% mutate(E_NvC_M = geneEnrichment$NLS_E_M / geneEnrichment$NES_E_M)
geneEnrichment = geneEnrichment %>% mutate(F_NvC_M = geneEnrichment$Nuc_F_M / geneEnrichment$Cyto_F_M)
geneEnrichment = geneEnrichment %>% mutate(E_NvC_S = geneEnrichment$NLS_E_S / geneEnrichment$NES_E_S)
geneEnrichment = geneEnrichment %>% mutate(F_NvC_S = geneEnrichment$Nuc_F_S / geneEnrichment$Cyto_F_S)

geneEnrichment = geneEnrichment %>% mutate(Input_SvM = geneEnrichment$I_S / geneEnrichment$I_M)
geneEnrichment = geneEnrichment %>% mutate(NLS_SvM = geneEnrichment$NLS_E_S / geneEnrichment$NLS_E_M)
geneEnrichment = geneEnrichment %>% mutate(NES_SvM = geneEnrichment$NES_E_S / geneEnrichment$NES_E_M)
geneEnrichment = geneEnrichment %>% mutate(G3BP_SvM = geneEnrichment$G3BP_E_S / geneEnrichment$G3BP_E_M)

####################################################################################################################

## RNASeq CountTable PRocessing
####################################################################################################################
## Read RNASeq count table and process it:
RNASeq_PATH = '/Users/soonyi/Desktop/Genomics/RNASeq/RNASeq_APEXCelllines/'
RNASeq_FILE = 'countTable.csv'

RNASeq = read_delim(paste0(RNASeq_PATH, RNASeq_FILE), show_col_types = FALSE)
geneIDs = getBM(attributes = c("ensembl_gene_id", "external_gene_name", "transcript_length", "transcript_tsl"),
                filters = "external_gene_name",
                values = RNASeq$gene[!grepl("^ENSG", RNASeq$gene)],
                mart = mart.hs)

RNASeq = cbind(RNASeq, ENSEMBL = geneIDs$ensembl_gene_id[match(RNASeq$gene, geneIDs$external_gene_name)])
RNASeq$ENSEMBL[grep("^ENSG", RNASeq$gene)] = RNASeq$gene[grep("^ENSG", RNASeq$gene)]
RNASeq = RNASeq[, c('ENSEMBL', 'gene', colnames(RNASeq)[2:17])]

geneIDs = geneIDs %>% group_by(ensembl_gene_id) %>% slice(which.max(transcript_length))
geneIDs$transcript_tsl = as.integer(str_extract(geneIDs$transcript_tsl, "(?<=tsl)\\d+"))

RNASeq = cbind(RNASeq, Length = geneIDs$transcript_length[match(RNASeq$ENSEMBL, geneIDs$ensembl_gene_id)])
RNASeq = RNASeq %>% filter(!is.na(Length))

## Calculate TPM:
RNASeq = RNASeq %>% mutate(across(all_of(colnames(RNASeq)[3:18]), ~ . / Length*1e6))

## Set up RNASeq sample table:
sampleTable = data.frame(sample = c('S1', 'S2', 'S3', 'S4', 
                                    'S5', 'S6', 'S7', 'S8', 
                                    'S9', 'S10', 'S11', 'S12', 
                                    'S13', 'S14', 'S15', 'S16'),
                         label = factor(c('Cyto', 'Cyto', 'Cyto', 'Cyto', 
                                          'ER', 'ER', 'ER', 'ER', 
                                          'Nuc', 'Nuc', 'Nuc', 'Nuc', 
                                          'SG', 'SG', 'SG', 'SG')),
                         condition = factor(c('WT', 'WT', 'Arsenite', 'Arsenite', 
                                              'WT', 'WT', 'Arsenite', 'Arsenite', 
                                              'WT', 'WT', 'Arsenite', 'Arsenite', 
                                              'WT', 'WT', 'Arsenite', 'Arsenite')))
####################################################################################################################