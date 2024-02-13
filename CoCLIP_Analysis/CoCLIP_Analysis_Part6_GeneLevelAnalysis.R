## CoCLIP Analysis: 
## Peak Processing for Gene Level Analysis:
## Written by Soon Yi
## Last Edit: 2024-01-31

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
library(rlang)

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
####################################################################################################################

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

## Make RowSum Table:
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
####################################################################################################################

## Make Enrichment Table:
####################################################################################################################
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
####################################################################################################################

## Gene Level Fold Change: Nuclear Fraction Mock vs Arsenite
####################################################################################################################
mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

Gene_Fr_Nuc_M = peakRowSum[, c('gene', 'external_gene_name', 'peak_names', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
Gene_Fr_Nuc_M = Gene_Fr_Nuc_M %>% filter(peak_names %in% Peak_F_Nuc_M$peak_names) %>% filter(!is.na(gene))
Gene_Fr_Nuc_M$peak_names = NULL
Gene_Fr_Nuc_M = Gene_Fr_Nuc_M %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_Fr_Nuc_M = Gene_Fr_Nuc_M[, c('gene', 'external_gene_name', 'Nuc_F_M')]

geneLength = getBM(attributes = c("ensembl_gene_id", "transcript_length"),
                   filters = "ensembl_gene_id",
                   values = Gene_Fr_Nuc_M$gene,
                   mart = mart.hs)
Gene_Fr_Nuc_M = cbind(Gene_Fr_Nuc_M, Length = geneLength$transcript_length[match(Gene_Fr_Nuc_M$gene, geneLength$ensembl_gene_id)])
Gene_Fr_Nuc_M = Gene_Fr_Nuc_M %>% mutate(across(all_of('Nuc_F_M'), ~ . / Length*1e6))
Gene_Fr_Nuc_M$Length = NULL

Gene_Fr_Nuc_S = peakRowSum[, c('gene', 'external_gene_name', 'peak_names', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
Gene_Fr_Nuc_S = Gene_Fr_Nuc_S %>% filter(peak_names %in% Peak_F_Nuc_S$peak_names) %>% filter(!is.na(gene))
Gene_Fr_Nuc_S$peak_names = NULL
Gene_Fr_Nuc_S = Gene_Fr_Nuc_S %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_Fr_Nuc_S = Gene_Fr_Nuc_S[, c('gene', 'external_gene_name', 'Nuc_F_S')]

geneLength = getBM(attributes = c("ensembl_gene_id", "transcript_length"),
                   filters = "ensembl_gene_id",
                   values = Gene_Fr_Nuc_S$gene,
                   mart = mart.hs)
Gene_Fr_Nuc_S = cbind(Gene_Fr_Nuc_S, Length = geneLength$transcript_length[match(Gene_Fr_Nuc_S$gene, geneLength$ensembl_gene_id)])
Gene_Fr_Nuc_S = Gene_Fr_Nuc_S %>% mutate(across(all_of('Nuc_F_S'), ~ . / Length*1e6))
Gene_Fr_Nuc_S$Length = NULL

geneOverlap = intersect(Gene_Fr_Nuc_M$gene, Gene_Fr_Nuc_S$gene)
Gene_Fr_Nuc_M_OL = Gene_Fr_Nuc_M %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Fr_Nuc_S_OL = Gene_Fr_Nuc_S %>% filter(gene %in% geneOverlap) %>% arrange(gene)

Gene_Fr_Nuc_SvM = Gene_Fr_Nuc_M_OL[, c('gene', 'external_gene_name')]
Gene_Fr_Nuc_SvM = Gene_Fr_Nuc_SvM %>% mutate(log2FoldChange = log2(Gene_Fr_Nuc_S_OL$Nuc_F_S/Gene_Fr_Nuc_M_OL$Nuc_F_M))

####################################################################################################################

## Gene Level Fold Change: Cyto Fraction Mock vs Arsenite
####################################################################################################################
Gene_Fr_Cyto_M = peakRowSum[, c('gene', 'external_gene_name', 'peak_names', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
Gene_Fr_Cyto_M = Gene_Fr_Cyto_M %>% filter(peak_names %in% Peak_F_Cyt_M$peak_names) %>% filter(!is.na(gene))
Gene_Fr_Cyto_M$peak_names = NULL
Gene_Fr_Cyto_M = Gene_Fr_Cyto_M %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_Fr_Cyto_M = Gene_Fr_Cyto_M[, c('gene', 'external_gene_name', 'Cyto_F_M')]

geneLength = getBM(attributes = c("ensembl_gene_id", "transcript_length"),
                   filters = "ensembl_gene_id",
                   values = Gene_Fr_Cyto_M$gene,
                   mart = mart.hs)
Gene_Fr_Cyto_M = cbind(Gene_Fr_Cyto_M, Length = geneLength$transcript_length[match(Gene_Fr_Cyto_M$gene, geneLength$ensembl_gene_id)])
Gene_Fr_Cyto_M = Gene_Fr_Cyto_M %>% mutate(across(all_of('Cyto_F_M'), ~ . / Length*1e6))
Gene_Fr_Cyto_M$Length = NULL

Gene_Fr_Cyto_S = peakRowSum[, c('gene', 'external_gene_name', 'peak_names', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
Gene_Fr_Cyto_S = Gene_Fr_Cyto_S %>% filter(peak_names %in% Peak_F_Cyt_S$peak_names) %>% filter(!is.na(gene))
Gene_Fr_Cyto_S$peak_names = NULL
Gene_Fr_Cyto_S = Gene_Fr_Cyto_S %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
Gene_Fr_Cyto_S = Gene_Fr_Cyto_S[, c('gene', 'external_gene_name', 'Cyto_F_S')]

geneLength = getBM(attributes = c("ensembl_gene_id", "transcript_length"),
                   filters = "ensembl_gene_id",
                   values = Gene_Fr_Cyto_S$gene,
                   mart = mart.hs)
Gene_Fr_Cyto_S = cbind(Gene_Fr_Cyto_S, Length = geneLength$transcript_length[match(Gene_Fr_Cyto_S$gene, geneLength$ensembl_gene_id)])
Gene_Fr_Cyto_S = Gene_Fr_Cyto_S %>% mutate(across(all_of('Cyto_F_S'), ~ . / Length*1e6))
Gene_Fr_Cyto_S$Length = NULL

geneOverlap = intersect(Gene_Fr_Cyto_M$gene, Gene_Fr_Cyto_S$gene)
Gene_Fr_Cyto_M_OL = Gene_Fr_Cyto_M %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Fr_Cyto_S_OL = Gene_Fr_Cyto_S %>% filter(gene %in% geneOverlap) %>% arrange(gene)

Gene_Fr_Cyto_SvM = Gene_Fr_Cyto_M_OL[, c('gene', 'external_gene_name')]
Gene_Fr_Cyto_SvM = Gene_Fr_Cyto_SvM %>% mutate(log2FoldChange = log2(Gene_Fr_Cyto_S_OL$Cyto_F_S/Gene_Fr_Cyto_M_OL$Cyto_F_M))

####################################################################################################################

## Gene Level Fold Change: Nuclear Mock vs Arsenite
####################################################################################################################
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

## Cytoplasm Mock --> G3BP Mock
geneOverlap = intersect(Gene_Co_NES_M$gene, Gene_Co_G3BP_M$gene)
Gene_Co_NES_M_OL = Gene_Co_NES_M %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Co_G3BP_M_OL = Gene_Co_G3BP_M %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Co_G3BP_MvNES_M = Gene_Co_NES_M_OL[, c('gene', 'external_gene_name')]
Gene_Co_G3BP_MvNES_M = Gene_Co_G3BP_MvNES_M %>% mutate(log2FoldChange = log2(Gene_Co_G3BP_M_OL$G3BP_E_M/Gene_Co_NES_M_OL$NES_E_M))

## Cytoplasm Arsenite --> G3BP Arsenite
geneOverlap = intersect(Gene_Co_NES_S$gene, Gene_Co_G3BP_S$gene)
Gene_Co_NES_S_OL = Gene_Co_NES_S %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Co_G3BP_S_OL = Gene_Co_G3BP_S %>% filter(gene %in% geneOverlap) %>% arrange(gene)
Gene_Co_G3BP_SvNES_S = Gene_Co_NES_S_OL[, c('gene', 'external_gene_name')]
Gene_Co_G3BP_SvNES_S = Gene_Co_G3BP_SvNES_S %>% mutate(log2FoldChange = log2(Gene_Co_G3BP_S_OL$G3BP_E_S/Gene_Co_NES_S_OL$NES_E_S))

####################################################################################################################

## RNASeq Fold Changes:
####################################################################################################################
## Read in DESeq results:
RNASeq_PATH = '/Users/soonyi/Desktop/Genomics/RNASeq/RNASeq_APEXCellLines/V43/'
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

## FIGURE5 CoCLIP vs RNASeq 
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
####################################################################################################################

## Supplement: RNASeq Internal Comparison
####################################################################################################################
## Read RNASeq count table and process it:
RNASeq_PATH = '/Users/soonyi/Desktop/Genomics/RNASeq/RNASeq_APEXCelllines/V43/'
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

RNASeq_Pairwise = RNASeq[, c('ENSEMBL', 'gene')]
RNASeq_Pairwise = RNASeq_Pairwise %>% mutate(Nuclear_Mock = (RNASeq$S9 + RNASeq$S10)/2)
RNASeq_Pairwise = RNASeq_Pairwise %>% mutate(Cytoplasm_Mock = (RNASeq$S1 + RNASeq$S2)/2)
RNASeq_Pairwise = RNASeq_Pairwise %>% mutate(SG_Mock = (RNASeq$S13 + RNASeq$S14)/2)
RNASeq_Pairwise = RNASeq_Pairwise %>% mutate(Nuclear_Arse = (RNASeq$S11 + RNASeq$S12)/2)
RNASeq_Pairwise = RNASeq_Pairwise %>% mutate(Cytoplasm_Arse = (RNASeq$S3 + RNASeq$S4)/2)
RNASeq_Pairwise = RNASeq_Pairwise %>% mutate(SG_Arse = (RNASeq$S15 + RNASeq$S16)/2)

ggplot(RNASeq_Pairwise %>% filter(Nuclear_Mock > median(Nuclear_Mock) & Nuclear_Arse > median(Nuclear_Arse)), aes(x = log10(Nuclear_Mock), y = log10(Nuclear_Arse))) +
  geom_point(pch = 16, size = 3, alpha = 0.25) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) +
  xlim(c(0, 10)) +
  ylim(c(0, 10)) +
  # coord_trans(y ='log10', x='log10') + 
  geom_abline(color = 'red')

cor((RNASeq_Pairwise %>% filter(Nuclear_Mock > median(Nuclear_Mock) & Nuclear_Arse > median(Nuclear_Arse)))$Nuclear_Mock, 
    (RNASeq_Pairwise %>% filter(Nuclear_Mock > median(Nuclear_Mock) & Nuclear_Arse > median(Nuclear_Arse)))$Nuclear_Arse) ^ 2

ggplot(RNASeq_Pairwise %>% filter(Cytoplasm_Mock > median(Cytoplasm_Mock) & Cytoplasm_Arse > median(Cytoplasm_Arse)), aes(x = log10(Cytoplasm_Mock), y = log10(Cytoplasm_Arse))) +
  geom_point(pch = 16, size = 3, alpha = 0.25) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) +
  xlim(c(0, 10)) +
  ylim(c(0, 10)) +
  # coord_trans(y ='log10', x='log10') + 
  geom_abline(color = 'red')

cor((RNASeq_Pairwise %>% filter(Cytoplasm_Mock > median(Cytoplasm_Mock) & Cytoplasm_Arse > median(Cytoplasm_Arse)))$Cytoplasm_Mock, 
    (RNASeq_Pairwise %>% filter(Cytoplasm_Mock > median(Cytoplasm_Mock) & Cytoplasm_Arse > median(Cytoplasm_Arse)))$Cytoplasm_Arse) ^ 2

ggplot(RNASeq_Pairwise %>% filter(SG_Mock > median(SG_Mock) & SG_Arse > median(SG_Arse)), aes(x = log10(SG_Mock), y = log10(SG_Arse))) +
  geom_point(pch = 16, size = 3, alpha = 0.25) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) +
  xlim(c(0, 10)) +
  ylim(c(0, 10)) +
  # coord_trans(y ='log10', x='log10') + 
  geom_abline(color = 'red')

cor((RNASeq_Pairwise %>% filter(SG_Mock > median(SG_Mock) & SG_Arse > median(SG_Arse)))$SG_Mock, 
    (RNASeq_Pairwise %>% filter(SG_Mock > median(SG_Mock) & SG_Arse > median(SG_Arse)))$SG_Arse) ^ 2

ggplot(RNASeq_Pairwise %>% filter(Nuclear_Arse > median(Nuclear_Arse) & Cytoplasm_Arse > median(Cytoplasm_Arse)), aes(x = log10(Nuclear_Arse), y = log10(Cytoplasm_Arse))) +
  geom_point(pch = 16, size = 3, alpha = 0.25) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) +
  xlim(c(0, 10)) +
  ylim(c(0, 10)) +
  # coord_trans(y ='log10', x='log10') + 
  geom_abline(color = 'red')

cor((RNASeq_Pairwise %>% filter(Nuclear_Arse > median(Nuclear_Arse) & Cytoplasm_Arse > median(Cytoplasm_Arse)))$Nuclear_Arse, 
    (RNASeq_Pairwise %>% filter(Nuclear_Arse > median(Nuclear_Arse) & Cytoplasm_Arse > median(Cytoplasm_Arse)))$Cytoplasm_Arse) ^ 2

ggplot(RNASeq_Pairwise %>% filter(Nuclear_Arse > median(Nuclear_Arse) & SG_Arse > median(SG_Arse)), aes(x = log10(Nuclear_Arse), y = log10(SG_Arse))) +
  geom_point(pch = 16, size = 3, alpha = 0.25) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) +
  xlim(c(0, 10)) +
  ylim(c(0, 10)) +
  # coord_trans(y ='log10', x='log10') + 
  geom_abline(color = 'red')

cor((RNASeq_Pairwise %>% filter(Nuclear_Arse > median(Nuclear_Arse) & SG_Arse > median(SG_Arse)))$Nuclear_Arse, 
    (RNASeq_Pairwise %>% filter(Nuclear_Arse > median(Nuclear_Arse) & SG_Arse > median(SG_Arse)))$SG_Arse) ^ 2

ggplot(RNASeq_Pairwise %>% filter(Cytoplasm_Arse > median(Cytoplasm_Arse) & SG_Arse > median(SG_Arse)), aes(x = log10(Cytoplasm_Arse), y = log10(SG_Arse))) +
  geom_point(pch = 16, size = 3, alpha = 0.25) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) +
  xlim(c(0, 10)) +
  ylim(c(0, 10)) +
  # coord_trans(y ='log10', x='log10') + 
  geom_abline(color = 'red')

cor((RNASeq_Pairwise %>% filter(Cytoplasm_Arse > median(Cytoplasm_Arse) & SG_Arse > median(SG_Arse)))$Cytoplasm_Arse, 
    (RNASeq_Pairwise %>% filter(Cytoplasm_Arse > median(Cytoplasm_Arse) & SG_Arse > median(SG_Arse)))$SG_Arse) ^ 2

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

## FIGURE5 Transcript length comparison 
####################################################################################################################
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
  
  mock_Length$transcript_tsl = as.integer(str_extract(mock_Length$transcript_tsl, "(?<=tsl)\\d+"))
  arse_Length$transcript_tsl = as.integer(str_extract(arse_Length$transcript_tsl, "(?<=tsl)\\d+"))
  
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
geneEnrichment = peakRowSum[, c('gene', 'external_gene_name', colnames(peakRowSum)[17:34], BC_columns, rowSum_columns)]
geneEnrichment = geneEnrichment %>% filter(!is.na(gene))
geneEnrichment = geneEnrichment %>% group_by(gene, external_gene_name) %>% summarise_all(sum) %>% ungroup()
####################################################################################################################

## FIGURE5 Location Specific GSEA:
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

## Calculate Score that combines both enrichment score and adjusted p_values:
fgsea_Gene_INP_M = fgsea_Gene_INP_M %>% mutate(newScore = NES * -log10(padj))
fgsea_Gene_INP_S = fgsea_Gene_INP_S %>% mutate(newScore = NES * -log10(padj))
fgsea_Gene_NLS_M = fgsea_Gene_NLS_M %>% mutate(newScore = NES * -log10(padj))
fgsea_Gene_NLS_S = fgsea_Gene_NLS_S %>% mutate(newScore = NES * -log10(padj))
fgsea_Gene_NES_M = fgsea_Gene_NES_M %>% mutate(newScore = NES * -log10(padj))
fgsea_Gene_NES_S = fgsea_Gene_NES_S %>% mutate(newScore = NES * -log10(padj))
fgsea_Gene_G3BP_M = fgsea_Gene_G3BP_M %>% mutate(newScore = NES * -log10(padj))
fgsea_Gene_G3BP_S = fgsea_Gene_G3BP_S %>% mutate(newScore = NES * -log10(padj))

fgsea_Gene_INP_M = fgsea_Gene_INP_M %>% arrange(pathway)
fgsea_Gene_INP_S = fgsea_Gene_INP_S %>% arrange(pathway)
fgsea_Gene_NLS_M = fgsea_Gene_NLS_M %>% arrange(pathway)
fgsea_Gene_NLS_S = fgsea_Gene_NLS_S %>% arrange(pathway)
fgsea_Gene_NES_M = fgsea_Gene_NES_M %>% arrange(pathway)
fgsea_Gene_NES_S = fgsea_Gene_NES_S %>% arrange(pathway)
fgsea_Gene_G3BP_M = fgsea_Gene_G3BP_M %>% arrange(pathway)
fgsea_Gene_G3BP_S = fgsea_Gene_G3BP_S %>% arrange(pathway)

## Filter by Go Terms:
fgsea_Gene_INP_M = fgsea_Gene_INP_M %>% mutate(ontology = unlist(lapply(str_split(fgsea_Gene_INP_M$pathway, '_'), function(x) (str_sub(x[1], 1, 2)))))  %>% filter(ontology == 'GO')
fgsea_Gene_INP_M = fgsea_Gene_INP_M %>% mutate(GO_Subset = unlist(lapply(str_split(fgsea_Gene_INP_M$pathway, '_'), function(x) (str_sub(x[1], 3, 4)))))
fgsea_Gene_INP_M_GO_BP = fgsea_Gene_INP_M %>% filter(GO_Subset == 'BP')
fgsea_Gene_INP_M_GO_CC = fgsea_Gene_INP_M %>% filter(GO_Subset == 'CC')
fgsea_Gene_INP_M_GO_MF = fgsea_Gene_INP_M %>% filter(GO_Subset == 'MF')

fgsea_Gene_INP_S = fgsea_Gene_INP_S %>% mutate(ontology = unlist(lapply(str_split(fgsea_Gene_INP_S$pathway, '_'), function(x) (str_sub(x[1], 1, 2)))))  %>% filter(ontology == 'GO')
fgsea_Gene_INP_S = fgsea_Gene_INP_S %>% mutate(GO_Subset = unlist(lapply(str_split(fgsea_Gene_INP_S$pathway, '_'), function(x) (str_sub(x[1], 3, 4)))))
fgsea_Gene_INP_S_GO_BP = fgsea_Gene_INP_S %>% filter(GO_Subset == 'BP')
fgsea_Gene_INP_S_GO_CC = fgsea_Gene_INP_S %>% filter(GO_Subset == 'CC')
fgsea_Gene_INP_S_GO_MF = fgsea_Gene_INP_S %>% filter(GO_Subset == 'MF')

fgsea_Gene_NLS_M = fgsea_Gene_NLS_M %>% mutate(ontology = unlist(lapply(str_split(fgsea_Gene_NLS_M$pathway, '_'), function(x) (str_sub(x[1], 1, 2)))))  %>% filter(ontology == 'GO')
fgsea_Gene_NLS_M = fgsea_Gene_NLS_M %>% mutate(GO_Subset = unlist(lapply(str_split(fgsea_Gene_NLS_M$pathway, '_'), function(x) (str_sub(x[1], 3, 4)))))
fgsea_Gene_NLS_M_GO_BP = fgsea_Gene_NLS_M %>% filter(GO_Subset == 'BP')
fgsea_Gene_NLS_M_GO_CC = fgsea_Gene_NLS_M %>% filter(GO_Subset == 'CC')
fgsea_Gene_NLS_M_GO_MF = fgsea_Gene_NLS_M %>% filter(GO_Subset == 'MF')

fgsea_Gene_NLS_S = fgsea_Gene_NLS_S %>% mutate(ontology = unlist(lapply(str_split(fgsea_Gene_NLS_S$pathway, '_'), function(x) (str_sub(x[1], 1, 2)))))  %>% filter(ontology == 'GO')
fgsea_Gene_NLS_S = fgsea_Gene_NLS_S %>% mutate(GO_Subset = unlist(lapply(str_split(fgsea_Gene_NLS_S$pathway, '_'), function(x) (str_sub(x[1], 3, 4)))))
fgsea_Gene_NLS_S_GO_BP = fgsea_Gene_NLS_S %>% filter(GO_Subset == 'BP')
fgsea_Gene_NLS_S_GO_CC = fgsea_Gene_NLS_S %>% filter(GO_Subset == 'CC')
fgsea_Gene_NLS_S_GO_MF = fgsea_Gene_NLS_S %>% filter(GO_Subset == 'MF')

fgsea_Gene_NES_M = fgsea_Gene_NES_M %>% mutate(ontology = unlist(lapply(str_split(fgsea_Gene_NES_M$pathway, '_'), function(x) (str_sub(x[1], 1, 2)))))  %>% filter(ontology == 'GO')
fgsea_Gene_NES_M = fgsea_Gene_NES_M %>% mutate(GO_Subset = unlist(lapply(str_split(fgsea_Gene_NES_M$pathway, '_'), function(x) (str_sub(x[1], 3, 4)))))
fgsea_Gene_NES_M_GO_BP = fgsea_Gene_NES_M %>% filter(GO_Subset == 'BP')
fgsea_Gene_NES_M_GO_CC = fgsea_Gene_NES_M %>% filter(GO_Subset == 'CC')
fgsea_Gene_NES_M_GO_MF = fgsea_Gene_NES_M %>% filter(GO_Subset == 'MF')

fgsea_Gene_NES_S = fgsea_Gene_NES_S %>% mutate(ontology = unlist(lapply(str_split(fgsea_Gene_NES_S$pathway, '_'), function(x) (str_sub(x[1], 1, 2)))))  %>% filter(ontology == 'GO')
fgsea_Gene_NES_S = fgsea_Gene_NES_S %>% mutate(GO_Subset = unlist(lapply(str_split(fgsea_Gene_NES_S$pathway, '_'), function(x) (str_sub(x[1], 3, 4)))))
fgsea_Gene_NES_S_GO_BP = fgsea_Gene_NES_S %>% filter(GO_Subset == 'BP')
fgsea_Gene_NES_S_GO_CC = fgsea_Gene_NES_S %>% filter(GO_Subset == 'CC')
fgsea_Gene_NES_S_GO_MF = fgsea_Gene_NES_S %>% filter(GO_Subset == 'MF')

fgsea_Gene_G3BP_M = fgsea_Gene_G3BP_M %>% mutate(ontology = unlist(lapply(str_split(fgsea_Gene_G3BP_M$pathway, '_'), function(x) (str_sub(x[1], 1, 2)))))  %>% filter(ontology == 'GO')
fgsea_Gene_G3BP_M = fgsea_Gene_G3BP_M %>% mutate(GO_Subset = unlist(lapply(str_split(fgsea_Gene_G3BP_M$pathway, '_'), function(x) (str_sub(x[1], 3, 4)))))
fgsea_Gene_G3BP_M_GO_BP = fgsea_Gene_G3BP_M %>% filter(GO_Subset == 'BP')
fgsea_Gene_G3BP_M_GO_CC = fgsea_Gene_G3BP_M %>% filter(GO_Subset == 'CC')
fgsea_Gene_G3BP_M_GO_MF = fgsea_Gene_G3BP_M %>% filter(GO_Subset == 'MF')

fgsea_Gene_G3BP_S = fgsea_Gene_G3BP_S %>% mutate(ontology = unlist(lapply(str_split(fgsea_Gene_G3BP_S$pathway, '_'), function(x) (str_sub(x[1], 1, 2)))))  %>% filter(ontology == 'GO')
fgsea_Gene_G3BP_S = fgsea_Gene_G3BP_S %>% mutate(GO_Subset = unlist(lapply(str_split(fgsea_Gene_G3BP_S$pathway, '_'), function(x) (str_sub(x[1], 3, 4)))))
fgsea_Gene_G3BP_S_GO_BP = fgsea_Gene_G3BP_S %>% filter(GO_Subset == 'BP')
fgsea_Gene_G3BP_S_GO_CC = fgsea_Gene_G3BP_S %>% filter(GO_Subset == 'CC')
fgsea_Gene_G3BP_S_GO_MF = fgsea_Gene_G3BP_S %>% filter(GO_Subset == 'MF')

## Set Threshold for p-value and normalized enrichment score:
padj_T = 0.00001
NES_T = 2

## All GO terms: Filter to statistically significant results:
fgsea_Gene_INP_M_GO_sig = fgsea_Gene_INP_M %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_INP_S_GO_sig = fgsea_Gene_INP_S %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)

fgsea_Gene_NLS_M_GO_sig = fgsea_Gene_NLS_M %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_NLS_S_GO_sig = fgsea_Gene_NLS_S %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_NES_M_GO_sig = fgsea_Gene_NES_M %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_NES_S_GO_sig = fgsea_Gene_NES_S %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_G3BP_M_GO_sig = fgsea_Gene_G3BP_M %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_G3BP_S_GO_sig = fgsea_Gene_G3BP_S %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)

## GO_BP: Filter to statistically significant results:
fgsea_Gene_INP_M_GO_BP_sig = fgsea_Gene_INP_M_GO_BP %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_INP_S_GO_BP_sig = fgsea_Gene_INP_S_GO_BP %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)

fgsea_Gene_NLS_M_GO_BP_sig = fgsea_Gene_NLS_M_GO_BP %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_NLS_S_GO_BP_sig = fgsea_Gene_NLS_S_GO_BP %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_NES_M_GO_BP_sig = fgsea_Gene_NES_M_GO_BP %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_NES_S_GO_BP_sig = fgsea_Gene_NES_S_GO_BP %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_G3BP_M_GO_BP_sig = fgsea_Gene_G3BP_M_GO_BP %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_G3BP_S_GO_BP_sig = fgsea_Gene_G3BP_S_GO_BP %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)

## GO_CC: Filter to statistically significant results:
fgsea_Gene_INP_M_GO_CC_sig = fgsea_Gene_INP_M_GO_CC %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_INP_S_GO_CC_sig = fgsea_Gene_INP_S_GO_CC %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)

fgsea_Gene_NLS_M_GO_CC_sig = fgsea_Gene_NLS_M_GO_CC %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_NLS_S_GO_CC_sig = fgsea_Gene_NLS_S_GO_CC %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_NES_M_GO_CC_sig = fgsea_Gene_NES_M_GO_CC %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_NES_S_GO_CC_sig = fgsea_Gene_NES_S_GO_CC %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_G3BP_M_GO_CC_sig = fgsea_Gene_G3BP_M_GO_CC %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_G3BP_S_GO_CC_sig = fgsea_Gene_G3BP_S_GO_CC %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)

## GO_MF: Filter to statistically significant results:
fgsea_Gene_INP_M_GO_MF_sig = fgsea_Gene_INP_M_GO_MF %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_INP_S_GO_MF_sig = fgsea_Gene_INP_S_GO_MF %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)

fgsea_Gene_NLS_M_GO_MF_sig = fgsea_Gene_NLS_M_GO_MF %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_NLS_S_GO_MF_sig = fgsea_Gene_NLS_S_GO_MF %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_NES_M_GO_MF_sig = fgsea_Gene_NES_M_GO_MF %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_NES_S_GO_MF_sig = fgsea_Gene_NES_S_GO_MF %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_G3BP_M_GO_MF_sig = fgsea_Gene_G3BP_M_GO_MF %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)
fgsea_Gene_G3BP_S_GO_MF_sig = fgsea_Gene_G3BP_S_GO_MF %>% filter(padj <= padj_T & NES >= NES_T) %>% arrange(padj)


fgsea_heatmap_colnames = c('pathway', 'INP_M', 'NLS_M', 'NES_M', 'G3BP_M', 'INP_S', 'NLS_S', 'NES_S', 'G3BP_S')
colSelection = c('padj')

## All GO terms: Create matrix for heatmap
Enriched_Terms = rbind(fgsea_Gene_INP_M_GO_sig, fgsea_Gene_INP_S_GO_sig, 
             fgsea_Gene_NLS_M_GO_sig, fgsea_Gene_NLS_S_GO_sig, 
             fgsea_Gene_NES_M_GO_sig, fgsea_Gene_NES_S_GO_sig, 
             fgsea_Gene_G3BP_M_GO_sig, fgsea_Gene_G3BP_S_GO_sig)
Enriched_Terms = unique(Enriched_Terms$pathway)

Enriched_GO_ALL = data.frame(matrix(NA, nrow = length(Enriched_Terms), ncol = 9))
colnames(Enriched_GO_ALL) = fgsea_heatmap_colnames
Enriched_GO_ALL = Enriched_GO_ALL %>% mutate(pathway = Enriched_Terms)

for (sample in fgsea_heatmap_colnames[2:9]) {
  temp = get(paste0('fgsea_Gene_', sample, '_GO_sig'))
  Enriched_GO_ALL[, sample] = temp[match(Enriched_GO_ALL$pathway, temp$pathway), ..colSelection]
}

## Make Dot Plot for GO Terms:
samplesOrder = fgsea_heatmap_colnames[2:9]

long_data = pivot_longer(Enriched_GO_ALL, cols = -everything('pathway'), names_to = "Type", values_to = "Value")
long_data$Type = factor(long_data$Type, levels = samplesOrder)
long_data$pathway = factor(long_data$pathway, levels = rev(Enriched_GO_ALL$pathway))

ggplot(long_data, aes(x = Type, y = pathway, size = -log(Value), color = Value)) + 
  geom_point() + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=8), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) +
  scale_color_gradientn(colours = rev(brewer.pal(9, "GnBu")), 
                        limits = c(1e-7, 1e-4), 
                        oob = oob_squish) +
  scale_size(range = c(3, 6)) 

####################################################################################################################



















