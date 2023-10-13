## CoCLIP Analysis: 
## Peak Processing for HuR Flow Analysis:
## Written by Soon Yi
## Last Edit: 2023-10-13

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

All_Annotation_List = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", 'ncRNA', "TE", "Other", "DS10K")
mRNA_List = c("5'UTR", "CDS", "3'UTR", "intron", 'CDS_RI', 'DS10K')
ncRNA_List = c('rRNA', 'miRNA', 'lncRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA', 'TE', 'Other', 'nC_RI')
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


## Create dataframes of annotation per gene per gene groups:
####################################################################################################################
dfNames = c('Nuclear_Only_M',
            'Cytoplasm_Only_M',
            'SGranule_Only_M',
            'Nuclear_Only_S',
            'Cytoplasm_Only_S',
            'SGranule_Only_S',
            'Nuclear_Both_MS',
            'Cytoplasm_Both_MS',
            'SGranule_Both_MS',
            'Nuclear_Only_M_Only',
            'Cytoplasm_Only_M_Only',
            'SGranule_Only_M_Only',
            'Nuclear_Only_S_Only',
            'Cytoplasm_Only_S_Only',
            'SGranule_Only_S_Only')

conditions = c('NLS_E_M',
               'NES_E_M',
               'G3BP_E_M',
               'NLS_E_S',
               'NES_E_S',
               'G3BP_E_S')

for (dfName in dfNames) {
  Gene_df_name = paste0('Gene_', dfName)
  Anno_df_name = paste0('Anno_', dfName)
  Gene_df = get(Gene_df_name)
  Peak_Subset = peaksMatrix %>% filter(gene %in% Gene_df$ENSEMBL) 
  
  AnnotationPerGene = Peak_Subset %>% 
    group_by(gene, grouped_annotation) %>% 
    summarise(count =n()) %>% 
    pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)
  
  Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
  columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
  for (column in columns_to_create) {
    column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
    colnames(column_to_bind) = column
    AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
  }
  AnnotationPerGene = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]
  
  assign(Anno_df_name, AnnotationPerGene)
}
####################################################################################################################

mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

## Persists in Nucleus between Mock and Stress
####################################################################################################################
overlapping_gene = Gene_Nuclear_Both_MS$ENSEMBL

## From Nuclear Mock Peaks.
## Find peaks that are in the genes that were in Nuclear during Mock, but also in Stress
## Count annotations per gene.

Peak_Subset = Peak_Co_NLS_M %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_Nuclear_Only_M_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]


## From Nuclear Stress Peaks.
## Find peaks that are in the genes that were in Nuclear during Stress, but also in Mock
## Count annotations per gene.

Peak_Subset = Peak_Co_NLS_S %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_Nuclear_Only_S_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]

write.table(Anno_Nuclear_Only_M_Only %>% arrange(gene), 'M_NLS_S_NLS_MGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
write.table(Anno_Nuclear_Only_S_Only %>% arrange(gene), 'M_NLS_S_NLS_SGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

## Add pseudo counts to the annotation counts per gene:
Anno_Nuclear_Only_M_Only[, Annotations] = Anno_Nuclear_Only_M_Only[, Annotations] + 1e-20
Anno_Nuclear_Only_S_Only[, Annotations] = Anno_Nuclear_Only_S_Only[, Annotations] + 1e-20

## Calcualte ln fold change for peak count per annotation
## 0 = no change
## > 0 greater in Stress
## < 0 greater in Mock
Anno_NucMock_2_NucStress = cbind(Anno_Nuclear_Only_S_Only$gene, log(Anno_Nuclear_Only_S_Only[, Annotations]/Anno_Nuclear_Only_M_Only[, Annotations]))
Anno_NucMock_2_NucStress = Anno_NucMock_2_NucStress %>% arrange(intron)
colnames(Anno_NucMock_2_NucStress) = c('gene', Annotations)

write.table(Anno_NucMock_2_NucStress %>% arrange(gene), 'M_NLS_S_NLS_FC.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

Anno_NucMock_2_NucStress_geneTable = data.frame(gene = Anno_NucMock_2_NucStress$gene)
geneSymbol = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = Anno_NucMock_2_NucStress_geneTable$gene,
                   mart = mart.hs)
Anno_NucMock_2_NucStress_geneTable = cbind(Anno_NucMock_2_NucStress_geneTable, symbol = geneSymbol$external_gene_name[match(Anno_NucMock_2_NucStress_geneTable$gene, geneSymbol$ensembl_gene_id)])
write.table(Anno_NucMock_2_NucStress_geneTable %>% arrange(gene), 'M_NLS_S_NLS_geneTable.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

####################################################################################################################

## Changes from Nuclear Mock to Cytoplasm Stress
####################################################################################################################
## First find overlapping genes in the two group:
overlapping_gene = intersect(Gene_Nuclear_Only_M_Only$ENSEMBL, Gene_Cytoplasm_Only_S_Only$ENSEMBL)

## From Nuclear Mock Peaks.
## Find peaks that are in the genes that were overlapping with Cytoplasm in Stress.
## Count annotations per gene.

Peak_Subset = Peak_Co_NLS_M %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_Nuclear_Only_M_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]


## From Cytoplasm Stress Peaks.
## Find peaks that are in the genes that were overlapping with Nuclear in Mock
## Count annotations per gene.

Peak_Subset = Peak_Co_NES_S %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_Cytoplasm_Only_S_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]

write.table(Anno_Nuclear_Only_M_Only %>% arrange(gene), 'M_NLS_S_NES_MGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
write.table(Anno_Cytoplasm_Only_S_Only %>% arrange(gene), 'M_NLS_S_NES_SGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

## Add pseudo counts to the annotation counts per gene:
Anno_Nuclear_Only_M_Only[, Annotations] = Anno_Nuclear_Only_M_Only[, Annotations] + 1e-20
Anno_Cytoplasm_Only_S_Only[, Annotations] = Anno_Cytoplasm_Only_S_Only[, Annotations] + 1e-20

## Calcualte ln fold change for peak count per annotation
## 0 = no change
## > 0 greater in Cytoplasm
## < 0 greater in Nucleus
Anno_NucMock_2_CytoStress = cbind(Anno_Cytoplasm_Only_S_Only$gene, log(Anno_Cytoplasm_Only_S_Only[, Annotations]/Anno_Nuclear_Only_M_Only[, Annotations]))
Anno_NucMock_2_CytoStress = Anno_NucMock_2_CytoStress %>% arrange(intron)
colnames(Anno_NucMock_2_CytoStress) = c('gene', Annotations)

write.table(Anno_NucMock_2_CytoStress %>% arrange(gene), 'M_NLS_S_NES_FC.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

####################################################################################################################

## Changes from Nuclear Mock to SG Stress
####################################################################################################################
## First find overlapping genes in the two group:
overlapping_gene = intersect(Gene_Nuclear_Only_M_Only$ENSEMBL, Gene_SGranule_Only_S_Only$ENSEMBL)

## From Nuclear Mock Peaks.
## Find peaks that are in the genes that were overlapping with SG in Stress.
## Count annotations per gene.

Peak_Subset = Peak_Co_NLS_M %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_Nuclear_Only_M_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]


## From G3BP Stress Peaks.
## Find peaks that are in the genes that were overlapping with Nuclear in Mock
## Count annotations per gene.

Peak_Subset = Peak_Co_G3BP_S %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_SGranule_Only_S_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]

write.table(Anno_Nuclear_Only_M_Only %>% arrange(gene), 'M_NLS_S_G3BP_MGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
write.table(Anno_SGranule_Only_S_Only %>% arrange(gene), 'M_NLS_S_G3BP_SGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

## Add pseudo counts to the annotation counts per gene:
Anno_Nuclear_Only_M_Only[, Annotations] = Anno_Nuclear_Only_M_Only[, Annotations] + 1e-20
Anno_SGranule_Only_S_Only[, Annotations] = Anno_SGranule_Only_S_Only[, Annotations] + 1e-20

## Calcualte ln fold change for peak count per annotation
## 0 = no change
## > 0 greater in SG
## < 0 greater in Nucleus
Anno_NucMock_2_GStress = cbind(Anno_SGranule_Only_S_Only$gene, log(Anno_SGranule_Only_S_Only[, Annotations]/Anno_Nuclear_Only_M_Only[, Annotations]))
Anno_NucMock_2_GStress = Anno_NucMock_2_GStress %>% arrange(intron)
colnames(Anno_NucMock_2_GStress) = c('gene', Annotations)

write.table(Anno_NucMock_2_GStress %>% arrange(gene), 'M_NLS_S_G3BP_FC.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

Anno_NucMock_2_GStress_geneTable = data.frame(gene = Anno_NucMock_2_GStress$gene)
geneSymbol = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = Anno_NucMock_2_GStress_geneTable$gene,
                   mart = mart.hs)
Anno_NucMock_2_GStress_geneTable = cbind(Anno_NucMock_2_GStress_geneTable, symbol = geneSymbol$external_gene_name[match(Anno_NucMock_2_GStress_geneTable$gene, geneSymbol$ensembl_gene_id)])
write.table(Anno_NucMock_2_GStress_geneTable %>% arrange(gene), 'M_NLS_S_G3BP_geneTable.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
####################################################################################################################




## Changes from Cytoplasm Mock to Nuclear Stress
####################################################################################################################
## First find overlapping genes in the two group:
overlapping_gene = intersect(Gene_Cytoplasm_Only_M_Only$ENSEMBL, Gene_Nuclear_Only_S_Only$ENSEMBL)

## From Cytoplasm Mock Peaks.
## Find peaks that are in the genes that were overlapping with Nuclear in Stress.
## Count annotations per gene.

Peak_Subset = Peak_Co_NES_M %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_Cytoplasm_Only_M_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]


## From Nuclear Stress Peaks.
## Find peaks that are in the genes that were overlapping with Cytoplasm in Mock
## Count annotations per gene.

Peak_Subset = Peak_Co_NLS_S %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_Nuclear_Only_S_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]

write.table(Anno_Cytoplasm_Only_M_Only %>% arrange(gene), 'M_NES_S_NLS_MGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
write.table(Anno_Nuclear_Only_S_Only %>% arrange(gene), 'M_NES_S_NLS_SGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

## Add pseudo counts to the annotation counts per gene:
Anno_Cytoplasm_Only_M_Only[, Annotations] = Anno_Cytoplasm_Only_M_Only[, Annotations] + 1e-20
Anno_Nuclear_Only_S_Only[, Annotations] = Anno_Nuclear_Only_S_Only[, Annotations] + 1e-20

## Calcualte ln fold change for peak count per annotation
## 0 = no change
## > 0 greater in Nucleus
## < 0 greater in Cytoplasm
Anno_CytoMock_2_NucStress = cbind(Anno_Nuclear_Only_S_Only$gene, log(Anno_Nuclear_Only_S_Only[, Annotations]/Anno_Cytoplasm_Only_M_Only[, Annotations]))
Anno_CytoMock_2_NucStress = Anno_CytoMock_2_NucStress %>% arrange(intron)
colnames(Anno_CytoMock_2_NucStress) = c('gene', Annotations)

write.table(Anno_CytoMock_2_NucStress %>% arrange(gene), 'M_NES_S_NLS_FC.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
####################################################################################################################

## Persists in Cytoplasm between Mock and Stress
####################################################################################################################
overlapping_gene = Gene_Cytoplasm_Both_MS$ENSEMBL

## From Cytoplasm Mock Peaks.
## Find peaks that are in the genes that were in Cytoplasm during Mock, but also in Stress
## Count annotations per gene.

Peak_Subset = Peak_Co_NES_M %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_Cytoplasm_Only_M_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]


## From Cytoplasm Stress Peaks.
## Find peaks that are in the genes that were in Cytoplasm during Stress, but also in Mock
## Count annotations per gene.

Peak_Subset = Peak_Co_NES_S %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_Cytoplasm_Only_S_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]

write.table(Anno_Cytoplasm_Only_M_Only %>% arrange(gene), 'M_NES_S_NES_MGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
write.table(Anno_Cytoplasm_Only_S_Only %>% arrange(gene), 'M_NES_S_NES_SGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

## Add pseudo counts to the annotation counts per gene:
Anno_Cytoplasm_Only_M_Only[, Annotations] = Anno_Cytoplasm_Only_M_Only[, Annotations] + 1e-20
Anno_Cytoplasm_Only_S_Only[, Annotations] = Anno_Cytoplasm_Only_S_Only[, Annotations] + 1e-20

## Calcualte ln fold change for peak count per annotation
## 0 = no change
## > 0 greater in Stress
## < 0 greater in Mock
Anno_CytoMock_2_CytoStress = cbind(Anno_Cytoplasm_Only_S_Only$gene, log(Anno_Cytoplasm_Only_S_Only[, Annotations]/Anno_Cytoplasm_Only_M_Only[, Annotations]))
Anno_CytoMock_2_CytoStress = Anno_CytoMock_2_CytoStress %>% arrange(intron)
colnames(Anno_CytoMock_2_CytoStress) = c('gene', Annotations)

write.table(Anno_CytoMock_2_CytoStress %>% arrange(gene), 'M_NES_S_NES_FC.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
####################################################################################################################

## Changes from Cytoplasm Mock to SG Stress
####################################################################################################################
## First find overlapping genes in the two group:
overlapping_gene = intersect(Gene_Cytoplasm_Only_M_Only$ENSEMBL, Gene_SGranule_Only_S_Only$ENSEMBL)

## From Cytoplasm Mock Peaks.
## Find peaks that are in the genes that were overlapping with SG in Stress.
## Count annotations per gene.

Peak_Subset = Peak_Co_NES_M %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_Cytoplasm_Only_M_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]


## From G3BP Stress Peaks.
## Find peaks that are in the genes that were overlapping with Cytoplasm in Mock
## Count annotations per gene.

Peak_Subset = Peak_Co_G3BP_S %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_SGranule_Only_S_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]

write.table(Anno_Cytoplasm_Only_M_Only %>% arrange(gene), 'M_NES_S_G3BP_MGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
write.table(Anno_SGranule_Only_S_Only %>% arrange(gene), 'M_NES_S_G3BP_SGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

## Add pseudo counts to the annotation counts per gene:
Anno_Cytoplasm_Only_M_Only[, Annotations] = Anno_Cytoplasm_Only_M_Only[, Annotations] + 1e-20
Anno_SGranule_Only_S_Only[, Annotations] = Anno_SGranule_Only_S_Only[, Annotations] + 1e-20

## Calcualte ln fold change for peak count per annotation
## 0 = no change
## > 0 greater in SG
## < 0 greater in Cytoplasm
Anno_CytoMock_2_GStress = cbind(Anno_SGranule_Only_S_Only$gene, log(Anno_SGranule_Only_S_Only[, Annotations]/Anno_Cytoplasm_Only_M_Only[, Annotations]))
Anno_CytoMock_2_GStress = Anno_CytoMock_2_GStress %>% arrange(intron)
colnames(Anno_CytoMock_2_GStress) = c('gene', Annotations)

write.table(Anno_CytoMock_2_GStress %>% arrange(gene), 'M_NES_S_G3BP_FC.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
####################################################################################################################




## Changes from SG Mock to Nuclear Stress
####################################################################################################################
## First find overlapping genes in the two group:
overlapping_gene = intersect(Gene_SGranule_Only_M_Only$ENSEMBL, Gene_Nuclear_Only_S_Only$ENSEMBL)

## From SG Mock Peaks.
## Find peaks that are in the genes that were overlapping with Nuclear in Stress.
## Count annotations per gene.

Peak_Subset = Peak_Co_G3BP_M %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_SGranule_Only_M_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]


## From Nuclear Stress Peaks.
## Find peaks that are in the genes that were overlapping with SG in Mock
## Count annotations per gene.

Peak_Subset = Peak_Co_NLS_S %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_Nuclear_Only_S_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]

write.table(Anno_SGranule_Only_M_Only %>% arrange(gene), 'M_G3BP_S_NLS_MGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
write.table(Anno_Nuclear_Only_S_Only %>% arrange(gene), 'M_G3BP_S_NLS_SGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

## Add pseudo counts to the annotation counts per gene:
Anno_SGranule_Only_M_Only[, Annotations] = Anno_SGranule_Only_M_Only[, Annotations] + 1e-20
Anno_Nuclear_Only_S_Only[, Annotations] = Anno_Nuclear_Only_S_Only[, Annotations] + 1e-20

## Calcualte ln fold change for peak count per annotation
## 0 = no change
## > 0 greater in Nucleus
## < 0 greater in SG
Anno_GMock_2_NucStress = cbind(Anno_Nuclear_Only_S_Only$gene, log(Anno_Nuclear_Only_S_Only[, Annotations]/Anno_SGranule_Only_M_Only[, Annotations]))
Anno_GMock_2_NucStress = Anno_GMock_2_NucStress %>% arrange(intron)
colnames(Anno_GMock_2_NucStress) = c('gene', Annotations)

write.table(Anno_GMock_2_NucStress %>% arrange(gene), 'M_G3BP_S_NLS_FC.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
####################################################################################################################

## Changes from SG Mock to Cytoplasm Stress
####################################################################################################################
## First find overlapping genes in the two group:
overlapping_gene = intersect(Gene_SGranule_Only_M_Only$ENSEMBL, Gene_Cytoplasm_Only_S_Only$ENSEMBL)

## From SGranule Mock Peaks.
## Find peaks that are in the genes that were overlapping with Cytoplasm in Stress.
## Count annotations per gene.

Peak_Subset = Peak_Co_G3BP_M %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_SGranule_Only_M_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]


## From Cytoplasm Stress Peaks.
## Find peaks that are in the genes that were overlapping with SGranule in Mock.
## Count annotations per gene.

Peak_Subset = Peak_Co_NES_S %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_Cytoplasm_Only_S_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]

write.table(Anno_SGranule_Only_M_Only %>% arrange(gene), 'M_G3BP_S_NES_MGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
write.table(Anno_Cytoplasm_Only_S_Only %>% arrange(gene), 'M_G3BP_S_NES_SGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

## Add pseudo counts to the annotation counts per gene:
Anno_SGranule_Only_M_Only[, Annotations] = Anno_SGranule_Only_M_Only[, Annotations] + 1e-20
Anno_Cytoplasm_Only_S_Only[, Annotations] = Anno_Cytoplasm_Only_S_Only[, Annotations] + 1e-20

## Calcualte ln fold change for peak count per annotation
## 0 = no change
## > 0 greater in Cytoplasm
## < 0 greater in Nucleus
Anno_GMock_2_CytoStress = cbind(Anno_Cytoplasm_Only_S_Only$gene, log(Anno_Cytoplasm_Only_S_Only[, Annotations]/Anno_SGranule_Only_M_Only[, Annotations]))
Anno_GMock_2_CytoStress = Anno_GMock_2_CytoStress %>% arrange(intron)
colnames(Anno_GMock_2_CytoStress) = c('gene', Annotations)

write.table(Anno_GMock_2_CytoStress %>% arrange(gene), 'M_G3BP_S_NES_FC.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
####################################################################################################################

## Persists in SG between Mock and Stress
####################################################################################################################
## First find overlapping genes in the two group:
overlapping_gene = Gene_SGranule_Both_MS$ENSEMBL

## From Cytoplasm Mock Peaks.
## Find peaks that are in the genes that were overlapping with SG in Stress.
## Count annotations per gene.

Peak_Subset = Peak_Co_G3BP_M %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_SGranule_Only_M_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]


## From G3BP Stress Peaks.
## Find peaks that are in the genes that were only in Stress Granule and only in Stress
## Count annotations per gene.

Peak_Subset = Peak_Co_G3BP_S %>% filter(gene %in% overlapping_gene) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_SGranule_Only_S_Only = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]

write.table(Anno_SGranule_Only_M_Only %>% arrange(gene), 'M_G3BP_S_G3BP_MGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
write.table(Anno_SGranule_Only_S_Only %>% arrange(gene), 'M_G3BP_S_G3BP_SGenes.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

## Add pseudo counts to the annotation counts per gene:
Anno_SGranule_Only_M_Only[, Annotations] = Anno_SGranule_Only_M_Only[, Annotations] + 1e-20
Anno_SGranule_Only_S_Only[, Annotations] = Anno_SGranule_Only_S_Only[, Annotations] + 1e-20

## Calcualte ln fold change for peak count per annotation
## 0 = no change
## > 0 greater in SG
## < 0 greater in Cytoplasm
Anno_GMock_2_GStress = cbind(Anno_SGranule_Only_S_Only$gene, log(Anno_SGranule_Only_S_Only[, Annotations]/Anno_SGranule_Only_M_Only[, Annotations]))
Anno_GMock_2_GStress = Anno_GMock_2_GStress %>% arrange(intron)
colnames(Anno_GMock_2_GStress) = c('gene', Annotations)

write.table(Anno_GMock_2_GStress %>% arrange(gene), 'M_G3BP_S_G3BP_FC.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
####################################################################################################################



Gene_Nuclear_S_deNovo = Gene_Nuclear_Only_S_Only %>% anti_join(Gene_Cytoplasm_Only_M_Only) %>% anti_join(Gene_SGranule_Only_M_Only)
Gene_Cytoplasm_S_deNovo = Gene_Cytoplasm_Only_S_Only %>% anti_join(Gene_Nuclear_Only_M_Only) %>% anti_join(Gene_SGranule_Only_M_Only)
Gene_SGranule_S_deNovo = Gene_SGranule_Only_S_Only %>% anti_join(Gene_Nuclear_Only_M_Only) %>% anti_join(Gene_Cytoplasm_Only_M_Only)

## De novo Nuclear binding in Stress
####################################################################################################################
Peak_Subset = Peak_Co_NLS_S %>% filter(gene %in% Gene_Nuclear_S_deNovo$ENSEMBL) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_Nuclear_S_denovo = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]
write.table(Anno_Nuclear_S_denovo %>% arrange(gene), 'Denovo_S_Nuclear.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

Anno_Nuclear_S_denovo_geneTable = data.frame(gene = Anno_Nuclear_S_denovo$gene)
geneSymbol = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = Anno_Nuclear_S_denovo_geneTable$gene,
                    mart = mart.hs)
Anno_Nuclear_S_denovo_geneTable = cbind(Anno_Nuclear_S_denovo_geneTable, symbol = geneSymbol$external_gene_name[match(Anno_Nuclear_S_denovo_geneTable$gene, geneSymbol$ensembl_gene_id)])
write.table(Anno_Nuclear_S_denovo_geneTable %>% arrange(gene), 'Denovo_S_Nuclear_geneTable.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

  
####################################################################################################################

## De novo Cytoplasm binding in Stress
####################################################################################################################
Peak_Subset = Peak_Co_NES_S %>% filter(gene %in% Gene_Cytoplasm_S_deNovo$ENSEMBL) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_Cytoplasm_S_denovo = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]
write.table(Anno_Cytoplasm_S_denovo %>% arrange(gene), 'Denovo_S_Cytoplasm.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

Anno_Cytoplasm_S_denovo_geneTable = data.frame(gene = Anno_Cytoplasm_S_denovo$gene)
geneSymbol = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = Anno_Cytoplasm_S_denovo_geneTable$gene,
                   mart = mart.hs)
Anno_Cytoplasm_S_denovo_geneTable = cbind(Anno_Cytoplasm_S_denovo_geneTable, symbol = geneSymbol$external_gene_name[match(Anno_Cytoplasm_S_denovo_geneTable$gene, geneSymbol$ensembl_gene_id)])
write.table(Anno_Cytoplasm_S_denovo_geneTable %>% arrange(gene), 'Denovo_S_Cytoplasm_geneTable.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
####################################################################################################################

## De novo SGranule binding in Stress
####################################################################################################################
Peak_Subset = Peak_Co_G3BP_S %>% filter(gene %in% Gene_SGranule_S_deNovo$ENSEMBL) %>% arrange(gene)
AnnotationPerGene = Peak_Subset %>% 
  group_by(gene, grouped_annotation) %>% 
  summarise(count =n()) %>% 
  pivot_wider(names_from = grouped_annotation, values_from = count, values_fill = 0)

Annotations = c("5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")
columns_to_create = Annotations[!Annotations %in% colnames(AnnotationPerGene)]
for (column in columns_to_create) {
  column_to_bind = data.frame(matrix(0, nrow = nrow(AnnotationPerGene), ncol = 1))
  colnames(column_to_bind) = column
  AnnotationPerGene = cbind(AnnotationPerGene, column_to_bind)
}
Anno_SGranule_S_denovo = AnnotationPerGene[, c("gene", "5'UTR", "CDS", "3'UTR", "intron", "DS10K", "snoRNA", "ncRNA", "TE", "Other", "UnAn")]
write.table(Anno_SGranule_S_denovo %>% arrange(gene), 'Denovo_S_SGranule.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

Anno_SGranule_S_denovo_geneTable = data.frame(gene = Anno_SGranule_S_denovo$gene)
geneSymbol = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = Anno_SGranule_S_denovo_geneTable$gene,
                   mart = mart.hs)
Anno_SGranule_S_denovo_geneTable = cbind(Anno_SGranule_S_denovo_geneTable, symbol = geneSymbol$external_gene_name[match(Anno_SGranule_S_denovo_geneTable$gene, geneSymbol$ensembl_gene_id)])
write.table(Anno_SGranule_S_denovo_geneTable %>% arrange(gene), 'Denovo_S_SGranule_geneTable.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
####################################################################################################################

## De novo binding Distribution
####################################################################################################################
denovo_AnnoCounts = data.frame(Annotation = Annotations[1:9])
denovo_AnnoCounts = denovo_AnnoCounts %>% mutate(Nuclear = colSums(Anno_Nuclear_S_denovo[, Annotations[1:9]]))
denovo_AnnoCounts = denovo_AnnoCounts %>% mutate(Cytoplasm = colSums(Anno_Cytoplasm_S_denovo[, Annotations[1:9]]))
denovo_AnnoCounts = denovo_AnnoCounts %>% mutate(SGranule = colSums(Anno_SGranule_S_denovo[, Annotations[1:9]]))
denovo_AnnoCounts = denovo_AnnoCounts %>% gather(key = "Source", value = "Freq", c('Nuclear', 'Cytoplasm', 'SGranule')) 
denovo_AnnoCounts$Source = factor(denovo_AnnoCounts$Source, levels = c('Nuclear', 'Cytoplasm', 'SGranule'))
denovo_AnnoCounts$Annotation = factor(denovo_AnnoCounts$Annotation, levels = Annotations[1:9])
plotStackedBar(denovo_AnnoCounts, c('Nuclear', 'Cytoplasm', 'SGranule'), c('NLS', 'NES', 'G3BP'), y_lim = c(0, 1500), 'CoCLIP Stress De Novo Binding') 

denovo_AnnoCounts = data.frame(Annotation = Annotations[1:9])
denovo_AnnoCounts = denovo_AnnoCounts %>% mutate(Nuclear = colSums(Anno_Nuclear_S_denovo[, Annotations[1:9]])/sum(colSums(Anno_Nuclear_S_denovo[, Annotations[1:9]])))
denovo_AnnoCounts = denovo_AnnoCounts %>% mutate(Cytoplasm = colSums(Anno_Cytoplasm_S_denovo[, Annotations[1:9]])/sum(colSums(Anno_Cytoplasm_S_denovo[, Annotations[1:9]])))
denovo_AnnoCounts = denovo_AnnoCounts %>% mutate(SGranule = colSums(Anno_SGranule_S_denovo[, Annotations[1:9]])/sum(colSums(Anno_SGranule_S_denovo[, Annotations[1:9]])))
denovo_AnnoCounts = denovo_AnnoCounts %>% gather(key = "Source", value = "Freq", c('Nuclear', 'Cytoplasm', 'SGranule')) 
denovo_AnnoCounts$Source = factor(denovo_AnnoCounts$Source, levels = c('Nuclear', 'Cytoplasm', 'SGranule'))
denovo_AnnoCounts$Annotation = factor(denovo_AnnoCounts$Annotation, levels = Annotations[1:9])
plotStackedBar(denovo_AnnoCounts, c('Nuclear', 'Cytoplasm', 'SGranule'), c('NLS', 'NES', 'G3BP'), y_lim = c(0, 1), 'CoCLIP Stress De Novo Binding') 

denovo_AnnoCounts = data.frame(Annotation = Annotations[1:5])
denovo_AnnoCounts = denovo_AnnoCounts %>% mutate(Nuclear = colSums(Anno_Nuclear_S_denovo[, Annotations[1:5]]))
denovo_AnnoCounts = denovo_AnnoCounts %>% mutate(Cytoplasm = colSums(Anno_Cytoplasm_S_denovo[, Annotations[1:5]]))
denovo_AnnoCounts = denovo_AnnoCounts %>% mutate(SGranule = colSums(Anno_SGranule_S_denovo[, Annotations[1:5]]))
denovo_AnnoCounts = denovo_AnnoCounts %>% gather(key = "Source", value = "Freq", c('Nuclear', 'Cytoplasm', 'SGranule')) 
denovo_AnnoCounts$Source = factor(denovo_AnnoCounts$Source, levels = c('Nuclear', 'Cytoplasm', 'SGranule'))
denovo_AnnoCounts$Annotation = factor(denovo_AnnoCounts$Annotation, levels = Annotations[1:9])
plotStackedBar(denovo_AnnoCounts, c('Nuclear', 'Cytoplasm', 'SGranule'), c('NLS', 'NES', 'G3BP'), y_lim = c(0, 1000), 'CoCLIP Stress De Novo Binding') 

denovo_AnnoCounts = data.frame(Annotation = Annotations[1:5])
denovo_AnnoCounts = denovo_AnnoCounts %>% mutate(Nuclear = colSums(Anno_Nuclear_S_denovo[, Annotations[1:5]])/sum(colSums(Anno_Nuclear_S_denovo[, Annotations[1:5]])))
denovo_AnnoCounts = denovo_AnnoCounts %>% mutate(Cytoplasm = colSums(Anno_Cytoplasm_S_denovo[, Annotations[1:5]])/sum(colSums(Anno_Cytoplasm_S_denovo[, Annotations[1:5]])))
denovo_AnnoCounts = denovo_AnnoCounts %>% mutate(SGranule = colSums(Anno_SGranule_S_denovo[, Annotations[1:5]])/sum(colSums(Anno_SGranule_S_denovo[, Annotations[1:5]])))
denovo_AnnoCounts = denovo_AnnoCounts %>% gather(key = "Source", value = "Freq", c('Nuclear', 'Cytoplasm', 'SGranule')) 
denovo_AnnoCounts$Source = factor(denovo_AnnoCounts$Source, levels = c('Nuclear', 'Cytoplasm', 'SGranule'))
denovo_AnnoCounts$Annotation = factor(denovo_AnnoCounts$Annotation, levels = Annotations[1:9])
plotStackedBar(denovo_AnnoCounts, c('Nuclear', 'Cytoplasm', 'SGranule'), c('NLS', 'NES', 'G3BP'), y_lim = c(0, 1), 'CoCLIP Stress De Novo Binding') 

####################################################################################################################


