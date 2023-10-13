## SubAnalysis

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


library(GenomicRanges)
library(IRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

## (extension + 10nt)*2 centered at peak
# extension = 15
extension = 0
# extension = 0

peaksGR = peaksMatrix[, c('chrom', 'start', 'end', 'strand', 'peak_names')]
peaksGR = peaksGR %>% mutate(chrom = ifelse(chrom == "chrMT", "chrM", chrom))
peaksGR$start = as.integer(peaksGR$start) + 1 - extension
peaksGR$end = as.integer(peaksGR$end) + extension
peaksGR = GRanges(peaksGR)

peaksGR_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38, peaksGR, as.character = TRUE)
peaksGR_seqs = cbind(peaksMatrix[, c('chrom', 'start', 'end', 'strand', 'peak_names', BC_columns, 'grouped_annotation', 'finalized_annotation')], data.frame(sequence = peaksGR_seqs))

ALL_M = peaksGR_seqs %>% filter((Nuc_F_M_BC + Cyto_F_M_BC + NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC + NLS_E_M_BC + NES_E_M_BC + G3BP_E_M_BC) >= 10)
I_M = peaksGR_seqs %>% filter((NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC) >= BC_Threshold_I)
NUC_M = peaksGR_seqs %>% filter(Nuc_F_M_BC >= BC_Threshold_F)
CYT_M = peaksGR_seqs %>% filter(Cyto_F_M_BC >= BC_Threshold_F)
NLS_M = peaksGR_seqs %>% filter(NLS_E_M_BC >= BC_Threshold_E)
NES_M = peaksGR_seqs %>% filter(NES_E_M_BC >= BC_Threshold_E)
G3BP_M = peaksGR_seqs %>% filter(G3BP_E_M_BC >= BC_Threshold_E)

ALL_S = peaksGR_seqs %>% filter((Nuc_F_S_BC + Cyto_F_S_BC + NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC + NLS_E_S_BC + NES_E_S_BC + G3BP_E_S_BC) >= 10)
I_S = peaksGR_seqs %>% filter((NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC) >= BC_Threshold_I)
NUC_S = peaksGR_seqs %>% filter(Nuc_F_S_BC >= BC_Threshold_F)
CYT_S = peaksGR_seqs %>% filter(Cyto_F_S_BC >= BC_Threshold_F)
NLS_S = peaksGR_seqs %>% filter(NLS_E_S_BC >= BC_Threshold_E)
NES_S = peaksGR_seqs %>% filter(NES_E_S_BC >= BC_Threshold_E)
G3BP_S = peaksGR_seqs %>% filter(G3BP_E_S_BC >= BC_Threshold_E)

Samples = c('ALL_M', 'I_M', 'NUC_M', 'CYT_M', 'NLS_M', 'NES_M', 'G3BP_M', 'ALL_S', 'I_S', 'NUC_S', 'CYT_S', 'NLS_S', 'NES_S', 'G3BP_S') 

nrow(CYT_S %>% filter(grepl("AAAAAAAAA", sequence)))
nrow(NUC_S %>% filter(grepl("AAAAAAAAA", sequence)))


View(G3BP_S %>% filter(grepl("AAAAAAAAAAAA", sequence)) %>% filter(finalized_annotation == "3'UTR") %>% arrange(G3BP_E_S_BC))





View(G3BP_S %>% filter(grepl("TAA", sequence)) %>% filter(grepl("TAG", sequence)) %>% filter(grepl("TGA", sequence)))







nrow(NLS_M %>% filter(grouped_annotation == "3'UTR") %>% filter(grepl("TAA", sequence)) %>% filter(grepl("TAG", sequence)) %>% filter(grepl("TGA", sequence))) / nrow(NLS_M %>% filter(grouped_annotation == "3'UTR")) * 100
nrow(NES_M %>% filter(grouped_annotation == "3'UTR") %>% filter(grepl("TAA", sequence)) %>% filter(grepl("TAG", sequence)) %>% filter(grepl("TGA", sequence))) / nrow(NES_M %>% filter(grouped_annotation == "3'UTR")) * 100
nrow(G3BP_M %>% filter(grouped_annotation == "3'UTR") %>% filter(grepl("TAA", sequence)) %>% filter(grepl("TAG", sequence)) %>% filter(grepl("TGA", sequence))) / nrow(G3BP_M %>% filter(grouped_annotation == "3'UTR")) * 100

nrow(NLS_S %>% filter(grouped_annotation == "3'UTR") %>% filter(grepl("TAA", sequence)) %>% filter(grepl("TAG", sequence)) %>% filter(grepl("TGA", sequence))) / nrow(NLS_S %>% filter(grouped_annotation == "3'UTR")) * 100
nrow(NES_S %>% filter(grouped_annotation == "3'UTR") %>% filter(grepl("TAA", sequence)) %>% filter(grepl("TAG", sequence)) %>% filter(grepl("TGA", sequence))) / nrow(NES_S %>% filter(grouped_annotation == "3'UTR")) * 100
nrow(G3BP_S %>% filter(grouped_annotation == "3'UTR") %>% filter(grepl("TAA", sequence)) %>% filter(grepl("TAG", sequence)) %>% filter(grepl("TGA", sequence))) / nrow(G3BP_S %>% filter(grouped_annotation == "3'UTR")) * 100

