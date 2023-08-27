## CoCLIP Analysis: 
## Peak Enrichment Calculation
## Written by Soon Yi
## Last Edit: 2023-08-25

library(stringr)
library(readr)
library(dplyr)

## Make initial peak enrichment table:
# peaksMatrix_PATH = 'L:/My Drive/CWRU/PhD/Luna Lab/1. coCLIP/Analysis/peaks/'  ## Use this for windows machine
peaksMatrix_PATH = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/'
peaksMatrix_FILE = 'Combined_peakCoverage_groomed_normalized_annotated.txt'

peakMatrix = read_delim(paste0(peaksMatrix_PATH, peaksMatrix_FILE), show_col_types = FALSE)
peakMatrix = peakMatrix %>% mutate_at('TOTAL_BC', as.numeric)

# peakCount_median = median(peakMatrix[, colnames(peakMatrix)[6:63]][peakMatrix[, colnames(peakMatrix)[6:63]] != 0], na.rm = TRUE)
# peakMatrix = peakMatrix %>% filter(if_any(all_of(comparison_vector), ~ . > peakCount_median ))

pseudoCount = min(peakMatrix[, colnames(peakMatrix)[6:63]][peakMatrix[, colnames(peakMatrix)[6:63]] != 0], na.rm = TRUE)
peakMatrix[, colnames(peakMatrix)[6:63]] = peakMatrix[, colnames(peakMatrix)[6:63]] + pseudoCount

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
write.table(peakEnrichment, paste0(peaksMatrix_PATH, str_replace(peaksMatrix_FILE, ".txt", "_enrichment.txt")), 
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t', na = "")



## snoRNA Comparison ##
snoRNA_Enrichment = peakEnrichment %>% filter(finalized_annotation == 'snoRNA')
snoRNA_Enrichment = snoRNA_Enrichment %>% filter(Nuc_F_M_BC >= 1, Cyto_F_M_BC >= 1, NLS_E_M_BC >= 1, NES_E_M_BC >= 1)

xlims = c(-5, 10)
ylims = xlims

snoRNA_F_M_NvC = log(snoRNA_Enrichment$Nuc_F_M / snoRNA_Enrichment$Cyto_F_M, 2)
snoRNA_E_M_NvC = log(snoRNA_Enrichment$NLS_E_M / snoRNA_Enrichment$NES_E_M, 2)
# snoRNA_F_M_NvC = (snoRNA_Enrichment$Nuc_F_M) / (snoRNA_Enrichment$Cyto_F_M)
# snoRNA_E_M_NvC = (snoRNA_Enrichment$NLS_E_M) / (snoRNA_Enrichment$NES_E_M)
plot(snoRNA_E_M_NvC, snoRNA_F_M_NvC, xlim = xlims, ylim = ylims, xlab = "LOG2(CoCLIP NLS/NES Mock)", ylab = "LOG2(FractionCLIP Nuc/Cyto Mock)", pch = 16, col = "blue",)
cor(snoRNA_E_M_NvC, snoRNA_F_M_NvC)

snoRNA_F_S_NvC = log(snoRNA_Enrichment$Nuc_F_S / snoRNA_Enrichment$Cyto_F_S, 2)
snoRNA_E_S_NvC = log(snoRNA_Enrichment$NLS_E_S / snoRNA_Enrichment$NES_E_S, 2)
# snoRNA_F_S_NvC = (snoRNA_Enrichment$Nuc_F_S) / (snoRNA_Enrichment$Cyto_F_S)
# snoRNA_E_S_NvC = (snoRNA_Enrichment$NLS_E_S) / (snoRNA_Enrichment$NES_E_S)
points(snoRNA_E_S_NvC, snoRNA_F_S_NvC, xlim = xlims, ylim = ylims, xlab = "LOG2(CoCLIP NLS/NES Stress)", ylab = "LOG2(FractionCLIP Nuc/Cyto Stress)", pch = 16, col = "red",)
cor(snoRNA_E_S_NvC, snoRNA_F_S_NvC)

## Nuclear Comparison 
NuclearPeaks_Filtered= peakEnrichment %>% filter(NLS_I_M_BC >= 1, NLS_I_S_BC >= 1, NLS_E_M_BC >= 2, NLS_E_S_BC >= 2)

# Input vs CoCLIP
Nu_M_IvE = NuclearPeaks_Filtered$NLS_E_M / NuclearPeaks_Filtered$NLS_I_M
Nu_S_IvE = NuclearPeaks_Filtered$NLS_E_S / NuclearPeaks_Filtered$NLS_I_S
plot(log(Nu_M_IvE, 2), log(Nu_S_IvE, 2), col = 'black', pch = 16,
     xlab = 'log2(mock Enriched/Input)',
     ylab = 'log2(stress Enriched/Input)',
     xlim = c(-4, 4), ylim = c(-4, 4))

# Mock vs Stress
Nu_I_MvS = NuclearPeaks_Filtered$NLS_I_M / NuclearPeaks_Filtered$NLS_I_S
Nu_E_MvS = NuclearPeaks_Filtered$NLS_E_M / NuclearPeaks_Filtered$NLS_E_S
plot(log(Nu_E_MvS, 2), log(Nu_I_MvS, 2), col = 'black', pch = 16,
     xlab = 'log2(CoCLIP Mock/Stress)',
     ylab = 'log2(Input Mock/Stress)',
     xlim = c(-4, 4), ylim = c(-4, 4))

## Cytoplasm Comparison 
CytopPeaks_Filtered= peakEnrichment %>% filter(NES_I_M_BC >= 1, NES_I_S_BC >= 1, NES_E_M_BC >= 2, NES_E_S_BC >= 2)

# Input vs CoCLIP
Cy_M_IvE = CytopPeaks_Filtered$NES_E_M / CytopPeaks_Filtered$NES_I_M
Cy_S_IvE = CytopPeaks_Filtered$NES_E_S / CytopPeaks_Filtered$NES_I_S
plot(log(Cy_M_IvE, 2), log(Cy_S_IvE, 2), col = 'black', pch = 16,
     xlab = 'log2(Mock Enriched/Input)',
     ylab = 'log2(Stress Enriched/Input)',
     xlim = c(-10, 10), ylim = c(-10, 10))
title('Cytoplasm HuR Peaks: Mock vs Stress of Enriched/Input')

# Mock vs Stress
Cy_I_MvS = CytopPeaks_Filtered$NES_I_M / CytopPeaks_Filtered$NES_I_S
Cy_E_MvS = CytopPeaks_Filtered$NES_E_M / CytopPeaks_Filtered$NES_E_S
plot(log(Cy_E_MvS, 2), log(Cy_I_MvS), col = 'black', pch = 16,
     xlab = 'log2(CoCLIP Mock/Stress)',
     ylab = 'log2(Input Mock/Stress)',
     xlim = c(-4, 4), ylim = c(-4, 4))
title('Cytoplasm HuR Peaks: CoCLIP vs Input of Mock/Stress')


Cy_M_IvE = peakEnrichment$NES_E_M / peakEnrichment$NES_I_M
Cy_S_IvE = peakEnrichment$NES_E_S / peakEnrichment$NES_I_S

SG_M_IvE = peakEnrichment$G3BP_E_M / peakEnrichment$G3BP_I_M
SG_S_IvE = peakEnrichment$G3BP_E_S / peakEnrichment$G3BP_I_S


plot(log(Cy_M_IvE, 2), log(Cy_S_IvE, 2), col = 'black', pch = 16)
plot(log(SG_M_IvE, 2), log(SG_S_IvE, 2), col = 'black', pch = 16)


## Mock vs Stress Comparison ## 


Cy_I_MvS = peakEnrichment$NES_I_M / peakEnrichment$NES_I_S
Cy_E_MvS = peakEnrichment$NES_E_M / peakEnrichment$NES_E_S

SG_I_MvS = peakEnrichment$G3BP_I_M / peakEnrichment$G3BP_I_S
SG_E_MvS = peakEnrichment$G3BP_E_M / peakEnrichment$G3BP_E_S

plot(log(Nu_E_MvS, 2), log(Nu_I_MvS, 2), col = 'black', pch = 16)
plot(log(Cy_E_MvS, 2), log(Cy_I_MvS, 2), col = 'black', pch = 16)
plot(log(SG_E_MvS, 2), log(SG_I_MvS, 2), col = 'black', pch = 16)



