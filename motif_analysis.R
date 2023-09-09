## CoCLIP Analysis: 
## Motif Analysis
## Written by Soon Yi
## Use Homer output for subsequent analysis.
## Last Edit: 2023-09-08

library(stringr)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(corrplot)
library(pheatmap)

baseDir = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/motifEnrichment/motifs/'
setwd(baseDir)
motifFiles = list.files(baseDir)
motifFiles = motifFiles[grep('MotifCounts.txt', motifFiles)]

raw_all_counts = list()
index = 1

for (motifFile in motifFiles) {
  motifCounts = read.delim(motifFile, header = T)
  sampleName = str_split(motifFile, '\\.')[[1]][1]
  
  counts = data.frame(motifCounts %>% group_by(Motif.Name) %>% summarise(Unique_Positions = n_distinct(PositionID)))
  rownames(counts) = counts$Motif.Name
  counts$Motif.Name = NULL
  colnames(counts)[1] = sampleName

  raw_all_counts[[index]] = counts
  index = index + 1
}

## Remove UUUUU and 35-WWWUAUUUAUUUW for analysis:
all_counts = data.frame(raw_all_counts)

idx2remove = which(rownames(all_counts) %in% c('0-UUUUU', '35-WWWUAUUUAUUUW'))
all_counts = all_counts[-idx2remove, ]

all_ranks = all_counts %>% mutate_all(~rank(-., ties.method = "min"))
# all_ranks = all_ranks[-idx2remove, ]

all_counts$PAR_CLIP = as.integer(sub("-.*", "", rownames(all_counts)))
all_counts$MOTIF = sub(".*-", "", rownames(all_counts))

all_counts = arrange(all_counts, PAR_CLIP)
all_counts = all_counts[, c('MOTIF', 'PAR_CLIP', 'All_Libraries', 
                          'CoCLIP_Input_Mock', 'CoCLIP_Input_Arsenite',
                          'CoCLIP_NLS_Mock', 'CoCLIP_NLS_Arsenite',
                          'CoCLIP_NES_Mock', 'CoCLIP_NES_Arsenite',
                          'CoCLIP_G3BP_Mock', 'CoCLIP_G3BP_Arsenite',
                          'FracCLIP_Nuclear_Mock', 'FracCLIP_Nuclear_Arsenite',
                          'FracCLIP_Cytoplasm_Mock', 'FracCLIP_Cytoplasm_Arsenite')]
rownames(all_counts) = all_counts$MOTIF

all_ranks$PAR_CLIP = as.integer(sub("-.*", "", rownames(all_ranks)))
all_ranks$MOTIF = sub(".*-", "", rownames(all_ranks))

all_ranks = arrange(all_ranks, PAR_CLIP)
all_ranks = all_ranks[, c('MOTIF', 'PAR_CLIP', 'All_Libraries', 
                                          'CoCLIP_Input_Mock', 'CoCLIP_Input_Arsenite',
                                          'CoCLIP_NLS_Mock', 'CoCLIP_NLS_Arsenite',
                                          'CoCLIP_NES_Mock', 'CoCLIP_NES_Arsenite',
                                          'CoCLIP_G3BP_Mock', 'CoCLIP_G3BP_Arsenite',
                                          'FracCLIP_Nuclear_Mock', 'FracCLIP_Nuclear_Arsenite',
                                          'FracCLIP_Cytoplasm_Mock', 'FracCLIP_Cytoplasm_Arsenite')]
rownames(all_ranks) = all_ranks$MOTIF


# write.table(all_counts, 'HuR_PAR_CLIP_MotifCounts.txt', row.names = T, col.names = T, quote = F, sep = '\t')
# write.table(all_ranks, 'HuR_PAR_CLIP_MotifRanks.txt', row.names = T, col.names = T, quote = F, sep = '\t')

## FIGURE3 CorrMatrix_PARCLIP_Motifs
CorrMatrix = cor(all_ranks[, colnames(all_ranks)[2:15]])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = 14)
colnames(CorrMatrix) = colnames(all_ranks)[2:15]
rownames(CorrMatrix) = colnames(all_ranks)[2:15]
pheatmap(CorrMatrix, cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1", 'indianred2'))(100)))

## Correlation Matrix With CoCLIP, and FracCLIP
CorrMatrix = cor(all_ranks[, colnames(all_ranks)[3:15]])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = 13)
colnames(CorrMatrix) = colnames(all_ranks)[3:15]
rownames(CorrMatrix) = colnames(all_ranks)[3:15]
pheatmap(CorrMatrix, cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

## FIGURE3 RankCorr_PARCLIP_Motifs
pheatmap(all_ranks[, colnames(all_ranks)[2:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))





####################
all_counts = data.frame(raw_all_counts)

idx2remove = which(rownames(all_counts) %in% c('35-WWWUAUUUAUUUW'))
all_counts = all_counts[-idx2remove, ]

all_ranks = all_counts %>% mutate_all(~rank(-., ties.method = "min"))
# all_ranks = all_ranks[-idx2remove, ]

all_counts$PAR_CLIP = as.integer(sub("-.*", "", rownames(all_counts)))
all_counts$MOTIF = sub(".*-", "", rownames(all_counts))

all_counts = arrange(all_counts, All_Libraries)
all_counts = all_counts[, c('MOTIF', 'PAR_CLIP', 'All_Libraries', 
                            'CoCLIP_Input_Mock', 'CoCLIP_Input_Arsenite',
                            'CoCLIP_NLS_Mock', 'CoCLIP_NLS_Arsenite',
                            'CoCLIP_NES_Mock', 'CoCLIP_NES_Arsenite',
                            'CoCLIP_G3BP_Mock', 'CoCLIP_G3BP_Arsenite',
                            'FracCLIP_Nuclear_Mock', 'FracCLIP_Nuclear_Arsenite',
                            'FracCLIP_Cytoplasm_Mock', 'FracCLIP_Cytoplasm_Arsenite')]
rownames(all_counts) = all_counts$MOTIF

all_ranks$PAR_CLIP = as.integer(sub("-.*", "", rownames(all_ranks)))
all_ranks$MOTIF = sub(".*-", "", rownames(all_ranks))

all_ranks = arrange(all_ranks, All_Libraries)
all_ranks = all_ranks[, c('MOTIF', 'PAR_CLIP', 'All_Libraries', 
                          'CoCLIP_Input_Mock', 'CoCLIP_Input_Arsenite',
                          'CoCLIP_NLS_Mock', 'CoCLIP_NLS_Arsenite',
                          'CoCLIP_NES_Mock', 'CoCLIP_NES_Arsenite',
                          'CoCLIP_G3BP_Mock', 'CoCLIP_G3BP_Arsenite',
                          'FracCLIP_Nuclear_Mock', 'FracCLIP_Nuclear_Arsenite',
                          'FracCLIP_Cytoplasm_Mock', 'FracCLIP_Cytoplasm_Arsenite')]
rownames(all_ranks) = all_ranks$MOTIF

write.table(all_counts, 'HuR_PAR_CLIP_MotifCounts.txt', row.names = T, col.names = T, quote = F, sep = '\t')
write.table(all_ranks, 'HuR_PAR_CLIP_MotifRanks.txt', row.names = T, col.names = T, quote = F, sep = '\t')
# 
# ## Correlation Matrix With PAR_CLIP, CoCLIP, and FracCLIP
# CorrMatrix = cor(all_ranks[, colnames(all_ranks)[2:15]])
# CorrMatrix = matrix(round(CorrMatrix,2), nrow = 14)
# colnames(CorrMatrix) = colnames(all_ranks)[2:15]
# rownames(CorrMatrix) = colnames(all_ranks)[2:15]
# pheatmap(CorrMatrix, cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1", 'indianred2'))(100)))

## FIGURE3 CorrMatrix_PARCLIP_Motifs_Co_and_Frac
CorrMatrix = cor(all_ranks[, colnames(all_ranks)[3:15]])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = 13)
colnames(CorrMatrix) = colnames(all_ranks)[3:15]
rownames(CorrMatrix) = colnames(all_ranks)[3:15]
pheatmap(CorrMatrix, cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

## FIGURE3 RankCorr_PARCLIP_Motifs_OrderedByCoCLIP
pheatmap(all_ranks[, colnames(all_ranks)[3:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))

