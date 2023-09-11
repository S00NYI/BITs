## CoCLIP Analysis: 
## Motif Analysis
## Written by Soon Yi
## Use Homer output for subsequent analysis.
## Last Edit: 2023-09-11

library(stringr)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(corrplot)
library(pheatmap)

## data processing for motif counts
################################################################################
baseDir = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/motifEnrichment/motifs/counts'
# baseDir = 'L:/.shortcut-targets-by-id/13hY9t_p6eUdnvP-c2OClHbzyeWMOruD_/CoCLIP_HuR_Paper/Data/Homer_Outputs/motifs/counts'
setwd(baseDir)
countFiles = list.files(baseDir)
countFiles = countFiles[grep('MotifCounts.txt', countFiles)]

peakCounts = data.frame(sample = sapply(str_split(countFiles, '\\.'), function(x) x[1]), 
                        peaks = c(287466, 190803, 10924, 9666, 81323, 48659, 3250, 694, 4725, 4699, 102500, 65923, 80733, 86254))

# All_libraries = 387880

unique_peaks_per_motif = list()
unique_peaks_per_motif_normalized = list()

raw_counts_per_motif = list()
raw_counts_per_motif_normalized = list()

index = 0

for (countFile in countFiles) {
  motifCounts = read.delim(countFile, header = T)
  sampleName = str_split(countFile, '\\.')[[1]][1]
  
  peaks_per_motif = data.frame(motifCounts %>% group_by(Motif.Name) %>% summarise(Unique_Positions = n_distinct(PositionID)))
  rownames(peaks_per_motif) = peaks_per_motif$Motif.Name
  peaks_per_motif$Motif.Name = NULL
  colnames(peaks_per_motif)[1] = sampleName
  
  index = index + 1
  unique_peaks_per_motif[[index]] = peaks_per_motif
  unique_peaks_per_motif_normalized[[index]] = peaks_per_motif/peakCounts$peaks[peakCounts$sample == sampleName]
  
  motif_appearance = data.frame(motifCounts %>% group_by(Motif.Name) %>% summarise(Count = n()))
  rownames(motif_appearance) = motif_appearance$Motif.Name
  motif_appearance$Motif.Name = NULL
  colnames(motif_appearance)[1] = sampleName
  
  raw_counts_per_motif[[index]] = motif_appearance
  raw_counts_per_motif_normalized[[index]] = motif_appearance/peakCounts$peaks[peakCounts$sample == sampleName]
}

all_unique_PPM = data.frame(unique_peaks_per_motif)
all_unique_PPM_normed = data.frame(unique_peaks_per_motif_normalized)
all_ranks = all_unique_PPM %>% mutate_all(~rank(-., ties.method = "min"))

all_raw_Motifs = data.frame(raw_counts_per_motif)
all_raw_Motifs_normed = data.frame(raw_counts_per_motif_normalized)

## Swap Order
swap1and2row = function(dataframe) {
  r_temp1 = row.names(dataframe)[1]
  r_temp2 = row.names(dataframe)[2]
  r_temp = dataframe[1, ]
  dataframe[1, ] = dataframe[2, ]
  dataframe[2, ] = r_temp
  row.names(dataframe)[2] = 'temp'
  row.names(dataframe)[1] = r_temp2
  row.names(dataframe)[2] = r_temp1
  
  return(dataframe)
}

all_unique_PPM = swap1and2row(all_unique_PPM)
all_unique_PPM_normed = swap1and2row(all_unique_PPM_normed)   ## Normalized by number of peaks
all_ranks = swap1and2row(all_ranks)
all_raw_Motifs = swap1and2row(all_raw_Motifs)
all_raw_Motifs_normed = swap1and2row(all_raw_Motifs_normed)   ## Normalized by number of peaks

new_column_names = c('MOTIF', 'PAR_CLIP', 'All_Mocks', 'All_Arsenites',
                     'CoCLIP_Input_Mock', 'CoCLIP_Input_Arsenite',
                     'CoCLIP_NLS_Mock', 'CoCLIP_NLS_Arsenite',
                     'CoCLIP_NES_Mock', 'CoCLIP_NES_Arsenite',
                     'CoCLIP_G3BP_Mock', 'CoCLIP_G3BP_Arsenite',
                     'FracCLIP_Nuclear_Mock', 'FracCLIP_Nuclear_Arsenite',
                     'FracCLIP_Cytoplasm_Mock', 'FracCLIP_Cytoplasm_Arsenite')

all_unique_PPM$PAR_CLIP = as.integer(sub("-.*", "", rownames(all_unique_PPM)))
all_unique_PPM$MOTIF = sub(".*-", "", rownames(all_unique_PPM))
all_unique_PPM = arrange(all_unique_PPM, PAR_CLIP)
all_unique_PPM = all_unique_PPM[, new_column_names]
rownames(all_unique_PPM) = all_unique_PPM$MOTIF

all_unique_PPM_normed$PAR_CLIP = as.integer(sub("-.*", "", rownames(all_unique_PPM_normed)))
all_unique_PPM_normed$MOTIF = sub(".*-", "", rownames(all_unique_PPM_normed))
all_unique_PPM_normed = arrange(all_unique_PPM_normed, PAR_CLIP)
all_unique_PPM_normed = all_unique_PPM_normed[, new_column_names]
rownames(all_unique_PPM_normed) = all_unique_PPM_normed$MOTIF

all_raw_Motifs$PAR_CLIP = as.integer(sub("-.*", "", rownames(all_raw_Motifs)))
all_raw_Motifs$MOTIF = sub(".*-", "", rownames(all_raw_Motifs))
all_raw_Motifs = arrange(all_raw_Motifs, PAR_CLIP)
all_raw_Motifs = all_raw_Motifs[, new_column_names]
rownames(all_raw_Motifs) = all_raw_Motifs$MOTIF

all_raw_Motifs_normed$PAR_CLIP = as.integer(sub("-.*", "", rownames(all_raw_Motifs_normed)))
all_raw_Motifs_normed$MOTIF = sub(".*-", "", rownames(all_raw_Motifs_normed))
all_raw_Motifs_normed = arrange(all_raw_Motifs_normed, PAR_CLIP)
all_raw_Motifs_normed = all_raw_Motifs_normed[, new_column_names]
rownames(all_raw_Motifs_normed) = all_raw_Motifs_normed$MOTIF

all_ranks$PAR_CLIP = as.integer(sub("-.*", "", rownames(all_ranks)))
all_ranks$MOTIF = sub(".*-", "", rownames(all_ranks))
all_ranks = arrange(all_ranks, PAR_CLIP)
all_ranks = all_ranks[, new_column_names]
rownames(all_ranks) = all_ranks$MOTIF

# write.table(all_unique_PPM, 'HuR_PAR_CLIP_MotifCounts.txt', row.names = T, col.names = T, quote = F, sep = '\t')
# write.table(all_ranks, 'HuR_PAR_CLIP_MotifRanks.txt', row.names = T, col.names = T, quote = F, sep = '\t')
################################################################################

## FIGURE3 Correlation Matrix with PAR-CLIP data and just the PAR-CLIP motifs
################################################################################
## With Ranks
row_selection = (all_ranks$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
col_selection = c('PAR_CLIP', colnames(all_ranks)[3:16])
col_selection = c('PAR_CLIP', 'All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock') # Mock Only
# col_selection = c('PAR_CLIP', 'All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

CorrMatrix = cor(all_ranks[row_selection, col_selection])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
colnames(CorrMatrix) = col_selection
rownames(CorrMatrix) = col_selection
pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1", "indianred2"))(100)))
################################################################################

## FIGURE3 Correlation Matrix without PAR-CLIP data and just the PAR-CLIP motifs
################################################################################
## With Ranks
row_selection = (all_ranks$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
col_selection = colnames(all_ranks)[3:16]
col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

CorrMatrix = cor(all_ranks[row_selection, col_selection])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
colnames(CorrMatrix) = col_selection
rownames(CorrMatrix) = col_selection
pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

# ## With unique peaks per motif counts
# row_selection = (all_unique_PPM$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
# col_selection = colnames(all_unique_PPM)[4:15]
# col_selection = c('CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
# col_selection = c('CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only
#
# CorrMatrix = cor(all_unique_PPM[row_selection, col_selection])
# CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
# colnames(CorrMatrix) = col_selection
# rownames(CorrMatrix) = col_selection
# pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

## With unique peaks per motif normalized counts
## This looks equal to un-normalized version
row_selection = (all_unique_PPM_normed$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
col_selection = colnames(all_unique_PPM_normed)[3:16]
col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

CorrMatrix = cor(all_unique_PPM_normed[row_selection, col_selection])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
colnames(CorrMatrix) = col_selection
rownames(CorrMatrix) = col_selection
pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

# ## With raw motif counts
# row_selection = (all_raw_Motifs$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
# col_selection = colnames(all_raw_Motifs)[4:15]
# col_selection = c('CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
# col_selection = c('CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only
#
# CorrMatrix = cor(all_raw_Motifs[row_selection, col_selection])
# CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
# colnames(CorrMatrix) = col_selection
# rownames(CorrMatrix) = col_selection
# pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

## With raw motif normalized counts
## This looks equal to un-normalized version
row_selection = (all_raw_Motifs_normed$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
col_selection = colnames(all_raw_Motifs_normed)[3:16]
col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

CorrMatrix = cor(all_raw_Motifs_normed[row_selection, col_selection])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
colnames(CorrMatrix) = col_selection
rownames(CorrMatrix) = col_selection
pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
################################################################################

## FIGURE3 Correlation Matrix without PAR-CLIP data and with UUUUU and AAAAA
################################################################################
## With Ranks
col_selection = colnames(all_ranks)[3:16]
col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

CorrMatrix = cor(all_ranks[, col_selection])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
colnames(CorrMatrix) = col_selection
rownames(CorrMatrix) = col_selection
pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

# ## With unique peaks per motif counts
# row_selection = (all_unique_PPM$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
# col_selection = colnames(all_unique_PPM)[4:15]
# col_selection = c('CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
# col_selection = c('CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only
# CorrMatrix = cor(all_unique_PPM[row_selection, col_selection])
# CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
# colnames(CorrMatrix) = col_selection
# rownames(CorrMatrix) = col_selection
# pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

## With unique peaks per motif normalized counts
## This looks equal to un-normalized version
col_selection = colnames(all_unique_PPM_normed)[3:16]
col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

CorrMatrix = cor(all_unique_PPM_normed[, col_selection])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
colnames(CorrMatrix) = col_selection
rownames(CorrMatrix) = col_selection
pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

# ## With raw motif counts
# row_selection = (all_raw_Motifs$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
# col_selection = colnames(all_raw_Motifs)[4:15]
# col_selection = c('CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
# col_selection = c('CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only
# CorrMatrix = cor(all_raw_Motifs[row_selection, col_selection])
# CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
# colnames(CorrMatrix) = col_selection
# rownames(CorrMatrix) = col_selection
# pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

## With raw motif normalized counts
## This looks equal to un-normalized version
col_selection = colnames(all_raw_Motifs_normed)[4:16]
col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

CorrMatrix = cor(all_raw_Motifs_normed[, col_selection])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
colnames(CorrMatrix) = col_selection
rownames(CorrMatrix) = col_selection
pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
################################################################################

## FIGURE3 Motif Score Heatmap with PAR-CLIP data and just the PAR-CLIP motifs
################################################################################
## with Ranks
row_selection = (all_ranks$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
col_selection = c('PAR_CLIP', colnames(all_ranks)[3:16])
col_selection = c('PAR_CLIP', 'All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock') # Mock Only
col_selection = c('PAR_CLIP', 'All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

pheatmap(all_ranks[row_selection, col_selection], cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))
################################################################################

## FIGURE3 Motif Score Heatmap without PAR-CLIP data and just the PAR-CLIP motifs
################################################################################
## With Ranks
row_selection = (all_ranks$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
col_selection = colnames(all_ranks)[3:16]
col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

pheatmap(all_ranks[row_selection, col_selection], cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))

# ## With unique peaks per motif counts
# ## Pointless to look at un-normalized ones, because scale varies too much
# row_selection = (all_unique_PPM$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
# col_selection = colnames(all_unique_PPM)[4:15]
# col_selection = c('CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
# col_selection = c('CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only
# 
# pheatmap(all_unique_PPM[row_selection, col_selection], cluster_rows=T, cluster_cols=T, color = (colorRampPalette(c("lemonchiffon1", "green4"))(100)))

## With unique peaks per motif normalized counts
row_selection = (all_unique_PPM_normed$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
col_selection = colnames(all_unique_PPM_normed)[3:15]
col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

pheatmap(all_unique_PPM_normed[row_selection, col_selection], cluster_rows=T, cluster_cols=T, color = (colorRampPalette(c("lemonchiffon1", "green4"))(100)))

# ## With raw motif counts
# ## Pointless to look at un-normalized ones, because scale varies too much
# row_selection = (all_raw_Motifs$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
# col_selection = colnames(all_raw_Motifs)[4:15]
# col_selection = c('CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
# col_selection = c('CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only
# 
# pheatmap(all_raw_Motifs[row_selection, col_selection], cluster_rows=T, cluster_cols=T, color = (colorRampPalette(c("lemonchiffon1", "green4"))(100)))

## With raw motif normalized counts
row_selection = (all_raw_Motifs_normed$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
col_selection = colnames(all_raw_Motifs_normed)[3:15]
col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

pheatmap(all_raw_Motifs_normed[row_selection, col_selection], cluster_rows=T, cluster_cols=T, color = (colorRampPalette(c("lemonchiffon1", "green4"))(100)))
################################################################################

## FIGURE3 Motif Score Heatmap without PAR-CLIP data and with UUUUU and AAAAA
################################################################################
## With Ranks
col_selection = colnames(all_ranks)[3:15]
col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

pheatmap(all_ranks[, col_selection], cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))

# ## With unique peaks per motif counts
# ## Pointless to look at un-normalized ones, because scale varies too much
# col_selection = colnames(all_unique_PPM)[4:15]
# col_selection = c('CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
# col_selection = c('CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only
# 
# pheatmap(all_unique_PPM[, col_selection], cluster_rows=T, cluster_cols=T, color = (colorRampPalette(c("lemonchiffon1", "green4"))(100)))

## With unique peaks per motif normalized counts
col_selection = colnames(all_unique_PPM_normed)[3:15]
col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

pheatmap(all_unique_PPM_normed[, col_selection], cluster_rows=T, cluster_cols=T, color = (colorRampPalette(c("lemonchiffon1", "green4"))(100)))

# ## With raw motif counts
# ## Pointless to look at un-normalized ones, because scale varies too much
# col_selection = colnames(all_raw_Motifs)[4:15]
# col_selection = c('CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
# col_selection = c('CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only
# 
# pheatmap(all_raw_Motifs[, col_selection], cluster_rows=T, cluster_cols=T, color = (colorRampPalette(c("lemonchiffon1", "green4"))(100)))

## With raw motif normalized counts
col_selection = colnames(all_raw_Motifs_normed)[3:15]
col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

pheatmap(all_raw_Motifs_normed[, col_selection], cluster_rows=T, cluster_cols=T, color = (colorRampPalette(c("lemonchiffon1", "green4"))(100)))
################################################################################

## data processing for motif density
################################################################################
# baseDir = 'L:/.shortcut-targets-by-id/13hY9t_p6eUdnvP-c2OClHbzyeWMOruD_/CoCLIP_HuR_Paper/Data/Homer_Outputs/motifs_density/50bp/'
baseDir = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/motifEnrichment/motifs/density'
setwd(baseDir)

densityFiles = list.files(paste0(baseDir))
densityFiles_CoCLIP = densityFiles[grep("^CoCLIP", densityFiles)]
densityFiles_FracCLIP  = densityFiles[grep("FracCLIP", densityFiles)]

return_Density = function(FILE, strand = NULL, normalize = NULL) {
  densityFile = read.delim(FILE)
  Localization = str_split(str_split(FILE , '\\.')[[1]][1], '_')[[1]][2]
  Condition = str_split(str_split(FILE , '\\.')[[1]][1], '_')[[1]][3]
  
  if (is.null(strand)) {
    densityFile = cbind(densityFile[, 1], densityFile[, c(2, seq(5, ncol(densityFile)-4, by = 3))])
  } else if (strand == '+') {
    densityFile = cbind(densityFile[, 1], densityFile[, c(3, seq(6, ncol(densityFile)-4, by = 3))])
  } else if (strand == '-') {
    densityFile = cbind(densityFile[, 1], densityFile[, c(4, seq(7, ncol(densityFile)-4, by = 3))])
  }
  
  colnames(densityFile) = c('position', sapply(str_split(colnames(densityFile)[2:ncol(densityFile)], '\\.'), function(x) x[2]))
  
  if (!is.null(normalize)) {
    densityFile = densityFile %>% mutate(across(-1, ~ . - min(.)))
  }
  
  return(densityFile)
}

plot_Density = function(density_data, motif_list, xaxis_lims = NULL, yaxis_lims = NULL, custom_colors = NULL, sampleName = NULL) {
  # Select the columns based on the motif_list
  plot_data = density_data %>% select(position, {{motif_list}})
  
  # Create a long-format dataframe for better legend handling
  plot_data_long = plot_data %>%
    pivot_longer(cols = {{motif_list}}, names_to = "Motif", values_to = "Density")
  
  # Define the order of levels for the Motif factor
  plot_data_long$Motif = factor(plot_data_long$Motif, levels = motif_list)
  
  # Create the ggplot object
  plot = ggplot(plot_data_long, aes(x = position, y = Density, color = Motif)) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    # ylim(yaxis_lims) +
    # xlim(xaxis_lims) +
    # labs(title = "Metagene Plot: Motif Density around Peaks",
    #      x = "Distance to peak center (nucleotides)",
    #      y = "Peak Density") +
    theme_minimal() +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = 'bold'))
  
  # Add smoothed lines
  plot = plot + geom_smooth(span = 0.2, se = FALSE)
  
  if (!is.null(custom_colors)) {
    plot = plot + scale_color_manual(values = custom_colors)
  }
  
  if (!is.null(xaxis_lims)) {
    plot = plot + xlim(xaxis_lims)
  }
  
  if (!is.null(yaxis_lims)) {
    plot = plot + ylim(yaxis_lims)
  }
  
  if (!is.null(sampleName)) {
    plot = plot + labs(title = paste0(sampleName, ": Motif Density around Peaks"),
                       x = "Distance to peak center (nucleotides)",
                       y = "Peak Density")
  } else {
    plot = plot + labs("Metagene Plot: Motif Density around Peaks",
                       x = "Distance to peak center (nucleotides)",
                       y = "Peak Density")
  }
  
  return(plot)
}
################################################################################

## FIGURE3 B Metagene Plot from All Libraries for Top Motifs UUUUU AAAAA 
################################################################################
## Mock
densityFile = return_Density(densityFiles[2], strand = '+')
plot_Density(densityFile, c('UUUUU', 'AAAAA'))

## Arsenite
densityFile = return_Density(densityFiles[1], strand = '+')
plot_Density(densityFile, c('UUUUU', 'AAAAA'))
################################################################################

## FIGURE3 C Metagene Plot from Input vs CoCLIP for Top Motifs UUUUU AAAAA 
################################################################################
## Mock
Mock_All = return_Density(densityFiles[2], strand = '+') 
Mock_Input = return_Density(densityFiles[6], strand = '+')
Mock_NES = return_Density(densityFiles[8], strand = '+')
Mock_NLS = return_Density(densityFiles[10], strand = '+')
Mock_G3BP = return_Density(densityFiles[4], strand = '+')
Mock_Nuclear = return_Density(densityFiles[14], strand = '+')
Mock_Cytoplasm = return_Density(densityFiles[12], strand = '+')

Mock_Combined = data.frame(cbind(Mock_All$position, 
                                 Mock_All$UUUUU, Mock_Input$UUUUU, Mock_NES$UUUUU, Mock_NLS$UUUUU, Mock_G3BP$UUUUU, Mock_Nuclear$UUUUU, Mock_Cytoplasm$UUUUU,
                                 Mock_All$AAAAA, Mock_Input$AAAAA, Mock_NES$AAAAA, Mock_NLS$AAAAA, Mock_G3BP$AAAAA, Mock_Nuclear$AAAAA, Mock_Cytoplasm$AAAAA))

colnames(Mock_Combined) = c('position', 
                            'UUUUU_All', 'UUUUU_Input', 'UUUUU_NES', 'UUUUU_NLS', 'UUUUU_G3BP', 'UUUUU_Nuclear', 'UUUUU_Cytoplasm',
                            'AAAAA_All', 'AAAAA_Input', 'AAAAA_NES', 'AAAAA_NLS', 'AAAAA_G3BP', 'AAAAA_Nuclear', 'AAAAA_Cytoplasm')

plot_Density(Mock_Combined, c('UUUUU_Input', 'UUUUU_NES', 'UUUUU_NLS', 'UUUUU_G3BP'), yaxis_lims = c(0, 0.15))
plot_Density(Mock_Combined, c('AAAAA_Input', 'AAAAA_NES', 'AAAAA_NLS', 'AAAAA_G3BP'), yaxis_lims = c(0, 0.15))

## Arsenite
Arsenite_All = return_Density(densityFiles[1], strand = '+')
Arsenite_Input = return_Density(densityFiles[5], strand = '+')
Arsenite_NES = return_Density(densityFiles[7], strand = '+')
Arsenite_NLS = return_Density(densityFiles[9], strand = '+')
Arsenite_G3BP = return_Density(densityFiles[3], strand = '+')
Arsenite_Nuclear = return_Density(densityFiles[13], strand = '+')
Arsenite_Cytoplasm = return_Density(densityFiles[11], strand = '+')

Arsenite_Combined = data.frame(cbind(Arsenite_All$position, 
                                     Arsenite_All$UUUUU, Arsenite_Input$UUUUU, Arsenite_NES$UUUUU, Arsenite_NLS$UUUUU, Arsenite_G3BP$UUUUU, Arsenite_Nuclear$UUUUU, Arsenite_Cytoplasm$UUUUU,
                                     Arsenite_All$AAAAA, Arsenite_Input$AAAAA, Arsenite_NES$AAAAA, Arsenite_NLS$AAAAA, Arsenite_G3BP$AAAAA, Arsenite_Nuclear$AAAAA, Arsenite_Cytoplasm$AAAAA))

colnames(Arsenite_Combined) = c('position', 
                                'UUUUU_All', 'UUUUU_Input', 'UUUUU_NES', 'UUUUU_NLS', 'UUUUU_G3BP', 'UUUUU_Nuclear', 'UUUUU_Cytoplasm',
                                'AAAAA_All', 'AAAAA_Input', 'AAAAA_NES', 'AAAAA_NLS', 'AAAAA_G3BP', 'AAAAA_Nuclear', 'AAAAA_Cytoplasm')

plot_Density(Arsenite_Combined, c('UUUUU_Input', 'UUUUU_NES', 'UUUUU_NLS', 'UUUUU_G3BP'), yaxis_lims = c(0, 0.15))
plot_Density(Arsenite_Combined, c('AAAAA_Input', 'AAAAA_NES', 'AAAAA_NLS', 'AAAAA_G3BP'), yaxis_lims = c(0, 0.15))
################################################################################

## FIGURE3 D
## Metagene Plot from Fractionation CLIP for Top Motifs UUUUU AAAAA
################################################################################
## Mock

plot_Density(Mock_Combined, c('UUUUU_Input', 'UUUUU_Nuclear', 'UUUUU_Cytoplasm'), yaxis_lims = c(0, 0.15))
plot_Density(Mock_Combined, c('AAAAA_Input', 'AAAAA_Nuclear', 'AAAAA_Cytoplasm'), yaxis_lims = c(0, 0.15))

plot_Density(Arsenite_Combined, c('UUUUU_Input', 'UUUUU_Nuclear', 'UUUUU_Cytoplasm'), yaxis_lims = c(0, 0.15))
plot_Density(Arsenite_Combined, c('AAAAA_Input', 'AAAAA_Nuclear', 'AAAAA_Cytoplasm'), yaxis_lims = c(0, 0.15))
################################################################################



