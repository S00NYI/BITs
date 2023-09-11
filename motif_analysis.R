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

# baseDir = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/motifEnrichment/motifs/'
baseDir = 'L:/.shortcut-targets-by-id/13hY9t_p6eUdnvP-c2OClHbzyeWMOruD_/CoCLIP_HuR_Paper/Data/Homer_Outputs/motifs/counts'
setwd(baseDir)
motifFiles = list.files(baseDir)
motifFiles = motifFiles[grep('MotifCounts.txt', motifFiles)]

peakCounts = data.frame(sample = sapply(str_split(motifFiles, '\\.'), function(x) x[1]), 
                        peaks = c(10924, 9666, 81323, 48659, 3250, 694, 4725, 4699, 102500, 65923, 80733, 86254))
# , peaks = c(387880, 10924, 9666, 81323, 48659, 3250, 694, 4725, 4699, 102500, 65923, 80733, 86254))

unique_peaks_per_motif = list()
unique_peaks_per_motif_normalized = list()

raw_counts_per_motif = list()
raw_counts_per_motif_normalized = list()

index = 1

for (motifFile in motifFiles) {
  motifCounts = read.delim(motifFile, header = T)
  sampleName = str_split(motifFile, '\\.')[[1]][1]
  
  peaks_per_motif = data.frame(motifCounts %>% group_by(Motif.Name) %>% summarise(Unique_Positions = n_distinct(PositionID)))
  rownames(peaks_per_motif) = peaks_per_motif$Motif.Name
  peaks_per_motif$Motif.Name = NULL
  colnames(peaks_per_motif)[1] = sampleName
  
  unique_peaks_per_motif[[index]] = peaks_per_motif
  unique_peaks_per_motif_normalized[[index]] = peaks_per_motif/peakCounts$peaks[peakCounts$sample == sampleName]
  
  motif_appearance = data.frame(motifCounts %>% group_by(Motif.Name) %>% summarise(Count = n()))
  rownames(motif_appearance) = motif_appearance$Motif.Name
  motif_appearance$Motif.Name = NULL
  colnames(motif_appearance)[1] = sampleName
  
  raw_counts_per_motif[[index]] = motif_appearance
  raw_counts_per_motif_normalized[[index]] = motif_appearance/peakCounts$peaks[peakCounts$sample == sampleName]
  
  index = index + 1
}


all_unique_PPM = data.frame(unique_peaks_per_motif)
all_unique_PPM_normed = data.frame(unique_peaks_per_motif_normalized)
all_raw_Motifs = data.frame(raw_counts_per_motif)
all_raw_Motifs_normed = data.frame(raw_counts_per_motif_normalized)

all_ranks = all_unique_PPM %>% mutate_all(~rank(-., ties.method = "min"))

new_column_names = c('MOTIF', 'PAR_CLIP', 
                     'CoCLIP_Input_Mock', 'CoCLIP_Input_Arsenite',
                     'CoCLIP_NLS_Mock', 'CoCLIP_NLS_Arsenite',
                     'CoCLIP_NES_Mock', 'CoCLIP_NES_Arsenite',
                     'CoCLIP_G3BP_Mock', 'CoCLIP_G3BP_Arsenite',
                     'FracCLIP_Nuclear_Mock', 'FracCLIP_Nuclear_Arsenite',
                     'FracCLIP_Cytoplasm_Mock', 'FracCLIP_Cytoplasm_Arsenite')

# new_column_names = c('MOTIF', 'PAR_CLIP', 'All_Libraries', 
#                      'CoCLIP_Input_Mock', 'CoCLIP_Input_Arsenite',
#                      'CoCLIP_NLS_Mock', 'CoCLIP_NLS_Arsenite',
#                      'CoCLIP_NES_Mock', 'CoCLIP_NES_Arsenite',
#                      'CoCLIP_G3BP_Mock', 'CoCLIP_G3BP_Arsenite',
#                      'FracCLIP_Nuclear_Mock', 'FracCLIP_Nuclear_Arsenite',
#                      'FracCLIP_Cytoplasm_Mock', 'FracCLIP_Cytoplasm_Arsenite')

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

####
#### FIGURE3 CorrMatrix_PARCLIP_Motifs
## With Ranks
CorrMatrix = cor(all_ranks[(all_ranks$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU'), colnames(all_ranks)[2:14]])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = 13)
colnames(CorrMatrix) = colnames(all_ranks)[2:14]
rownames(CorrMatrix) = colnames(all_ranks)[2:14]
pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "indianred2"))(100)))
pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1", "indianred2"))(100)))
pheatmap(CorrMatrix, cluster_rows=F, cluster_cols=F, color = rev(colorRampPalette(c("green4", "lemonchiffon1", "indianred2"))(100)))

# ## With unique peaks per motif counts
# CorrMatrix = cor(all_unique_PPM[(all_unique_PPM$MOTIF != 'AAAAA' & all_unique_PPM$MOTIF != 'UUUUU'), colnames(all_unique_PPM)[2:15]])
# CorrMatrix = matrix(round(CorrMatrix,2), nrow = 14)
# colnames(CorrMatrix) = colnames(all_unique_PPM)[2:15]
# rownames(CorrMatrix) = colnames(all_unique_PPM)[2:15]
# pheatmap(CorrMatrix, cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1", 'indianred2'))(100)))
# 
# ## With unique peaks per motif normalized counts
# CorrMatrix = cor(all_unique_PPM_normed[(all_unique_PPM_normed$MOTIF != 'AAAAA' & all_unique_PPM_normed$MOTIF != 'UUUUU'), colnames(all_unique_PPM_normed)[2:15]])
# CorrMatrix = matrix(round(CorrMatrix,2), nrow = 14)
# colnames(CorrMatrix) = colnames(all_unique_PPM_normed)[2:15]
# rownames(CorrMatrix) = colnames(all_unique_PPM_normed)[2:15]
# pheatmap(CorrMatrix, cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1", 'indianred2'))(100)))
# 
# ## With raw motif counts
# CorrMatrix = cor(all_raw_Motifs[(all_raw_Motifs$MOTIF != 'AAAAA' & all_raw_Motifs$MOTIF != 'UUUUU'), colnames(all_raw_Motifs)[2:15]])
# CorrMatrix = matrix(round(CorrMatrix,2), nrow = 14)
# colnames(CorrMatrix) = colnames(all_raw_Motifs)[2:15]
# rownames(CorrMatrix) = colnames(all_raw_Motifs)[2:15]
# pheatmap(CorrMatrix, cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1", 'indianred2'))(100)))
# 
# ## With raw motif normalized counts
# CorrMatrix = cor(all_raw_Motifs_normed[(all_raw_Motifs_normed$MOTIF != 'AAAAA' & all_raw_Motifs_normed$MOTIF != 'UUUUU'), colnames(all_raw_Motifs_normed)[2:15]])
# CorrMatrix = matrix(round(CorrMatrix,2), nrow = 14)
# colnames(CorrMatrix) = colnames(all_raw_Motifs_normed)[2:15]
# rownames(CorrMatrix) = colnames(all_raw_Motifs_normed)[2:15]
# pheatmap(CorrMatrix, cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1", 'indianred2'))(100)))

####
#### Correlation Matrix Without PAR CLIP
## With Ranks
CorrMatrix = cor(all_ranks[, colnames(all_ranks)[3:14]])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = 12)
colnames(CorrMatrix) = colnames(all_ranks)[3:14]
rownames(CorrMatrix) = colnames(all_ranks)[3:14]
pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
pheatmap(CorrMatrix, cluster_rows=F, cluster_cols=F, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

# ## With unique peaks per motif counts
# CorrMatrix = cor(all_unique_PPM[, colnames(all_unique_PPM)[3:14]])
# CorrMatrix = matrix(round(CorrMatrix,2), nrow = 12)
# colnames(CorrMatrix) = colnames(all_unique_PPM)[3:14]
# rownames(CorrMatrix) = colnames(all_unique_PPM)[3:14]
# pheatmap(CorrMatrix, cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

# ## With unique peaks per motif normalized counts
# CorrMatrix = cor(all_unique_PPM_normed[, colnames(all_unique_PPM_normed)[3:14]])
# CorrMatrix = matrix(round(CorrMatrix,2), nrow = 12)
# colnames(CorrMatrix) = colnames(all_unique_PPM_normed)[3:14]
# rownames(CorrMatrix) = colnames(all_unique_PPM_normed)[3:14]
# pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

# ## With raw motif counts
# CorrMatrix = cor(all_raw_Motifs[, colnames(all_raw_Motifs)[3:14]])
# CorrMatrix = matrix(round(CorrMatrix,2), nrow = 12)
# colnames(CorrMatrix) = colnames(all_raw_Motifs)[3:14]
# rownames(CorrMatrix) = colnames(all_raw_Motifs)[3:14]
# pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

## With raw motif normalized counts
CorrMatrix = cor(all_raw_Motifs_normed[, colnames(all_raw_Motifs_normed)[3:14]])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = 12)
colnames(CorrMatrix) = colnames(all_raw_Motifs_normed)[3:14]
rownames(CorrMatrix) = colnames(all_raw_Motifs_normed)[3:14]
pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1", "indianred2"))(100)))

## With raw motif normalized counts with JUST PAR-CLIP Motifs
CorrMatrix = cor(all_raw_Motifs_normed[(all_raw_Motifs_normed$MOTIF != 'AAAAA' & all_raw_Motifs_normed$MOTIF != 'UUUUU'), colnames(all_raw_Motifs_normed)[3:14]])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = 12)
colnames(CorrMatrix) = colnames(all_raw_Motifs_normed)[3:14]
rownames(CorrMatrix) = colnames(all_raw_Motifs_normed)[3:14]
pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1", "indianred2"))(100)))

####
#### FIGURE3 Corr_PARCLIP_Motifs
#### with PAR-CLIP
## with Ranks
pheatmap(all_ranks[(all_ranks$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU'), colnames(all_ranks)[2:14]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))
# pheatmap(all_ranks[(all_ranks$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_ranks)[2:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))

# ## With unique peaks per motif counts
# pheatmap(all_unique_PPM[(all_unique_PPM$MOTIF != 'WWWUAUUUAUUUW') & (all_unique_PPM$MOTIF != 'UUUUU'), colnames(all_unique_PPM)[2:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
# pheatmap(all_unique_PPM[(all_unique_PPM$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_unique_PPM)[2:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
# 
# ## With unique peaks per motif normalized counts
# pheatmap(all_unique_PPM_normed[(all_unique_PPM_normed$MOTIF != 'WWWUAUUUAUUUW') & (all_unique_PPM_normed$MOTIF != 'UUUUU'), colnames(all_unique_PPM_normed)[2:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
# pheatmap(all_unique_PPM_normed[(all_unique_PPM_normed$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_unique_PPM_normed)[2:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
# 
# ## With raw motif counts
# pheatmap(all_raw_Motifs[(all_raw_Motifs$MOTIF != 'WWWUAUUUAUUUW') & (all_raw_Motifs$MOTIF != 'UUUUU'), colnames(all_raw_Motifs)[2:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
# pheatmap(all_raw_Motifs[(all_raw_Motifs$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_raw_Motifs)[2:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
# 
# ## With raw motif normalized counts
# pheatmap(all_raw_Motifs_normed[(all_raw_Motifs_normed$MOTIF != 'WWWUAUUUAUUUW') & (all_raw_Motifs_normed$MOTIF != 'UUUUU'), colnames(all_raw_Motifs_normed)[2:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
# pheatmap(all_raw_Motifs_normed[(all_raw_Motifs_normed$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_raw_Motifs_normed)[2:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

#### without PAR-CLIP
## with Ranks
pheatmap(all_ranks[, colnames(all_ranks)[3:14]], cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))
pheatmap(all_ranks[, colnames(all_ranks)[3:14]], cluster_rows=F, cluster_cols=F, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))
# pheatmap(all_ranks[(all_ranks$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_ranks)[3:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))

# ## With unique peaks per motif counts
# pheatmap(all_unique_PPM[(all_unique_PPM$MOTIF != 'WWWUAUUUAUUUW') & (all_unique_PPM$MOTIF != 'UUUUU'), colnames(all_unique_PPM)[3:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
# pheatmap(all_unique_PPM[(all_unique_PPM$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_unique_PPM)[3:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

## With unique peaks per motif normalized counts
pheatmap(all_unique_PPM_normed[, colnames(all_unique_PPM_normed)[3:14]], cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
pheatmap(all_unique_PPM_normed[, colnames(all_unique_PPM_normed)[3:14]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
# pheatmap(all_unique_PPM_normed[, colnames(all_unique_PPM_normed)[3:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

# ## With raw motif counts
# pheatmap(all_raw_Motifs[(all_raw_Motifs$MOTIF != 'WWWUAUUUAUUUW') & (all_raw_Motifs$MOTIF != 'UUUUU'), colnames(all_raw_Motifs)[3:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
# pheatmap(all_raw_Motifs[(all_raw_Motifs$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_raw_Motifs)[3:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

## With raw motif normalized counts
pheatmap(all_raw_Motifs_normed[, colnames(all_raw_Motifs_normed)[3:14]], cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
pheatmap(all_raw_Motifs_normed[, colnames(all_raw_Motifs_normed)[3:14]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
# pheatmap(all_raw_Motifs_normed[(all_raw_Motifs_normed$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_raw_Motifs_normed)[3:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

####################
baseDir = 'L:/.shortcut-targets-by-id/13hY9t_p6eUdnvP-c2OClHbzyeWMOruD_/CoCLIP_HuR_Paper/Data/Homer_Outputs/motifs_density/50bp/'
# baseDir = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/motifEnrichment/sampled/motifs_density/50bp/'
setwd(baseDir)


countsFile = list.files(paste0(baseDir))

CoCLIP_Files = countsFile[grep("^CoCLIP", countsFile)]
FracCLIP_Files = countsFile[grep("FracCLIP", countsFile)]

## FIGURE3 B

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

## FIGURE3 C
## Metagene plots centered on peaks across peaks but for top motif only: Input vs CoCLIP

mock_all = list()
arsenite_all = list()

# for (CoCLIP_File in CoCLIP_Files) {
#   count = read.delim(CoCLIP_File)
#   Localization = str_split(str_split(CoCLIP_File , '\\.')[[1]][1], '_')[2]
#   Condition = str_split(str_split(CoCLIP_File , '\\.')[[1]][1], '_')[2]
#   
#   count = cbind(count[, 1], count[, c(2, seq(5, ncol(count)-4, by = 3))])
#   # count = cbind(count[, 1], count[, c(3, seq(6, ncol(count)-4, by = 3))])     # + strand
#   # count = cbind(count[, 1], count[, c(4, seq(7, ncol(count)-4, by = 3))])     # - strand
#   colnames(count) = c('position', 'WUUUA', 'UUUUU', 'YUUUA', 'AUUUY', 'UWUAA', 'WAAAA', 'UYAAA', 'UAAAW', 'WAAAU', 'WUAAA')
# }

return_counts = function(FILE) {
  count = read.delim(FILE)
  Localization = str_split(str_split(FILE , '\\.')[[1]][1], '_')[2]
  Condition = str_split(str_split(FILE , '\\.')[[1]][1], '_')[2]
  
  count = cbind(count[, 1], count[, c(2, seq(5, ncol(count)-4, by = 3))])
  # count = cbind(count[, 1], count[, c(3, seq(6, ncol(count)-4, by = 3))])     # + strand
  # count = cbind(count[, 1], count[, c(4, seq(7, ncol(count)-4, by = 3))])     # - strand
  colnames(count) = c('position', 'WUUUA', 'UUUUU', 'YUUUA', 'AUUUY', 'UWUAA', 'WAAAA', 'UYAAA', 'UAAAW', 'WAAAU', 'WUAAA')
  
  return(count)
}

return_normed_counts = function(FILE) {
  count = read.delim(FILE)
  Localization = str_split(str_split(FILE , '\\.')[[1]][1], '_')[2]
  Condition = str_split(str_split(FILE , '\\.')[[1]][1], '_')[2]
  
  count = cbind(count[, 1], count[, c(2, seq(5, ncol(count)-4, by = 3))])
  # count = cbind(count[, 1], count[, c(3, seq(6, ncol(count)-4, by = 3))])     # + strand
  # count = cbind(count[, 1], count[, c(4, seq(7, ncol(count)-4, by = 3))])     # - strand
  colnames(count) = c('position', 'WUUUA', 'UUUUU', 'YUUUA', 'AUUUY', 'UWUAA', 'WAAAA', 'UYAAA', 'UAAAW', 'WAAAU', 'WUAAA')
  normalized_count = count %>% mutate(across(-1, ~ . - min(.)))
  
  return(normalized_count)
}

plot_density = function(density_data, motif_list, xaxis_lims = NULL, yaxis_lims = NULL, custom_colors = NULL, sampleName = NULL) {
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

# motif_list = c('WUUUA', 'UUUUU', 'YUUUA', 'AUUUY', 'UWUAA', 'WAAAA', 'UYAAA', 'UAAAW', 'WAAAU', 'WUAAA')
motif_list = c('UUUUU', 'WAAAA')

CoCLIP_File = CoCLIP_Files[8]
plot_density(return_counts(CoCLIP_File), motif_list, yaxis_lims = c(0, 1), sampleName = CoCLIP_File)


## FIGURE3 D
## Metagene plots centered on peaks across peaks but for top motif only: Fractionation CLIP

for (FracCLIP_File in FracCLIP_Files) {
  count = read.delim(FracCLIP_File)
  
  count = cbind(count[, 1], count[, c(3, seq(6, ncol(count)-4, by = 3))])
  colnames(count) = c('position', 'WUUUA', 'UUUUU', 'YUUUA', 'AUUUY', 'UWUAA', 'WAAAA', 'UYAAA', 'UAAAW', 'WAAAU', 'WUAAA')
  
  normalized_count = count %>% mutate(across(-1, ~ . - min(.)))
  
}




motif_list = c('WAAAA', 'UUUUU')
plot_density(normalized_count, motif_list)




