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
baseDir = 'L:/.shortcut-targets-by-id/13hY9t_p6eUdnvP-c2OClHbzyeWMOruD_/CoCLIP_HuR_Paper/Data/Homer_Outputs/motifs'
setwd(baseDir)
motifFiles = list.files(baseDir)
motifFiles = motifFiles[grep('MotifCounts.txt', motifFiles)]

peakCounts = data.frame(sample = sapply(str_split(motifFiles, '\\.'), function(x) x[1]), 
                        peaks = c(387880, 10924, 9666, 81323, 48659, 3250, 694, 4725, 4699, 102500, 65923, 80733, 86254))

raw_all_counts = list()
normalized_all_counts = list()
index = 1

for (motifFile in motifFiles) {
  motifCounts = read.delim(motifFile, header = T)
  sampleName = str_split(motifFile, '\\.')[[1]][1]
  
  peaks = 
  counts = data.frame(motifCounts %>% group_by(Motif.Name) %>% summarise(Unique_Positions = n_distinct(PositionID)))
  rownames(counts) = counts$Motif.Name
  counts$Motif.Name = NULL
  colnames(counts)[1] = sampleName

  raw_all_counts[[index]] = counts
  normalized_all_counts[[index]] = counts/peakCounts$peaks[peakCounts$sample == sampleName]
  index = index + 1
}

## Remove UUUUU and 35-WWWUAUUUAUUUW for analysis:
all_counts = data.frame(raw_all_counts)
all_ncounts = data.frame(normalized_all_counts)

# idx2remove = which(rownames(all_counts) %in% c('0-UUUUU', '35-WWWUAUUUAUUUW'))
# all_counts = all_counts[-idx2remove, ]
# all_ncounts = all_ncounts[-idx2remove, ]

all_ranks = all_counts %>% mutate_all(~rank(-., ties.method = "min"))

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

all_ncounts$PAR_CLIP = as.integer(sub("-.*", "", rownames(all_ncounts)))
all_ncounts$MOTIF = sub(".*-", "", rownames(all_ncounts))
all_ncounts = arrange(all_ncounts, PAR_CLIP)
all_ncounts = all_ncounts[, c('MOTIF', 'PAR_CLIP', 'All_Libraries', 
                            'CoCLIP_Input_Mock', 'CoCLIP_Input_Arsenite',
                            'CoCLIP_NLS_Mock', 'CoCLIP_NLS_Arsenite',
                            'CoCLIP_NES_Mock', 'CoCLIP_NES_Arsenite',
                            'CoCLIP_G3BP_Mock', 'CoCLIP_G3BP_Arsenite',
                            'FracCLIP_Nuclear_Mock', 'FracCLIP_Nuclear_Arsenite',
                            'FracCLIP_Cytoplasm_Mock', 'FracCLIP_Cytoplasm_Arsenite')]
rownames(all_ncounts) = all_ncounts$MOTIF

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

####
#### FIGURE3 CorrMatrix_PARCLIP_Motifs
## With Ranks
CorrMatrix = cor(all_ranks[(all_ranks$MOTIF != 'WWWUAUUUAUUUW') & all_ranks$MOTIF != 'UUUUU', colnames(all_ranks)[2:15]])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = 14)
colnames(CorrMatrix) = colnames(all_ranks)[2:15]
rownames(CorrMatrix) = colnames(all_ranks)[2:15]
pheatmap(CorrMatrix, cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1", 'indianred2'))(100)))

## With Counts
CorrMatrix = cor(all_counts[(all_counts$MOTIF != 'WWWUAUUUAUUUW') & all_counts$MOTIF != 'UUUUU', colnames(all_counts)[2:15]])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = 14)
colnames(CorrMatrix) = colnames(all_counts)[2:15]
rownames(CorrMatrix) = colnames(all_counts)[2:15]
pheatmap(CorrMatrix, cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1", 'indianred2'))(100)))

## With normalized counts
CorrMatrix = cor(all_ncounts[(all_ncounts$MOTIF != 'WWWUAUUUAUUUW') & all_ncounts$MOTIF != 'UUUUU', colnames(all_ncounts)[2:15]])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = 14)
colnames(CorrMatrix) = colnames(all_ncounts)[2:15]
rownames(CorrMatrix) = colnames(all_ncounts)[2:15]
pheatmap(CorrMatrix, cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1", 'indianred2'))(100)))

####
#### Correlation Matrix Without PAR CLIP
## With Ranks
CorrMatrix = cor(all_ranks[(all_ranks$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_ranks)[3:15]])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = 13)
colnames(CorrMatrix) = colnames(all_ranks)[3:15]
rownames(CorrMatrix) = colnames(all_ranks)[3:15]
pheatmap(CorrMatrix, cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

## With Counts
CorrMatrix = cor(all_counts[(all_counts$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_counts)[3:15]])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = 13)
colnames(CorrMatrix) = colnames(all_counts)[3:15]
rownames(CorrMatrix) = colnames(all_counts)[3:15]
pheatmap(CorrMatrix, cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

## With normalized counts
CorrMatrix = cor(all_ncounts[(all_ncounts$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_ncounts)[3:15]])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = 13)
colnames(CorrMatrix) = colnames(all_ncounts)[3:15]
rownames(CorrMatrix) = colnames(all_ncounts)[3:15]
pheatmap(CorrMatrix, cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

####
#### FIGURE3 Corr_PARCLIP_Motifs
#### with PAR-CLIP
## with Ranks
pheatmap(all_ranks[(all_ranks$MOTIF != 'WWWUAUUUAUUUW') & (all_ranks$MOTIF != 'UUUUU'), colnames(all_ranks)[2:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))
pheatmap(all_ranks[(all_ranks$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_ranks)[2:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))

## with counts
pheatmap(all_counts[(all_counts$MOTIF != 'WWWUAUUUAUUUW') & (all_counts$MOTIF != 'UUUUU'), colnames(all_counts)[2:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
pheatmap(all_counts[(all_counts$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_counts)[2:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

## with ncounts
pheatmap(all_ncounts[(all_ncounts$MOTIF != 'WWWUAUUUAUUUW') & (all_ncounts$MOTIF != 'UUUUU'), colnames(all_ncounts)[2:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
pheatmap(all_ncounts[(all_ncounts$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_ncounts)[2:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

#### without PAR-CLIP
## without Ranks
pheatmap(all_ranks[(all_ranks$MOTIF != 'WWWUAUUUAUUUW') & (all_ranks$MOTIF != 'UUUUU'), colnames(all_ranks)[3:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))
pheatmap(all_ranks[(all_ranks$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_ranks)[3:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))

## with counts
pheatmap(all_counts[(all_counts$MOTIF != 'WWWUAUUUAUUUW') & (all_counts$MOTIF != 'UUUUU'), colnames(all_counts)[3:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
pheatmap(all_counts[(all_counts$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_counts)[3:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

## with ncounts
pheatmap(all_ncounts[(all_ncounts$MOTIF != 'WWWUAUUUAUUUW') & (all_ncounts$MOTIF != 'UUUUU'), colnames(all_ncounts)[3:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
pheatmap(all_ncounts[(all_ncounts$MOTIF != 'WWWUAUUUAUUUW'), colnames(all_ncounts)[3:15]], cluster_rows=FALSE, cluster_cols=FALSE, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))


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




