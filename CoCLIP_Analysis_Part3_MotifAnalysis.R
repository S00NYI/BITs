## CoCLIP Analysis: 
## Motif Analysis
## Written by Soon Yi
## Use Homer output for subsequent analysis.
## Last Edit: 2023-09-15

library(stringr)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(corrplot)
library(pheatmap)
library(scales)
library(RColorBrewer)

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

## Count motif occurrence:
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

## Return density table that fits my format:
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

## Plot density plot:
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
####################################################################################################################

## Build peak sequence table and motif enrichment ranks:
####################################################################################################################
library(GenomicRanges)
library(IRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd('/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/motifEnrichment/peaks/R_Ready')
peakFiles = list.files(getwd())

motifs = c("UUUUU", "AAAAA", 
           "WUUUA", "YUUUA", "AUUUY", "WUUAA", "UUAAW", "UWUAA", 
           "YUUAA", "UYUAA", "AWUUA", "WAAUU", "AUWUA", "WUAAU", 
           "UAAUW", "AAUWU", "AAUUY", "AUUYA", "AAUYU", "AYUUA", 
           "WUAAA", "UAAAW", "UAUYA", "WAAAU", "YAUUA", "AUUAY", 
           "AYUAU", "UYAAA", "YUUCA", "ACUUY", "AUAUY", "AAAYU", 
           "AAAUW", "AAUWA", "WAAAA", "AACWU")

## extension + 10nt around peak center
extension = 65

motif_ranks = data.frame(PAR_CLIP = c(0, 0, 1:34))
rownames(motif_ranks) = motifs

for (peakFile in peakFiles) {
  print(paste0('processing ', peakFile))
  file = read.delim(peakFile, header = F)
  file = file[, c('V2', 'V3', 'V4', 'V1', 'V5')]
  colnames(file) = c('chrom', 'start', 'end', 'name', 'strand')
  file = file %>% filter(chrom %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
  peaksGR = file[, c('chrom', 'start', 'end', 'strand', 'name')]
  peaksGR$start = as.integer(peaksGR$start) + 1 - extension
  peaksGR$end = as.integer(peaksGR$end) + extension
  peaksGR = GRanges(peaksGR)
  peaksGR_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38, peaksGR, as.character = TRUE)
  peaksGR_seqs = cbind(file, data.frame(sequence = peaksGR_seqs))
  
  for (motif in motifs) {
    print(paste0('processing ', motif))
    if (sum(unlist(strsplit(motif, split = "")) %in% c('W', 'Y'))) {
      original_motif = motif
      NonStandard_nt = unlist(strsplit(motif, split = ""))[unlist(strsplit(motif, split = "")) %in% c('W', 'Y')]
      if (NonStandard_nt == 'W') {
        temp = sub('W', 'A', motif)
        motif = c(temp, sub('W', 'U', motif))
        motif_counts = motifCounts(peaksGR_seqs, motif)
        motif_counts = data.frame(rowSums(motif_counts))
        colnames(motif_counts) = original_motif
        peaksGR_seqs = cbind(peaksGR_seqs, motif_counts)
      } else if (NonStandard_nt == 'Y') {
        temp = sub('Y', 'C', motif)
        motif = c(temp, sub('Y', 'U', motif))
        motif_counts = motifCounts(peaksGR_seqs, motif)
        motif_counts = data.frame(rowSums(motif_counts))
        colnames(motif_counts) = original_motif
        peaksGR_seqs = cbind(peaksGR_seqs, motif_counts)
      }
    } else {
      motif_counts = motifCounts(peaksGR_seqs, motif)
      peaksGR_seqs = cbind(peaksGR_seqs, motif_counts)
    }
  }
  
  motif_rank = data.frame(rank(-colSums(peaksGR_seqs[, motifs])))
  colnames(motif_rank) = str_split(peakFile, '\\.')[[1]][1]
  
  motif_ranks = cbind(motif_ranks, motif_rank)
  
}

eCLIP = read.delim('/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/motifEnrichment/eCLIP_bed/ELAVL1_eCLIP_ENCFF566LNK.bed', header = F)
eCLIP = eCLIP[, c('V1', 'V2', 'V3', 'V4', 'V6')]
colnames(eCLIP) = c('chrom', 'start', 'end', 'name', 'strand')
eCLIP = eCLIP %>% filter(chrom %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
peaksGR = eCLIP[, c('chrom', 'start', 'end', 'strand', 'name')]
peaksGR$start = as.integer(peaksGR$start) + 1 - extension
peaksGR$end = as.integer(peaksGR$end) + extension
peaksGR = GRanges(peaksGR)
peaksGR_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38, peaksGR, as.character = TRUE)
peaksGR_seqs = cbind(eCLIP, data.frame(sequence = peaksGR_seqs))

for (motif in motifs) {
  print(paste0('processing ', motif))
  if (sum(unlist(strsplit(motif, split = "")) %in% c('W', 'Y'))) {
    original_motif = motif
    NonStandard_nt = unlist(strsplit(motif, split = ""))[unlist(strsplit(motif, split = "")) %in% c('W', 'Y')]
    if (NonStandard_nt == 'W') {
      temp = sub('W', 'A', motif)
      motif = c(temp, sub('W', 'U', motif))
      motif_counts = motifCounts(peaksGR_seqs, motif)
      motif_counts = data.frame(rowSums(motif_counts))
      colnames(motif_counts) = original_motif
      peaksGR_seqs = cbind(peaksGR_seqs, motif_counts)
    } else if (NonStandard_nt == 'Y') {
      temp = sub('Y', 'C', motif)
      motif = c(temp, sub('Y', 'U', motif))
      motif_counts = motifCounts(peaksGR_seqs, motif)
      motif_counts = data.frame(rowSums(motif_counts))
      colnames(motif_counts) = original_motif
      peaksGR_seqs = cbind(peaksGR_seqs, motif_counts)
    }
  } else {
    motif_counts = motifCounts(peaksGR_seqs, motif)
    peaksGR_seqs = cbind(peaksGR_seqs, motif_counts)
  }
}

motif_rank = data.frame(rank(-colSums(peaksGR_seqs[, motifs])))
colnames(motif_rank) = 'eCLIP'

motif_ranks = cbind(motif_ranks, motif_rank)
motif_ranks = cbind(motifs, motif_ranks)

####################################################################################################################

## Make heatmap:
####################################################################################################################
col_selection = c('ALL_M', 'I_M', 'NUC_M', 'CYT_M', 'NLS_M', 'NES_M', 'G3BP_M', 
                  'ALL_S', 'I_S', 'NUC_S', 'CYT_S', 'NLS_S', 'NES_S', 'G3BP_S')

CorrMatrix = cor(motif_ranks[, col_selection])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
colnames(CorrMatrix) = col_selection
rownames(CorrMatrix) = col_selection

pheatmap(CorrMatrix, cluster_rows=F, cluster_cols=F, color = colorRampPalette(brewer.pal(9, "GnBu"))(100))
# pheatmap(CorrMatrix, cluster_rows=F, cluster_cols=F, color = (colorRampPalette(c("lemonchiffon1", "skyblue", "skyblue4"))(100)))
# pheatmap(CorrMatrix, clustering_method = 'ward.D2', cluster_rows=T, cluster_cols=T, color = colorRampPalette(brewer.pal(9, "GnBu"))(100))


col_selection = c('PAR_CLIP', 'eCLIP', 
                  'ALL_M', 'I_M', 'NUC_M', 'CYT_M', 'NLS_M', 'NES_M', 'G3BP_M', 
                  'ALL_S', 'I_S', 'NUC_S', 'CYT_S', 'NLS_S', 'NES_S', 'G3BP_S')

# pheatmap(motif_ranks[, col_selection], cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))

pheatmap(motif_ranks[, col_selection], cluster_rows=F, cluster_cols=F, color = rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))

## Figure 5
col_selection = c('CoCLIP_Input_Mock', 'CoCLIP_G3BP_Mock', 
                  'CoCLIP_Input_Arsenite', 'CoCLIP_G3BP_Arsenite')

pheatmap((all_ranks %>% arrange(CoCLIP_G3BP_Arsenite))[, col_selection], cluster_rows=F, cluster_cols=F, color = rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))
pheatmap((all_ranks %>% arrange(CoCLIP_Input_Mock))[, col_selection], cluster_rows=F, cluster_cols=F, color = rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))
pheatmap(all_ranks[, col_selection], cluster_rows=F, cluster_cols=F, color = rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))
####################################################################################################################



## Data processing for motif density
################################################################################
# baseDir = 'L:/.shortcut-targets-by-id/13hY9t_p6eUdnvP-c2OClHbzyeWMOruD_/CoCLIP_HuR_Paper/Data/Homer_Outputs/motifs_density/50bp/'
baseDir = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/motifEnrichment/motifs/density'
setwd(baseDir)

densityFiles = list.files(paste0(baseDir))
densityFiles_CoCLIP = densityFiles[grep("^CoCLIP", densityFiles)]
densityFiles_FracCLIP  = densityFiles[grep("FracCLIP", densityFiles)]
################################################################################

## Metagene Plot from All Libraries for Top Motifs UUUUU AAAAA 
################################################################################
## Mock
densityFile = return_Density(densityFiles[2], strand = '+', normalize = T)
plot_Density(densityFile, c('UUUUU', 'AAAAA'), yaxis_lims = c(0, 0.06))

## Arsenite
densityFile = return_Density(densityFiles[1], strand = '+', normalize = T)
plot_Density(densityFile, c('UUUUU', 'AAAAA'), yaxis_lims = c(0, 0.06))
################################################################################

## FIGURE3/5 C Metagene Plot from Input vs CoCLIP for Top Motifs UUUUU AAAAA 
################################################################################
## Mock
Mock_All = return_Density(densityFiles[2], strand = '+', normalize = T) 
Mock_Input = return_Density(densityFiles[6], strand = '+', normalize = T)
Mock_NES = return_Density(densityFiles[8], strand = '+', normalize = T)
Mock_NLS = return_Density(densityFiles[10], strand = '+', normalize = T)
Mock_G3BP = return_Density(densityFiles[4], strand = '+', normalize = T)
Mock_Nuclear = return_Density(densityFiles[14], strand = '+', normalize = T)
Mock_Cytoplasm = return_Density(densityFiles[12], strand = '+', normalize = T)

Mock_Combined = data.frame(cbind(Mock_All$position, 
                                 Mock_All$UUUUU, Mock_Input$UUUUU, Mock_NES$UUUUU, Mock_NLS$UUUUU, Mock_G3BP$UUUUU, Mock_Nuclear$UUUUU, Mock_Cytoplasm$UUUUU,
                                 Mock_All$AAAAA, Mock_Input$AAAAA, Mock_NES$AAAAA, Mock_NLS$AAAAA, Mock_G3BP$AAAAA, Mock_Nuclear$AAAAA, Mock_Cytoplasm$AAAAA))

colnames(Mock_Combined) = c('position', 
                            'UUUUU_All', 'UUUUU_Inp', 'UUUUU_NES', 'UUUUU_NLS', 'UUUUU_G3BP', 'UUUUU_Nuc', 'UUUUU_Cyt',
                            'AAAAA_All', 'AAAAA_Inp', 'AAAAA_NES', 'AAAAA_NLS', 'AAAAA_G3BP', 'AAAAA_Nuc', 'AAAAA_Cyt')

plot_Density(Mock_Combined, c('UUUUU_Inp', 'UUUUU_NES', 'UUUUU_NLS'), yaxis_lims = c(0, 0.15))
plot_Density(Mock_Combined, c('AAAAA_Inp', 'AAAAA_NES', 'AAAAA_NLS'), yaxis_lims = c(0, 0.15))

## Arsenite
Arsenite_All = return_Density(densityFiles[1], strand = '+', normalize = T) 
Arsenite_Input = return_Density(densityFiles[5], strand = '+', normalize = T) 
Arsenite_NES = return_Density(densityFiles[7], strand = '+', normalize = T) 
Arsenite_NLS = return_Density(densityFiles[9], strand = '+', normalize = T) 
Arsenite_G3BP = return_Density(densityFiles[3], strand = '+', normalize = T) 
Arsenite_Nuclear = return_Density(densityFiles[13], strand = '+', normalize = T) 
Arsenite_Cytoplasm = return_Density(densityFiles[11], strand = '+', normalize = T) 

Arsenite_Combined = data.frame(cbind(Arsenite_All$position, 
                                     Arsenite_All$UUUUU, Arsenite_Input$UUUUU, Arsenite_NES$UUUUU, Arsenite_NLS$UUUUU, Arsenite_G3BP$UUUUU, Arsenite_Nuclear$UUUUU, Arsenite_Cytoplasm$UUUUU,
                                     Arsenite_All$AAAAA, Arsenite_Input$AAAAA, Arsenite_NES$AAAAA, Arsenite_NLS$AAAAA, Arsenite_G3BP$AAAAA, Arsenite_Nuclear$AAAAA, Arsenite_Cytoplasm$AAAAA))

colnames(Arsenite_Combined) = c('position', 
                                'UUUUU_All', 'UUUUU_Inp', 'UUUUU_NES', 'UUUUU_NLS', 'UUUUU_G3BP', 'UUUUU_Nuc', 'UUUUU_Cyt',
                                'AAAAA_All', 'AAAAA_Inp', 'AAAAA_NES', 'AAAAA_NLS', 'AAAAA_G3BP', 'AAAAA_Nuc', 'AAAAA_Cyt')

plot_Density(Arsenite_Combined, c('UUUUU_Inp', 'UUUUU_NES', 'UUUUU_NLS'), yaxis_lims = c(0, 0.15))
plot_Density(Arsenite_Combined, c('AAAAA_Inp', 'AAAAA_NES', 'AAAAA_NLS'), yaxis_lims = c(0, 0.15))

## Figure 5
plot_Density(Mock_Combined, c('UUUUU_Inp', 'UUUUU_G3BP'), yaxis_lims = c(0, 0.15))
plot_Density(Mock_Combined, c('AAAAA_Inp', 'AAAAA_G3BP'), yaxis_lims = c(0, 0.15))
plot_Density(Arsenite_Combined, c('UUUUU_Inp', 'UUUUU_G3BP'), yaxis_lims = c(0, 0.15))
plot_Density(Arsenite_Combined, c('AAAAA_Inp', 'AAAAA_G3BP'), yaxis_lims = c(0, 0.15))

################################################################################

## FIGURE3 D
## Metagene Plot from Fractionation CLIP for Top Motifs UUUUU AAAAA
################################################################################
## Mock

plot_Density(Mock_Combined, c('UUUUU_Inp', 'UUUUU_Nuc', 'UUUUU_Cyt'), yaxis_lims = c(0, 0.15))
plot_Density(Mock_Combined, c('AAAAA_Inp', 'AAAAA_Nuc', 'AAAAA_Cyt'), yaxis_lims = c(0, 0.15))

plot_Density(Arsenite_Combined, c('UUUUU_Inp', 'UUUUU_Nuc', 'UUUUU_Cyt'), yaxis_lims = c(0, 0.15))
plot_Density(Arsenite_Combined, c('AAAAA_Inp', 'AAAAA_Nuc', 'AAAAA_Cyt'), yaxis_lims = c(0, 0.15))
################################################################################









################################################################################
###################################DEPRECATED###################################
################################################################################








## data processing for motif counts
################################################################################
baseDir = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/motifEnrichment/motifs/counts'
# baseDir = 'L:/.shortcut-targets-by-id/13hY9t_p6eUdnvP-c2OClHbzyeWMOruD_/CoCLIP_HuR_Paper/Data/Homer_Outputs/motifs/counts'
setwd(baseDir)
countFiles = list.files(baseDir)
countFiles = countFiles[grep('MotifCounts.txt', countFiles)]

peakCounts = data.frame(sample = sapply(str_split(countFiles, '\\.'), function(x) x[1]), 
                        peaks = c(287466, 190803, 10924, 9666, 81323, 48659, 3250, 694, 4725, 4699, 8823, 102500, 65923, 80733, 86254))

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

# new_column_names = c('MOTIF', 'PAR_CLIP', 'eCLIP', 
#                      'All_Mocks', 'All_Arsenites',
#                      'CoCLIP_Input_Mock', 'CoCLIP_Input_Arsenite',
#                      'CoCLIP_NLS_Mock', 'CoCLIP_NLS_Arsenite',
#                      'CoCLIP_NES_Mock', 'CoCLIP_NES_Arsenite',
#                      'CoCLIP_G3BP_Mock', 'CoCLIP_G3BP_Arsenite',
#                      'FracCLIP_Nuclear_Mock', 'FracCLIP_Nuclear_Arsenite',
#                      'FracCLIP_Cytoplasm_Mock', 'FracCLIP_Cytoplasm_Arsenite')

new_column_names = c('MOTIF', 'PAR_CLIP', 'eCLIP', 
                     'All_Mocks', 
                     'CoCLIP_Input_Mock', 
                     'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock', 
                     'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'CoCLIP_G3BP_Mock', 
                     'All_Arsenites',
                     'CoCLIP_Input_Arsenite', 
                     'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite',
                     'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'CoCLIP_G3BP_Arsenite')

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

# all_ranks$eCLIP = all_ranks$eCLIP - 1
# all_ranks$eCLIP[1] = all_ranks$eCLIP[1] = 0
# all_ranks$eCLIP[2] = all_ranks$eCLIP[2] = 0
# write.table(all_unique_PPM, 'HuR_PAR_CLIP_MotifCounts.txt', row.names = T, col.names = T, quote = F, sep = '\t')
# write.table(all_ranks, 'HuR_PAR_CLIP_MotifRanks.txt', row.names = T, col.names = T, quote = F, sep = '\t')
################################################################################

## Correlation Matrix with PAR-CLIP data and just the PAR-CLIP motifs
################################################################################
## With Ranks
row_selection = (all_ranks$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
col_selection = colnames(all_ranks)[2:17]
col_selection = c('PAR_CLIP', 'eCLIP',  'All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock') # Mock Only
col_selection = c('PAR_CLIP', 'All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

col_selection = c('PAR_CLIP', 'eCLIP', 
                  'All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock', 
                  'All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite')

CorrMatrix = cor(all_ranks[row_selection, col_selection])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
colnames(CorrMatrix) = col_selection
rownames(CorrMatrix) = col_selection
# pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1", "indianred2"))(100)))

pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = colorRampPalette(brewer.pal(9, "YlGnBu"))(100))

################################################################################

## Correlation Matrix without PAR-CLIP data and just the PAR-CLIP motifs
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

## FIGURE3/5 Correlation Matrix without PAR-CLIP data and with UUUUU and AAAAA
################################################################################
## With Ranks
# col_selection = colnames(all_ranks)[3:16]
# col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
# col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 
                  'All_Arsenites', 'CoCLIP_Input_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite')

CorrMatrix = cor(all_ranks[, col_selection])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
colnames(CorrMatrix) = col_selection
rownames(CorrMatrix) = col_selection

pheatmap(CorrMatrix, clustering_method = 'ward.D2', cluster_rows=F, cluster_cols=F, color = colorRampPalette(brewer.pal(9, "GnBu"))(100))

## Figure 5
col_selection = c('CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'CoCLIP_G3BP_Mock', 
                  'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'CoCLIP_G3BP_Arsenite')

col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'CoCLIP_G3BP_Mock', 
                  'All_Arsenites', 'CoCLIP_Input_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'CoCLIP_G3BP_Arsenite')

CorrMatrix = cor(all_ranks[, col_selection])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
colnames(CorrMatrix) = col_selection
rownames(CorrMatrix) = col_selection

pheatmap(CorrMatrix, clustering_method = 'ward.D2', cluster_rows=F, cluster_cols=F, color = colorRampPalette(brewer.pal(9, "GnBu"))(100))

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

# ## With unique peaks per motif normalized counts
# ## This looks equal to un-normalized version
# col_selection = colnames(all_unique_PPM_normed)[3:16]
# col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
# col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only
# 
# CorrMatrix = cor(all_unique_PPM_normed[, col_selection])
# CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
# colnames(CorrMatrix) = col_selection
# rownames(CorrMatrix) = col_selection
# pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))

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
# 
# ## With raw motif normalized counts
# ## This looks equal to un-normalized version
# col_selection = colnames(all_raw_Motifs_normed)[4:16]
# col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
# col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only
# 
# CorrMatrix = cor(all_raw_Motifs_normed[, col_selection])
# CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
# colnames(CorrMatrix) = col_selection
# rownames(CorrMatrix) = col_selection
# pheatmap(CorrMatrix, cluster_rows=T, cluster_cols=T, color = rev(colorRampPalette(c("green4", "lemonchiffon1"))(100)))
################################################################################

## FIGURE3/5 Motif Score Heatmap with PAR-CLIP data and just the PAR-CLIP motifs
################################################################################
## with Ranks
# row_selection = (all_ranks$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
# col_selection = colnames(all_ranks)[2:17]
# col_selection = c('PAR_CLIP', 'eCLIP', 'All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock') # Mock Only
# col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

col_selection = c('PAR_CLIP', 'eCLIP', 
                  'All_Mocks', 'CoCLIP_Input_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 
                  'All_Arsenites', 'CoCLIP_Input_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite')

# pheatmap(all_ranks[row_selection, col_selection], cluster_rows=F, cluster_cols=T, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))

pheatmap(all_ranks[, col_selection], cluster_rows=F, cluster_cols=F, color = rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))

## Figure 5
col_selection = c('CoCLIP_Input_Mock', 'CoCLIP_G3BP_Mock', 
                  'CoCLIP_Input_Arsenite', 'CoCLIP_G3BP_Arsenite')

pheatmap((all_ranks %>% arrange(CoCLIP_G3BP_Arsenite))[, col_selection], cluster_rows=F, cluster_cols=F, color = rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))
pheatmap((all_ranks %>% arrange(CoCLIP_Input_Mock))[, col_selection], cluster_rows=F, cluster_cols=F, color = rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))
pheatmap(all_ranks[, col_selection], cluster_rows=F, cluster_cols=F, color = rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))


################################################################################

## Motif Score Heatmap without PAR-CLIP data and just the PAR-CLIP motifs
################################################################################
## With Ranks
row_selection = (all_ranks$MOTIF != 'AAAAA' & all_ranks$MOTIF != 'UUUUU')
col_selection = colnames(all_ranks)[3:16]
col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

pheatmap(all_ranks[row_selection, col_selection], cluster_rows=F, cluster_cols=T, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))

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

## Motif Score Heatmap without PAR-CLIP data and with UUUUU and AAAAA
################################################################################
## With Ranks
col_selection = colnames(all_ranks)[3:15]
col_selection = c('All_Mocks', 'CoCLIP_Input_Mock', 'CoCLIP_NLS_Mock', 'CoCLIP_NES_Mock', 'FracCLIP_Nuclear_Mock', 'FracCLIP_Cytoplasm_Mock')
col_selection = c('All_Arsenites', 'CoCLIP_Input_Arsenite', 'CoCLIP_NLS_Arsenite', 'CoCLIP_NES_Arsenite', 'FracCLIP_Nuclear_Arsenite', 'FracCLIP_Cytoplasm_Arsenite') # Arsenite Only

pheatmap(all_ranks[, col_selection], cluster_rows=F, cluster_cols=T, color = rev(colorRampPalette(c("lemonchiffon1", "green4"))(100)))

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

pheatmap(all_raw_Motifs_normed[, col_selection], cluster_rows=F, cluster_cols=T, color = (colorRampPalette(c("lemonchiffon1", "green4"))(100)))
################################################################################
