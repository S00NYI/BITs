## CoCLIP Analysis: 
## Motif Analysis
## Written by Soon Yi
## Use Homer output for subsequent analysis.
## Last Edit: 2023-10-12

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
library(GenomicRanges)
library(IRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

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
plot_Density = function(density_data, columns_list, xaxis_lims = NULL, yaxis_lims = NULL, custom_colors = NULL, densityType = NULL, sampleName = NULL, featureName = NULL, smoothing = NULL) {
  # Select the columns based on the columns_list
  # plot_data = density_data %>% select(position, {{columns_list}})
  plot_data = density_data[, c('position', columns_list)]
  
  # Create a long-format dataframe for better legend handling
  plot_data_long = plot_data %>%
    pivot_longer(cols = {{columns_list}}, names_to = "Data", values_to = "Density")
  
  # Define the order of levels for the Data factor
  plot_data_long$Data = factor(plot_data_long$Data, levels = columns_list)
  
  # Create the ggplot object
  plot = ggplot(plot_data_long, aes(x = position, y = Density, color = Data)) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = 'bold')) +
    labs(y = "Peak Density")
  
  # Add smoothed lines
  if (!is.null(smoothing)) {
    plot = plot + geom_smooth(span = smoothing, se = FALSE)
  } else {
    plot = plot + geom_line()
  }
  
  if (!is.null(custom_colors)) {
    plot = plot + scale_color_manual(values = custom_colors)
  }
  
  if (!is.null(xaxis_lims)) {
    plot = plot + xlim(xaxis_lims)
  }
  
  if (!is.null(yaxis_lims)) {
    plot = plot + ylim(yaxis_lims)
  }
  
  if (is.null(densityType)) {
    if (!is.null(sampleName)) {
      plot = plot + labs(title = paste0('Density for ', sampleName),
                         x = "Distance from center (nucleotides)")
      
    } else {
      plot = plot + labs(title = 'Density Plot',
                         x = "Distance from center (nucleotides)")
    }
    
  } else if (densityType == 'motif_density') {
    if (!is.null(sampleName)) {
      plot = plot + labs(title = paste0('Motif Density around Peaks for ', sampleName),
                         x = "Distance from peak center (nucleotides)")
    } else {
      plot = plot + labs(title = 'Motif Density around Peaks',
                         x = "Distance from peak center (nucleotides)")
    }
  } else if (densityType == 'feature_metagene') {
    if (!is.null(sampleName)) {
      if (!is.null(featureName)) {
        plot = plot + labs(title = paste0('Peaks density around ', featureName, ' for ', sampleName),
                           x = paste0("Distance from ", featureName, " (nucleotides)"))
      } else {
        plot = plot + labs(title = 'Peaks density around feature',
                           x = "Distance from feature (nucleotides)")
      }
    } else {
      if (!is.null(featureName)) {
        plot = plot + labs(title = paste0('Peaks density around ', featureName),
                           x = paste0("Distance from ", featureName, " (nucleotides)"))
      } else {
        plot = plot + labs(title = 'Peaks density around feature',
                           x = "Distance from feature (nucleotides)")
      }
    }
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

## Build peak sequence table and motif enrichment ranks based on peaksMatrix:
####################################################################################################################
## (extension + 10nt)*2 centered at peak
# extension = 15
extension = 65
# extension = 0

peaksGR = peaksMatrix[, c('chrom', 'start', 'end', 'strand', 'peak_names')]
peaksGR = peaksGR %>% mutate(chrom = ifelse(chrom == "chrMT", "chrM", chrom))
peaksGR$start = as.integer(peaksGR$start) + 1 - extension
peaksGR$end = as.integer(peaksGR$end) + extension
peaksGR = GRanges(peaksGR)

peaksGR_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38, peaksGR, as.character = TRUE)
peaksGR_seqs = cbind(peaksMatrix[, c('chrom', 'start', 'end', 'strand', 'peak_names', BC_columns, 'grouped_annotation', 'finalized_annotation')], data.frame(sequence = peaksGR_seqs))

motifs = c("UUUUU", "AAAAA", 
           "WUUUA", "YUUUA", "AUUUY", "WUUAA", "UUAAW", "UWUAA", 
           "YUUAA", "UYUAA", "AWUUA", "WAAUU", "AUWUA", "WUAAU", 
           "UAAUW", "AAUWU", "AAUUY", "AUUYA", "AAUYU", "AYUUA", 
           "WUAAA", "UAAAW", "UAUYA", "WAAAU", "YAUUA", "AUUAY", 
           "AYUAU", "UYAAA", "YUUCA", "ACUUY", "AUAUY", "AAAYU", 
           "AAAUW", "AAUWA", "WAAAA", "AACWU")


motif_ranks = data.frame(PAR_CLIP = c(0, 0, 1:34))
rownames(motif_ranks) = motifs

ALL_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_F_Nuc_M$peak_names, Peak_F_Cyt_M$peak_names, Peak_Co_Input_M$peak_names, Peak_Co_NLS_M$peak_names, Peak_Co_NES_M$peak_names, Peak_Co_G3BP_M$peak_names))
I_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_Input_M$peak_names))
NUC_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_F_Nuc_M$peak_names))
CYT_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_F_Cyt_M$peak_names))
NLS_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NLS_M$peak_names))
NES_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NES_M$peak_names))
G3BP_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_G3BP_M$peak_names))

ALL_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_F_Nuc_S$peak_names, Peak_F_Cyt_S$peak_names, Peak_Co_Input_S$peak_names, Peak_Co_NLS_S$peak_names, Peak_Co_NES_S$peak_names, Peak_Co_G3BP_S$peak_names))
I_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_Input_S$peak_names))
NUC_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_F_Nuc_S$peak_names))
CYT_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_F_Cyt_S$peak_names))
NLS_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NLS_S$peak_names))
NES_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NES_S$peak_names))
G3BP_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_G3BP_S$peak_names))

Samples = c('ALL_M', 'I_M', 'NUC_M', 'CYT_M', 'NLS_M', 'NES_M', 'G3BP_M', 'ALL_S', 'I_S', 'NUC_S', 'CYT_S', 'NLS_S', 'NES_S', 'G3BP_S') 

for (sample in Samples) {
  INPUT = get(sample)
  for (motif in motifs) {
    if (sum(unlist(strsplit(motif, split = "")) %in% c('W', 'Y'))) {
      original_motif = motif
      NonStandard_nt = unlist(strsplit(motif, split = ""))[unlist(strsplit(motif, split = "")) %in% c('W', 'Y')]
      if (NonStandard_nt == 'W') {
        temp = sub('W', 'A', motif)
        motif = c(temp, sub('W', 'U', motif))
        motif_counts = motifCounts(INPUT, motif)
        motif_counts = data.frame(rowSums(motif_counts))
        colnames(motif_counts) = original_motif
        INPUT = cbind(INPUT, motif_counts)
      } else if (NonStandard_nt == 'Y') {
        temp = sub('Y', 'C', motif)
        motif = c(temp, sub('Y', 'U', motif))
        motif_counts = motifCounts(INPUT, motif)
        motif_counts = data.frame(rowSums(motif_counts))
        colnames(motif_counts) = original_motif
        INPUT = cbind(INPUT, motif_counts)
      }
    } else {
      motif_counts = motifCounts(INPUT, motif)
      INPUT = cbind(INPUT, motif_counts)
    }
  }
  motif_rank = data.frame(rank(-colSums(INPUT[, motifs])))
  colnames(motif_rank) = sample
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

## FIGURE4 Make heatmap:
####################################################################################################################
col_selection = c('ALL_M', 'I_M', 'NUC_M', 'CYT_M', 'NLS_M', 'NES_M', 'G3BP_M', 
                  'ALL_S', 'I_S', 'NUC_S', 'CYT_S', 'NLS_S', 'NES_S', 'G3BP_S')

CorrMatrix = cor(motif_ranks[, col_selection])
CorrMatrix = matrix(round(CorrMatrix,2), nrow = length(col_selection))
colnames(CorrMatrix) = col_selection
rownames(CorrMatrix) = col_selection

pheatmap(CorrMatrix, cluster_rows=F, cluster_cols=F, color = colorRampPalette(brewer.pal(9, "GnBu"))(100))

col_selection = c('PAR_CLIP', 'eCLIP', 
                  'ALL_M', 'I_M', 'NUC_M', 'CYT_M', 'NLS_M', 'NES_M', 'G3BP_M', 
                  'ALL_S', 'I_S', 'NUC_S', 'CYT_S', 'NLS_S', 'NES_S', 'G3BP_S')


pheatmap(motif_ranks[, col_selection], cluster_rows=F, cluster_cols=F, color = rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))
####################################################################################################################

## Read motif density calculations:
## motif density calculated using Homer (see Homer_calls.sh)
################################################################################
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
plot_Density(densityFile, c('UUUUU', 'AAAAA'), yaxis_lims = c(0, 0.06), densityType = 'motif_density', sampleName = 'All Mocks')

## Arsenite
densityFile = return_Density(densityFiles[1], strand = '+', normalize = T)
plot_Density(densityFile, c('UUUUU', 'AAAAA'), yaxis_lims = c(0, 0.06), densityType = 'motif_density', sampleName = 'All Mocks')
################################################################################

## FIGURE4 Metagene Plot from Input vs CoCLIP for Top Motifs UUUUU AAAAA 
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

plot_Density(Mock_Combined, c('UUUUU_Inp', 'UUUUU_NLS', 'UUUUU_NES', 'UUUUU_G3BP'), yaxis_lims = c(0, 0.15), custom_colors = c('ivory4', 'skyblue', 'darkseagreen2', 'salmon'), densityType = 'motif_density', sampleName = 'UUUUU')
plot_Density(Mock_Combined, c('AAAAA_Inp', 'AAAAA_NLS', 'AAAAA_NES', 'AAAAA_G3BP'), yaxis_lims = c(0, 0.15), custom_colors = c('ivory4', 'skyblue', 'darkseagreen2', 'salmon'), densityType = 'motif_density', sampleName = 'AAAAA')

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

plot_Density(Arsenite_Combined, c('UUUUU_Inp', 'UUUUU_NLS', 'UUUUU_NES', 'UUUUU_G3BP'), yaxis_lims = c(0, 0.15), custom_colors = c('ivory4', 'skyblue', 'darkseagreen2', 'salmon'), densityType = 'motif_density', sampleName = 'UUUUU')
plot_Density(Arsenite_Combined, c('AAAAA_Inp', 'AAAAA_NLS', 'AAAAA_NES', 'AAAAA_G3BP'), yaxis_lims = c(0, 0.15), custom_colors = c('ivory4', 'skyblue', 'darkseagreen2', 'salmon'), densityType = 'motif_density', sampleName = 'AAAAA')
################################################################################

## FIGURE4 Motif Cumulative Distribution Analysis
####################################################################################################################
peaksGR = peaksMatrix[, c('chrom', 'start', 'end', 'strand', 'peak_names')]
peaksGR = peaksGR %>% mutate(chrom = ifelse(chrom == "chrMT", "chrM", chrom))
peaksGR$start = as.integer(peaksGR$start) + 1 - 15
peaksGR$end = as.integer(peaksGR$end) + 15
peaksGR = GRanges(peaksGR)

peaksGR_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38, peaksGR, as.character = TRUE)
peaksGR_seqs = cbind(peaksMatrix[, c('chrom', 'start', 'end', 'strand', 'peak_names', BC_columns, 'grouped_annotation', 'finalized_annotation')], data.frame(sequence = peaksGR_seqs))

motifs = c('AAAAA', 'UUUUU')
motif_counts = motifCounts(peaksGR_seqs, motifs)
peaksGR_seqs = cbind(peaksGR_seqs, motif_counts)

peaksGR_seqs_org = peaksGR_seqs

peaksGR_seqs = peaksGR_seqs_org
peaksGR_seqs = peaksGR_seqs %>% filter(grouped_annotation != 'UnAn')

ALL_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_F_Nuc_M$peak_names, Peak_F_Cyt_M$peak_names, Peak_Co_Input_M$peak_names, Peak_Co_NLS_M$peak_names, Peak_Co_NES_M$peak_names, Peak_Co_G3BP_M$peak_names))
I_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_Input_M$peak_names))
NUC_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_F_Nuc_M$peak_names))
CYT_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_F_Cyt_M$peak_names))
NLS_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NLS_M$peak_names))
NES_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NES_M$peak_names))
G3BP_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_G3BP_M$peak_names))

ALL_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_F_Nuc_S$peak_names, Peak_F_Cyt_S$peak_names, Peak_Co_Input_S$peak_names, Peak_Co_NLS_S$peak_names, Peak_Co_NES_S$peak_names, Peak_Co_G3BP_S$peak_names))
I_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_Input_S$peak_names))
NUC_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_F_Nuc_S$peak_names))
CYT_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_F_Cyt_S$peak_names))
NLS_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NLS_S$peak_names))
NES_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NES_S$peak_names))
G3BP_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_G3BP_S$peak_names))

peaksGR_Co_Input_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_Input_M$peak_names))
peaksGR_Co_Input_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_Input_S$peak_names))

peaksGR_Co_NLS_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NLS_M$peak_names))
peaksGR_Co_NLS_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NLS_S$peak_names))

peaksGR_Co_NES_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NES_M$peak_names))
peaksGR_Co_NES_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NES_S$peak_names))

peaksGR_Co_G3BP_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_G3BP_M$peak_names))
peaksGR_Co_G3BP_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_G3BP_S$peak_names))


## Filter peaks to specific genomic locus
freq_table = data.frame(counts = 0:50)

freq_table = freq_table %>% mutate(I_M_AAAAA = rep(NA, 51))
freq_table = freq_table %>% mutate(I_M_UUUUU = rep(NA, 51))
freq_table = freq_table %>% mutate(NLS_M_AAAAA = rep(NA, 51))
freq_table = freq_table %>% mutate(NLS_M_UUUUU = rep(NA, 51))
freq_table = freq_table %>% mutate(NES_M_AAAAA = rep(NA, 51))
freq_table = freq_table %>% mutate(NES_M_UUUUU = rep(NA, 51))
freq_table = freq_table %>% mutate(G3BP_M_AAAAA = rep(NA, 51))
freq_table = freq_table %>% mutate(G3BP_M_UUUUU = rep(NA, 51))

freq_table = freq_table %>% mutate(I_S_AAAAA = rep(NA, 51))
freq_table = freq_table %>% mutate(I_S_UUUUU = rep(NA, 51))
freq_table = freq_table %>% mutate(NLS_S_AAAAA = rep(NA, 51))
freq_table = freq_table %>% mutate(NLS_S_UUUUU = rep(NA, 51))
freq_table = freq_table %>% mutate(NES_S_AAAAA = rep(NA, 51))
freq_table = freq_table %>% mutate(NES_S_UUUUU = rep(NA, 51))
freq_table = freq_table %>% mutate(G3BP_S_AAAAA = rep(NA, 51))
freq_table = freq_table %>% mutate(G3BP_S_UUUUU = rep(NA, 51))

freq_table$I_M_AAAAA[freq_table$counts %in% data.frame(table(peaksGR_Co_Input_M$AAAAA))$Var1] = data.frame(table(peaksGR_Co_Input_M$AAAAA))$Freq
freq_table$I_M_UUUUU[freq_table$counts %in% data.frame(table(peaksGR_Co_Input_M$UUUUU))$Var1] = data.frame(table(peaksGR_Co_Input_M$UUUUU))$Freq
freq_table$NLS_M_AAAAA[freq_table$counts %in% data.frame(table(peaksGR_Co_NLS_M$AAAAA))$Var1] = data.frame(table(peaksGR_Co_NLS_M$AAAAA))$Freq
freq_table$NLS_M_UUUUU[freq_table$counts %in% data.frame(table(peaksGR_Co_NLS_M$UUUUU))$Var1] = data.frame(table(peaksGR_Co_NLS_M$UUUUU))$Freq
freq_table$NES_M_AAAAA[freq_table$counts %in% data.frame(table(peaksGR_Co_NES_M$AAAAA))$Var1] = data.frame(table(peaksGR_Co_NES_M$AAAAA))$Freq
freq_table$NES_M_UUUUU[freq_table$counts %in% data.frame(table(peaksGR_Co_NES_M$UUUUU))$Var1] = data.frame(table(peaksGR_Co_NES_M$UUUUU))$Freq
freq_table$G3BP_M_AAAAA[freq_table$counts %in% data.frame(table(peaksGR_Co_G3BP_M$AAAAA))$Var1] = data.frame(table(peaksGR_Co_G3BP_M$AAAAA))$Freq
freq_table$G3BP_M_UUUUU[freq_table$counts %in% data.frame(table(peaksGR_Co_G3BP_M$UUUUU))$Var1] = data.frame(table(peaksGR_Co_G3BP_M$UUUUU))$Freq

freq_table$I_S_AAAAA[freq_table$counts %in% data.frame(table(peaksGR_Co_Input_S$AAAAA))$Var1] = data.frame(table(peaksGR_Co_Input_S$AAAAA))$Freq
freq_table$I_S_UUUUU[freq_table$counts %in% data.frame(table(peaksGR_Co_Input_S$UUUUU))$Var1] = data.frame(table(peaksGR_Co_Input_S$UUUUU))$Freq
freq_table$NLS_S_AAAAA[freq_table$counts %in% data.frame(table(peaksGR_Co_NLS_S$AAAAA))$Var1] = data.frame(table(peaksGR_Co_NLS_S$AAAAA))$Freq
freq_table$NLS_S_UUUUU[freq_table$counts %in% data.frame(table(peaksGR_Co_NLS_S$UUUUU))$Var1] = data.frame(table(peaksGR_Co_NLS_S$UUUUU))$Freq
freq_table$NES_S_AAAAA[freq_table$counts %in% data.frame(table(peaksGR_Co_NES_S$AAAAA))$Var1] = data.frame(table(peaksGR_Co_NES_S$AAAAA))$Freq
freq_table$NES_S_UUUUU[freq_table$counts %in% data.frame(table(peaksGR_Co_NES_S$UUUUU))$Var1] = data.frame(table(peaksGR_Co_NES_S$UUUUU))$Freq
freq_table$G3BP_S_AAAAA[freq_table$counts %in% data.frame(table(peaksGR_Co_G3BP_S$AAAAA))$Var1] = data.frame(table(peaksGR_Co_G3BP_S$AAAAA))$Freq
freq_table$G3BP_S_UUUUU[freq_table$counts %in% data.frame(table(peaksGR_Co_G3BP_S$UUUUU))$Var1] = data.frame(table(peaksGR_Co_G3BP_S$UUUUU))$Freq

freq_table[is.na(freq_table)] = 0

cumulative_table = cbind(freq_table$counts, cumsum(freq_table[, colnames(freq_table)[2:ncol(freq_table)]]))
colnames(cumulative_table)[1] = 'counts'

normed_table = cbind(freq_table$counts, data.frame(sapply(cumsum(freq_table[, colnames(freq_table)[2:ncol(freq_table)]]), rescale)))
colnames(normed_table)[1] = 'counts'

pVal_Threshold = 0.01

# CD plots and accompanying KS-test
ks.test(freq_table$I_M_AAAAA, freq_table$I_S_AAAAA)$p.value < pVal_Threshold
ks.test(freq_table$I_M_UUUUU, freq_table$I_S_UUUUU)$p.value < pVal_Threshold

ks.test(freq_table$NLS_M_AAAAA, freq_table$NLS_S_AAAAA)$p.value < pVal_Threshold
ks.test(freq_table$NLS_M_UUUUU, freq_table$NLS_S_UUUUU)$p.value < pVal_Threshold

ks.test(freq_table$NES_M_AAAAA, freq_table$NES_S_AAAAA)$p.value < pVal_Threshold
ks.test(freq_table$NES_M_UUUUU, freq_table$NES_S_UUUUU)$p.value < pVal_Threshold

ks.test(freq_table$G3BP_M_AAAAA, freq_table$G3BP_S_AAAAA)$p.value < pVal_Threshold
ks.test(freq_table$G3BP_M_UUUUU, freq_table$G3BP_S_UUUUU)$p.value < pVal_Threshold

plot_CD(normed_table, y_data = c('I_M_AAAAA', 'NLS_M_AAAAA', 'NES_M_AAAAA', 'G3BP_M_AAAAA'), colormaps = c('ivory4', 'skyblue', 'darkseagreen2', 'salmon'), linetypes = c(1, 1, 1, 1), linesize = 2)
plot_CD(normed_table, y_data = c('I_S_AAAAA', 'NLS_S_AAAAA', 'NES_S_AAAAA', 'G3BP_S_AAAAA'), colormaps = c('ivory4', 'skyblue', 'darkseagreen2', 'salmon'), linetypes = c(1, 1, 1, 1), linesize = 2)
plot_CD(normed_table, y_data = c('I_M_UUUUU', 'NLS_M_UUUUU', 'NES_M_UUUUU', 'G3BP_M_UUUUU'), colormaps = c('ivory4', 'skyblue', 'darkseagreen2', 'salmon'), linetypes = c(1, 1, 1, 1), linesize = 2)
plot_CD(normed_table, y_data = c('I_S_UUUUU', 'NLS_S_UUUUU', 'NES_S_UUUUU', 'G3BP_S_UUUUU'), colormaps = c('ivory4', 'skyblue', 'darkseagreen2', 'salmon'), linetypes = c(1, 1, 1, 1), linesize = 2)
####################################################################################################################









