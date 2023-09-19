## CoCLIP Analysis:
## Metagene Analysis
## Written by Soon Yi
## Last edit: 2023-09-18

library(dplyr)
library(tidyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(IRanges)
library(GenomicRanges)
library(biomaRt)
library(data.table)
library(rtracklayer)
library(stringr)
library(GenomicAlignments)
library(ggplot2)

## Custom Functions
####################################################################################################################
# # Calculate densities for each feature. Can playwith window width in the function
# calculate_density = function(feature_windows, feature_starts, peaks_gr) {
#   window_def = 500
#   window_width = 7
#   
#   start_positions = seq(-window_def, window_def - window_width, by=window_width)
#   end_positions = seq(-window_def + window_width, window_def, by=window_width)
# 
#   densities = numeric(length(start_positions))
# 
#   for (i in 1:length(start_positions)) {
#     window <- GRanges(seqnames = seqnames(feature_windows),
#                       ranges = IRanges(start = feature_starts + start_positions[i],
#                                      end = feature_starts + end_positions[i]))
#     densities[i] <- sum(countOverlaps(window, peaks_gr))
#   }
# 
#   densities <- densities / (window_width * length(feature_windows))
#   return(data.frame(midpoint=(start_positions + end_positions) / 2, density=densities))
# }

# Normalize densities by base-correcting.
normalize_density = function(density_data) {
  density_data$density = density_data$density - min(density_data$density)
  return(density_data)
}

# Function to calculate densities for metagene:
metaDensity = function(feature_df, peaks_df, feature_center, window_width = 20, window_definition = 500) {
  feature_df_sense = feature_df %>% filter(strand == '+')
  feature_df_antisense = feature_df %>% filter(strand == '-')
  
  peaks_df_sense = peaks_df %>% filter(strand == '+')
  peaks_df_antisense = peaks_df %>% filter(strand == '-')
  
  if (feature_center == 'start') {
    feature_center_sense = 'start'
    feature_center_antisense = 'end'
  } else if (feature_center == 'end') {
    feature_center_sense = 'end'
    feature_center_antisense = 'start'
  }
  
  ## Make bins for density calcualtion:  
  bins = list()
  bin_counts = floor(window_definition / window_width)
  for (i in -(bin_counts):(bin_counts-1)) {
    start = 0 + (i * window_width)
    end = start + window_width - 1
    bins[[i + bin_counts + 1]] = c(start, end)
  }
  
  ## For Sense Strands:
  chroms = feature_df_sense$seqid
  centers = unlist(feature_df_sense[, feature_center_sense])
  strands = feature_df_sense$strand
  
  peaksGR = GRanges(seqnames = peaks_df_sense$chrom, 
                    ranges = IRanges(start = peaks_df_sense$start, 
                                     end = peaks_df_sense$end), 
                    strand = peaks_df_sense$strand)
  
  counts_sense = matrix(0, nrow = 1, ncol = length(bins))
  
  for (bin_id in 1:length(bins)) {
    bin = bins[bin_id]
    bin_Ranges = matrix(0, nrow = nrow(feature_df_sense), ncol = 4)
    for (center_id in 1:length(centers)) {
      bin_Ranges[center_id, 1] = chroms[center_id]
      bin_Ranges[center_id, 2] = centers[center_id] + bin[[1]][1]
      bin_Ranges[center_id, 3] = centers[center_id] + bin[[1]][2]
      bin_Ranges[center_id, 4] = strands[center_id]
    }
    bin_Ranges = data.frame(bin_Ranges)
    colnames(bin_Ranges) = c('chrom', 'start', 'end', 'strand')
    bin_Ranges_GR = GRanges(bin_Ranges)
    
    counts_sense[1, bin_id] = sum(countOverlaps(peaksGR, bin_Ranges_GR, minoverlap = 10))
  }
  
  
  ## For AntiSense Strands:
  chroms = feature_df_antisense$seqid
  centers = unlist(feature_df_antisense[, feature_center_antisense])
  strands = feature_df_antisense$strand
  
  peaksGR = GRanges(seqnames = peaks_df_antisense$chrom, 
                    ranges = IRanges(start = peaks_df_antisense$start, 
                                     end = peaks_df_antisense$end), 
                    strand = peaks_df_antisense$strand)
  
  counts_antisense = matrix(0, nrow = 1, ncol = length(bins))
  
  bins_antisense = rev(lapply(bins, function(bin) -(bin + 1)))
  
  for (bin_id in 1:length(bins)) {
    bin = bins_antisense[bin_id]
    bin_Ranges = matrix(0, nrow = nrow(feature_df_antisense), ncol = 4)
    for (center_id in 1:length(centers)) {
      bin_Ranges[center_id, 1] = chroms[center_id]
      bin_Ranges[center_id, 2] = centers[center_id] - bin[[1]][1]
      bin_Ranges[center_id, 3] = centers[center_id] - bin[[1]][2]
      bin_Ranges[center_id, 4] = strands[center_id]
    }
    bin_Ranges = data.frame(bin_Ranges)
    colnames(bin_Ranges) = c('chrom', 'start', 'end', 'strand')
    bin_Ranges_GR = GRanges(bin_Ranges)
    
    counts_antisense[1, bin_id] = sum(countOverlaps(peaksGR, bin_Ranges_GR, minoverlap = 10))
    
  }
  
  counts = data.frame(sense = t(data.frame(counts_sense)), 
                      antisense = t(data.frame(counts_antisense)))
  
  counts$both = counts$sense + counts$antisense
  
  counts$density = counts$both / nrow(peaks_df)
  counts$midpoint = unlist(lapply(bins, function(bin) mean(bin)))
  
  counts = counts[, c('midpoint', 'sense', 'antisense', 'both', 'density')]
  
  return(counts)
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
    plot = plot + geom_line(linewidth = 1)
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
####################################################################################################################

## Load Annotation:
####################################################################################################################
# Read in GTF file
gtfFile = '~/Desktop/Genomics/Annotations/Homo_sapiens.GRCh38.110.gtf'
gtf_raw = readGFF(gtfFile)

# Add chr to chromosome names and filter out genes (transcripts only), focus on protein coding 
gtf_raw$seqid = ifelse(substr(gtf_raw$seqid, 1, 3) != "chr", paste0("chr", gtf_raw$seqid), gtf_raw$seqid)
gtf_raw = gtf_raw[gtf_raw$type != "gene", ]

gtf = subset(gtf_raw, transcript_biotype == "protein_coding")
####################################################################################################################

## Define RNA Features:
## Ones that can be directly defined based on annotations: 
## TSS (Transcription Start Sites)
## CDS (Coding Start sites; a.k.a., translation start sites)
## TLS (Translation Stop sites)
## TTS (Transcription Termination Sites)
####################################################################################################################
# Transcription Start Sites (TSS)
TSS = subset(gtf, type == "five_prime_utr")

# Translation Start Sites (aka CDS starts)
CDS = subset(gtf, type == "start_codon")

# Translation Stop Sites (aka 3'UTR starts) 
# Transcription Termination Sites (aka 3'UTR ends)
UTR3 = subset(gtf, type == "three_prime_utr")
####################################################################################################################

## Splice sites need little bit more customization:
## For sanity sake, consider sense and antisense separately:
####################################################################################################################
# Subset exons
exons = subset(gtf, type == "exon")
# exons = exons %>% filter(abs(end-start) >= 500)

# Group by transcript and arrange by start position
exons$exon_number = as.numeric(exons$exon_number)
exons = exons %>% arrange(transcript_id, exon_number) %>% group_by(transcript_id)

# Calculate intron boundaries on sense strand
exons_sense = exons %>% filter(strand == '+')
introns_sense = exons_sense %>%
  mutate(next_start = lead(start), 
         next_end = lead(end)) %>%
  rowwise() %>%
  filter(!is.na(next_start)) %>%
  summarize(intron_start = end + 1,
            intron_end = next_start - 1,
            seqid = seqid[1],
            strand = strand[1],
            intron_number = exon_number) %>%
  ungroup()

# Calculate intron boundaries on anti-sense strand
exons_antisense = exons %>% filter(strand == '-')
introns_antisense = exons_antisense %>%
  mutate(next_start = lead(start), 
         next_end = lead(end)) %>%
  rowwise() %>%
  filter(!is.na(next_start)) %>%
  summarize(intron_start = next_end + 1,
            intron_end = start - 1,
            seqid = seqid[1],
            strand = strand[1],
            intron_number = exon_number) %>%
  ungroup()

INTRONs = rbind(introns_sense, introns_antisense)
INTRONs = INTRONs %>% arrange(transcript_id, intron_number)
colnames(INTRONs) = c('transcript_id', 'start', 'end', 'seqid', 'strand', 'intron_number')
####################################################################################################################

# Read in the peaks and filters to be used:
####################################################################################################################
peakFile = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/Combined_peakCoverage_groomed_annotated.txt'
peaks_org = read.delim(peakFile, header=TRUE, sep="\t")

Nuc_F_M = c('Nuc_F_M_1', 'Nuc_F_M_2', 'Nuc_F_M_3')
Cyto_F_M = c('Cyto_F_M_1', 'Cyto_F_M_2', 'Cyto_F_M_3')
Nuc_F_S = c('Nuc_F_S_1', 'Nuc_F_S_2', 'Nuc_F_S_3')
Cyto_F_S = c('Cyto_F_S_1', 'Cyto_F_S_2', 'Cyto_F_S_3')
NLS_I_M = c('NLS_I_M_1', 'NLS_I_M_2')
NES_I_M = c('NES_I_M_1', 'NES_I_M_2')
G3BP_I_M = c('G3BP_I_M_1', 'G3BP_I_M_2', 'G3BP_I_M_3', 'G3BP_I_M_4')
NLS_I_S = c('NLS_I_S_1', 'NLS_I_S_2')
NES_I_S = c('NES_I_S_1', 'NES_I_S_2')
G3BP_I_S = c('G3BP_I_S_1', 'G3BP_I_S_2', 'G3BP_I_S_3', 'G3BP_I_S_4', 'G3BP_I_S_5')
NLS_E_M = c('NLS_E_M_1', 'NLS_E_M_2', 'NLS_E_M_3', 'NLS_E_M_4')
NES_E_M = c('NES_E_M_1', 'NES_E_M_2', 'NES_E_M_3', 'NES_E_M_4')
G3BP_E_M = c('G3BP_E_M_1', 'G3BP_E_M_2', 'G3BP_E_M_3', 'G3BP_E_M_4', 'G3BP_E_M_5', 'G3BP_E_M_6')
NLS_E_S = c('NLS_E_S_1', 'NLS_E_S_2', 'NLS_E_S_3', 'NLS_E_S_4')
NES_E_S = c('NES_E_S_1', 'NES_E_S_2', 'NES_E_S_3', 'NES_E_S_4')
G3BP_E_S = c('G3BP_E_S_1', 'G3BP_E_S_2', 'G3BP_E_S_3', 'G3BP_E_S_4', 'G3BP_E_S_5', 'G3BP_E_S_6', 'G3BP_E_S_7')

BC_columns = c("Nuc_F_M_BC", "Nuc_F_S_BC", "Cyto_F_M_BC", "Cyto_F_S_BC",
               "NLS_I_M_BC", "NLS_I_S_BC", "NES_I_M_BC", "NES_I_S_BC", "G3BP_I_M_BC", "G3BP_I_S_BC",
               "NLS_E_M_BC", "NLS_E_S_BC", "NES_E_M_BC", "NES_E_S_BC", "G3BP_E_M_BC", "G3BP_E_S_BC",
               "TOTAL_BC")

peaks = peaks_org
####################################################################################################################

## Filter Criteria:
####################################################################################################################
BC_Threshold_F = 2
BC_Threshold_I = 4
BC_Threshold_E = 1
BC_Threshold_E_SG = 2

rowSum_Multiplier_F = 4
rowSum_Multiplier_I = 1
rowSum_Multiplier_E = 1
####################################################################################################################

## Filter peaks as necessary:
####################################################################################################################
## Filter for Input Mock
I_M_peaks = peaks_org %>% filter(NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC >= BC_Threshold_I)
I_M_peaks = I_M_peaks %>% filter(rowSums(I_M_peaks[, c(NLS_I_M, NES_I_M, G3BP_I_M)]) >= median(rowSums(I_M_peaks[, c(NLS_I_M, NES_I_M, G3BP_I_M)]) * rowSum_Multiplier_I))
I_M_peaks = I_M_peaks %>% mutate(selectRowSum = rowSums(I_M_peaks[, c(NLS_I_M, NES_I_M, G3BP_I_M)])) %>% uncount(selectRowSum)

## Filter for Input Stress
I_S_peaks = peaks_org %>% filter(NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC >= BC_Threshold_I)
I_S_peaks = I_S_peaks %>% filter(rowSums(I_S_peaks[, c(NLS_I_S, NES_I_S, G3BP_I_S)]) >= median(rowSums(I_S_peaks[, c(NLS_I_S, NES_I_S, G3BP_I_S)]) * rowSum_Multiplier_I))
I_S_peaks = I_S_peaks %>% mutate(selectRowSum = rowSums(I_S_peaks[, c(NLS_I_S, NES_I_S, G3BP_I_S)])) %>% uncount(selectRowSum)

## Filter for NLS Mock
NLS_M_peaks = peaks_org %>% filter(NLS_E_M_BC >= BC_Threshold_E)
NLS_M_peaks = NLS_M_peaks %>% filter(rowSums(NLS_M_peaks[, c(NLS_E_M)]) >= median(rowSums(NLS_M_peaks[, c(NLS_E_M)])) * rowSum_Multiplier_E)
NLS_M_peaks = NLS_M_peaks %>% mutate(selectRowSum = rowSums(NLS_M_peaks[, c(NLS_E_M)])) %>% uncount(selectRowSum)

## Filter for NLS Stress
NLS_S_peaks = peaks_org %>% filter(NLS_E_S_BC >= BC_Threshold_E)
NLS_S_peaks = NLS_S_peaks %>% filter(rowSums(NLS_S_peaks[, c(NLS_E_S)]) >= median(rowSums(NLS_S_peaks[, c(NLS_E_S)])) * rowSum_Multiplier_E)
NLS_S_peaks = NLS_S_peaks %>% mutate(selectRowSum = rowSums(NLS_S_peaks[, c(NLS_E_S)])) %>% uncount(selectRowSum)

## Filter for NES Mock
NES_M_peaks = peaks_org %>% filter(NES_E_M_BC >= BC_Threshold_E)
NES_M_peaks = NES_M_peaks %>% filter(rowSums(NES_M_peaks[, c(NES_E_M)]) >= median(rowSums(NES_M_peaks[, c(NES_E_M)])) * rowSum_Multiplier_E)
NES_M_peaks = NES_M_peaks %>% mutate(selectRowSum = rowSums(NES_M_peaks[, c(NES_E_M)])) %>% uncount(selectRowSum)

## Filter for NES Stress
NES_S_peaks = peaks_org %>% filter(NES_E_S_BC >= BC_Threshold_E)
NES_S_peaks = NES_S_peaks %>% filter(rowSums(NES_S_peaks[, c(NES_E_S)]) >= median(rowSums(NES_S_peaks[, c(NES_E_S)])) * rowSum_Multiplier_E)
NES_S_peaks = NES_S_peaks %>% mutate(selectRowSum = rowSums(NES_S_peaks[, c(NES_E_S)])) %>% uncount(selectRowSum)

## Filter for G3BP Mock
G3BP_M_peaks = peaks_org %>% filter(G3BP_E_M_BC >= BC_Threshold_E_SG)
G3BP_M_peaks = G3BP_M_peaks %>% filter(rowSums(G3BP_M_peaks[, c(G3BP_E_M)]) >= median(rowSums(G3BP_M_peaks[, c(G3BP_E_M)])) * rowSum_Multiplier_E)
G3BP_M_peaks = G3BP_M_peaks %>% mutate(selectRowSum = rowSums(G3BP_M_peaks[, c(G3BP_E_M)])) %>% uncount(selectRowSum)

## Filter for G3BP Stress
G3BP_S_peaks = peaks_org %>% filter(G3BP_E_S_BC >= BC_Threshold_E_SG)
G3BP_S_peaks = G3BP_S_peaks %>% filter(rowSums(G3BP_S_peaks[, c(G3BP_E_S)]) >= median(rowSums(G3BP_S_peaks[, c(G3BP_E_S)])) * rowSum_Multiplier_E)
G3BP_S_peaks = G3BP_S_peaks %>% mutate(selectRowSum = rowSums(G3BP_S_peaks[, c(G3BP_E_S)])) %>% uncount(selectRowSum)

## Filter for Nuclear Fraction Mock
Nuc_M_peaks = peaks_org %>% filter(Nuc_F_M_BC >= BC_Threshold_F)
Nuc_M_peaks = Nuc_M_peaks %>% filter(rowSums(Nuc_M_peaks[, c(Nuc_F_M)]) >= median(rowSums(Nuc_M_peaks[, c(Nuc_F_M)])) * rowSum_Multiplier_F)
Nuc_M_peaks = Nuc_M_peaks %>% mutate(selectRowSum = rowSums(Nuc_M_peaks[, c(Nuc_F_M)])) %>% uncount(selectRowSum)

## Filter for Nuclear Fraction Stress
Nuc_S_peaks = peaks_org %>% filter(Nuc_F_S_BC >= BC_Threshold_F)
Nuc_S_peaks = Nuc_S_peaks %>% filter(rowSums(Nuc_S_peaks[, c(Nuc_F_S)]) >= median(rowSums(Nuc_S_peaks[, c(Nuc_F_S)])) * rowSum_Multiplier_F)
Nuc_S_peaks = Nuc_S_peaks %>% mutate(selectRowSum = rowSums(Nuc_S_peaks[, c(Nuc_F_S)])) %>% uncount(selectRowSum)

## Filter for Cytoplasm Fraction Mock
Cyto_M_peaks = peaks_org %>% filter(Cyto_F_M_BC >= BC_Threshold_F)
Cyto_M_peaks = Cyto_M_peaks %>% filter(rowSums(Cyto_M_peaks[, c(Cyto_F_M)]) >= median(rowSums(Cyto_M_peaks[, c(Cyto_F_M)])) * rowSum_Multiplier_F)
Cyto_M_peaks = Cyto_M_peaks %>% mutate(selectRowSum = rowSums(Cyto_M_peaks[, c(Cyto_F_M)])) %>% uncount(selectRowSum)

## Filter for Cytoplasm Fraction Stress
Cyto_S_peaks = peaks_org %>% filter(Cyto_F_S_BC >= BC_Threshold_F)
Cyto_S_peaks = Cyto_S_peaks %>% filter(rowSums(Cyto_S_peaks[, c(Cyto_F_S)]) >= median(rowSums(Cyto_S_peaks[, c(Cyto_F_S)])) * rowSum_Multiplier_F)
Cyto_S_peaks = Cyto_S_peaks %>% mutate(selectRowSum = rowSums(Cyto_S_peaks[, c(Cyto_F_S)])) %>% uncount(selectRowSum)
####################################################################################################################

## Transcription Start Site (TSS):
#########################################################################################################################
I_M_metaDensity = metaDensity(feature_df = TSS, peaks_df = I_M_peaks, feature_center = 'start')
NLS_M_metaDensity = metaDensity(feature_df = TSS, peaks_df = NLS_M_peaks, feature_center = 'start')
NES_M_metaDensity = metaDensity(feature_df = TSS, peaks_df = NES_M_peaks, feature_center = 'start')
G3BP_M_metaDensity = metaDensity(feature_df = TSS, peaks_df = G3BP_M_peaks, feature_center = 'start')

I_S_metaDensity = metaDensity(feature_df = TSS, peaks_df = I_S_peaks, feature_center = 'start')
NLS_S_metaDensity = metaDensity(feature_df = TSS, peaks_df = NLS_S_peaks, feature_center = 'start')
NES_S_metaDensity = metaDensity(feature_df = TSS, peaks_df = NES_S_peaks, feature_center = 'start')
G3BP_S_metaDensity = metaDensity(feature_df = TSS, peaks_df = G3BP_S_peaks, feature_center = 'start')

TSS_densities = data.frame(position = NLS_M_metaDensity$midpoint,
                           I_M = I_M_metaDensity$density,
                           NLS_M = NLS_M_metaDensity$density,
                           NES_M = NES_M_metaDensity$density,
                           G3BP_M = G3BP_M_metaDensity$density,
                           I_S = I_S_metaDensity$density,
                           NLS_S = NLS_S_metaDensity$density,
                           NES_S = NES_S_metaDensity$density,
                           G3BP_S = G3BP_S_metaDensity$density)
#########################################################################################################################

## Translation Start (Coding Start) Site (CDS):
#########################################################################################################################
I_M_metaDensity = metaDensity(feature_df = CDS, peaks_df = I_M_peaks, feature_center = 'start')
NLS_M_metaDensity = metaDensity(feature_df = CDS, peaks_df = NLS_M_peaks, feature_center = 'start')
NES_M_metaDensity = metaDensity(feature_df = CDS, peaks_df = NES_M_peaks, feature_center = 'start')
G3BP_M_metaDensity = metaDensity(feature_df = CDS, peaks_df = G3BP_M_peaks, feature_center = 'start')

I_S_metaDensity = metaDensity(feature_df = CDS, peaks_df = I_S_peaks, feature_center = 'start')
NLS_S_metaDensity = metaDensity(feature_df = CDS, peaks_df = NLS_S_peaks, feature_center = 'start')
NES_S_metaDensity = metaDensity(feature_df = CDS, peaks_df = NES_S_peaks, feature_center = 'start')
G3BP_S_metaDensity = metaDensity(feature_df = CDS, peaks_df = G3BP_S_peaks, feature_center = 'start')

CDS_densities = data.frame(position = NLS_M_metaDensity$midpoint,
                           I_M = I_M_metaDensity$density,
                           NLS_M = NLS_M_metaDensity$density,
                           NES_M = NES_M_metaDensity$density,
                           G3BP_M = G3BP_M_metaDensity$density,
                           I_S = I_S_metaDensity$density,
                           NLS_S = NLS_S_metaDensity$density,
                           NES_S = NES_S_metaDensity$density,
                           G3BP_S = G3BP_S_metaDensity$density)
#########################################################################################################################

## 5' Splice Site (SS5):
#########################################################################################################################
I_M_metaDensity = metaDensity(feature_df = INTRONs, peaks_df = I_M_peaks, feature_center = 'start')
NLS_M_metaDensity = metaDensity(feature_df = INTRONs, peaks_df = NLS_M_peaks, feature_center = 'start')
NES_M_metaDensity = metaDensity(feature_df = INTRONs, peaks_df = NES_M_peaks, feature_center = 'start')
G3BP_M_metaDensity = metaDensity(feature_df = INTRONs, peaks_df = G3BP_M_peaks, feature_center = 'start')

I_S_metaDensity = metaDensity(feature_df = INTRONs, peaks_df = I_S_peaks, feature_center = 'start')
NLS_S_metaDensity = metaDensity(feature_df = INTRONs, peaks_df = NLS_S_peaks, feature_center = 'start')
NES_S_metaDensity = metaDensity(feature_df = INTRONs, peaks_df = NES_S_peaks, feature_center = 'start')
G3BP_S_metaDensity = metaDensity(feature_df = INTRONs, peaks_df = G3BP_S_peaks, feature_center = 'start')

SS5_densities = data.frame(position = NLS_M_metaDensity$midpoint,
                           I_M = I_M_metaDensity$density,
                           NLS_M = NLS_M_metaDensity$density,
                           NES_M = NES_M_metaDensity$density,
                           G3BP_M = G3BP_M_metaDensity$density,
                           I_S = I_S_metaDensity$density,
                           NLS_S = NLS_S_metaDensity$density,
                           NES_S = NES_S_metaDensity$density,
                           G3BP_S = G3BP_S_metaDensity$density)
#########################################################################################################################

## 3' Splice Site (SS3):
#########################################################################################################################
INTRONs_noLast = INTRONs %>% group_by(transcript_id) %>% filter(intron_number != max(intron_number)) %>% ungroup()

I_M_metaDensity = metaDensity(feature_df = INTRONs_noLast, peaks_df = I_M_peaks, feature_center = 'end')
NLS_M_metaDensity = metaDensity(feature_df = INTRONs_noLast, peaks_df = NLS_M_peaks, feature_center = 'end')
NES_M_metaDensity = metaDensity(feature_df = INTRONs_noLast, peaks_df = NES_M_peaks, feature_center = 'end')
G3BP_M_metaDensity = metaDensity(feature_df = INTRONs_noLast, peaks_df = G3BP_M_peaks, feature_center = 'end')

I_S_metaDensity = metaDensity(feature_df = INTRONs_noLast, peaks_df = I_S_peaks, feature_center = 'end')
NLS_S_metaDensity = metaDensity(feature_df = INTRONs_noLast, peaks_df = NLS_S_peaks, feature_center = 'end')
NES_S_metaDensity = metaDensity(feature_df = INTRONs_noLast, peaks_df = NES_S_peaks, feature_center = 'end')
G3BP_S_metaDensity = metaDensity(feature_df = INTRONs_noLast, peaks_df = G3BP_S_peaks, feature_center = 'end')

SS3_densities = data.frame(position = NLS_M_metaDensity$midpoint,
                           I_M = I_M_metaDensity$density,
                           NLS_M = NLS_M_metaDensity$density,
                           NES_M = NES_M_metaDensity$density,
                           G3BP_M = G3BP_M_metaDensity$density,
                           I_S = I_S_metaDensity$density,
                           NLS_S = NLS_S_metaDensity$density,
                           NES_S = NES_S_metaDensity$density,
                           G3BP_S = G3BP_S_metaDensity$density)
#########################################################################################################################

## Translation Stop Site (TLS):
#########################################################################################################################
I_M_metaDensity = metaDensity(feature_df = UTR3, peaks_df = I_M_peaks, feature_center = 'start')
NLS_M_metaDensity = metaDensity(feature_df = UTR3, peaks_df = NLS_M_peaks, feature_center = 'start')
NES_M_metaDensity = metaDensity(feature_df = UTR3, peaks_df = NES_M_peaks, feature_center = 'start')
G3BP_M_metaDensity = metaDensity(feature_df = UTR3, peaks_df = G3BP_M_peaks, feature_center = 'start')

I_S_metaDensity = metaDensity(feature_df = UTR3, peaks_df = I_S_peaks, feature_center = 'start')
NLS_S_metaDensity = metaDensity(feature_df = UTR3, peaks_df = NLS_S_peaks, feature_center = 'start')
NES_S_metaDensity = metaDensity(feature_df = UTR3, peaks_df = NES_S_peaks, feature_center = 'start')
G3BP_S_metaDensity = metaDensity(feature_df = UTR3, peaks_df = G3BP_S_peaks, feature_center = 'start')

TLS_densities = data.frame(position = NLS_M_metaDensity$midpoint,
                           I_M = I_M_metaDensity$density,
                           NLS_M = NLS_M_metaDensity$density,
                           NES_M = NES_M_metaDensity$density,
                           G3BP_M = G3BP_M_metaDensity$density,
                           I_S = I_S_metaDensity$density,
                           NLS_S = NLS_S_metaDensity$density,
                           NES_S = NES_S_metaDensity$density,
                           G3BP_S = G3BP_S_metaDensity$density)
#########################################################################################################################

## Transcription Termination Site (TTS):
#########################################################################################################################
I_M_metaDensity = metaDensity(feature_df = UTR3, peaks_df = I_M_peaks, feature_center = 'end')
NLS_M_metaDensity = metaDensity(feature_df = UTR3, peaks_df = NLS_M_peaks, feature_center = 'end')
NES_M_metaDensity = metaDensity(feature_df = UTR3, peaks_df = NES_M_peaks, feature_center = 'end')
G3BP_M_metaDensity = metaDensity(feature_df = UTR3, peaks_df = G3BP_M_peaks, feature_center = 'end')

I_S_metaDensity = metaDensity(feature_df = UTR3, peaks_df = I_S_peaks, feature_center = 'end')
NLS_S_metaDensity = metaDensity(feature_df = UTR3, peaks_df = NLS_S_peaks, feature_center = 'end')
NES_S_metaDensity = metaDensity(feature_df = UTR3, peaks_df = NES_S_peaks, feature_center = 'end')
G3BP_S_metaDensity = metaDensity(feature_df = UTR3, peaks_df = G3BP_S_peaks, feature_center = 'end')

TTS_densities = data.frame(position = NLS_M_metaDensity$midpoint,
                           I_M = I_M_metaDensity$density,
                           NLS_M = NLS_M_metaDensity$density,
                           NES_M = NES_M_metaDensity$density,
                           G3BP_M = G3BP_M_metaDensity$density,
                           I_S = I_S_metaDensity$density,
                           NLS_S = NLS_S_metaDensity$density,
                           NES_S = NES_S_metaDensity$density,
                           G3BP_S = G3BP_S_metaDensity$density)

#########################################################################################################################

## Samples
#########################################################################################################################

# TSS_densities_amplified = TSS_densities
# CDS_densities_amplified = CDS_densities
# SS5_densities_amplified = SS5_densities
# SS3_densities_amplified = SS3_densities
# TLS_densities_amplified = TLS_densities
# TTS_densities_amplified = TTS_densities

# TSS_densities = TSS_densities_amplified
# CDS_densities = CDS_densities_amplified
# SS5_densities = SS5_densities_amplified
# SS3_densities = SS3_densities_amplified
# TLS_densities = TLS_densities_amplified
# TTS_densities = TTS_densities_amplified

# TSS_densities_peakLevel = TSS_densities
# CDS_densities_peakLevel = CDS_densities
# SS5_densities_peakLevel = SS5_densities
# SS3_densities_peakLevel = SS3_densities
# TLS_densities_peakLevel = TLS_densities
# TTS_densities_peakLevel = TTS_densities

TSS_densities = TSS_densities_peakLevel
CDS_densities = CDS_densities_peakLevel
SS5_densities = SS5_densities_peakLevel
SS3_densities = SS3_densities_peakLevel
TLS_densities = TLS_densities_peakLevel
TTS_densities = TTS_densities_peakLevel

#########################################################################################################################

## Inputs 
#########################################################################################################################
y_lim = c(0, 0.1)
smoothing = 0.2

plot_Density(density_data = TSS_densities,
             columns_list = c('I_S', 'I_M'),
             custom_colors = c('grey3', 'grey'),
             densityType = 'feature_metagene',
             featureName = "Transcription Start Sites",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Input Mocks and Stress')

plot_Density(density_data = CDS_densities,
             columns_list = c('I_S', 'I_M'),
             custom_colors = c('grey3', 'grey'),
             densityType = 'feature_metagene',
             featureName = "Translation Start Sites",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Input Mocks and Stress')

plot_Density(density_data = SS5_densities,
             columns_list = c('I_S', 'I_M'),
             custom_colors = c('grey3', 'grey'),
             densityType = 'feature_metagene',
             featureName = "5' Splice Site",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Input Mocks and Stress')

plot_Density(density_data = SS3_densities,
             columns_list = c('I_S', 'I_M'),
             custom_colors = c('grey3', 'grey'),
             densityType = 'feature_metagene',
             featureName = "3' Splice Site",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Input Mocks and Stress')

plot_Density(density_data = TLS_densities,
             columns_list = c('I_S', 'I_M'),
             custom_colors = c('grey3', 'grey'),
             densityType = 'feature_metagene',
             featureName = "Translation Stop Sites",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Input Mocks and Stress')

plot_Density(density_data = TTS_densities,
             columns_list = c('I_S', 'I_M'),
             custom_colors = c('grey3', 'grey'),
             densityType = 'feature_metagene',
             featureName = "Transcription Termination Sites",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Input Mocks and Stress')
#########################################################################################################################

## Mocks: Use y_lim = c(0, 0.05) and smoothing = 0.2 for figures
#########################################################################################################################
y_lim = c(0, 0.05)
smoothing = 0.2

plot_Density(density_data = TSS_densities,
             columns_list = c('G3BP_M', 'NES_M', 'NLS_M', 'I_M'),
             custom_colors = c('salmon', 'darkseagreen2', 'skyblue','grey'),
             densityType = 'feature_metagene',
             featureName = "Transcription Start Sites",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Mocks')

plot_Density(density_data = CDS_densities,
             columns_list = c('G3BP_M', 'NES_M', 'NLS_M', 'I_M'),
             custom_colors = c('salmon', 'darkseagreen2', 'skyblue','grey'),
             densityType = 'feature_metagene',
             featureName = "Translation Start Sites",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Mocks')

plot_Density(density_data = SS5_densities,
             columns_list = c('G3BP_M', 'NES_M', 'NLS_M', 'I_M'),
             custom_colors = c('salmon', 'darkseagreen2', 'skyblue','grey'),
             densityType = 'feature_metagene',
             featureName = "5' Splice Site",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Mocks')

plot_Density(density_data = SS3_densities,
             columns_list = c('G3BP_M', 'NES_M', 'NLS_M', 'I_M'),
             custom_colors = c('salmon', 'darkseagreen2', 'skyblue','grey'),
             densityType = 'feature_metagene',
             featureName = "3' Splice Site",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Mocks')

plot_Density(density_data = TLS_densities,
             columns_list = c('G3BP_M', 'NES_M', 'NLS_M', 'I_M'),
             custom_colors = c('salmon', 'darkseagreen2', 'skyblue','grey'),
             densityType = 'feature_metagene',
             featureName = "Translation Stop Sites",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Mocks')

plot_Density(density_data = TTS_densities,
             columns_list = c('G3BP_M', 'NES_M', 'NLS_M', 'I_M'),
             custom_colors = c('salmon', 'darkseagreen2', 'skyblue','grey'),
             densityType = 'feature_metagene',
             featureName = "Transcription Termination Sites",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Mocks')
#########################################################################################################################

## Arsenites
#########################################################################################################################
y_lim = c(0, 0.05)
smoothing = 0.2

plot_Density(density_data = TSS_densities,
             columns_list = c('G3BP_S', 'NES_S', 'NLS_S', 'I_S'),
             custom_colors = c('salmon', 'darkseagreen2', 'skyblue', 'grey'),
             densityType = 'feature_metagene',
             featureName = "Transcription Start Sites",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Arsenite')

plot_Density(density_data = CDS_densities,
             columns_list = c('G3BP_S', 'NES_S', 'NLS_S', 'I_S'),
             custom_colors = c('salmon', 'darkseagreen2', 'skyblue', 'grey'),
             densityType = 'feature_metagene',
             featureName = "Translation Start Sites",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Arsenite')

plot_Density(density_data = SS5_densities,
             columns_list = c('G3BP_S', 'NES_S', 'NLS_S', 'I_S'),
             custom_colors = c('salmon', 'darkseagreen2', 'skyblue', 'grey'),
             densityType = 'feature_metagene',
             featureName = "5' Splice Site",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Arsenite')

plot_Density(density_data = SS3_densities,
             columns_list = c('G3BP_S', 'NES_S', 'NLS_S', 'I_S'),
             custom_colors = c('salmon', 'darkseagreen2', 'skyblue', 'grey'),
             densityType = 'feature_metagene',
             featureName = "3' Splice Site",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Arsenite')

plot_Density(density_data = TLS_densities,
             columns_list = c('G3BP_S', 'NES_S', 'NLS_S', 'I_S'),
             custom_colors = c('salmon', 'darkseagreen2', 'skyblue', 'grey'),
             densityType = 'feature_metagene',
             featureName = "Translation Stop Sites",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Arsenite')

plot_Density(density_data = TTS_densities,
             columns_list = c('G3BP_S', 'NES_S', 'NLS_S', 'I_S'),
             custom_colors = c('salmon', 'darkseagreen2', 'skyblue', 'grey'),
             densityType = 'feature_metagene',
             featureName = "Transcription Termination Sites",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Arsenite')
#########################################################################################################################

## Fold Change
#########################################################################################################################
TSS_densities_FC = data.frame(position = TSS_densities$position,
                              Inp = log2(TSS_densities$I_S / TSS_densities$I_M),
                              NLS = log2(TSS_densities$NLS_S / TSS_densities$NLS_M),
                              NES = log2(TSS_densities$NES_S / TSS_densities$NES_M),
                              G3BP = log2(TSS_densities$G3BP_S / TSS_densities$G3BP_M)) 

CDS_densities_FC = data.frame(position = CDS_densities$position,
                              Inp = log2(CDS_densities$I_S / CDS_densities$I_M),
                              NLS = log2(CDS_densities$NLS_S / CDS_densities$NLS_M),
                              NES = log2(CDS_densities$NES_S / CDS_densities$NES_M),
                              G3BP = log2(CDS_densities$G3BP_S / CDS_densities$G3BP_M)) 

SS5_densities_FC = data.frame(position = SS5_densities$position,
                              Inp = log2(SS5_densities$I_S / SS5_densities$I_M),
                              NLS = log2(SS5_densities$NLS_S / SS5_densities$NLS_M),
                              NES = log2(SS5_densities$NES_S / SS5_densities$NES_M),
                              G3BP = log2(SS5_densities$G3BP_S / SS5_densities$G3BP_M)) 

SS3_densities_FC = data.frame(position = SS3_densities$position,
                              Inp = log2(SS3_densities$I_S / SS3_densities$I_M),
                              NLS = log2(SS3_densities$NLS_S / SS3_densities$NLS_M),
                              NES = log2(SS3_densities$NES_S / SS3_densities$NES_M),
                              G3BP = log2(SS3_densities$G3BP_S / SS3_densities$G3BP_M)) 

TLS_densities_FC = data.frame(position = TLS_densities$position,
                              Inp = log2(TLS_densities$I_S / TLS_densities$I_M),
                              NLS = log2(TLS_densities$NLS_S / TLS_densities$NLS_M),
                              NES = log2(TLS_densities$NES_S / TLS_densities$NES_M),
                              G3BP = log2(TLS_densities$G3BP_S / TLS_densities$G3BP_M)) 

TTS_densities_FC = data.frame(position = TTS_densities$position,
                              Inp = log2(TTS_densities$I_S / TTS_densities$I_M),
                              NLS = log2(TTS_densities$NLS_S / TTS_densities$NLS_M),
                              NES = log2(TTS_densities$NES_S / TTS_densities$NES_M),
                              G3BP = log2(TTS_densities$G3BP_S / TTS_densities$G3BP_M)) 

y_lim = c(-8, 8)
smoothing = 0.2

plot_Density(density_data = TSS_densities_FC,
             columns_list = c('Inp', 'NLS', 'NES', 'G3BP'),
             custom_colors = c('grey', 'skyblue', 'darkseagreen2', 'salmon'),
             densityType = 'feature_metagene',
             featureName = "Transcription Start Sites",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Arsenite/Mock')

plot_Density(density_data = CDS_densities_FC,
             columns_list = c('Inp', 'NLS', 'NES', 'G3BP'),
             custom_colors = c('grey', 'skyblue', 'darkseagreen2', 'salmon'),
             densityType = 'feature_metagene',
             featureName = "Translation Start Sites",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Arsenite/Mock')

plot_Density(density_data = SS5_densities_FC,
             columns_list = c('Inp', 'NLS', 'NES', 'G3BP'),
             custom_colors = c('grey', 'skyblue', 'darkseagreen2', 'salmon'),
             densityType = 'feature_metagene',
             featureName = "5' Splice Site",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Arsenite/Mock')

plot_Density(density_data = SS3_densities_FC,
             columns_list = c('Inp', 'NLS', 'NES', 'G3BP'),
             custom_colors = c('grey', 'skyblue', 'darkseagreen2', 'salmon'),
             densityType = 'feature_metagene',
             featureName = "3' Splice Site",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Arsenite/Mock')

plot_Density(density_data = TLS_densities_FC,
             columns_list = c('Inp', 'NLS', 'NES', 'G3BP'),
             custom_colors = c('grey', 'skyblue', 'darkseagreen2', 'salmon'),
             densityType = 'feature_metagene',
             featureName = "Translation Stop Sites",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Arsenite/Mock')

plot_Density(density_data = TTS_densities_FC,
             columns_list = c('Inp', 'NLS', 'NES', 'G3BP'),
             custom_colors = c('grey', 'skyblue', 'darkseagreen2', 'salmon'),
             densityType = 'feature_metagene',
             featureName = "Transcription Termination Sites",
             yaxis_lims = y_lim,
             smoothing = smoothing,
             sampleName = 'Arsenite/Mock')
#########################################################################################################################





#########################################################################################################################
#######################################################DEPRECATED########################################################
#########################################################################################################################

# Convert peaks to GRanges objects
####################################################################################################################
# I_M_peaks_gr = GRanges(seqnames=I_M_peaks$chrom, ranges=IRanges(start=I_M_peaks$start, end=I_M_peaks$end), strand = I_M_peaks$strand)
# NLS_M_peaks_gr = GRanges(seqnames=NLS_M_peaks$chrom, ranges=IRanges(start=NLS_M_peaks$start, end=NLS_M_peaks$end), strand = NLS_M_peaks$strand)
# NES_M_peaks_gr = GRanges(seqnames=NES_M_peaks$chrom, ranges=IRanges(start=NES_M_peaks$start, end=NES_M_peaks$end), strand = NES_M_peaks$strand)
# G3BP_M_peaks_gr = GRanges(seqnames=G3BP_M_peaks$chrom, ranges=IRanges(start=G3BP_M_peaks$start, end=G3BP_M_peaks$end), strand = G3BP_M_peaks$strand)
# 
# I_S_peaks_gr = GRanges(seqnames=I_S_peaks$chrom, ranges=IRanges(start=I_S_peaks$start, end=I_S_peaks$end), strand = I_S_peaks$strand)
# NLS_S_peaks_gr = GRanges(seqnames=NLS_S_peaks$chrom, ranges=IRanges(start=NLS_S_peaks$start, end=NLS_S_peaks$end), strand = NLS_S_peaks$strand)
# NES_S_peaks_gr = GRanges(seqnames=NES_S_peaks$chrom, ranges=IRanges(start=NES_S_peaks$start, end=NES_S_peaks$end), strand = NES_S_peaks$strand)
# G3BP_S_peaks_gr = GRanges(seqnames=G3BP_S_peaks$chrom, ranges=IRanges(start=G3BP_S_peaks$start, end=G3BP_S_peaks$end), strand = G3BP_S_peaks$strand)

#########################################################################################################################

## Calculate Density:
####################################################################################################################
# samples = c('I_M_peaks_gr', 'NLS_M_peaks_gr', 'NES_M_peaks_gr', 'G3BP_M_peaks_gr', 
#             'I_S_peaks_gr', 'NLS_S_peaks_gr', 'NES_S_peaks_gr', 'G3BP_S_peaks_gr')
# 
# density_lists = list()
# density_lists_sense = list()
# density_lists_antisense = list()
# 
# for (sample in samples) {
#   peaks_gr = get(sample)
#   tss_density_data = normalize_density(calculate_density(tss_windows, tss_starts$start, peaks_gr))
#   cds_density_data = normalize_density(calculate_density(cds_windows, cds_starts$start, peaks_gr))
#   five_prime_splice_density_data = normalize_density(calculate_density(five_prime_splice_windows, five_prime_splice$start, peaks_gr))
#   three_prime_splice_density_data = normalize_density(calculate_density(three_prime_splice_windows, three_prime_splice$end, peaks_gr))
#   translation_stop_density_data = normalize_density(calculate_density(translation_stop_windows, utr3_starts$start, peaks_gr))
#   tts_density_data = normalize_density(calculate_density(tts_windows, tts_starts$end, peaks_gr))
#   
#   density_lists[[str_replace(sample, '_peaks_gr', '')]] = list(TSS = tss_density_data$density,
#                                                                      CDS = cds_density_data$density,
#                                                                      SS5 = five_prime_splice_density_data$density,
#                                                                      SS3 = three_prime_splice_density_data$density,
#                                                                      TLS = translation_stop_density_data$density,
#                                                                      TTS = tts_density_data$density)
#   density_lists[['midpoint']] = list(tss_density_data$midpoint)
#   
#   strandSelection = '+'
#   tss_density_data = normalize_density(calculate_density(tss_windows[strand(tss_windows) == strandSelection], (tss_starts %>% filter(strand == strandSelection))$start, peaks_gr))
#   cds_density_data = normalize_density(calculate_density(cds_windows[strand(cds_windows) == strandSelection], (cds_starts %>% filter(strand == strandSelection))$start, peaks_gr))
#   five_prime_splice_density_data = normalize_density(calculate_density(five_prime_splice_windows[strand(five_prime_splice_windows) == strandSelection], (five_prime_splice %>% filter(strand == strandSelection))$start, peaks_gr))
#   three_prime_splice_density_data = normalize_density(calculate_density(three_prime_splice_windows[strand(three_prime_splice_windows) == strandSelection], (three_prime_splice %>% filter(strand == strandSelection))$end, peaks_gr))
#   translation_stop_density_data = normalize_density(calculate_density(translation_stop_windows[strand(translation_stop_windows) == strandSelection], (utr3_starts %>% filter(strand == strandSelection))$start, peaks_gr))
#   tts_density_data = normalize_density(calculate_density(tts_windows[strand(tts_windows) == strandSelection], (tts_starts %>% filter(strand == strandSelection))$end, peaks_gr))
#   
#   density_lists_sense[[str_replace(sample, '_peaks_gr', '')]] = list(TSS = tss_density_data$density,
#                                                                CDS = cds_density_data$density,
#                                                                SS5 = five_prime_splice_density_data$density,
#                                                                SS3 = three_prime_splice_density_data$density,
#                                                                TLS = translation_stop_density_data$density,
#                                                                TTS = tts_density_data$density)
# density_lists_sense[['midpoint']] = list(tss_density_data$midpoint)
#   
#   strandSelection = '-'
#   tss_density_data = normalize_density(calculate_density(tss_windows[strand(tss_windows) == strandSelection], (tss_starts %>% filter(strand == strandSelection))$end, peaks_gr))
#   cds_density_data = normalize_density(calculate_density(cds_windows[strand(cds_windows) == strandSelection], (cds_starts %>% filter(strand == strandSelection))$end, peaks_gr))
#   five_prime_splice_density_data = normalize_density(calculate_density(five_prime_splice_windows[strand(five_prime_splice_windows) == strandSelection], (five_prime_splice %>% filter(strand == strandSelection))$end, peaks_gr))
#   three_prime_splice_density_data = normalize_density(calculate_density(three_prime_splice_windows[strand(three_prime_splice_windows) == strandSelection], (three_prime_splice %>% filter(strand == strandSelection))$start, peaks_gr))
#   translation_stop_density_data = normalize_density(calculate_density(translation_stop_windows[strand(translation_stop_windows) == strandSelection], (utr3_starts %>% filter(strand == strandSelection))$end, peaks_gr))
#   tts_density_data = normalize_density(calculate_density(tts_windows[strand(tts_windows) == strandSelection], (tts_starts %>% filter(strand == strandSelection))$start, peaks_gr))
#   
#   density_lists_antisense[[str_replace(sample, '_peaks_gr', '')]] = list(TSS = tss_density_data$density,
#                                                                      CDS = cds_density_data$density,
#                                                                      SS5 = five_prime_splice_density_data$density,
#                                                                      SS3 = three_prime_splice_density_data$density,
#                                                                      TLS = translation_stop_density_data$density,
#                                                                      TTS = tts_density_data$density)
# density_lists_antisense[['midpoint']] = list(tss_density_data$midpoint)
# }
# 
# featuresList = c('TSS', 'CDS', 'SS5', 'SS3', 'TLS', 'TTS')
# InputsList = c('I_M')
# mocksList = c('NLS_M', 'NES_M', 'G3BP_M')
# arsenitesList = c('NLS_S', 'NES_S', 'G3BP_S')
# 
# featureSelection = 4
# plot_Density(density_data = data.frame(position = density_lists_sense$midpoint[[1]],
#                                        I_M = density_lists_sense$I_M[featuresList[featureSelection]][[1]] + rev(density_lists_antisense$I_M[featuresList[featureSelection]][[1]]),
#                                        NLS_M = density_lists_sense$NLS_M[featuresList[featureSelection]][[1]] + rev(density_lists_antisense$NLS_M[featuresList[featureSelection]][[1]]),
#                                        NES_M = density_lists_sense$NES_M[featuresList[featureSelection]][[1]] + rev(density_lists_antisense$NES_M[featuresList[featureSelection]][[1]]),
#                                        G3BP_M = density_lists_sense$G3BP_M[featuresList[featureSelection]][[1]] + rev(density_lists_antisense$G3BP_M[featuresList[featureSelection]][[1]]),
#                                        I_S = density_lists_sense$I_S[featuresList[featureSelection]][[1]] + rev(density_lists_antisense$I_S[featuresList[featureSelection]][[1]]),
#                                        NLS_S = density_lists_sense$NLS_S[featuresList[featureSelection]][[1]] + rev(density_lists_antisense$NLS_S[featuresList[featureSelection]][[1]]),
#                                        NES_S = density_lists_sense$NES_S[featuresList[featureSelection]][[1]] + rev(density_lists_antisense$NES_S[featuresList[featureSelection]][[1]]),
#                                        G3BP_S = density_lists_sense$G3BP_S[featuresList[featureSelection]][[1]] + rev(density_lists_antisense$G3BP_S[featuresList[featureSelection]][[1]])),
#              columns_list = mocksList,
#              custom_colors = c('skyblue', 'darkseagreen2', 'salmon'),
#              densityType = 'feature_metagene',
#              featureName = featuresList[featureSelection],
#              yaxis_lims = c(0, 0.04),
#              smoothing = 0.2)
# 
# plot_Density(density_data = data.frame(position = density_lists_sense$midpoint[[1]],
#                                        I_M = density_lists_sense$I_M[featuresList[featureSelection]][[1]] + rev(density_lists_antisense$I_M[featuresList[featureSelection]][[1]]),
#                                        NLS_M = density_lists_sense$NLS_M[featuresList[featureSelection]][[1]] + rev(density_lists_antisense$NLS_M[featuresList[featureSelection]][[1]]),
#                                        NES_M = density_lists_sense$NES_M[featuresList[featureSelection]][[1]] + rev(density_lists_antisense$NES_M[featuresList[featureSelection]][[1]]),
#                                        G3BP_M = density_lists_sense$G3BP_M[featuresList[featureSelection]][[1]] + rev(density_lists_antisense$G3BP_M[featuresList[featureSelection]][[1]]),
#                                        I_S = density_lists_sense$I_S[featuresList[featureSelection]][[1]] + rev(density_lists_antisense$I_S[featuresList[featureSelection]][[1]]),
#                                        NLS_S = density_lists_sense$NLS_S[featuresList[featureSelection]][[1]] + rev(density_lists_antisense$NLS_S[featuresList[featureSelection]][[1]]),
#                                        NES_S = density_lists_sense$NES_S[featuresList[featureSelection]][[1]] + rev(density_lists_antisense$NES_S[featuresList[featureSelection]][[1]]),
#                                        G3BP_S = density_lists_sense$G3BP_S[featuresList[featureSelection]][[1]] + rev(density_lists_antisense$G3BP_S[featuresList[featureSelection]][[1]])),
#              columns_list = arsenitesList,
#              custom_colors = c('skyblue', 'darkseagreen2', 'salmon'),
#              densityType = 'feature_metagene',
#              featureName = featuresList[featureSelection],
#              yaxis_lims = c(0, 0.04),
#              smoothing = 0.2)
# 
# featureSelection = 6
# plot_Density(density_data = data.frame(position = density_lists$midpoint[[1]],
#                                        I_M = density_lists$I_M[featuresList[featureSelection]][[1]],
#                                        NLS_M = density_lists$NLS_M[featuresList[featureSelection]][[1]],
#                                        NES_M = density_lists$NES_M[featuresList[featureSelection]][[1]],
#                                        G3BP_M = density_lists$G3BP_M[featuresList[featureSelection]][[1]],
#                                        I_S = density_lists$I_S[featuresList[featureSelection]][[1]],
#                                        NLS_S = density_lists$NLS_S[featuresList[featureSelection]][[1]],
#                                        NES_S = density_lists$NES_S[featuresList[featureSelection]][[1]],
#                                        G3BP_S = density_lists$G3BP_S[featuresList[featureSelection]][[1]]),
#              columns_list = c('G3BP_M', 'G3BP_S'),
#              custom_colors = c('skyblue', 'darkseagreen2', 'salmon'),
#              densityType = 'feature_metagene',
#              featureName = featuresList[featureSelection],
#              yaxis_lims = c(0, 0.01),
#              smoothing = 0.3)
# 
# 
# 
# sampleID = 'Input'
# condition = 'Mock'
# y_lim = c(0, 0.005)
# 
# # Plot for each feature
# plot_density(tss_density_data, "TSS", sampleID, condition, y_lim, 'salmon')
# plot_density(cds_density_data, "Translation Start Sites", sampleID, condition, y_lim, 'salmon')
# plot_density(five_prime_splice_density_data, "5' Splice Sites", sampleID, condition, y_lim, 'salmon')
# plot_density(three_prime_splice_density_data, "3' Splice Sites", sampleID, condition, y_lim, 'salmon')
# plot_density(translation_stop_density_data,"Translation Stop Sites",sampleID, condition,  y_lim, 'salmon')
# plot_density(tts_density_data, "TTS", sampleID, condition, y_lim, 'salmon')

#########################################################################################################################

####################################################################################################################
# scale_values <- function(density){
#   (density - min(density)) / (max(density) - min(density))
#   }
# 
# # Calculate densities for each feature
# I_M_tss_density_data <- calculate_density(tss_windows, tss_starts$start, I_M_peaks_gr)
# I_M_cds_density_data <- calculate_density(cds_windows, cds_starts$start, I_M_peaks_gr)
# I_M_five_prime_splice_density_data <- calculate_density(five_prime_splice_windows, five_prime_splice$intron_start, I_M_peaks_gr)
# I_M_three_prime_splice_density_data <- calculate_density(three_prime_splice_windows, three_prime_splice$intron_end, I_M_peaks_gr)
# I_M_translation_stop_density_data <- calculate_density(translation_stop_windows, utr3_starts$start, I_M_peaks_gr)
# I_M_tts_density_data <- calculate_density(tts_windows, tts_starts$end, I_M_peaks_gr)
# 
# NLS_M_tss_density_data <- calculate_density(tss_windows, tss_starts$start, NLS_M_peaks_gr)
# NLS_M_cds_density_data <- calculate_density(cds_windows, cds_starts$start, NLS_M_peaks_gr)
# NLS_M_five_prime_splice_density_data <- calculate_density(five_prime_splice_windows, five_prime_splice$intron_start, NLS_M_peaks_gr)
# NLS_M_three_prime_splice_density_data <- calculate_density(three_prime_splice_windows, three_prime_splice$intron_end, NLS_M_peaks_gr)
# NLS_M_translation_stop_density_data <- calculate_density(translation_stop_windows, utr3_starts$start, NLS_M_peaks_gr)
# NLS_M_tts_density_data <- calculate_density(tts_windows, tts_starts$end, NLS_M_peaks_gr)
# 
# NES_M_tss_density_data <- calculate_density(tss_windows, tss_starts$start, NES_M_peaks_gr)
# NES_M_cds_density_data <- calculate_density(cds_windows, cds_starts$start, NES_M_peaks_gr)
# NES_M_five_prime_splice_density_data <- calculate_density(five_prime_splice_windows, five_prime_splice$intron_start, NES_M_peaks_gr)
# NES_M_three_prime_splice_density_data <- calculate_density(three_prime_splice_windows, three_prime_splice$intron_end, NES_M_peaks_gr)
# NES_M_translation_stop_density_data <- calculate_density(translation_stop_windows, utr3_starts$start, NES_M_peaks_gr)
# NES_M_tts_density_data <- calculate_density(tts_windows, tts_starts$end, NES_M_peaks_gr)
# 
# tss_density_data = data.frame(cbind(I_M_tss_density_data$midpoint, scale_values(I_M_tss_density_data$density), scale_values(NLS_M_tss_density_data$density), scale_values(NES_M_tss_density_data$density)))
# cds_density_data = cbind(I_M_cds_density_data$midpoint, I_M_cds_density_data$density, NLS_M_cds_density_data$density, NES_M_cds_density_data$density)
# five_prime_splice_density_data = cbind(I_M_five_prime_splice_density_data$midpoint, I_M_five_prime_splice_density_data$density, NLS_M_five_prime_splice_density_data$density, NES_M_five_prime_splice_density_data$density)
# three_prime_splice_density_data = cbind(I_M_three_prime_splice_density_data$midpoint, I_M_three_prime_splice_density_data$density, NLS_M_three_prime_splice_density_data$density, NES_M_three_prime_splice_density_data$density)
# translation_stop_density_data = cbind(I_M_translation_stop_density_data$midpoint, I_M_translation_stop_density_data$density, NLS_M_translation_stop_density_data$density, NES_M_translation_stop_density_data$density)
# tts_density_data = cbind(I_M_tts_density_data$midpoint, I_M_tts_density_data$density, NLS_M_tts_density_data$density, NES_M_tts_density_data$density)
# 
# colnames(tss_density_data) = c('midpoint', 'Input', 'NLS', 'NES')
# colnames(cds_density_data) = c('midpoint', 'Input', 'NLS', 'NES')
# colnames(five_prime_splice_density_data) = c('midpoint', 'Input', 'NLS', 'NES')
# colnames(three_prime_splice_density_data) = c('midpoint', 'Input', 'NLS', 'NES')
# colnames(translation_stop_density_data) = c('midpoint', 'Input', 'NLS', 'NES')
# colnames(tts_density_data) = c('midpoint', 'Input', 'NLS', 'NES')
#
# # Function to plot the densities
# plot_density <- function(density_data1, feature_name) {
#   ggplot(density_data1) + 
#     geom_smooth(aes(x=midpoint, y=Input), span = 0.15, colour = 'grey', se = FALSE) +
#     geom_smooth(aes(x=midpoint, y=NLS), span = 0.15, colour = 'salmon', se = FALSE) +
#     geom_smooth(aes(x=midpoint, y=NES), span = 0.15, colour = 'blue', se = FALSE) +
#     # geom_line(color="blue") +
#     geom_vline(xintercept=0, color="red", linetype="dashed") +
#     ylim(0, 2) +
#     labs(title=paste("Metagene Plot: Peak Density around", feature_name),
#          x=paste("Distance to", feature_name, "(nucleotides)"), y="Peak Density") +
#     theme_minimal() +
#     theme_bw() + 
#     theme(axis.text = element_text(size=14), 
#           axis.title = element_text(size=14, face = 'bold'), 
#           legend.text = element_text(size=14))
# }



## Peak Level Plotting Info:
## Input Mock   BC: 3   ylim: 0.001 or 0.002
## Input Ars    BC: 3   ylim: 0.002 (0.001 for splice sites?)

## NLS Mock     BC: 1   ylim: 0.001 (0.0005 for splice sites?)
## NLS Ars      BC: 1   ylim: 0.001 (0.0005 for splice sites?)
## NES Mock     BC: 1   ylim: 0.0005 (nothing on splice sites)
## NES Ars      BC: 1   ylim: 0.001  (nothing on splice sites)
## SG Mock      BC: 2   ylim: 0.0005
## SG Ars       BC: 2   ylim: 0.001 

## Nuc Mock     BC: 2   ylim: 0.002 
## Nuc Ars      BC: 2   ylim: 0.003 (0.001 for splice sites) 
## Cyto Mock    BC: 2   ylim: 0.003 (nothing on splice sites)
## Cyto Ars     BC: 2   ylim: 0.003 (some splice site signal)
###Above is good for metagene plots. Now need to create for each main category (Input/Total, Nuc, NLS,, G3BP, Frac etc)

#########################################################################################################################

## Based on ChipSeeker
####################################################################################################################
library(ChIPpeakAnno)

I_M_peaks_gr <- GRanges(seqnames=I_M_peaks$chrom, ranges=IRanges(start=I_M_peaks$start, end=I_M_peaks$end))
NLS_M_peaks_gr <- GRanges(seqnames=NLS_M_peaks$chrom, ranges=IRanges(start=NLS_M_peaks$start, end=NLS_M_peaks$end))
NES_M_peaks_gr <- GRanges(seqnames=NES_M_peaks$chrom, ranges=IRanges(start=NES_M_peaks$start, end=NES_M_peaks$end))
G3BP_M_peaks_gr <- GRanges(seqnames=G3BP_M_peaks$chrom, ranges=IRanges(start=G3BP_M_peaks$start, end=G3BP_M_peaks$end))

I_S_peaks_gr <- GRanges(seqnames=I_S_peaks$chrom, ranges=IRanges(start=I_S_peaks$start, end=I_S_peaks$end))
NLS_S_peaks_gr <- GRanges(seqnames=NLS_S_peaks$chrom, ranges=IRanges(start=NLS_S_peaks$start, end=NLS_S_peaks$end))
NES_S_peaks_gr <- GRanges(seqnames=NES_S_peaks$chrom, ranges=IRanges(start=NES_S_peaks$start, end=NES_S_peaks$end))
G3BP_S_peaks_gr <- GRanges(seqnames=G3BP_S_peaks$chrom, ranges=IRanges(start=G3BP_S_peaks$start, end=G3BP_S_peaks$end))

peaks_gr = I_M_peaks_gr

density_TSS = metagenePlot(peaks_gr, tss_windows, PeakLocForDistance = 'middle', FeatureLocForDistance = 'middle', upstream = 5e2, downstream = 5e2)
density_CDS = metagenePlot(peaks_gr, cds_windows, PeakLocForDistance = 'middle', FeatureLocForDistance = 'middle', upstream = 5e2, downstream = 5e2)
density_5SS = metagenePlot(peaks_gr, five_prime_splice_windows, PeakLocForDistance = 'middle', FeatureLocForDistance = 'middle', upstream = 5e2, downstream = 5e2)
density_3SS = metagenePlot(peaks_gr, three_prime_splice_windows, PeakLocForDistance = 'middle', FeatureLocForDistance = 'middle', upstream = 5e2, downstream = 5e2)
density_STC = metagenePlot(peaks_gr, translation_stop_windows, PeakLocForDistance = 'middle', FeatureLocForDistance = 'middle', upstream = 5e2, downstream = 5e2)
density_TTS = metagenePlot(peaks_gr, tts_windows, PeakLocForDistance = 'middle', FeatureLocForDistance = 'middle', upstream = 5e2, downstream = 5e2)


density_TSS = as.data.frame(table(density_TSS$data$distance))
colnames(density_TSS) = c('position', 'Freq')
density_TSS = density_TSS %>% mutate(density = Freq/sum(Freq))
ggplot(density_TSS, aes(x = position, y = density)) + geom_line()
  geom_smooth()
plot(density_TSS$position, density_TSS$density, type = "l")

density_CDS = as.data.frame(table(density_CDS$data$distance))
colnames(density_CDS) = c('position', 'Freq')
density_CDS = density_CDS %>% mutate(density = Freq/sum(Freq))
plot(density_CDS$position, density_CDS$density, type = "l")

density_5SS = as.data.frame(table(density_5SS$data$distance))
colnames(density_5SS) = c('position', 'Freq')
density_5SS = density_5SS %>% mutate(density = Freq/sum(Freq))
plot(density_5SS$position, density_5SS$density, type = "l")

density_3SS = as.data.frame(table(density_3SS$data$distance))
colnames(density_3SS) = c('position', 'Freq')
density_3SS = density_3SS %>% mutate(density = Freq/sum(Freq))
plot(density_3SS$position, density_3SS$density)

density_STC = as.data.frame(table(density_STC$data$distance))
colnames(density_STC) = c('position', 'Freq')
density_STC = density_STC %>% mutate(density = Freq/sum(Freq))
plot(density_STC$position, density_STC$density)

density_TTS = as.data.frame(table(density_TTS$data$distance))
colnames(density_TTS) = c('position', 'Freq')
density_TTS = density_TTS %>% mutate(density = Freq/sum(Freq))
plot(density_TTS$position, density_TTS$density)

#########################################################################################################################

## How about at bed file (Reads) level
####################################################################################################################
# bedFile = '/Users/soonyi/Desktop/collapsed_bed/filtered_sorted_combined_collapsed_bed/Input_Mock.sorted.peaks.bed'
# bedFile = '/Users/soonyi/Desktop/collapsed_bed/filtered_sorted_combined_collapsed_bed/Input_Stress.sorted.peaks.bed'

# bedFile = '/Users/soonyi/Desktop/collapsed_bed/filtered_sorted_combined_collapsed_bed/NLS_Mock.sorted.peaks.bed'
bedFile = '/Users/soonyi/Desktop/collapsed_bed/filtered_sorted_combined_collapsed_bed/NLS_Stress.sorted.peaks.bed'
# bedFile = '/Users/soonyi/Desktop/collapsed_bed/filtered_sorted_combined_collapsed_bed/NES_Mock.sorted.peaks.bed'
# bedFile = '/Users/soonyi/Desktop/collapsed_bed/filtered_sorted_combined_collapsed_bed/NES_Stress.sorted.peaks.bed'
# bedFile = '/Users/soonyi/Desktop/collapsed_bed/filtered_sorted_combined_collapsed_bed/G3BP_Mock.sorted.peaks.bed'
# bedFile = '/Users/soonyi/Desktop/collapsed_bed/filtered_sorted_combined_collapsed_bed/G3BP_Stress.sorted.peaks.bed'

# bedFile = '/Users/soonyi/Desktop/collapsed_bed/filtered_sorted_combined_collapsed_bed/Nuc_Mock.sorted.peaks.bed'
# bedFile = '/Users/soonyi/Desktop/collapsed_bed/filtered_sorted_combined_collapsed_bed/Nuc_Stress.sorted.peaks.bed'

# bedFile = '/Users/soonyi/Desktop/collapsed_bed/filtered_sorted_combined_collapsed_bed/Cyto_Mock.sorted.peaks.bed'
# bedFile = '/Users/soonyi/Desktop/collapsed_bed/filtered_sorted_combined_collapsed_bed/Cyto_Stress.sorted.peaks.bed'


reads_org = read.delim(bedFile, header=F, sep="\t")
colnames(reads_org) = c("chrom", "start", "end", "name", "score", "strand" )
reads = reads_org
reads = reads %>% filter(score > mean(score))

# Convert reads and various RNA landmarks  to GRanges objects
reads_gr <- GRanges(seqnames=reads$chrom, ranges=IRanges(start=reads$start, end=reads$end))
tss_windows <- GRanges(seqnames=tss_starts$seqid, 
                       ranges=IRanges(start=tss_starts$window_start, end=tss_starts$window_end))
cds_windows <- GRanges(seqnames=cds_starts$seqid, 
                       ranges=IRanges(start=cds_starts$window_start, end=cds_starts$window_end))
five_prime_splice_windows <- GRanges(seqnames=five_prime_splice$seqid, 
                                     ranges=IRanges(start=five_prime_splice$window_start, end=five_prime_splice$window_end))
three_prime_splice_windows <- GRanges(seqnames=three_prime_splice$seqid, 
                                      ranges=IRanges(start=three_prime_splice$window_start, end=three_prime_splice$window_end))
translation_stop_windows <- GRanges(seqnames=utr3_starts$seqid, 
                                    ranges=IRanges(start=utr3_starts$window_start, end=utr3_starts$window_end))
tts_windows <- GRanges(seqnames=tts_starts$seqid, 
                       ranges=IRanges(start=tts_starts$window_start, end=tts_starts$window_end))


# Calculate densities for each feature, function. Can playwith window width in the function
calculate_density <- function(feature_windows, feature_starts, reads_gr) {
  window_width <- 50
  window_def <- 500
  start_positions <- seq(-window_def, window_def - window_width, by=window_width)
  end_positions <- seq(-window_def + window_width, window_def, by=window_width)
  densities <- numeric(length(start_positions))
  
  for (i in 1:length(start_positions)) {
    window <- GRanges(seqnames=seqnames(feature_windows), 
                      ranges=IRanges(start=feature_starts + start_positions[i], 
                                     end=feature_starts + end_positions[i]))
    densities[i] <- sum(countOverlaps(window, reads_gr))
  }
  
  densities <- densities / (window_width * length(feature_windows))
  return(data.frame(midpoint=(start_positions + end_positions) / 2, density=densities))
}

# Calculate densities for each feature
tss_density_data <- calculate_density(tss_windows, tss_starts$start, reads_gr)
cds_density_data <- calculate_density(cds_windows, cds_starts$start, reads_gr)
five_prime_splice_density_data <- calculate_density(five_prime_splice_windows, five_prime_splice$intron_start, reads_gr)
three_prime_splice_density_data <- calculate_density(three_prime_splice_windows, three_prime_splice$intron_end, reads_gr)
translation_stop_density_data <- calculate_density(translation_stop_windows, utr3_starts$start, reads_gr)
tts_density_data <- calculate_density(tts_windows, tts_starts$end, reads_gr)  # Note: using end for TTS

# Function to plot the densities
plot_density <- function(density_data, feature_name) {
  ggplot(density_data, aes(x=midpoint, y=density)) +
    geom_line(color="blue") +
    geom_vline(xintercept=0, color="red", linetype="dashed") +
    ylim(0, 0.002) +
    labs(title=paste("NLS Mock Metagene Plot: Peak Density around", feature_name),
         x=paste("Distance to", feature_name, "(nucleotides)"), y="Peak Density") +
    theme_minimal()
}

# Plot for each feature
plot_density(tss_density_data, "TSS")
plot_density(cds_density_data, "Translation Start Sites")
plot_density(five_prime_splice_density_data, "5' Splice Sites")
plot_density(three_prime_splice_density_data, "3' Splice Sites")
plot_density(translation_stop_density_data, "Translation Stop Sites")
plot_density(tts_density_data, "TTS")
#########################################################################################################################

####################################################################################################################


#####Intron as percentage
###Need to convert from a peak based to read based analysis, easiest way is to amplify total Tagcount as number of rows per peak.

library(dplyr)
library(GenomicRanges)

head(peaks)
table(peaks$TOTAL_BC)
table(peaks$TOTAL_TagCount)
peaks_for_intron <- subset(peaks, TOTAL_BC >= 10 & TOTAL_TagCount >= 100)
table(peaks_for_intron$finalized_annotation)

# Convert introns and peaks to GRanges objects
introns_gr <- GRanges(seqnames=introns$seqid, ranges=IRanges(start=introns$intron_start, end=introns$intron_end))
peaks_gr <- GRanges(seqnames=peaks_for_intron$chrom, ranges=IRanges(start=peaks_for_intron$start, end=peaks_for_intron$end))

# Number of bins for metagene profile
num_bins <- 100

# Tile introns into bins
tiled_introns <- tile(introns_gr, width=width(introns_gr)/num_bins)

# Count overlaps between binned introns and peaks
overlap_counts <- countOverlaps(tiled_introns, peaks_gr)

# Assign a bin number to each overlap count
bin_assignments <- rep(1:num_bins, times=length(introns_gr))

# Ensure lengths match
overlap_counts <- overlap_counts[1:length(bin_assignments)]

# Grouping by bin and calculating the mean density
df <- data.frame(bin = bin_assignments, density = overlap_counts) %>%
  group_by(bin) %>%
  summarize(density = sum(density, na.rm = TRUE) / length(density), .groups = "drop")

# Checking the first few rows of df
head(df)

# Plot the metagene profile
library(ggplot2)
ggplot(df, aes(x=bin, y=density)) + 
  ylim(0, 0.005) +
  geom_line(color="blue") +
  labs(title="Metagene Profile across Intron Bodies", x="Intron (binned from 0% to 100%)", y="Average Peak Density") +
  theme_minimal()
#########################################################################################################################

####################################################################################################################
# Function to plot the densities
# plot_density <- function(density_data, feature_name, sampleName, condition, y_lim, line_color) {
#   ggplot(density_data, aes(x=midpoint, y=density)) +
#     # geom_line(color="blue") +
#     geom_vline(xintercept=0, color="red", linetype="dashed") +
#     ylim(y_lim) +
#     labs(title=paste0(sampleName, " ", condition, " Metagene Plot: Peak Density around ", feature_name),
#          x=paste("Distance to", feature_name, "(nucleotides)"), y="Peak Density") +
#     theme_minimal() +
#     geom_smooth(span = 0.1, se = F, color = line_color)
# }

#########################################################################################################################