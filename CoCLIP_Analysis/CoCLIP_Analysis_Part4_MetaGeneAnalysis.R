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



