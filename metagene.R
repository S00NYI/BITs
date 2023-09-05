## CoCLIP Analysis:
## Metagene Analysis
## Original script by Dr. Joseph M. Luna.
## Edited by Soon Yi
## Last edit: 2023-09-04

library(dplyr)
library(tidyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(IRanges)
library(GenomicRanges)
library(biomaRt)
library(data.table)
library(rtracklayer)
# library(Guitar)
library(stringr)
library(GenomicAlignments)
library(ggplot2)

# Read in GTF file
gtfFile = '~/Desktop/Genomics/Annotations/Homo_sapiens.GRCh38.110.gtf'
gtf_raw = readGFF(gtfFile)

# Add chr to chromosome names and filter out genes (transcripts only), focus on protein coding 
gtf_raw$seqid = ifelse(substr(gtf_raw$seqid, 1, 3) != "chr", paste0("chr", gtf_raw$seqid), gtf_raw$seqid)
gtf_raw = gtf_raw[gtf_raw$type != "gene", ]

gtf = subset(gtf_raw, transcript_biotype == "protein_coding")

# head(gtf)
# table(gtf$transcript_biotype)
# table(gtf$transcript_support_level)
# table(gtf$type)

#Define window on either side of feature
window_def = 500

# Define windows across RNA features
#Transcription Start Sites (TSS)
tss_starts <- subset(gtf, type == "five_prime_utr")
tss_starts$window_start <- ifelse(tss_starts$strand == "+", tss_starts$start - window_def, tss_starts$end - window_def)
tss_starts$window_end <- ifelse(tss_starts$strand == "+", tss_starts$start + window_def, tss_starts$end + window_def)

# Translation Start Sites (aka CDS starts)
cds_starts <- subset(gtf, type == "start_codon")
cds_starts$window_start <- ifelse(cds_starts$strand == "+", cds_starts$start - window_def, cds_starts$end - window_def)
cds_starts$window_end <- ifelse(cds_starts$strand == "+", cds_starts$start + window_def, cds_starts$end + window_def)

#### Splice acceptor sites
# Subset exons
exons <- subset(gtf, type == "exon")

# Group by transcript and arrange by start position
exons <- exons %>% 
  arrange(transcript_id, start) %>%
  group_by(transcript_id)
# head(exons)

# Calculate intron boundaries
introns <- exons %>%
  mutate(next_start = lead(start), 
         next_end = lead(end)) %>%
  rowwise() %>%
  filter(!is.na(next_start)) %>%
  summarize(intron_start = end + 1,
            intron_end = next_start - 1,
            seqid = seqid[1],
            strand = strand[1]) %>%
  ungroup()

# Display the intron boundaries
# head(introns)

# 5' Splice Sites (centered around intron starts)
five_prime_splice <- introns
five_prime_splice$window_start <- ifelse(five_prime_splice$strand == "+", five_prime_splice$intron_start - window_def, five_prime_splice$intron_start - window_def)
five_prime_splice$window_end <- ifelse(five_prime_splice$strand == "+", five_prime_splice$intron_start + window_def, five_prime_splice$intron_start + window_def)

# 3' Splice Sites (centered around intron ends)
three_prime_splice <- introns
three_prime_splice$window_start <- ifelse(three_prime_splice$strand == "+", three_prime_splice$intron_end - window_def, three_prime_splice$intron_end - window_def)
three_prime_splice$window_end <- ifelse(three_prime_splice$strand == "+", three_prime_splice$intron_end + window_def, three_prime_splice$intron_end + window_def)

# Confirm the number of entries in five_prime_splice and three_prime_splice, should be equivalent
cat("Number of entries in five_prime_splice:", nrow(five_prime_splice), "\n")
cat("Number of entries in three_prime_splice:", nrow(three_prime_splice), "\n")

# Translation Stop Sites (aka 3'UTR starts)
utr3_starts <- subset(gtf, type == "three_prime_utr")
utr3_starts$window_start <- ifelse(utr3_starts$strand == "+", utr3_starts$start - window_def, utr3_starts$end - window_def)
utr3_starts$window_end <- ifelse(utr3_starts$strand == "+", utr3_starts$start + window_def, utr3_starts$end + window_def)

# Transcription Stop Sites (TTS)
tts_starts <- subset(gtf, type == "three_prime_utr")
tts_starts$window_start <- ifelse(tts_starts$strand == "+", tts_starts$end - window_def, tts_starts$start - window_def)
tts_starts$window_end <- ifelse(tts_starts$strand == "+", tts_starts$end + window_def, tts_starts$start + window_def)


# Read in the peaks and filter as necessary:
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

## Filter for Inputs
peaks = peaks_org %>% filter(NLS_I_M_BC + NES_I_M_BC + G3BP_I_M_BC >= 3)
peaks = peaks_org %>% filter(NLS_I_S_BC + NES_I_S_BC + G3BP_I_S_BC >= 3)

## Filter for CoCLIPs
peaks = peaks_org %>% filter(NLS_E_M_BC >= 1)
peaks = peaks_org %>% filter(NLS_E_S_BC >= 1)

peaks = peaks_org %>% filter(NES_E_M_BC >= 1)
peaks = peaks_org %>% filter(NES_E_S_BC >= 1)

peaks = peaks_org %>% filter(G3BP_E_M_BC >= 2)
peaks = peaks_org %>% filter(G3BP_E_S_BC >= 2)

## Filter for Fractionation CLIPs
peaks = peaks_org %>% filter(Nuc_F_M_BC >= 2)
peaks = peaks_org %>% filter(Nuc_F_S_BC >= 2)
peaks = peaks_org %>% filter(Cyto_F_M_BC >= 2)
peaks = peaks_org %>% filter(Cyto_F_S_BC >= 2)



# Convert peaks and various RNA landmarks  to GRanges objects
peaks_gr <- GRanges(seqnames=peaks$chrom, ranges=IRanges(start=peaks$start, end=peaks$end))
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
calculate_density <- function(feature_windows, feature_starts, peaks_gr) {
  window_width <- 10
  window_def <- 500
  start_positions <- seq(-window_def, window_def - window_width, by=window_width)
  end_positions <- seq(-window_def + window_width, window_def, by=window_width)
  densities <- numeric(length(start_positions))
  
  for (i in 1:length(start_positions)) {
    window <- GRanges(seqnames=seqnames(feature_windows), 
                      ranges=IRanges(start=feature_starts + start_positions[i], 
                                     end=feature_starts + end_positions[i]))
    densities[i] <- sum(countOverlaps(window, peaks_gr))
  }
  
  densities <- densities / (window_width * length(feature_windows))
  return(data.frame(midpoint=(start_positions + end_positions) / 2, density=densities))
}

# Calculate densities for each feature
tss_density_data <- calculate_density(tss_windows, tss_starts$start, peaks_gr)
cds_density_data <- calculate_density(cds_windows, cds_starts$start, peaks_gr)
five_prime_splice_density_data <- calculate_density(five_prime_splice_windows, five_prime_splice$intron_start, peaks_gr)
three_prime_splice_density_data <- calculate_density(three_prime_splice_windows, three_prime_splice$intron_end, peaks_gr)
translation_stop_density_data <- calculate_density(translation_stop_windows, utr3_starts$start, peaks_gr)
tts_density_data <- calculate_density(tts_windows, tts_starts$end, peaks_gr)  # Note: using end for TTS

# 7. Generate the Metagene Plot
library(ggplot2)

# Function to plot the densities
plot_density <- function(density_data, feature_name) {
  ggplot(density_data, aes(x=midpoint, y=density)) +
    geom_line(color="blue") +
    geom_vline(xintercept=0, color="red", linetype="dashed") +
    ylim(0, 0.005) + 
    labs(title=paste("Cyto Ars Metagene Plot: Peak Density around", feature_name),
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

#############################################
## How about at bed file (Reads) level

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


############################################################


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


