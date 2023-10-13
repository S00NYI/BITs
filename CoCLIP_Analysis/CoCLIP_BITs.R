## Bioinformatics Tools and Scripts
## Collection of custom R functions for analysis.
## Written by Soon Yi
## Last Edit: 2023-10-13

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
filterPeaksByAnnotation = function(peak_matrix, annotation_column, annotation_list) {
  temp = peak_matrix %>% filter({{ annotation_column }} %in% annotation_list)
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

## Plot Stacked bar:
plotStackedBar = function(annotation_counts, sample_list, sample_label, title, y_lim = NULL, y_tick = NULL) {
  plot = ggplot(annotation_counts %>% filter(Source %in% sample_list), aes(fill = Annotation, y=Freq, x=Source)) + 
    geom_bar(position='stack', stat='identity') +
    scale_x_discrete(labels = sample_label) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_brewer(palette = "Set3") +
    theme_bw() + 
    theme(axis.text = element_text(size=14), 
          axis.title = element_text(size=14, face = 'bold'), 
          legend.text = element_text(size=14))
  
  if(!is.null(y_lim)) {
    plot = plot + ylim(y_lim)
  }
  
  if (!is.null(y_tick)) {
    plot = plot + scale_y_continuous(breaks = seq(0, y_lim[2], by=y_tick), limits=c(0, y_lim[2]))
  }
  
  return(plot)
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


####################################################################################################################