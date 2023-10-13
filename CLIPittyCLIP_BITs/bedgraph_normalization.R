## CoCLIP Analysis: 
## BedGraph Normalization by Mapped Read Depth
## Written by Soon Yi
## Last Edit: 2023-09-27
## Use this after running normalized_bedgraph.sh

library(tidyr)
library(dplyr)
library(stringr)
library(data.table)

baseDir = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/normalize_BEDGRAPH/'
bgDir = paste0(baseDir, 'normalized_combined_bedgraph')
setwd(bgDir)
bgFiles = list.files(bgDir)

for (bgFile in bgFiles) {
  out_bg = str_replace(bgFile, '.bedgraph', '.normalized.bedGraph')
  
  expType = str_split(str_split(bgFile, "\\.")[[1]][1], "_")[[1]][1]
  
  compartment = str_split(str_split(bgFile, "\\.")[[1]][1], "_")[[1]][2]
  if (compartment == 'NLS') {compartment = 'Nuc'} else if (compartment == 'NES') {compartment = 'Cyto'} else if (compartment == 'G3BP') {compartment = 'SG'}
  
  condition = str_split(str_split(bgFile, "\\.")[[1]][1], "_")[[1]][3]
  if (condition == 'Ars') {condition = 'Stress'} else {condition = 'Mock'}
  
  strand = str_split(bgFile, "\\.")[[1]][4]
  if (strand == 'pos') {strand = '(+)'} else {strand = '(-)'}
  
  bg = read.delim(bgFile, header = F)
  bg = bg[2:nrow(bg), ]
  bg = bg %>% mutate(final = (rowSums(bg[, colnames(bg)[4:length(colnames(bg))]]) / length(4:length(colnames(bg)))))
  # bg = bg %>% filter(final != 0)
  bg = bg[, c('V1', 'V2', 'V3', 'final')]
  
  bg_header = paste0('track type=bedGraph name=', expType, '_', compartment, '_', condition, '_', strand)
  
  writeLines(bg_header, paste0(baseDir, '/finished_bedgraph/', out_bg))
  write.table(bg, paste0(baseDir, '/finished_bedgraph/', out_bg), row.names = F, col.names = F, quote = F, append = T)
}




# exp_IDs = c('JL0361_Input', 'JL0361_Enrich', 'JL0380_Fraction', 'JL0388_Input', 'JL0388_Enrich', 'JL1024_Pool')
# 
# mapped_read_depth = t(read.table('/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/Combined_peakCoverage_mappedDepth.txt', header = T))
# colnames(mapped_read_depth) = 'depth'
# 
# for (exp_ID in exp_IDs) {
#   bgDir = paste0(baseDir, str_split(exp_ID, '_')[[1]][1], '/', exp_ID, '/', exp_ID, '_mapped_reads/BEDGRAPH/')
#   setwd(bgDir)
#   bgFiles = list.files()
#   bgNames = unique(sapply(strsplit(bgFiles, '\\.'), '[', 1))
#   
#   for (bgName in bgNames) {
#     depth = mapped_read_depth[bgName, ]
#     pos_bgDir = paste0(bgDir, bgName, '.sorted.collapsed.pos.bedgraph')
#     rev_bgDir = paste0(bgDir, bgName, '.sorted.collapsed.rev.bedgraph')
#     
#     pos_bgDir = read.table(pos_bgDir)
#     rev_bgDir = read.table(rev_bgDir)
#     
#     pos_bgDir$V4 = as.integer(round(pos_bgDir$V4/depth*1e6))
#     rev_bgDir$V4 = as.integer(round(rev_bgDir$V4/depth*1e6))
#     
#     write.table(pos_bgDir, paste0(baseDir, 'Analysis/normalized_BEDGRAPH/', bgName, '.sorted.collapsed.pos.normalized.bedgraph'), row.names = F, col.names = F, quote = F)
#     write.table(rev_bgDir, paste0(baseDir, 'Analysis/normalized_BEDGRAPH/', bgName, '.sorted.collapsed.rev.normalized.bedgraph'), row.names = F, col.names = F, quote = F)
#   }
# }
