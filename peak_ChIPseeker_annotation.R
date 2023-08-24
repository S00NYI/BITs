# BiocManager::install('ChIPseeker')

library(ChIPseeker)
library(GenomicFeatures)
library(dplyr)

baseDir = "/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/"
# peakFile = "Combined_peakCoverage_biotyped.bed"
peakFile = "JL1024_Pool_peakCoverage_biotyped.bed"

gtfFile = "/Users/soonyi/Desktop/Genomics/Annotations/gencode.v44.primary_assembly.annotation.gtf"  #custom gtf or gff
TxDb = makeTxDbFromGFF(gtfFile)

peakPath = paste0(baseDir, peakFile)
TagMatrix = read.table(peakPath, sep = '\t')
TagMatrix = TagMatrix %>% select(-c(ncol(TagMatrix)-6, ncol(TagMatrix)-5, ncol(TagMatrix)-4, ncol(TagMatrix)-1, ncol(TagMatrix)))

allanno = annotatePeak(peakPath,
                        TxDb = TxDb,
                        sameStrand = FALSE,
                        genomicAnnotationPriority = c("3UTR", "5UTR", "Exon", "Intron", "Intergenic"),
                        overlap = "all")

annotations = allanno@anno$annotation
annotations[grepl("^Intron", annotations)] = "Intron"
annotations[grepl("^Exon", annotations)] = "Exon"
TagMatrix$annotation = annotations

peakFile_org = sub('_biotyped.bed', '.txt', peakFile)
TagMatrix_columns = c(colnames(read.table(paste0(baseDir, peakFile_org), sep = '\t', header = TRUE)), 'gene_name', 'gene_biotype', 'annotation')
colnames(TagMatrix) = TagMatrix_columns

write.table(TagMatrix, sub("\\.bed$", "_annotated.txt", peakPath), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')


