## CoCLIP Analysis:
## Custom Peak Annotation
## Original script by Dr. Joseph M. Luna.
## Edited by Soon Yi
## Last edit: 2023-08-24

## Load necessary libraries
library(stringr)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(IRanges)
library(GenomicRanges)
library(biomaRt)
library(data.table)
library(rtracklayer)

## diffloop is deprecated in higher version of bioconductor, so will need to install from sources, including some of the dependencies.
# BiocManager::install(c('foreach', 'readr', 'pbapply', 'statmod', 'zoo'))
# install.packages('https://www.bioconductor.org/packages/3.8/bioc/src/contrib/Sushi_1.20.0.tar.gz', repos=NULL, type="source")
# install.packages('https://bioconductor.riken.jp/packages/3.8/bioc/src/contrib/diffloop_1.10.0.tar.gz', repos=NULL, type="source")
library(diffloop)

## Set working directory and name of GTF file to be used for annotation.
## Import and make a dataframe out of the GTF:
setwd("~/Desktop/Genomics/Annotations")
gtfFile = 'Homo_sapiens.GRCh38.110.gtf'
gtf = import(gtfFile)
gtf_df = as.data.frame(gtf)

## Subset the data_frame by gene/transcript types:
protCoding = subset(gtf_df, gene_biotype == "protein_coding")
miR = subset(gtf_df, transcript_biotype == "miRNA" & type == "transcript")
lncRNA = subset(gtf_df, transcript_biotype == "lncRNA" & type == "transcript")
rRNA = subset(gtf_df, transcript_biotype == "rRNA" & type == "transcript")
snoRNA = subset(gtf_df, transcript_biotype == "snoRNA" & type == "transcript")
scaRNA = subset(gtf_df, transcript_biotype == "scaRNA" & type == "transcript")
snRNA = subset(gtf_df, transcript_biotype == "snRNA" & type == "transcript")
miscRNA = subset(gtf_df, transcript_biotype == "misc_RNA" & type == "transcript")
Prot_retained_int = subset(gtf_df, transcript_biotype == "retained_intron" & gene_biotype == "protein_coding" & type == "transcript")
nc_retained_int = subset(gtf_df, transcript_biotype == "retained_intron" & gene_biotype != "protein_coding" & type == "transcript")

## Make Granges for GTF derived categories above:
miR.gr = GRanges(seqnames=miR$seqnames, ranges=IRanges(start=miR$start, end=miR$end, names=miR$gene_name), strand=miR$strand)
lncRNA.gr = GRanges(seqnames=lncRNA$seqnames, ranges=IRanges(start=lncRNA$start, end=lncRNA$end, names=lncRNA$gene_name), strand=lncRNA$strand)
rRNA.gr = GRanges(seqnames=rRNA$seqnames, ranges=IRanges(start=rRNA$start, end=rRNA$end, names=rRNA$gene_name), strand=rRNA$strand)
snoRNA.gr = GRanges(seqnames=snoRNA$seqnames, ranges=IRanges(start=snoRNA$start, end=snoRNA$end, names=snoRNA$gene_name), strand=snoRNA$strand)
scaRNA.gr = GRanges(seqnames=scaRNA$seqnames, ranges=IRanges(start=scaRNA$start, end=scaRNA$end, names=scaRNA$gene_name), strand=scaRNA$strand)
snRNA.gr = GRanges(seqnames=snRNA$seqnames, ranges=IRanges(start=snRNA$start, end=snRNA$end, names=snRNA$gene_name), strand=snRNA$strand)
miscRNA.gr = GRanges(seqnames=miscRNA$seqnames, ranges=IRanges(start=miscRNA$start, end=miscRNA$end, names=miscRNA$gene_name), strand=miscRNA$strand)
Prot_retained_int.gr = GRanges(seqnames=Prot_retained_int$seqnames, ranges=IRanges(start=Prot_retained_int$start, end=Prot_retained_int$end, names=Prot_retained_int$gene_name), strand=Prot_retained_int$strand)
nc_retained_int.gr = GRanges(seqnames=nc_retained_int$seqnames, ranges=IRanges(start=nc_retained_int$start, end=nc_retained_int$end, names=nc_retained_int$gene_name), strand=nc_retained_int$strand)

## Use the same GTF to make txdb
txdb = makeTxDbFromGFF(gtfFile, format="gtf")
hg38.tx = transcriptsBy(txdb, "gene")

## Extract features of interest from hg38 TxDb:
exons = unique(unlist(exonsBy(txdb, "tx", use.names=T)))    
CDS = unique(unlist(cdsBy(txdb, "tx", use.names=T)))
introns = unique(unlist(intronsByTranscript(txdb, use.names=T)))
fiveUTRs = unique(unlist(fiveUTRsByTranscript(txdb,use.names=T)))
threeUTRs = unique(unlist(threeUTRsByTranscript(txdb,use.names=T)))

## We often consider the region 10K downstream of a gene potential extended 3'UTR. 
## I will create a separate category corresponding to these regions, as it is not an option in GenomicFeatures.
## First, use elements of threeUTRs GRange object to create a BED-like variable.  This is a very useful command!:
threeUTRs.bed = data.frame(chr=seqnames(threeUTRs), start=start(threeUTRs)-1, end=end(threeUTRs), name=names(threeUTRs), score=0, strand=strand(threeUTRs))

## Here, it is necessary to separate "+" and "-" strands, because I am only extending in one direction (downstream):
UTRs.bed.plus = threeUTRs.bed[threeUTRs.bed$strand=="+",]
UTRs.bed.plus$new.start = UTRs.bed.plus$end
UTRs.bed.plus$new.end = UTRs.bed.plus$end + 10000
DNS10K.plus = UTRs.bed.plus[,c("chr", "new.start", "new.end", "name", "score", "strand")]
UTRs.bed.minus = threeUTRs.bed[threeUTRs.bed$strand=="-",]
UTRs.bed.minus$new.start = UTRs.bed.minus$start - 10000
UTRs.bed.minus$new.end = UTRs.bed.minus$start
DNS10K.minus = UTRs.bed.minus[,c("chr", "new.start", "new.end", "name", "score", "strand")]

## Re-combine plus and minus strands
DNS10K.bed = rbind(DNS10K.plus, DNS10K.minus)
colnames(DNS10K.bed) = c("chr", "start", "end", "name", "score", "strand")

## Now, re-creating GRange object from BED. Note, this is a useful command!  
DNS10K = GRanges(seqnames=DNS10K.bed$chr, ranges=IRanges(start=DNS10K.bed$start, end=DNS10K.bed$end, names=DNS10K.bed$name), strand=DNS10K.bed$strand)

## We are also interested in defining deep intergenic sites, which we define as >10000 bp away from the nearest gene in either direction
## Extract trancripts and create BED:
transcripts = transcripts(txdb)
transcripts.bed = data.frame(chr=seqnames(transcripts), start=start(transcripts)-1, end=end(transcripts), name=transcripts$tx_name, score=0, strand=strand(transcripts))

## Note, it isn't necessary to bother with strand here because we're extending symetrically:
transcripts.bed$new.start = transcripts.bed$start - 10000
transcripts.bed$new.end = transcripts.bed$end + 10000
tx_ext.bed = transcripts.bed[,c("chr", "new.start", "new.end", "name", "score", "strand")]
colnames(tx_ext.bed) = c("chr", "start", "end", "name", "score", "strand")

## Creating GRange correpsonding to transcript coordinates extended in both directions
tx_ext.gr = GRanges(seqnames=tx_ext.bed$chr, ranges=IRanges(start=tx_ext.bed$start, end=tx_ext.bed$end, names=tx_ext.bed$name), strand=tx_ext.bed$strand)

## We are also interested in identifying sites overlapping from repmasked objects
repMask_Main = read.delim('Rmsk_hg38_2023Aug.txt')
repMask_Main = repMask_Main[, c('genoName', 'genoStart', 'genoEnd', 'strand', 'repName', 'repClass', 'repFamily')]
colnames(repMask_Main) = c('chr', 'start', 'end', 'strand', 'name', 'repClass', 'repFamily')
repMask_Main = repMask_Main %>% mutate(chr = ifelse(chr == "chrM", "chrMT", chr))

tRNA.bed = subset(repMask_Main, repClass == 'tRNA')
LINE.bed = subset(repMask_Main, repClass == 'LINE')
LC_SR.bed = subset(repMask_Main, (repClass == 'Low_complexity' | repClass == 'Simple_repeat'))
LTR.bed = subset(repMask_Main, repClass == 'LTR')
Satellite.bed = subset(repMask_Main, repClass == 'Satellite')
SINE.bed = subset(repMask_Main, repClass == 'SINE')

tRNA.gr = GRanges(seqnames=tRNA.bed$chr, ranges=IRanges(start=tRNA.bed$start, end=tRNA.bed$end, names=tRNA.bed$name), strand=tRNA.bed$strand)
LINE.gr = GRanges(seqnames=LINE.bed$chr, ranges=IRanges(start=LINE.bed$start, end=LINE.bed$end, names=LINE.bed$name), strand=LINE.bed$strand)
LC_SR.gr = GRanges(seqnames=LC_SR.bed$chr, ranges=IRanges(start=LC_SR.bed$start, end=LC_SR.bed$end, names=LC_SR.bed$name), strand=LC_SR.bed$strand)
LTR.gr = GRanges(seqnames=LTR.bed$chr, ranges=IRanges(start=LTR.bed$start, end=LTR.bed$end, names=LTR.bed$name), strand=LTR.bed$strand)
Satellite.gr = GRanges(seqnames=Satellite.bed$chr, ranges=IRanges(start=Satellite.bed$start, end=Satellite.bed$end, names=Satellite.bed$name), strand=Satellite.bed$strand)
SINE.gr = GRanges(seqnames=SINE.bed$chr, ranges=IRanges(start=SINE.bed$start, end=SINE.bed$end, names=SINE.bed$name), strand=SINE.bed$strand)

tRNA.gr = rmchr(tRNA.gr)
LINE.gr = rmchr(LINE.gr)
LC_SR.gr = rmchr(LC_SR.gr)
LTR.gr = rmchr(LTR.gr)
Satellite.gr = rmchr(Satellite.gr)
SINE.gr = rmchr(SINE.gr)

## Extract gene regions for annotation of genes:
genes = genes(txdb)
genes.bed = data.frame(chr=seqnames(genes), start=start(genes)-1, end=end(genes), name=unlist(genes$gene_id), score=0, strand=strand(genes))
genes.bed.plus = genes.bed[genes.bed$strand=="+",]
genes.bed.minus = genes.bed[genes.bed$strand=="-",]
genes.bed.plus$new.start = genes.bed.plus$start
genes.bed.plus$new.end = genes.bed.plus$end + 10000
genes.bed.minus$new.start = genes.bed.minus$start - 10000
genes.bed.minus$new.end = genes.bed.minus$end
genes.bed.ext = rbind(genes.bed.plus, genes.bed.minus)
genes.bed.ext = genes.bed.ext[,c("chr", "new.start", "new.end", "name", "score", "strand")]
colnames(genes.bed.ext) = c("chr", "start", "end", "name", "score", "strand")
genes.ext.gr = GRanges(seqnames=genes.bed.ext$chr, ranges=IRanges(start=genes.bed.ext$start, end=genes.bed.ext$end, names=genes.bed.ext$name), strand=genes.bed.ext$strand) 

## Load peak file:
peaksMatrix_PATH = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/'
peaksMatrix_FILE = 'Combined_peakCoverage_groomed.txt'
peaksMatrix = read.delim(paste0(peaksMatrix_PATH, peaksMatrix_FILE), header = TRUE, sep = "\t")
peaksMatrix = peaksMatrix %>% mutate(chr = ifelse(chr == "chrM", "chrMT", chr))
names(peaksMatrix)[names(peaksMatrix) == 'chr'] = 'chrom'
peaksMatrix = data.frame(peaksMatrix, row.names = 4)

## Create GRanges object for CLIP peak regions
gr.CombinedCoCLIP = GRanges(seqnames=peaksMatrix$chrom, ranges=IRanges(start=peaksMatrix$start, end=peaksMatrix$end, names=row.names(peaksMatrix)), strand=peaksMatrix$strand)
gr.CombinedCoCLIP = rmchr(gr.CombinedCoCLIP)

## Use "findOverlaps" function to identify features that overlap our peaks:
fiveUTRs.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=fiveUTRs,  minoverlap=1, select="first"))
threeUTRs.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=threeUTRs,  minoverlap=1, select="first"))
CDS.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=CDS,  minoverlap=1, select="first"))
DNS10K.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=DNS10K,  minoverlap=1, select="first"))
introns.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=introns,  minoverlap=1, select="first"))

## This command looks for overlap with transcripts extended in both directions.
## It should in principle be inclusive of all of the above.
## This is later used to define "deep intergenic."
tx_ext.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=tx_ext.gr,  minoverlap=1, select="first"))
genes.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=genes,  minoverlap=1, select="first"))
tRNA.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=tRNA.gr,  minoverlap=1, select="first"))
LINE.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=LINE.gr,  minoverlap=1, select="first"))
LC_SR.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=LC_SR.gr,  minoverlap=1, select="first"))
LTR.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=LTR.gr,  minoverlap=1, select="first"))
Satellite.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=Satellite.gr,  minoverlap=1, select="first"))
SINE.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=SINE.gr,  minoverlap=1, select="first"))
miR.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=miR.gr,  minoverlap=1, select="first"))
lncRNA.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=lncRNA.gr,  minoverlap=1, select="first"))
rRNA.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=rRNA.gr,  minoverlap=1, select="first"))
snoRNA.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=snoRNA.gr,  minoverlap=1, select="first"))
scaRNA.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=scaRNA.gr,  minoverlap=1, select="first"))
snRNA.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=snRNA.gr,  minoverlap=1, select="first"))
miscRNA.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=miscRNA.gr,  minoverlap=1, select="first"))
Prot_retained_int.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=Prot_retained_int.gr,  minoverlap=1, select="first"))
nc_retained_int.CombinedCoCLIP = as.data.frame(findOverlaps(query=gr.CombinedCoCLIP, subject=nc_retained_int.gr,  minoverlap=1, select="first"))

## Extract index information for each element and add to peak table.
## Note, converting x to a GR kept everything in the same order, so we can just construct the table by adding new columns.
peaksMatrix$fiveUTRs = fiveUTRs.CombinedCoCLIP[,1]
peaksMatrix$threeUTRs = threeUTRs.CombinedCoCLIP[,1]
peaksMatrix$CDS = CDS.CombinedCoCLIP[,1]
peaksMatrix$DNS10K = DNS10K.CombinedCoCLIP[,1]
peaksMatrix$introns = introns.CombinedCoCLIP[,1]
peaksMatrix$tx_ext = tx_ext.CombinedCoCLIP[,1]
peaksMatrix$gene = genes.bed[genes.CombinedCoCLIP[,1], "name"]

## For the above regions, I am only interested in if a peak overlaps ANY 3'UTR (or intron or CDS, etc).  

## For the below queries, I am interested in knowing WHICH genes they overlap, so I use the returned indexes to extract that info from the appropriate BED files.
## Note that we could have written similarly structured commands above if you were interested in retrieving WHICH features (e.g. 3'UTRs) your peaks overlapped, rather than IF they simply overlapped 3'UTRs.
peaksMatrix$tRNA = tRNA.bed[tRNA.CombinedCoCLIP[,1], "name"]
peaksMatrix$tRNA1 = tRNA.CombinedCoCLIP[,1]
peaksMatrix$LINE = LINE.bed[LINE.CombinedCoCLIP[,1], "name"]
peaksMatrix$LINE1 = LINE.CombinedCoCLIP[,1]
peaksMatrix$LTR = LTR.bed[LTR.CombinedCoCLIP[,1], "name"]
peaksMatrix$LTR1 = LTR.CombinedCoCLIP[,1]
peaksMatrix$LC_SR = LC_SR.bed[LC_SR.CombinedCoCLIP[,1], "name"]
peaksMatrix$LC_SR1 = LC_SR.CombinedCoCLIP[,1]
peaksMatrix$Satellite = Satellite.bed[Satellite.CombinedCoCLIP[,1], "name"]
peaksMatrix$Satellite1 = Satellite.CombinedCoCLIP[,1]
peaksMatrix$SINE = SINE.bed[SINE.CombinedCoCLIP[,1], "name"]
peaksMatrix$SINE1 = SINE.CombinedCoCLIP[,1]
peaksMatrix$miR = miR[miR.CombinedCoCLIP[,1], "gene_name"]
peaksMatrix$miR1 = miR.CombinedCoCLIP[,1]
peaksMatrix$lncRNA = lncRNA[lncRNA.CombinedCoCLIP[,1], "gene_name"]
peaksMatrix$lncRNA1 = lncRNA.CombinedCoCLIP[,1]
peaksMatrix$rRNA = rRNA[rRNA.CombinedCoCLIP[,1], "gene_name"]
peaksMatrix$rRNA1 = rRNA.CombinedCoCLIP[,1]
peaksMatrix$snoRNA = snoRNA[snoRNA.CombinedCoCLIP[,1], "gene_name"]
peaksMatrix$snoRNA1 = snoRNA.CombinedCoCLIP[,1]
peaksMatrix$scaRNA = scaRNA[scaRNA.CombinedCoCLIP[,1], "gene_name"]
peaksMatrix$scaRNA1 = scaRNA.CombinedCoCLIP[,1]
peaksMatrix$snRNA = snRNA[snRNA.CombinedCoCLIP[,1], "gene_name"]
peaksMatrix$snRNA1 = snRNA.CombinedCoCLIP[,1]
peaksMatrix$miscRNA = miscRNA[miscRNA.CombinedCoCLIP[,1], "gene_name"]
peaksMatrix$miscRNA1 = miscRNA.CombinedCoCLIP[,1]
peaksMatrix$Prot_retained_int = Prot_retained_int[Prot_retained_int.CombinedCoCLIP[,1], "transcript_name"]
peaksMatrix$Prot_retained_int1 = Prot_retained_int.CombinedCoCLIP[,1]
peaksMatrix$nc_retained_int = nc_retained_int[nc_retained_int.CombinedCoCLIP[,1], "transcript_name"]
peaksMatrix$nc_retained_int1 = nc_retained_int.CombinedCoCLIP[,1]

## Now, a series of ifelse statments to make the output readable.  
## Here, if a NA value was returned by findOverlaps, we keep it.  
## If a numeric index was returned (indicating our peak overlapped a feature of interest), we re-assign a name.  
## Note that these commands erase the index information.
peaksMatrix$fiveUTRs = ifelse(is.na(peaksMatrix$fiveUTRs), NA, "5'UTR")
peaksMatrix$threeUTRs = ifelse(is.na(peaksMatrix$threeUTRs), NA, "3'UTR")
peaksMatrix$CDS = ifelse(is.na(peaksMatrix$CDS), NA, "CDS")
peaksMatrix$DNS10K = ifelse(is.na(peaksMatrix$DNS10K), NA, "downstream 10K")
peaksMatrix$introns = ifelse(is.na(peaksMatrix$introns), NA, "intron")
peaksMatrix$intergenic = ifelse(is.na(peaksMatrix$tx_ext), "unannotated", NA)
peaksMatrix$tRNA1 = ifelse(is.na(peaksMatrix$tRNA1), NA, "tRNA")
peaksMatrix$LINE1 = ifelse(is.na(peaksMatrix$LINE1), NA, "TE")
peaksMatrix$LTR1 = ifelse(is.na(peaksMatrix$LTR1), NA, "Other")
peaksMatrix$LC_SR1 = ifelse(is.na(peaksMatrix$LC_SR1), NA, "Other")
peaksMatrix$Satellite1 = ifelse(is.na(peaksMatrix$Satellite1), NA, "Other")
peaksMatrix$SINE1 = ifelse(is.na(peaksMatrix$SINE1), NA, "TE")
peaksMatrix$miR1 = ifelse(is.na(peaksMatrix$miR1), NA, "miRNA")
peaksMatrix$lncRNA1 = ifelse(is.na(peaksMatrix$lncRNA1), NA, "lncRNA")
peaksMatrix$rRNA1 = ifelse(is.na(peaksMatrix$rRNA1), NA, "rRNA")
peaksMatrix$snoRNA1 = ifelse(is.na(peaksMatrix$snoRNA1), NA, "snoRNA")
peaksMatrix$scaRNA1 = ifelse(is.na(peaksMatrix$scaRNA1 ), NA, "scaRNA")
peaksMatrix$snRNA1 = ifelse(is.na(peaksMatrix$snRNA1), NA, "snRNA")
peaksMatrix$miscRNA1 = ifelse(is.na(peaksMatrix$miscRNA1 ), NA, "Other")
peaksMatrix$Prot_retained_int1  = ifelse(is.na(peaksMatrix$Prot_retained_int1), NA, "CDS_Retained_intron")
peaksMatrix$nc_retained_int1  = ifelse(is.na(peaksMatrix$nc_retained_int1), NA, "ncRNA_Retained_intron")

## Now, a series of commands to format the output as I like it, and to impose my preferred hierarchy of annotation.  
## Many downstream 10K also overlap annotated 3'UTRs. 
## Actual 3'UTRs are given priority-- this statement removes "downstream 10K" annotation where an annotated 3'UTR exists
peaksMatrix$DNS10K = apply(peaksMatrix, 1, function(i) ifelse(is.na(i["DNS10K"]), NA, ifelse(is.na(i["threeUTRs"]), "downstream 10K", NA)))

## Collapse info into a single column giving region annotation of peak.  Note this uses a customized application of the paste function within an apply loop.
peaksMatrix$annotation = apply(peaksMatrix, 1, function(i) paste(i[c("fiveUTRs","threeUTRs","CDS","DNS10K","introns","tRNA1","LTR1","LINE1","SINE1", "Satellite1","LC_SR1","intergenic", "miR1", "lncRNA1", "rRNA1", "snoRNA1", "scaRNA1", "snRNA1", "miscRNA1", "Prot_retained_int1", "nc_retained_int1")][which(!is.na(c(i["fiveUTRs"], i["threeUTRs"], i["CDS"], i["DNS10K"], i["introns"], i["tRNA1"], i["LTR1"], i["LINE1"], i["SINE1"],  i["Satellite1"], i["LC_SR1"], i["intergenic"], i["miR1"], i["lncRNA1"], i["rRNA1"], i["snoRNA1"], i["scaRNA1"], i["snRNA1"], i["miscRNA1"], i["Prot_retained_int1"], i["nc_retained_int1"])))], collapse="|"))


#########################################################################################################
# Priority list
priority_list = c("3'UTR", "5'UTR", "CDS",
                  "miRNA", "lncRNA",  "rRNA", "snoRNA",  "scaRNA",  "snRNA",  "miscRNA",  "tRNA",  "TE",  "Other",  
                  "CDS_Retained_intron", "ncRNA_Retained_intron", 
                  "intron", "downstream 10K", "unannotated")

# Extract the highest priority term for finalized_annotation
peaksMatrix = peaksMatrix %>%
  mutate(finalized_annotation = sapply(strsplit(annotation, "\\|"), function(terms) {
    for (term in priority_list) {
      if (term %in% terms) {
        return(term)
      }
    }
    return(NA)
  }))

# Map terms to grouped_annotation
peaksMatrix = peaksMatrix %>%
  mutate(grouped_annotation = ifelse(finalized_annotation %in% c("5'UTR", "CDS", "3'UTR", "TE", "Other", "downstream 10K", "unannotated", "snoRNA"), finalized_annotation,
                                     ifelse(finalized_annotation %in% c("miRNA", "lncRNA", "rRNA", "scaRNA", "snRNA", "miscRNA", "tRNA"), "ncRNA",
                                            ifelse(finalized_annotation %in% c("CDS_Retained_intron", "ncRNA_Retained_intron", "intron"), "intron", "unannotated"))),
  annotation_count = sapply(strsplit(annotation, "\\|"), length))  # Count elements after splitting

#########################################################################################################
peaksMatrix[peaksMatrix == '' | peaksMatrix == 'NA'] = NA

drops = c('fiveUTRs', 'threeUTRs', 'CDS', 'DNS10K', 'introns', 'tx_ext',
          'tRNA', 'tRNA1', 'LINE', 'LINE1', 'LTR', 'LTR1', 'LC_SR', 'LC_SR1', 'Satellite', 'Satellite1', 'SINE', 'SINE1',
          'miR', 'miR1', 'lncRNA', 'lncRNA1', 'rRNA', 'rRNA1', 'snoRNA', 'snoRNA1', 'scaRNA', 'scaRNA1', 'snRNA', 'snRNA1', 'miscRNA', 'miscRNA1',
          'Prot_retained_int', 'Prot_retained_int1', 'nc_retained_int', 'nc_retained_int1', 'intergenic')
peaksMatrix = peaksMatrix[ , -which(names(peaksMatrix) %in% drops)]

peaksMatrix$peak_names = row.names(peaksMatrix)

## Use biomaRt to import ensembl necessary attributes:
mart.hs = useMart("ensembl", host = "https://useast.ensembl.org", dataset="hsapiens_gene_ensembl")
gene_names = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = peaksMatrix$gene,
                    mart = mart.hs)

peaksMatrix = peaksMatrix %>% left_join(gene_names, by = c("gene" = "ensembl_gene_id"), relationship = "many-to-many")
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(is.na(finalized_annotation), 'unannotated', finalized_annotation))

write.table(apply(peaksMatrix,2,as.character), paste0(peaksMatrix_PATH, str_replace(peaksMatrix_FILE, ".txt", "_annotated.txt")), 
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t', na = "")
