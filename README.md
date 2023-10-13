# BITs
Collection of custom <ins>**B**</ins>io<ins>**I**</ins>nformatics <ins>**T**</ins>ools and <ins>**S**</ins>cripts.

Mostly for Post-CLIPittyClip pipeline analysis, but not always. 

Last Update: 2023-10-12

----
## CoCLIP_Analysis
CoCLIP_Analysis scripts contain many custom functions. See *BITs.R* for more details.<br>
All CoCLIP_Analysis scripts except for the part 0 uses annotated peak matrix as input.
    
#### CoCLIP_Analysis_Part0_PeakAnnotation.R
  - R script for peak annotation.
  - Gene features are extracted from:
    - Gene annotation (GTF file either from Ensembl or Gencode)
    - RepeatMasked region (from UCSC Table Browser)
  -  Annotations will be given priority as defined in "priority_list".
  -  Outputs peak matrix with new annotation columns appended.

#### CoCLIP_Analysis_Part1_StackedBar.R
  - R script to create stacked bar graphs for peaks distribution.

#### CoCLIP_Analysis_Part2_Comparison_to_Fractionation.R
  - R script comparing fractionation CLIP to coCLIP.

#### CoCLIP_Analysis_Part3_MotifAnalysis.R
  - R script for all motif analysis.
  - In addition to peak matrix, following files are utilized:
    - HuR eCLIP peak file from ENCODE consortium.
    - Motif density calculated using HOMER (see Homer_Calls.sh) 

#### CoCLIP_Analysis_Part4_MetaGeneAnalysis.R
  - R script for custom metagene plot generation.

#### CoCLIP_Analysis_Part5_StressGranuleSpecificAnalysis.R
  - R script for SG specific analysis.

#### CoCLIP_Analysis_Part6_GeneLevelAnalysis.R
  - R script for analysis at gene levels (GSEA, transcript length, and RNASeq comparison)

#### CoCLIP_Analysis_Part7_HuRFlowAnalysis.R
  - R script for context specific target gene analysis.

----
## CLIPittyCLIP_BITs
Scripts associated with CLIPittyClip but not necessarily with CoCLIP Analysis.

#### CLIPittyClip_Calls.sh
  - Command lines for CLIPittyClip calls.
    
#### check_barcode_mismatches.py:
  - Python script to check how many mismatches are allowed for barcode demultiplexing.

#### peak_grooming.py:
  - Python script to organize peak matrix from CLIPittyClip pipeline.
  - Renames and reorganizes columns to make the matrix more readable.
  - Calculates biological complexities and total peak heights.
  - Normalizes raw tag counts with total mapped reads (from separate file).

#### Homer_Calls.sh
  - Command lines for HOMER calls.

#### bedgraph_normalization.R and .sh
  - Command line and R script for making mapped depth normalized bedGraphs.
  - Requires table containing mapped depth information.

#### demultiplex_FASTQ_Calls.sh
  - Command line calls for CTK and fastX for demultiplexing initial fastq files.

----
## Random_BITs
Scripts that are for testing or some other random stuffs.

