Last Update: 2023-10-12

## CoCLIP Specific Scripts:
    
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

#### CoCLIP_BITs.R
  - R script for all custom functions.
