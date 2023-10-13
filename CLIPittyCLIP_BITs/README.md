Last Update: 2023-10-12

## CLIPittyCLIP_BITs

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
