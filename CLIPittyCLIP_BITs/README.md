Last Update: 2024-01-06

## CLIPittyCLIP_BITs
    
#### AllowedBarcodeMismatches.py:
  - Python script to check how many mismatches are allowed for barcode demultiplexing.

#### peak_grooming.py:
  - Python script to organize peak matrix from CLIPittyClip pipeline.
  - Renames and reorganizes columns to make the matrix more readable.
  - Calculates biological complexities and total peak heights.
  - Normalizes raw tag counts with total mapped reads (from separate file).

#### bedgraph_normalization.R and .sh
  - Command line and R script for making mapped depth normalized bedGraphs.
  - Requires table containing mapped depth information.
    
#### fastx_barcode_splitter_custom.pl
  - Modified version of fastx_barcode_splitter.pl for demultiplexing.
  - This script was used to generate GEO submitted demultiplexed fastq.gz files.

