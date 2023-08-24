# BITs
Collection of custom BioInformatics Tools and scripts.

Last Update: 2023-08-24

----
check_barcode_mismatches.py:
  - Python script to check how many mismatches are allowed for barcode demultiplexing.

peak_grooming.py:
  - Python script to organize peak matrix from CLIPittyClip pipeline.
  - Renames and reorganizes columns to make the matrix more readable.
  - Calculates biological complexities and total peak heights.
  - Normalizes raw tag counts with total mapped reads (from separate file).

peak_custom_annotation.R:
  - R script to annotate peak matrix from CLIPittyClip pipeline.
  - Requires gtf files.

