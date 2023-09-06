## shell script to generate normalized and combined bedgraph for sample types.
## need separate tsv file for mapped depth information.
## tsv file contains normalization factor, which is 1e6/mapped_depth. The value will be used as scaling factor for genomecov.
## also need to run bedgraph_normalization.R to finish generating the bedgraph files.

depth_file="Combined_peakCoverage_mappedDepth.txt"

while IFS=$'\t' read -r NAME DEPTH; do
    for id in "/Users/soonyi/Desktop/Genomics/CoCLIP/combined_collapsedBED/collapsedBED/${NAME}.sorted.collapsed.bed"; do
        echo "Running intersect to filter to peaks..."
        bedtools intersect -s -wa -a "${id}" -b peaks_Sorted_Filtered.bed > ./filtered_bed/$(basename "$id" .bed).filtered.bed
        echo "Running genomecov for ./filtered_bed/$(basename "$id" .bed).filtered.bed with normalization factor ${DEPTH}..."
        bedtools genomecov -i ./filtered_bed/$(basename "$id" .bed).filtered.bed -g /Users/soonyi/Desktop/Genomics/Annotations/GRCh38.primary_assembly.genome.fa.fai -bg -strand + -scale ${DEPTH} > ./normalized_bedgraph/$(basename "$id" .bed).filtered.pos.bedgraph
        bedtools genomecov -i ./filtered_bed/$(basename "$id" .bed).filtered.bed -g /Users/soonyi/Desktop/Genomics/Annotations/GRCh38.primary_assembly.genome.fa.fai -bg -strand - -scale ${DEPTH} > ./normalized_bedgraph/$(basename "$id" .bed).filtered.rev.bedgraph
    done
done < "$depth_file"

# depth_file="Combined_peakCoverage_mappedDepth.txt"

# while IFS=$'\t' read -r NAME DEPTH; do
#     for id in "/Users/soonyi/Desktop/Genomics/CoCLIP/combined_collapsedBED/collapsedBED/${NAME}.sorted.collapsed.bed"; do
#         echo "Running intersect to filter to peaks..."
#         bedtools intersect -s -wa -a "${id}" -b peaks_Sorted_Filtered.bed > ./filtered_bed/$(basename "$id" .bed).filtered.bed
#         echo "Running genomecov for ./filtered_bed/$(basename "$id" .bed).filtered.bed with normalization factor ${DEPTH}..."
#         bedtools genomecov -i ./filtered_bed/$(basename "$id" .bed).filtered.bed -g /Users/soonyi/Desktop/Genomics/Annotations/GRCh38.primary_assembly.genome.fa.fai -bg -strand + > ./normalized_bedgraph/$(basename "$id" .bed).filtered.pos.bedgraph
#         bedtools genomecov -i ./filtered_bed/$(basename "$id" .bed).filtered.bed -g /Users/soonyi/Desktop/Genomics/Annotations/GRCh38.primary_assembly.genome.fa.fai -bg -strand - > ./normalized_bedgraph/$(basename "$id" .bed).filtered.rev.bedgraph
#     done
# done < "$depth_file"


# Input Mock
bedtools unionbedg -i ./normalized_bedgraph/JL0361_Input_G3BP_Input_Mock*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL0388_Input_NLS_Input_Mock*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL0388_Input_NES_Input_Mock*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL0388_Input_SG_Input_Mock*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_NLS_I_Mock*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_NES_I_Mock*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_G3BP_I_Mock*.sorted.collapsed.filtered.pos.bedgraph -header > ./normalized_combined_bedgraph/CoCLIP_Input_Mock.combined.filtered.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0361_Input_G3BP_Input_Mock*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL0388_Input_NLS_Input_Mock*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL0388_Input_NES_Input_Mock*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL0388_Input_SG_Input_Mock*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_NLS_I_Mock*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_NES_I_Mock*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_G3BP_I_Mock*.sorted.collapsed.filtered.rev.bedgraph -header > ./normalized_combined_bedgraph/CoCLIP_Input_Mock.combined.filtered.rev.bedgraph

# Input Ars
bedtools unionbedg -i ./normalized_bedgraph/JL0361_Input_G3BP_Input_Ars*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL0388_Input_NLS_Input_Ars*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL0388_Input_NES_Input_Ars*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL0388_Input_SG_Input_Ars*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_NLS_I_Ars*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_NES_I_Ars*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_G3BP_I_Ars*.sorted.collapsed.filtered.pos.bedgraph -header > ./normalized_combined_bedgraph/CoCLIP_Input_Ars.combined.filtered.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0361_Input_G3BP_Input_Ars*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL0388_Input_NLS_Input_Ars*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL0388_Input_NES_Input_Ars*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL0388_Input_SG_Input_Ars*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_NLS_I_Ars*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_NES_I_Ars*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_G3BP_I_Ars*.sorted.collapsed.filtered.rev.bedgraph -header > ./normalized_combined_bedgraph/CoCLIP_Input_Ars.combined.filtered.rev.bedgraph

# NLS Mock
bedtools unionbedg -i ./normalized_bedgraph/JL0388_Enrich_NLS_Enrich_Mock*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_NLS_E_Mock*.sorted.collapsed.filtered.pos.bedgraph -header > ./normalized_combined_bedgraph/CoCLIP_NLS_Mock.combined.filtered.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0388_Enrich_NLS_Enrich_Mock*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_NLS_E_Mock*.sorted.collapsed.filtered.rev.bedgraph -header > ./normalized_combined_bedgraph/CoCLIP_NLS_Mock.combined.filtered.rev.bedgraph

# NLS Ars
bedtools unionbedg -i ./normalized_bedgraph/JL0388_Enrich_NLS_Enrich_Ars*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_NLS_E_Ars*.sorted.collapsed.filtered.pos.bedgraph -header > ./normalized_combined_bedgraph/CoCLIP_NLS_Ars.combined.filtered.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0388_Enrich_NLS_Enrich_Ars*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_NLS_E_Ars*.sorted.collapsed.filtered.rev.bedgraph -header > ./normalized_combined_bedgraph/CoCLIP_NLS_Ars.combined.filtered.rev.bedgraph

# NES Mock
bedtools unionbedg -i ./normalized_bedgraph/JL0388_Enrich_NES_Enrich_Mock*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_NES_E_Mock*.sorted.collapsed.filtered.pos.bedgraph -header > ./normalized_combined_bedgraph/CoCLIP_NES_Mock.combined.filtered.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0388_Enrich_NES_Enrich_Mock*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_NES_E_Mock*.sorted.collapsed.filtered.rev.bedgraph -header > ./normalized_combined_bedgraph/CoCLIP_NES_Mock.combined.filtered.rev.bedgraph

# NES Ars
bedtools unionbedg -i ./normalized_bedgraph/JL0388_Enrich_NES_Enrich_Ars*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_NES_E_Ars*.sorted.collapsed.filtered.pos.bedgraph -header > ./normalized_combined_bedgraph/CoCLIP_NES_Ars.combined.filtered.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0388_Enrich_NES_Enrich_Ars*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_NES_E_Ars*.sorted.collapsed.filtered.rev.bedgraph -header > ./normalized_combined_bedgraph/CoCLIP_NES_Ars.combined.filtered.rev.bedgraph

# G3BP Mock
bedtools unionbedg -i ./normalized_bedgraph/JL0361_Enrich_G3BP_Enrich_Mock*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL0388_Enrich_SG_Enrich_Mock*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_G3BP_E_Mock*.sorted.collapsed.filtered.pos.bedgraph -header > ./normalized_combined_bedgraph/CoCLIP_G3BP_Mock.combined.filtered.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0361_Enrich_G3BP_Enrich_Mock*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL0388_Enrich_SG_Enrich_Mock*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_G3BP_E_Mock*.sorted.collapsed.filtered.rev.bedgraph -header > ./normalized_combined_bedgraph/CoCLIP_G3BP_Mock.combined.filtered.rev.bedgraph

# G3BP Ars
bedtools unionbedg -i ./normalized_bedgraph/JL0361_Enrich_G3BP_Enrich_Ars*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL0388_Enrich_SG_Enrich_Ars*.sorted.collapsed.filtered.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_G3BP_E_Ars*.sorted.collapsed.filtered.pos.bedgraph -header > ./normalized_combined_bedgraph/CoCLIP_G3BP_Ars.combined.filtered.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0361_Enrich_G3BP_Enrich_Ars*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL0388_Enrich_SG_Enrich_Ars*.sorted.collapsed.filtered.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_G3BP_E_Ars*.sorted.collapsed.filtered.rev.bedgraph -header > ./normalized_combined_bedgraph/CoCLIP_G3BP_Ars.combined.filtered.rev.bedgraph

# Nuc Fraction Mock
bedtools unionbedg -i ./normalized_bedgraph/JL0380_Fraction_Nuc_Fraction_Mock*.sorted.collapsed.filtered.pos.bedgraph -header > ./normalized_combined_bedgraph/FracCLIP_Nuc_Mock.combined.filtered.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0380_Fraction_Nuc_Fraction_Mock*.sorted.collapsed.filtered.rev.bedgraph -header > ./normalized_combined_bedgraph/FracCLIP_Nuc_Mock.combined.filtered.rev.bedgraph
# Nuc Fraction Ars
bedtools unionbedg -i ./normalized_bedgraph/JL0380_Fraction_Nuc_Fraction_Ars*.sorted.collapsed.filtered.pos.bedgraph -header > ./normalized_combined_bedgraph/FracCLIP_Nuc_Ars.combined.filtered.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0380_Fraction_Nuc_Fraction_Ars*.sorted.collapsed.filtered.rev.bedgraph -header > ./normalized_combined_bedgraph/FracCLIP_Nuc_Ars.combined.filtered.rev.bedgraph

# Cyt Fraction Mock
bedtools unionbedg -i ./normalized_bedgraph/JL0380_Fraction_Cyto_Fraction_Mock*.sorted.collapsed.filtered.pos.bedgraph -header > ./normalized_combined_bedgraph/FracCLIP_Cyt_Mock.combined.filtered.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0380_Fraction_Cyto_Fraction_Mock*.sorted.collapsed.filtered.rev.bedgraph -header > ./normalized_combined_bedgraph/FracCLIP_Cyt_Mock.combined.filtered.rev.bedgraph
# Cyt Fraction Ars
bedtools unionbedg -i ./normalized_bedgraph/JL0380_Fraction_Cyto_Fraction_Ars*.sorted.collapsed.filtered.pos.bedgraph -header > ./normalized_combined_bedgraph/FracCLIP_Cyt_Ars.combined.filtered.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0380_Fraction_Cyto_Fraction_Ars*.sorted.collapsed.filtered.rev.bedgraph -header > ./normalized_combined_bedgraph/FracCLIP_Cyt_Ars.combined.filtered.rev.bedgraph

## Proceed to bedgraph_normalization.R



## Deprecated
# cat JL0361_Input_G3BP_Input_Mock*.bed JL0388_Input_SG_Input_Mock*.bed JL1024_Pool_G3BP_I_Mock*.bed JL0388_Input_NES_Input_Mock*.bed JL1024_Pool_NES_I_Mock*.bed JL0388_Input_NLS_Input_Mock*.bed JL1024_Pool_NLS_I_Mock*.bed > CoCLIP_Input_Mock.bed
# cat JL0361_Input_G3BP_Input_Ars*.bed JL0388_Input_SG_Input_Ars*.bed JL1024_Pool_G3BP_I_Ars*.bed JL0388_Input_NES_Input_Ars*.bed JL1024_Pool_NES_I_Ars*.bed JL0388_Input_NLS_Input_Ars*.bed JL1024_Pool_NLS_I_Ars*.bed > CoCLIP_Input_Arsenite.bed

# cat JL0388_Enrich_NLS_Enrich_Mock*.bed JL1024_Pool_NLS_E_Mock*.bed > CoCLIP_NLS_Mock.bed
# cat JL0388_Enrich_NLS_Enrich_Ars*.bed JL1024_Pool_NLS_E_Ars*.bed > CoCLIP_NLS_Arsenite.bed

# cat JL0388_Enrich_NES_Enrich_Mock*.bed JL1024_Pool_NES_E_Mock*.bed > CoCLIP_NES_Mock.bed
# cat JL0388_Enrich_NES_Enrich_Ars*.bed JL1024_Pool_NES_E_Ars*.bed > CoCLIP_NES_Arsenite.bed

# cat JL0361_Enrich_G3BP_Enrich_Mock*.bed JL0388_Enrich_SG_Enrich_Mock*.bed JL1024_Pool_G3BP_E_Mock*.bed > CoCLIP_G3BP_Mock.bed
# cat JL0361_Enrich_G3BP_Enrich_Ars*.bed JL0388_Enrich_SG_Enrich_Ars*.bed JL1024_Pool_G3BP_E_Ars*.bed > CoCLIP_G3BP_Arsenite.bed

# cat JL0380_Fraction_Nuc_Fraction_Mock*.bed > FracCLIP_Nuclear_Mock.bed
# cat JL0380_Fraction_Nuc_Fraction_Ars*.bed > FracCLIP_Nuclear_Arsenite.bed
# cat JL0380_Fraction_Cyto_Fraction_Mock*.bed > FracCLIP_Cytoplasm_Mock.bed
# cat JL0380_Fraction_Cyto_Fraction_Ars*.bed > FracCLIP_Cytoplasm_Arsenite.bed

# sort -k 1,1 -k2,2n CoCLIP_Input_Mock.bed > CoCLIP_Input_Mock.sorted.bed
# sort -k 1,1 -k2,2n CoCLIP_Input_Arsenite.bed > CoCLIP_Input_Arsenite.sorted.bed
# sort -k 1,1 -k2,2n CoCLIP_NLS_Mock.bed > CoCLIP_NLS_Mock.sorted.bed
# sort -k 1,1 -k2,2n CoCLIP_NLS_Arsenite.bed > CoCLIP_NLS_Arsenite.sorted.bed
# sort -k 1,1 -k2,2n CoCLIP_NES_Mock.bed > CoCLIP_NES_Mock.sorted.bed
# sort -k 1,1 -k2,2n CoCLIP_NES_Arsenite.bed > CoCLIP_NES_Arsenite.sorted.bed
# sort -k 1,1 -k2,2n CoCLIP_G3BP_Mock.bed > CoCLIP_G3BP_Mock.sorted.bed
# sort -k 1,1 -k2,2n CoCLIP_G3BP_Arsenite.bed > CoCLIP_G3BP_Arsenite.sorted.bed
# sort -k 1,1 -k2,2n FracCLIP_Nuclear_Mock.bed > FracCLIP_Nuclear_Mock.sorted.bed
# sort -k 1,1 -k2,2n FracCLIP_Nuclear_Arsenite.bed > FracCLIP_Nuclear_Arsenite.sorted.bed
# sort -k 1,1 -k2,2n FracCLIP_Cytoplasm_Mock.bed > FracCLIP_Cytoplasm_Mock.sorted.bed
# sort -k 1,1 -k2,2n FracCLIP_Cytoplasm_Arsenite.bed > FracCLIP_Cytoplasm_Arsenite.sorted.bed
