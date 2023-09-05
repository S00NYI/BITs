
depth_file="Combined_peakCoverage_mappedDepth.txt"

while IFS=$'\t' read -r NAME DEPTH; do
    for id in "/Users/soonyi/Desktop/Genomics/CoCLIP/combined_collapsedBED/collapsedBED/${NAME}.sorted.collapsed.bed"; do
        echo "running genomecov for ${id} with depth ${DEPTH}"
        bedtools genomecov -i ${id} -g /Users/soonyi/Desktop/Genomics/Annotations/GRCh38.primary_assembly.genome.fa.fai -bg -strand + -scale ${DEPTH} > ./normalized_bedgraph/$(basename "$id" .bed).pos.bedgraph
        bedtools genomecov -i ${id} -g /Users/soonyi/Desktop/Genomics/Annotations/GRCh38.primary_assembly.genome.fa.fai -bg -strand - -scale ${DEPTH} > ./normalized_bedgraph/$(basename "$id" .bed).rev.bedgraph
    done
done < "$depth_file"


# Input Mock
bedtools unionbedg -i ./normalized_bedgraph/JL0361_Input_G3BP_Input_Mock*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL0388_Input_NLS_Input_Mock*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL0388_Input_NES_Input_Mock*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL0388_Input_SG_Input_Mock*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_NLS_I_Mock*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_NES_I_Mock*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_G3BP_I_Mock*.sorted.collapsed.pos.bedgraph -header > ./combined_bedgraph/CoCLIP_Input_Mock.combined.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0361_Input_G3BP_Input_Mock*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL0388_Input_NLS_Input_Mock*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL0388_Input_NES_Input_Mock*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL0388_Input_SG_Input_Mock*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_NLS_I_Mock*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_NES_I_Mock*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_G3BP_I_Mock*.sorted.collapsed.rev.bedgraph -header > ./combined_bedgraph/CoCLIP_Input_Mock.combined.rev.bedgraph

# Input Ars
bedtools unionbedg -i ./normalized_bedgraph/JL0361_Input_G3BP_Input_Ars*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL0388_Input_NLS_Input_Ars*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL0388_Input_NES_Input_Ars*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL0388_Input_SG_Input_Ars*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_NLS_I_Ars*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_NES_I_Ars*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_G3BP_I_Ars*.sorted.collapsed.pos.bedgraph -header > ./combined_bedgraph/CoCLIP_Input_Ars.combined.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0361_Input_G3BP_Input_Ars*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL0388_Input_NLS_Input_Ars*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL0388_Input_NES_Input_Ars*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL0388_Input_SG_Input_Ars*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_NLS_I_Ars*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_NES_I_Ars*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_G3BP_I_Ars*.sorted.collapsed.rev.bedgraph -header > ./combined_bedgraph/CoCLIP_Input_Ars.combined.rev.bedgraph

# NLS Mock
bedtools unionbedg -i ./normalized_bedgraph/JL0388_Enrich_NLS_Enrich_Mock*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_NLS_E_Mock*.sorted.collapsed.pos.bedgraph -header > ./combined_bedgraph/CoCLIP_NLS_Mock.combined.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0388_Enrich_NLS_Enrich_Mock*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_NLS_E_Mock*.sorted.collapsed.rev.bedgraph -header > ./combined_bedgraph/CoCLIP_NLS_Mock.combined.rev.bedgraph

# NLS Ars
bedtools unionbedg -i ./normalized_bedgraph/JL0388_Enrich_NLS_Enrich_Ars*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_NLS_E_Ars*.sorted.collapsed.pos.bedgraph -header > ./combined_bedgraph/CoCLIP_NLS_Ars.combined.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0388_Enrich_NLS_Enrich_Ars*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_NLS_E_Ars*.sorted.collapsed.rev.bedgraph -header > ./combined_bedgraph/CoCLIP_NLS_Ars.combined.rev.bedgraph

# NES Mock
bedtools unionbedg -i ./normalized_bedgraph/JL0388_Enrich_NES_Enrich_Mock*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_NES_E_Mock*.sorted.collapsed.pos.bedgraph -header > ./combined_bedgraph/CoCLIP_NES_Mock.combined.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0388_Enrich_NES_Enrich_Mock*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_NES_E_Mock*.sorted.collapsed.rev.bedgraph -header > ./combined_bedgraph/CoCLIP_NES_Mock.combined.rev.bedgraph

# NES Ars
bedtools unionbedg -i ./normalized_bedgraph/JL0388_Enrich_NES_Enrich_Ars*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_NES_E_Ars*.sorted.collapsed.pos.bedgraph -header > ./combined_bedgraph/CoCLIP_NES_Ars.combined.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0388_Enrich_NES_Enrich_Ars*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_NES_E_Ars*.sorted.collapsed.rev.bedgraph -header > ./combined_bedgraph/CoCLIP_NES_Ars.combined.rev.bedgraph

# G3BP Mock
bedtools unionbedg -i ./normalized_bedgraph/JL0361_Enrich_G3BP_Enrich_Mock*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL0388_Enrich_SG_Enrich_Mock*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_G3BP_E_Mock*.sorted.collapsed.pos.bedgraph -header > ./combined_bedgraph/CoCLIP_G3BP_Mock.combined.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0361_Enrich_G3BP_Enrich_Mock*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL0388_Enrich_SG_Enrich_Mock*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_G3BP_E_Mock*.sorted.collapsed.rev.bedgraph -header > ./combined_bedgraph/CoCLIP_G3BP_Mock.combined.rev.bedgraph

# G3BP Ars
bedtools unionbedg -i ./normalized_bedgraph/JL0361_Enrich_G3BP_Enrich_Ars*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL0388_Enrich_SG_Enrich_Ars*.sorted.collapsed.pos.bedgraph ./normalized_bedgraph/JL1024_Pool_G3BP_E_Ars*.sorted.collapsed.pos.bedgraph -header > ./combined_bedgraph/CoCLIP_G3BP_Ars.combined.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0361_Enrich_G3BP_Enrich_Ars*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL0388_Enrich_SG_Enrich_Ars*.sorted.collapsed.rev.bedgraph ./normalized_bedgraph/JL1024_Pool_G3BP_E_Ars*.sorted.collapsed.rev.bedgraph -header > ./combined_bedgraph/CoCLIP_G3BP_Ars.combined.rev.bedgraph

# Nuc Fraction Mock
bedtools unionbedg -i ./normalized_bedgraph/JL0380_Fraction_Nuc_Fraction_Mock*.sorted.collapsed.pos.bedgraph -header > ./combined_bedgraph/FracCLIP_Nuc_Mock.combined.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0380_Fraction_Nuc_Fraction_Mock*.sorted.collapsed.rev.bedgraph -header > ./combined_bedgraph/FracCLIP_Nuc_Mock.combined.rev.bedgraph
# Nuc Fraction Ars
bedtools unionbedg -i ./normalized_bedgraph/JL0380_Fraction_Nuc_Fraction_Ars*.sorted.collapsed.pos.bedgraph -header > ./combined_bedgraph/FracCLIP_Nuc_Ars.combined.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0380_Fraction_Nuc_Fraction_Ars*.sorted.collapsed.rev.bedgraph -header > ./combined_bedgraph/FracCLIP_Nuc_Ars.combined.rev.bedgraph

# Cyt Fraction Mock
bedtools unionbedg -i ./normalized_bedgraph/JL0380_Fraction_Cyto_Fraction_Mock*.sorted.collapsed.pos.bedgraph -header > ./combined_bedgraph/FracCLIP_Cyt_Mock.combined.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0380_Fraction_Cyto_Fraction_Mock*.sorted.collapsed.rev.bedgraph -header > ./combined_bedgraph/racCLIP_Cyt_Mock.combined.rev.bedgraph
# Cyt Fraction Ars
bedtools unionbedg -i ./normalized_bedgraph/JL0380_Fraction_Cyto_Fraction_Ars*.sorted.collapsed.pos.bedgraph -header > ./combined_bedgraph/FracCLIP_Cyt_Ars.combined.pos.bedgraph
bedtools unionbedg -i ./normalized_bedgraph/JL0380_Fraction_Cyto_Fraction_Ars*.sorted.collapsed.rev.bedgraph -header > ./combined_bedgraph/FracCLIP_Cyt_Ars.combined.rev.bedgraph

## Proceed to bedgraph_normalization.R