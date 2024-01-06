## command lines for HOMER program used for CoCLIP analysis.
## Use this as a reference for other HOMER usage.
## Written by Soon Yi
## Last Update: 2024-01-06

## Use findMotifsGenome.pl for de novo motif finding: 
for peakFILE in *.txt; do
    FILE_NAME="${peakFILE%.txt}"
    findMotifsGenome.pl ${FILE_NAME}.txt hg38 ./motif/${FILE_NAME} -rna -len 7 -S 10 -noknown
done

## FOR FIGURE 4:
## Put all collapsedBed files from CLIPittyClip into a folder:
cat JL0361_Input_G3BP_Input_Mock*.bed JL0388_Input_SG_Input_Mock*.bed JL1024_Pool_G3BP_I_Mock*.bed JL0388_Input_NES_Input_Mock*.bed JL1024_Pool_NES_I_Mock*.bed JL0388_Input_NLS_Input_Mock*.bed JL1024_Pool_NLS_I_Mock*.bed > CoCLIP_Input_Mock.bed
cat JL0361_Input_G3BP_Input_Ars*.bed JL0388_Input_SG_Input_Ars*.bed JL1024_Pool_G3BP_I_Ars*.bed JL0388_Input_NES_Input_Ars*.bed JL1024_Pool_NES_I_Ars*.bed JL0388_Input_NLS_Input_Ars*.bed JL1024_Pool_NLS_I_Ars*.bed > CoCLIP_Input_Arsenite.bed
cat JL0388_Enrich_NLS_Enrich_Mock*.bed JL1024_Pool_NLS_E_Mock*.bed > CoCLIP_NLS_Mock.bed
cat JL0388_Enrich_NLS_Enrich_Ars*.bed JL1024_Pool_NLS_E_Ars*.bed > CoCLIP_NLS_Arsenite.bed
cat JL0388_Enrich_NES_Enrich_Mock*.bed JL1024_Pool_NES_E_Mock*.bed > CoCLIP_NES_Mock.bed
cat JL0388_Enrich_NES_Enrich_Ars*.bed JL1024_Pool_NES_E_Ars*.bed > CoCLIP_NES_Arsenite.bed
cat JL0361_Enrich_G3BP_Enrich_Mock*.bed JL0388_Enrich_SG_Enrich_Mock*.bed JL1024_Pool_G3BP_E_Mock*.bed > CoCLIP_G3BP_Mock.bed
cat JL0361_Enrich_G3BP_Enrich_Ars*.bed JL0388_Enrich_SG_Enrich_Ars*.bed JL1024_Pool_G3BP_E_Ars*.bed > CoCLIP_G3BP_Arsenite.bed
cat JL0380_Fraction_Nuc_Fraction_Mock*.bed > FracCLIP_Nuclear_Mock.bed
cat JL0380_Fraction_Nuc_Fraction_Ars*.bed > FracCLIP_Nuclear_Arsenite.bed
cat JL0380_Fraction_Cyto_Fraction_Mock*.bed > FracCLIP_Cytoplasm_Mock.bed
cat JL0380_Fraction_Cyto_Fraction_Ars*.bed > FracCLIP_Cytoplasm_Arsenite.bed

sort -k 1,1 -k2,2n CoCLIP_Input_Mock.bed > CoCLIP_Input_Mock.sorted.bed
sort -k 1,1 -k2,2n CoCLIP_Input_Arsenite.bed > CoCLIP_Input_Arsenite.sorted.bed
sort -k 1,1 -k2,2n CoCLIP_NLS_Mock.bed > CoCLIP_NLS_Mock.sorted.bed
sort -k 1,1 -k2,2n CoCLIP_NLS_Arsenite.bed > CoCLIP_NLS_Arsenite.sorted.bed
sort -k 1,1 -k2,2n CoCLIP_NES_Mock.bed > CoCLIP_NES_Mock.sorted.bed
sort -k 1,1 -k2,2n CoCLIP_NES_Arsenite.bed > CoCLIP_NES_Arsenite.sorted.bed
sort -k 1,1 -k2,2n CoCLIP_G3BP_Mock.bed > CoCLIP_G3BP_Mock.sorted.bed
sort -k 1,1 -k2,2n CoCLIP_G3BP_Arsenite.bed > CoCLIP_G3BP_Arsenite.sorted.bed
sort -k 1,1 -k2,2n FracCLIP_Nuclear_Mock.bed > FracCLIP_Nuclear_Mock.sorted.bed
sort -k 1,1 -k2,2n FracCLIP_Nuclear_Arsenite.bed > FracCLIP_Nuclear_Arsenite.sorted.bed
sort -k 1,1 -k2,2n FracCLIP_Cytoplasm_Mock.bed > FracCLIP_Cytoplasm_Mock.sorted.bed
sort -k 1,1 -k2,2n FracCLIP_Cytoplasm_Arsenite.bed > FracCLIP_Cytoplasm_Arsenite.sorted.bed

cat *.sorted.bed > All_Libraries.bed
cat *Mock.sorted.bed > All_Mocks.bed
cat *Arsenite.sorted.bed > All_Arsenite.bed

## Call peaks for each combined bed file:
## -strand separate option will get peaks on both ends.
for bedFILE in ./combined_bed/*.bed; do
    FILE_NAME="${bedFILE##*/}"
    FILE_NAME="${FILE_NAME%.bed}"
    makeTagDirectory ./peaks/${FILE_NAME}_combined_TagDir/ ./combined_bed/${FILE_NAME}.bed -single -format bed
    findPeaks ./peaks/${FILE_NAME}_combined_TagDir/ -o auto -style factor -L 2 -localSize 10000 -strand separate -minDist 50 -size 20 -fragLength 25
done

## Search for de novo motifs --> Supplement
## -norevopp only search for motifs in the sense strand relative to the peaks. So only searches strand specified by the peak.
for bedFILE in ./combined_bed/*.bed; do
    FILE_NAME="${bedFILE##*/}"
    FILE_NAME="${FILE_NAME%.bed}"
    findMotifsGenome.pl ./peaks/${FILE_NAME}_combined_TagDir/peaks.txt hg38 ./motifs/denovo/${FILE_NAME}/7mer/ -rna -len 7 -size 50 -S 25 -noknown -norevopp -p 8 -mask -mis 2 -bits
done

## annotatePeaks.pl for motif density 
## Positive strand here means relative to the original peak strand.
## When making density plot, use the one from the "positive" only.
for bedFILE in ./combined_bed/*.bed; do
    FILE_NAME="${bedFILE##*/}"
    FILE_NAME="${FILE_NAME%.bed}"
    annotatePeaks.pl ./peaks/${FILE_NAME}_combined_TagDir/peaks.txt hg38 -cpu 10 -size 1000 -hist 50 -m ./motifs/HuR.subsets.motifs > ./motifs/density/${FILE_NAME}.txt
done
