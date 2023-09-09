## annotatePeaks.pl for peak annotation 
annotatePeaks.pl /Users/soonyi/Desktop/Genomics/CoCLIP/combined_collapsedBED/Combined_peaks/peaks.bed hg38 > /Users/soonyi/Desktop/Genomics/CoCLIP/combined_collapsedBED/Combined_peaks/peaks.annotated.bed

## findMotifsGenome.pl 
for peakFILE in *.txt; do
    FILE_NAME="${peakFILE%.txt}"
    findMotifsGenome.pl ${FILE_NAME}.txt hg38 ./motif/${FILE_NAME} -rna -len 7 -S 10 -noknown
done

## makeMetaGeneProfile.pl 
for peakFILE in *.txt; do
    FILE_NAME="${peakFILE%.txt}"
    makeMetaGeneProfile.pl rna hg38 -size 500 -bin 50 -gRatio 2 -gbin 20 -p ${FILE_NAME}.txt > ./metaGene/${FILE_NAME}_metaGene.txt
done

## annotatePeaks.pl for motif density 
for bedFILE in *.bed; do
    FILE_NAME="${peakFILE%.txt}"
    annotatePeaks.pl ./6mer_motif/${FILE_NAME}_combined_TagDir/peaks.txt hg38 -size 1000 -hist 20 -m ./6mer_motif/${FILE_NAME}_combined_TagDir/homerMotifs.motifs7 > ./motifDensity/${FILE_NAME}.txt
done


## in the collapsedBed directory:
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

## Homer Motif Enrichment with Bed file after filtering to peaks:
# for id in *.bed; do
#     bedtools intersect -s -wa -a "${id}" -b peaks_Sorted_Filtered.bed > "./filtered/$(basename "$id" .bed).filtered.bed"
# done

## Call peaks for each combined bed file:

for bedFILE in ./combined_bed/*.bed; do
    FILE_NAME="${bedFILE##*/}"
    FILE_NAME="${FILE_NAME%.bed}"
    makeTagDirectory ./peaks/${FILE_NAME}_combined_TagDir/ ./combined_bed/${FILE_NAME}.bed -single -format bed
    findPeaks ./peaks/${FILE_NAME}_combined_TagDir/ -o auto -style factor -L 2 -localSize 10000 -strand separate -minDist 50 -size 20 -fragLength 25
done

## Search for PAR-CLIP motifs:

for bedFILE in ./combined_bed/*.bed; do
    FILE_NAME="${bedFILE##*/}"
    FILE_NAME="${FILE_NAME%.bed}"
    findMotifsGenome.pl ./peaks/${FILE_NAME}_combined_TagDir/peaks.txt hg38 ./motifs/ -rna -len 5 -size 50 -find ./motifs/HuR_PAR_CLIP.motifs > ./motifs/${FILE_NAME}.HuR_PAR_CLIP.MotifCounts.txt
done

## Search for enriched motifs:

for bedFILE in ./combined_bed/*.bed; do
    FILE_NAME="${bedFILE##*/}"
    FILE_NAME="${FILE_NAME%.bed}"
    findMotifsGenome.pl ./peaks/${FILE_NAME}_combined_TagDir/peaks.txt hg38 ./motifs/${FILE_NAME}/5mer/ -rna -len 5 -size 50 -S 25 -noknown -norevopp -p 8 -mask -mis 2 -bits
    findMotifsGenome.pl ./peaks/${FILE_NAME}_combined_TagDir/peaks.txt hg38 ./motifs/${FILE_NAME}/6mer/ -rna -len 6 -size 50 -S 25 -noknown -norevopp -p 8 -mask -mis 2 -bits
    findMotifsGenome.pl ./peaks/${FILE_NAME}_combined_TagDir/peaks.txt hg38 ./motifs/${FILE_NAME}/7mer/ -rna -len 7 -size 50 -S 25 -noknown -norevopp -p 8 -mask -mis 2 -bits
done

## annotatePeaks.pl for motif density 

for bedFILE in ./combined_bed/*.bed; do
    FILE_NAME="${bedFILE##*/}"
    FILE_NAME="${FILE_NAME%.bed}"
    annotatePeaks.pl ./peaks/${FILE_NAME}_combined_TagDir/peaks.txt hg38 -size 1000 -hist 50 -m ./motifs/HuR_PAR_CLIP.motifs -cpu 7> ./motifs_density/${FILE_NAME}.txt
done
    
