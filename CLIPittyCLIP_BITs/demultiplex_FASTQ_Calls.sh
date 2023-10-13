cd ..
cd JL0361
stripBarcode.pl -len 7 -format fastq -v JL0361_Input.fastq JL0361_Input_rmBC.fastq
cat JL0361_Input_rmBC.fastq | fastx_barcode_splitter.pl --bcfile JL0361_Input_BC.txt --bol --mismatches 2 --prefix "JL0361_Input_" --suffix ".fastq"
stripBarcode.pl -len 7 -format fastq -v JL0361_Enrich.fastq JL0361_Enrich_rmBC.fastq
cat JL0361_Enrich_rmBC.fastq | fastx_barcode_splitter.pl --bcfile JL0361_Enrich_BC.txt --bol --mismatches 2 --prefix "JL0361_Enrich_" --suffix ".fastq"

cd ..
cd JL0380
stripBarcode.pl -len 7 -format fastq -v JL0380_Fraction.fastq JL0380_Fraction_rmBC.fastq
cat JL0380_Fraction_rmBC.fastq | fastx_barcode_splitter.pl --bcfile JL0380_Fraction_BC.txt --bol --mismatches 2 --prefix "JL0380_Fraction_" --suffix ".fastq"

cd ..
cd JL0388
stripBarcode.pl -len 7 -format fastq -v JL0388_Input.fastq JL0388_Input_rmBC.fastq
cat JL0388_Input_rmBC.fastq | fastx_barcode_splitter.pl --bcfile JL0388_Input_BC.txt --bol --mismatches 2 --prefix "JL0388_Input_" --suffix ".fastq"
stripBarcode.pl -len 7 -format fastq -v JL0388_Enrich.fastq JL0388_Enrich_rmBC.fastq
cat JL0388_Enrich_rmBC.fastq | fastx_barcode_splitter.pl --bcfile JL0388_Enrich_BC.txt --bol --mismatches 2 --prefix "JL0388_Enrich_" --suffix ".fastq"

cd ..
cd JL1024
stripBarcode.pl -len 7 -format fastq -v JL1024_Pool.fastq JL1024_Pool_rmBC.fastq
cat JL1024_Pool_rmBC.fastq | fastx_barcode_splitter.pl --bcfile JL1024_Pool_BC.txt --bol --mismatches 2 --prefix "JL1024_Pool_" --suffix ".fastq"