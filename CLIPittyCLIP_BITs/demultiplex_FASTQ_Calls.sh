cd JL0361
gunzip -k *.gz
cat JL0361_Input.fastq | perl fastx_barcode_splitter_custom.pl --bcfile JL0361_Input_BC.txt --bol --mismatches 2 --prefix "JL0361_Input_" --suffix ".fastq"
cat JL0361_Enrich.fastq | perl fastx_barcode_splitter_custom.pl --bcfile JL0361_Enrich_BC.txt --bol --mismatches 2 --prefix "JL0361_Enrich_" --suffix ".fastq"
rm JL0361_Input.fastq
rm JL0361_Enrich.fastq
mkdir demultiplexed
mv *.fastq demultiplexed
cd demultiplexed
gzip -k *.fastq
mkdir fastq_gz
mkdir fastq_raw
mv *.gz fastq_gz
mv *.fastq fastq_raw
cd ..

cd ..
cd JL0380
gunzip -k *.gz
cat JL0380_Fraction.fastq | perl fastx_barcode_splitter_custom.pl --bcfile JL0380_Fraction_BC.txt --bol --mismatches 2 --prefix "JL0380_Fraction_" --suffix ".fastq"
rm JL0380_Fraction.fastq
mkdir demultiplexed
mv *.fastq demultiplexed
cd demultiplexed
gzip -k *.fastq
mkdir fastq_gz
mkdir fastq_raw
mv *.gz fastq_gz
mv *.fastq fastq_raw
cd ..

cd ..
cd JL0388
gunzip -k *.gz
cat JL0388_Input.fastq | perl fastx_barcode_splitter_custom.pl --bcfile JL0388_Input_BC.txt --bol --mismatches 2 --prefix "JL0388_Input_" --suffix ".fastq"
cat JL0388_Enrich.fastq | perl fastx_barcode_splitter_custom.pl --bcfile JL0388_Enrich_BC.txt --bol --mismatches 2 --prefix "JL0388_Enrich_" --suffix ".fastq"
rm JL0388_Input.fastq
rm JL0388_Enrich.fastq
mkdir demultiplexed
mv *.fastq demultiplexed
cd demultiplexed
gzip -k *.fastq
mkdir fastq_gz
mkdir fastq_raw
mv *.gz fastq_gz
mv *.fastq fastq_raw
cd ..

cd ..
cd JL1024
gunzip -k *.gz
cat JL1024_Pool.fastq | perl fastx_barcode_splitter_custom.pl --bcfile JL1024_Pool_BC.txt --bol --mismatches 2 --prefix "JL1024_Pool_" --suffix ".fastq"
rm JL1024_Pool.fastq
mkdir demultiplexed
mv *.fastq demultiplexed
cd demultiplexed
gzip -k *.fastq
mkdir fastq_gz
mkdir fastq_raw
mv *.gz fastq_gz
mv *.fastq fastq_raw
cd ..
cd ..