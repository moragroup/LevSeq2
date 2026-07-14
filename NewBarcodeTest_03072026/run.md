python -m levseq.cmd plasmida /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/ /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/ref_new_barcodes_03072026.csv --barcodes /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/barcodes.fasta


python -m levseq.cmd old_run fastq_pass  ref_oligo.csv --oligopool

python -m levseq.cmd running_example fastq_pass ref_oligo.csv --barcodes barcodes.fasta


python -m levseq.cmd plasmidaV2 /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/JQW22M /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/ref_new_barcodes_03072026.csv --barcodes /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/barcodes.fasta

python -m levseq.cmd plasmid_old /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/CF69T3_fastq/old /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/ref_new_barcodes_12072026.csv

python -m levseq.cmd plasmid_new3_rb01 /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/CF69T3_fastq/new /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/ref_new_barcodes_12072026.csv --barcodes /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/barcodes.fasta

python -m levseq.cmd oligoFH /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/CF69T3_fastq/new /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/plasmid_new_plate3_oligopool.csv --oligopool --barcodes /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/barcodes.fasta


# Running with the new barcodes
python -m levseq.cmd plasmid_new3_rb01_inhouse /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/fastq_pass_levseq_13072026 /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/ref_new_barcodes_12072026.csv --barcodes /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/barcodes.fasta

python -m levseq.cmd plasmid_old1_rb01_inhouse /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/fastq_pass_levseq_13072026 /mnt/storage01/home/amora/code/LevSeq2/NewBarcodeTest_03072026/ref_old_barcodes_12072026.csv