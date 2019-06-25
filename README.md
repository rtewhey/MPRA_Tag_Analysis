# MPRA_Tag_Analysis



*(1) Extract barcode from the Illumina fastq. Run this seperately for each replicate. Remove A/H flags to increase stringency and require a strong match to the flanking adapter sequence.*

`gzip -dc [Illumina_R1.fastq.gz] | matchadapter_TagRead.pl  -A -H 20 TCTAGAGGTTCGTCG [Output_Filename]`



*(2) Create a temporary file. Used only if you need to collapse similar oligo sequences to a single feature. Typically you just pass an empty file*

`touch tmp.out`


*(3) Assign oligos to barcodes*

`cat [matchadapter_output.match] | perl associate_tags.pl stdin [Oligo-Barcode_LookupFile.barcode.ct.parsed] tmp.out  > output_file.tag`


*(4) Build the count matrix*

`perl compile_bc.pl samples.txt [output_file.out]`

samples.txt has the sampleID in column 1 and the .tag in column 2 (see example in the repo)
