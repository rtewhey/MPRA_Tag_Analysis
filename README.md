# MPRA_Tag_Analysis



*Extract barcode from the Illumina fastq. Run this seperately for each replicate. Remove A/H flags to increase stringency and require a strong match to the flanking adapter sequence.*
`gzip -dc <Illumina_R1.fastq.gz> | matchadapter_TagRead.pl  -A -H 20 TCTAGAGGTTCGTCG <Output_Filename>`


*Create a temporary file. Used only if you need to collapse similar oligo sequences to a single feature. Typically you just pass an empty file*
`touch tmp.out`

*Assign oligos to barcodes*

`cat <matchadapter_output.match> | perl associate_tags.pl stdin <Oligo-Barcode_LookupFile.barcode.ct.parsed> tmp.out`
*Build the count matrix*
``perl compile_bc.pl samples.txt <output_file.out>``
