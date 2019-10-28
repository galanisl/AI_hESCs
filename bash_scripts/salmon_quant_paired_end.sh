#!/bin/bash/

for fq in fastq_raw/*_R1.fq.gz
do 
	sample=`basename $fq _R1.fq.gz`
	r1=$sample"_R1.fq.gz"
	r2=$sample"_R2.fq.gz"
	~/Documents/salmon-0.11.3-linux_x86_64/bin/salmon quant -i ../ref_transcriptome/salmon_index -l A -1 fastq_raw/$r1 -2 fastq_raw/$r2 --seqBias --gcBias --posBias -o tx_quantification/$sample
done

 
