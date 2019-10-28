#!/bin/bash/

for fq in fastq_trimmed/*_trimmed.fq.gz
do 
	sample=`basename $fq _trimmed.fq.gz`
	r1=$sample"_trimmed.fq.gz"
	~/Documents/salmon-0.11.3-linux_x86_64/bin/salmon quant -i ../ref_transcriptome/salmon_index -l A -r fastq_trimmed/$r1 --seqBias --gcBias --posBias -o tx_quantification/$sample
done

 
