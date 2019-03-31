DNA copy number analysis of shallow whole-genome sequencing samples
================
Gregorio Alanis-Lobato

Intro
-----

This [R Markdown](http://rmarkdown.rstudio.com) Notebook presents a step-by-step guide to perform copy number analysis from low-pass whole-genome sequencing.

Data processing and quality control
-----------------------------------

This report covers the analysis for the following sample:

| ID       | Sample  | Files                                |
|:---------|:--------|:-------------------------------------|
| NIC216A8 | H9\_p75 | NIC216A8\_S7\_L001\_R1\_001.fastq.gz |
|          |         | NIC216A8\_S7\_L001\_R2\_001.fastq.gz |

We performed quality control with `FastQC v0.11.8` and adapter and quality trimming with `Trim Galore! v0.6.0`. The code below assumes that the `FastQ` files are stored in the folder `H9_p75` and that there is a directory called `qc_report` to store the results.

``` bash
fastqc H9_p75/*.fastq.gz --outdir qc_report
trim_galore --paired -o H9_p75/ H9_p75/NIC216A8_S7_L001_R1_001.fastq.gz H9_p75/NIC216A8_S7_L001_R2_001.fastq.gz
mv H9_p75/NIC216A8_S7_L001_R1_001_val1.fq H9_p75/H9_p75_R1.fastq.gz
mv H9_p75/NIC216A8_S7_L001_R2_001_val2.fq H9_p75/H9_p75_R2.fastq.gz
```

Then, we aligned the samples to the `UCSC hg19` reference genome the two-pass alignment pipeline of `STAR v2.7.0d`. The code below shows the process for sample `H9_p75`.

``` bash
# Genome index generation
STAR --runMode genomeGenerate --genomeDir idx_dir --genomeFastaFiles ucsc.hg19.fa --runThreadN 8

# Alignment (1st pass)
STAR --genomeDir idx_dir --readFilesIn H9_p75/H9_p75_R1.fastq.gz H9_p75/H9_p75_R2.fastq.gz --runThreadN 8

# 2nd index
STAR --runMode genomeGenerate --genomeDir idx2_dir --genomeFastaFiles ucsc.hg19.fa --sjdbFileChrStartEnd SJ.out.tab --sjdbOverhang 75 --runThreadN 8

# Alignment (2nd pass)
STAR --genomeDir idx2_dir --readFilesIn H9_p75/H9_p75_R1.fastq.gz H9_p75/H9_p75_R2.fastq.gz --runThreadN 8
```

Finally, the resulting SAM file are parsed into BAM file and sorted with `samtools v1.9`:

``` bash
# Convert to BAM
samtools view -S -b Aligned.out.sam > H9_p75.bam

# Sort
samtools sort H9_p75.bam -o H9_p75_sorted.bam
```

Copy number analysis
--------------------

We do copy number analysis using the `R` package `QDNAseq v1.18.0`. Briefly, `QDNAseq` performs quantitative analysis of chromosomal aberrations by dividing the genome into non-overlapping fixed-sized bins, counting the number of sequence reads in each, adjusting with a simultaneous two-dimensional loess correction for sequence mappability and GC content, and filtering of spurious regions in the genome.

First of all, we obtain annotations for the hg19 reference genome in bins of 1000kb in size.

``` r
library(QDNAseq)
bins <- getBinAnnotations(binSize = 1000)
```

### H9\_p75

Now, we read and bin the sample:

``` r
readCounts <- binReadCounts(bins, bamfiles = "data/H9_p75_sorted.bam")
```

The rest is applying filters, estimating the correction for GC content and mappability, and applying this correction.

``` r
readCountsFiltered <- applyFilters(readCounts, residual = TRUE, 
                                   blacklist = TRUE, chromosomes = c("MT"))
readCountsFiltered <- estimateCorrection(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)
```

The resulting `QDNAseqCopyNumbers` object is then normalised and outliers are smoothed

``` r
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
```

Finally, we perform segmentation and plot the copy number profile.

``` r
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
#plot(copyNumbersSegmented)
```
