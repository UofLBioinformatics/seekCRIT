# seekCRIT
seek for circular RNA in transcriptome (identifies deferentially expressed circRNAs between two samples)

Version: 

Last Modified: 2017-##-##

Authors: 


## Prerequisites

### Software / Package
[STAR](https://github.com/alexdobin/STAR) Aligner: v2.5.2b



#### Others

## Installation





## Usage

```bash 
usage: CRIT.py [-h] -s1 S1 -s2 S2 -gtf GTF -o OUTDIR -t {SE,PE} --aligner
               ALIGNER --genomeIndex GENOMEINDEX -fa FASTA -ref REFSEQ
               [--threadNumber numThreads]
               [--deltaPSI DELTAPSI] [--highConfidence HIGHCONFIDENCE]
               [--libType {fr-unstranded,fr-firststrand,fr-secondstrand}]
               [--keepTemp {Y,N}]

Identifying and Characterizing Differentially Spliced circular RNAs between
two samples

Required arguments:
=================== 
  -s1 S1, --sample1 S1  fastq files for sample_1. Replicates are separated by
                        comma. Paired-end reads are separated by colon.
                        e.g.,s1-1.fastq,s1-2.fastq for single-end read. s1-1.R
                        1.fastq:s1-1.R2.fastq,s1-2.R1.fastq:s1-2.R2.fastq for
                        single-end read
                        
  -s2 S2, --sample2 S2  fastq files for sample_2. Replicates are separated by
                        comma. Paired-end reads are separated by colon.
                        e.g.,s2-1.fastq,s2-2.fastq for single-end read. s2-1.R
                        1.fastq:s2-1.R2.fastq,s2-2.R1.fastq:s2-2.R2.fastq for
                        single-end read
                        
  -gtf GTF, --gtf GTF   The gtf annotation file. e.g., hg38.gtf
  
  -o OUTDIR, --output OUTDIR
                        Output directory
                        
  -t {SE,PE}, --readType {SE,PE}
                        Read type. SE for Single-end read, PE for Paired-end read
                        
  --aligner ALIGNER     Aligner to use. e.g. STAR, Tophat 
  
  --genomeIndex GENOMEINDEX
                        Genome indexes for the aligner
                        
  -fa FASTA, --fasta FASTA
                        Genome sequence. e.g., hg38.fa
                        
  -ref REFSEQ, --refseq REFSEQ
                        Transcriptome in refseq format. e.g., hg38.ref.txt
                        
 optional arguments:
====================

   -h, --help            show this help message and exit
   --threadNumber numberOfThreadsk
                        Number of threads for multi-threading feature [default = 4]

  --deltaPSI DELTAPSI   Delta PSI cutoff. i.e., significant event must show
                        bigger deltaPSI than this cutoff [default = 0.05]
                        
  --highConfidence HIGHCONFIDENCE
                        Minimum number of circular junction counts required [default = 1]
                        
  --libType  {fr-unstranded,fr-firststrand,fr-secondstrand}
                        library type used by Tophat aligner [default ='fr-unstranded']
                        
  --keepTemp {Y,N}      Keep temp files or not  [default='Y']

```

### Example


```bash
python3 CRIT.py -o PEtest -t PE --aligner STAR -fa fa/hg19.fa -ref ref/hg19.ref.txt --genomeIndex /media/bio/data/STARIndex/hg19 -s1 testData/231ESRP.25K.rep-1.R1.fastq:testData/231ESRP.25K.rep-1.R2.fastq,testData/231ESRP.25K.rep-2.R1.fastq:testData/231ESRP.25K.rep-2.R2.fastq -s2 testData/231EV.25K.rep-1.R1.fastq:testData/231EV.25K.rep-1.R2.fastq,testData/231EV.25K.rep-2.R1.fastq:testData/231EV.25K.rep-2.R2.fastq -gtf testData/test.gtf --threadNumber 12 

```


### Note

* ref.txt is in the format ([Gene Predictions and RefSeq Genes with Gene Names](https://genome.ucsc.edu/FAQ/FAQformat.html#format9)) below (see details in [the example file]())

| Field       | Description                   |
| :---------: | :---------------------------- |
| geneName    | Name of gene                  |
| isoformName | Name of isoform               |
| chrom       | Reference sequence            |
| strand      | + or - for strand             |
| txStart     | Transcription start position  |
| txEnd       | Transcription end position    |
| cdsStart    | Coding region start           |
| cdsEnd      | Coding region end             |
| exonCount   | Number of exons               |
| exonStarts  | Exon start positions          |
| exonEnds    | Exon end positions            |

* hg19.fa is genome sequence in FASTA format.

## Output

See details in [the example file]()

| Field       | Description                           |
| :---------: | :------------------------------------ |
| chrom       | Chromosome                            |



## License

Copyright (C) 2017 .  See the [LICENSE](https://github.com/UofLBioinformatics/seekCRIT/blob/master/LICENSE)
file for license rights and limitations (MIT).
