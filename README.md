**NOTE: UNDER CONSTRUCTION, the actual software is NOT released yet**

# seekCRIT
seek for circular RNA in transcriptome (identifies deferentially expressed circRNAs between two samples)

Version: (no official version yet)

Last Modified: 2017-06-02

Authors: [Bioinformatics Lab](http://bioinformatics.louisville.edu/lab/index.php), University of Louisville, [Kentucky Biomedical Research Infrastructure Network (KBRIN)](http://louisville.edu/research/kbrin/)


## Prerequisites

### Software / Package
- [STAR](https://github.com/alexdobin/STAR) Aligner: v2.5.2b

#### Others
- [pysam](https://github.com/pysam-developers/pysam) >=0.9.1.4

- [numpy](https://github.com/numpy/numpy) >=1.11.2

- [scipy](https://docs.scipy.org/doc/scipy/reference/dev/index.html)

- [fisher](https://pypi.python.org/pypi/fisher/) (Fisher’s exact test)

- [mne](https://github.com/mne-tools/mne-python/tree/master/mne) (FDR calculation)


## Installation
1 Download seekCRIT
```bash
git clone https://github.com/UofLBioinformatics/seekCRIT.git

cd seekCRIT
```

2 install required packages
```bash
pip install -r Prerequisites.txt
```



## Usage

```bash 
usage: seekCRIT.py [-h] -s1 S1 -s2 S2 -gtf GTF -o OUTDIR -t {SE,PE} --aligner
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
#### Paired-end reads
```bash
python3 seekCRIT.py -o PEtest -t PE --aligner STAR -fa fa/hg19.fa -ref ref/hg19.ref.txt --genomeIndex /media/bio/data/STARIndex/hg19 -s1 testData/231ESRP.25K.rep-1.R1.fastq:testData/231ESRP.25K.rep-1.R2.fastq,testData/231ESRP.25K.rep-2.R1.fastq:testData/231ESRP.25K.rep-2.R2.fastq -s2 testData/231EV.25K.rep-1.R1.fastq:testData/231EV.25K.rep-1.R2.fastq,testData/231EV.25K.rep-2.R1.fastq:testData/231EV.25K.rep-2.R2.fastq -gtf testData/test.gtf --threadNumber 12 
```
#### Single-end reads

```bash
python3 seekCRIT.py -o SEtest -t SE --aligner STAR -fa fa/hg19.fa -ref ref/hg19.ref.txt --genomeIndex /media/bio/data/STARIndex/hg19 -s1 testData/231ESRP.25K.rep-1.R1.fastq,testData/231ESRP.25K.rep-1.R2.fastq,testData/231ESRP.25K.rep-2.R1.fastq,testData/231ESRP.25K.rep-2.R2.fastq -s2 testData/231EV.25K.rep-1.R1.fastq,testData/231EV.25K.rep-1.R2.fastq,testData/231EV.25K.rep-2.R1.fastq,testData/231EV.25K.rep-2.R2.fastq -gtf testData/test.gtf --threadNumber 12 
```

## Note
- hg19.fasta is genome sequence in FASTA format and can be downloaded from [UCSC Genome Browser](https://genome.ucsc.edu/).
- Transcriptome should be in refseq format below (see more details in the [example](https://github.com/UofLBioinformatics/seekCRIT/blob/master/example/hg19.ref.txt) ):

| Field                           | Description                                                                  |
| :------------------------------:|:---------------------------------------------------------------------------- |
| geneName                        | Name of gene                                                                 |
| isoform_name                    | name of isoform                                                              |
| chrom                           | chromosme                                                                    |
| strand                          |  strand  (+|-)                                                             |
| txStart                         | Transcription start position                                                 |
| txEnd                           | Transcription end position                                                   |
| cdsStart                        |Coding region end   		                                                       |
| exonCount						            | 		    Number of exons 		                    						                 |
| exonStarts					            | 		 Exon start positions       								                             |
| exonEnds						            | 		    Exon end positions           						                             |

- It is not obligatory to provide REFSEQ file, we made script (GTFtoREFSEQ) to convert from gtf to refseq that is used in the main code if no refseq file is provided.


## Output

See details in [the example file](https://github.com/UofLBioinformatics/seekCRIT/blob/master/example/circRNAs.pVal.FDR.txt)

| Field                           | Description                                                                  |
| :------------------------------:|:---------------------------------------------------------------------------- |
| chrom                           | chromosome                                                                   |
| circRNA_start                   | circular RNA 5' end position                                                 |
| circRNA_end                     | circular RNA 3' end position                                                 |
| strand                          | DNA strand  (+|-)                                                            |
| exonCount                       | number of exons included in the circular RNA transcript                      |
| exonSizes                       | size of exons included in the circular RNA transcript                        |
| exonOffsets                     | offsets of exons included in the circular RNA transcript  					         |  
| circType					          	  | 		        						                    		                             |
| geneName						            | 		        		                    						                             |
| isoformName					            | 		                             								                             |
| exonIndexOrIntronIndex		      |                     		        								                             |
| FlankingIntrons			        	  | 		                             								                             |
| CircularJunctionCount_Sample_1  | read count of the circular junction in sample # 1                            |
| LinearJunctionCount_Sample_1	  | 	                        											                             |
| CircularJunctionCount_Sample_2  | 	                        											                             |
| LinearJunctionCount_Sample_2	  | 		                        										                             |
| PBI_Sample_1					          | 		                        										                             |
| PBI_Sample_2				         	  | 				                        								                             |
| deltaPBI(PBI_1-PBI_2)			      | 			                            							                             |
| pValue						              | 				                        								                             |
| FDR							                |                          												                             |


## License

Copyright (C) 2017 .  See the [LICENSE](https://github.com/UofLBioinformatics/seekCRIT/blob/master/LICENSE)
file for license rights and limitations (MIT).
