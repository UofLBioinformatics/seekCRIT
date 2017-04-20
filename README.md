# seekCRIT
seek for circular RNA in transcriptome (identifies deferentially expressed circRNAs between two samples)

Version: 

Last Modified: 2017-##-##

Authors: 


## Prerequisites

### Software / Package

#### TopHat or STAR

* TopHat & TopHat-Fusion
    + [TopHat](http://ccb.jhu.edu/software/tophat/index.shtml) 2.0.9
    + [TopHat-Fusion](http://ccb.jhu.edu/software/tophat/fusion_index.html) included in TopHat 2.0.9
* STAR (optional)
    + [STAR](https://github.com/alexdobin/STAR) 2.4.0j

#### Others

## Installation





## Usage


Usage: python3 CRIT.py [options]

Options:
    -h --help                      Show this screen.
    --version                      Show version.
    -f FUSION --fusion=FUSION      TopHat-Fusion fusion BAM file. (used in TopHat-Fusion mapping)
    -j JUNC --junc=JUNC            STAR Chimeric junction file. (used in STAR mapping)
    -g GENOME --genome=GENOME      Genome FASTA file.
    -r REF --ref=REF               Gene annotation.
    -o PREFIX --output=PREFIX      Output prefix [default: ].
    --tmp                          Keep temporary files.
    --no-fix                       No-fix mode (useful for species with poor gene annotations)
```

### Example


```bash
python3 
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

* You could use fetch_ucsc.py script to download relevant ref.txt (Known Genes, RefSeq or Ensembl) and the genome fasta file for hg19, hg38 or mm10 from UCSC.

```bash
fetch_ucsc.py hg19/hg38/mm10 ref/kg/ens/fa out
```

Example (download hg19 RefSeq gene annotation file):

```bash
fetch_ucsc.py hg19 ref ref.txt
```

## Output

See details in [the example file]()

| Field       | Description                           |
| :---------: | :------------------------------------ |
| chrom       | Chromosome                            |
| start       | Start of junction                     |
| end         | End of junction                       |
| name        | Circular RNA/Junction reads           |
| score       | Flag to indicate realignment of fusion junctions      |
| strand      | + or - for strand                     |
| thickStart  | No meaning                            |
| thickEnd    | No meaning                            |
| itemRgb     | 0,0,0                                 |
| exonCount   | Number of exons                       |
| exonSizes   | Exon sizes                            |
| exonOffsets | Exon offsets                          |
| readNumber  | Number of junction reads              |
| circType    | 'Yes' for ciRNA, and 'No' for circRNA (before 1.1.0); 'circRNA' or 'ciRNA' (after 1.1.1)|
| geneName    | Name of gene                          |
| isoformName | Name of isoform                       |
| exonIndex/intronIndex | Index (start from 1) of exon (for circRNA) or intron (for ciRNA) in given isoform (newly added in 1.1.6) |
| flankIntron | Left intron/Right intron              |

***Note: The first 12 columns are in [BED12 format](http://genome.ucsc.edu/FAQ/FAQformat.html#format1).***

## Citation


## License

Copyright (C) 2017 .  See the [LICENSE](https://github.com/UofLBioinformatics/seekCRIT/blob/master/LICENSE)
file for license rights and limitations (MIT).
