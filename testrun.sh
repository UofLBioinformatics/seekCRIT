#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "seekCRIT Setup Tool usage: testrun.sh /path/to/gtf.file /path/to/genome.fasta.file  "
    exit 1
fi




if ! [ -x "$(command -v STAR)" ]; then
  echo ' STAR is not installed.' >&2
  echo 'install STAR first from https://github.com/alexdobin/STAR '
  exit 1
fi

cd seekCRIT
mkdir genomeindex

echo "Building genome index for star... "

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir genomeindex/ --genomeFastaFiles $2 --sjdbGTFfile $1 --sjdbOverhang 100

echo "Starting testing seekCRIT... "


python3 seekCRIT.py -o SEtest -t SE --aligner STAR -fa $2 -ref ref/rn6.ref.txt --genomeIndex genomeindex -s1 testData/CRTL12.fastq -s2 testData/IR12.fastq -gtf $1 --threadNumber 12
echo "seekCRIT tested with success "


