

import pysam
from collections import defaultdict
import os
import argparse


parser = argparse.ArgumentParser(description='Detection circRNAs from chimeric file')
parser.add_argument('-g', '--genome', dest='genome', required=True, help='The UCSC genome build fasta file: hg38.fasta, hg19.fasta, mm10.fasta ..')
parser.add_argument('-o', '--output', dest='output', required=True, help='output directory')
parser.add_argument('-j', '--junction', dest='junction', required=True, help='The chimeric file')
parser.add_argument('-r', '--ref', dest='ref_file', help='ref file ')


def extract_jucnctions(junc_file):
    jFile = open(junc_file, 'r');
    junctions= {};
    juctions_chr = defaultdict(list)
    for line in jFile:  ## process each line

        ele = line.strip().split('\t');
        chr1=(ele[0]);
        site1=int(ele[1]);
        strand1=(ele[2])
        chr2=(ele[3])
        site2=int(ele[4])
        strand2=(ele[5])
        flag=int(ele[6])
        if flag>=0 and strand1==strand2 and chr1==chr2:
            if strand1 == '+':
                start = int(site2)
                end = int(site1) - 1
            else:
                start = int(site1)
                end = int(site2) - 1
            if start<end:
                key = ':'.join([chr1,str(start),str(end)]);
                if key not in junctions:
                    juctions_chr[chr1].append([start, end,strand1, key])
                    value = [1, '\t'.join([chr1,str(start),str(end)])];
                    junctions[key]=value;
                else:
                    junctions[key][0]=junctions[key][0]+1;

    return (juctions_chr,junctions)



def extract_genes(ref):
    rFile = open(ref, 'r');

    iso_chr = defaultdict(list)
    for line in rFile:  ## process each line

        ele = line.strip().split('\t');
        gene_id = (ele[0]);
        iso_id = (ele[1]);
        chr = (ele[2])
        strand = (ele[3])
        exons_start=(ele[9])
        exons_ends = (ele[10])
        start_points=list(eval(exons_start))
        end_points = list(eval(exons_ends))
        key = ':'.join([gene_id,iso_id,chr,strand]);
        iso_chr[chr].append([start_points[0],end_points[-1],strand, key,start_points,end_points])
    return (iso_chr)


def check_circ(isoform,junc,t,chr,genome_fa):

    junc_start=junc[0]
    junc_end=junc[1]
    strand=junc[2]
    start_points=isoform[4]
    end_points=isoform[5]
    start_exon=-1
    end_exon = -1
    circtype=None
    starts=[];
    ends=[];
    coord_flank_intron1='None';
    coord_flank_intron2='None'
    flank_intron='';
    reallignment=False
    start_index_intron=-2
    end_index_intron=-2

    for i in range (len(start_points)):
        if junc_start>=start_points[i]-t and junc_start<=start_points[i]+t:
            start_exon=i;


    for i in range (len(end_points)-1):
        if  start_exon==-1:
            if strand =='+' and junc_start >= end_points[i] - t and junc_start <= end_points[i] + t and junc_end <start_points[i+1]:

                starts = [junc_start]
                ends = [junc_end]
                if junc_start !=  end_points[i]:
                    if realign_ciRNA(junc_start, end_points[i], junc_end, genome_fa,chr):
                        reallignment = True
                        starts = [end_points[i]]
                        ends = [junc_end+end_points[i]-junc_start]

                start_exon=i;
                circtype='ciRNA'
                coord_flank_intron='-'.join([str(end_points[i]),str(start_points[i+1])]);
                flank_intron=':'.join([chr,coord_flank_intron]);


    for i in range (len(end_points)):
        if junc_end>=end_points[i]-t and junc_end<=end_points[i]+t:
            end_exon=i;


    for i in range (len(start_points)-1):
        if end_exon==-1:
            if strand =='-' and junc_end >= start_points[i+1] - t and junc_end <= start_points[i+1] + t and junc_start >end_points[i]:
                starts=[junc_start]
                ends=[junc_end]
                if junc_end !=  start_points[i+1]:
                    if realign_ciRNA(junc_end, start_points[i+1], junc_start, genome_fa,chr):
                        reallignment = True
                        starts = [junc_start+start_points[i+1]-junc_end]
                        ends = [start_points[i+1] ]
                start_exon = i;
                circtype = 'ciRNA'
                coord_flank_intron = '-'.join([str(end_points[i]), str(start_points[i + 1])]);

                flank_intron = ':'.join([chr, coord_flank_intron]);


    if start_exon>=0 and end_exon>=0:
        if realign_circRNA(junc_start,start_points[start_exon],junc_end,end_points[end_exon],genome_fa,chr):
            reallignment=True
            starts =start_points[start_exon:end_exon+1]
            ends = end_points[start_exon:end_exon+1]
        else:
            starts=[junc_start]+start_points[start_exon+1:end_exon+1]
            ends=end_points[start_exon:end_exon]+[junc_end]
        circtype = 'circRNA'
        if start_exon>0:
            coord_flank_intron1 = '-'.join([str(end_points[start_exon-1]), str(start_points[start_exon])]);
        if end_exon<len(end_points)-1:
            coord_flank_intron2 = '-'.join([str(end_points[end_exon]), str(start_points[end_exon + 1])]);
        intron='|'.join([coord_flank_intron1,coord_flank_intron2]);
        flank_intron = ':'.join([chr,intron ]);
    if circtype==None:
        if junc_start >= start_points[0] - t and junc_start <= start_points[0] + t:
            start_index_intron=-1;
        if junc_end >= end_points[len(start_points)-1] - t and junc_end <= end_points[len(start_points)-1] + t:
            end_index_intron=len(start_points)-1

        for i in range (len(start_points)-1):
            if junc_start>end_points[i] -t and junc_start<start_points[i+1]+t:
                start_index_intron=i
            if junc_end>end_points[i]-t and junc_end<start_points[i+1] +t:
                end_index_intron=i


        if start_index_intron>0 and end_index_intron>0:
            if end_index_intron>start_index_intron:
                starts=start_points[start_index_intron+1:end_index_intron+1]
                ends=end_points[start_index_intron+1:end_index_intron+1]
                circtype = 'ccRNA'
                if start_index_intron>-1:
                    coord_flank_intron1 = '-'.join([str(end_points[start_index_intron ]), str(start_points[start_index_intron+1])]);
                if end_index_intron<len(start_points)-1:
                    coord_flank_intron2 = '-'.join([str(end_points[end_index_intron]), str(start_points[end_index_intron + 1])]);
               
                intron = '|'.join([coord_flank_intron1, coord_flank_intron2]);
                flank_intron = ':'.join([chr, intron]);
                reallignment=False
                start_exon=start_index_intron+1;
                end_exon=end_index_intron;




    return (starts,ends,circtype,flank_intron,reallignment,start_exon,end_exon)

def realign_ciRNA(junc_start, exon_end, junc_end, genome_fa,chr):
    if junc_start < exon_end:
        seq1 = genome_fa.fetch(chr, junc_start, exon_end)
        seq2 = genome_fa.fetch(chr, junc_end, junc_end+exon_end-junc_start)
    else:
        seq1 = genome_fa.fetch(chr, exon_end, junc_start)
        seq2 = genome_fa.fetch(chr, junc_end-junc_start+exon_end, junc_end)
    return seq1.upper() == seq2.upper()



def realign_circRNA(junc_start,first_exon_start,junc_end,last_exon_end,genome_fa,chr):

    if junc_start-first_exon_start == junc_end-last_exon_end and junc_start-first_exon_start>0 :
        if junc_start<first_exon_start:
            seq1=genome_fa.fetch(chr, junc_start, first_exon_start)
            seq2 = genome_fa.fetch(chr, junc_end, last_exon_end)
        else:
            seq1 = genome_fa.fetch(chr, first_exon_start, junc_start)
            seq2 = genome_fa.fetch(chr, last_exon_end, junc_end)
        return seq1.upper() == seq2.upper()

    else:
        return False


def find_cirRNAS(ref, junc_file, genome_fa,outPath,t):
    (juctions_chr, junctions)=extract_jucnctions(junc_file)

    (iso_chr) = extract_genes(ref)
    circular_RNAs = {};
    output = open(outPath, 'w');
    for chr in iso_chr:
        for isoform in iso_chr[chr]:
            for junc in juctions_chr[chr]:

                if  junc[0]>= isoform[0]-t and junc[1]<=isoform[1]+t:
                    (starts, ends, circtype, flank_intron,reallignment,start_exon,end_exon)=check_circ(isoform,junc,t,chr,genome_fa)

                    if circtype!=None and circtype!='ccRNA':

                        strand = isoform[2]
                        num_exons = len(starts)
                        key = junc[3]
                        reads = junctions[key][0]
                        iso_key = isoform[3]
                        gene_id = iso_key.split(':')[0]
                        iso_id = iso_key.split(':')[1]
                        flag = '0'
                        exons_index=str(start_exon+1)
                        if reallignment:
                            flag = '1'
                        if circtype=='circRNA' or circtype=='ccRNA':
                            for i in range(end_exon-start_exon):
                                exons_index=exons_index+','+str(2+start_exon+i)

                        sizes=str(ends[0]-starts[0])+','
                        offsets = str(0)+','

                        for i in range(len(starts)-1):
                            sizes=sizes+str(ends[i+1]-starts[i+1])+','
                            offsets = offsets + str(starts[i + 1]-starts[0]) + ','

                        key = ':'.join([chr, strand, str(starts[0]), str(ends[-1]), sizes, offsets]);
                        donor=genome_fa.fetch(chr, ends[-1], ends[-1]+2)
                        acceptor=genome_fa.fetch(chr, starts[0]-2, starts[0])
                        if key not in circular_RNAs:
                            output.write(
                                chr + '\t' + str(starts[0]) + '\t' + str(ends[-1]) + '\t' + 'circular_RNA/' + str(
                                    reads) + '\t' + flag + '\t' + strand + '\t' + str(starts[0]) + '\t' +
                                str(starts[0]) + '\t' + '0,0,0' + '\t' + str(num_exons) + '\t' +
                                sizes + '\t' + offsets + '\t' + str(
                                    reads) + '\t' + circtype + '\t' + gene_id + '\t' + iso_id + '\t' + exons_index + '\t' + flank_intron + '\t' + donor+'-'+acceptor + '\t' + '\n');

                            value = [1];
                            circular_RNAs[key] = value;



def main():

    args = parser.parse_args()
    outDir = args.output;
    ref = args.ref_file;
    junc_file = args.junction;
    genome = args.genome;
    genome_fa = pysam.FastaFile(genome)
    outPath = os.path.abspath(outDir);  ## absolute output path

    find_cirRNAS(ref, junc_file, genome_fa, outPath, 10)


if __name__ == '__main__':
    main()
