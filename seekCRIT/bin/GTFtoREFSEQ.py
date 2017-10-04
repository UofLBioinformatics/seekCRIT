import sys
from collections import defaultdict
import os


def main():

    rFile = open(sys.argv[1], 'r');
    outPath = os.path.abspath(sys.argv[2]);
    output = open(outPath, 'w');
    exons = defaultdict(list)
    CDS = defaultdict(list)

    for line in rFile:  ## process each line
        if line[0]!='#':
            ele = line.strip().split('\t');
            if len(ele)>7 :

                if ele[2]=='exon':
                    id=ele[8]
                    chr = (ele[0])
                    strand = (ele[6])
                    start = int(ele[3])
                    end = int(ele[4])
                    id_split = id.strip().split(';');
                    for j in range(len(id_split)):
                        part=id_split[j]
                        part_id = part.strip().split(' ');
                        if part_id[0]=='gene_id':
                            gene_id = part_id[1].replace('"','');
                        elif part_id[0] == 'transcript_id':
                            iso_id = part_id[1].replace('"', '');
                        elif  part_id[0] == 'gene_name':
                            gene_id = part_id[1].replace('"', '');

                    key = ':'.join([gene_id,iso_id,chr,strand]);
                    exons[key].append([start,end])
                if ele[2] == 'CDS':
                    id = ele[8]
                    chr = (ele[0])
                    strand = (ele[6])
                    start = int(ele[3])
                    end = int(ele[4])
                    id_split = id.strip().split(';');
                    for j in range(len(id_split)):
                        part = id_split[j]
                        part_id = part.strip().split(' ');
                        if part_id[0] == 'gene_id':
                            gene_id = part_id[1].replace('"', '');
                        elif part_id[0] == 'transcript_id':
                            iso_id = part_id[1].replace('"', '');
                        elif part_id[0] == 'gene_name':
                            gene_id = part_id[1].replace('"', '');
                    key = ':'.join([gene_id, iso_id, chr, strand]);
                    CDS[key].append([start, end])

    for key in exons:

        start_points=[];
        end_points=[]
        start_points_CDS = [];
        end_points_CDS = []
        list_exons=exons[key]
        exons_start='';
        exons_end =  '';
        for i in range(len(list_exons)):
            start_points.append(list_exons[i][0]-1)
            end_points.append(list_exons[i][1])
        start_points=sorted(start_points)
        end_points=sorted(end_points)
        for i in range(len(list_exons)):
            exons_start=exons_start+str(start_points[i])+',';
            exons_end=exons_end+str(end_points[i])+',';

        list_CDS = CDS[key]
        if list_CDS==[]:
            list_CDS=[[start_points[0],end_points[-1]]]
        for i in range(len(list_CDS)):
            start_points_CDS.append(list_CDS[i][0]-1)
            end_points_CDS.append(list_CDS[i][1])
        start_points_CDS = sorted(start_points_CDS)
        end_points_CDS = sorted(end_points_CDS)


        gene_id = key.strip().split(':')[0]
        iso_id = key.strip().split(':')[1]
        chr = key.strip().split(':')[2]
        strand = key.strip().split(':')[3]

        output.write(gene_id + '\t' +iso_id+'\t'+ chr+'\t'+strand+'\t'+ str(start_points[0]) + '\t' + str(end_points[-1]) + '\t'+
                     str(start_points_CDS[0]) +
                 '\t' + str(end_points_CDS[-1])+ '\t' +  str(len(list_exons)) + '\t' + exons_start + '\t' +exons_end+'\t'+ '\n');



if __name__ == '__main__':
    main()


