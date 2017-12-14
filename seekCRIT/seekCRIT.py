#
## this program indentifies and characterizes differentially spliced circRNAs 
#

### import necessary libraries
import re,os,sys,logging,time,datetime,scipy,numpy,argparse;
#import random, 
import subprocess;
#import fisher,mne;  ## for p-value and FDR calculation
import pysam; ## use pysam package to access bam/sam files
import shutil;   ## to find a path of executables
from distutils.version import LooseVersion

seekCRIT_ver = "0.0.2";

### checking out the number of arguments
parser = argparse.ArgumentParser(description='Identifying and Characterizing Differentially Spliced circular RNAs between two samples');

parser.add_argument('-s1', '--sample1', dest='s1', required=True, help='fastq files for sample_1. Replicates are separated by comma. Paired-end reads are separated by colon. e.g.,s1-1.fastq,s1-2.fastq for single-end read. s1-1.R1.fastq:s1-1.R2.fastq,s1-2.R1.fastq:s1-2.R2.fastq for single-end read');
parser.add_argument('-s2', '--sample2', dest='s2', required=True, help='fastq files for sample_2. Replicates are separated by comma. Paired-end reads are separated by colon. e.g.,s2-1.fastq,s2-2.fastq for single-end read. s2-1.R1.fastq:s2-1.R2.fastq,s2-2.R1.fastq:s2-2.R2.fastq for single-end read');
parser.add_argument('-gtf', '--gtf', dest='gtf', required=True, help='The gtf annotation file. e.g., hg38.gtf');
parser.add_argument('-o', '--output', dest='outDir', required=True, help='Output directory');
parser.add_argument('-t', '--readType', dest='readType', required=True, choices=['SE','PE'], help='Read type. SE for Single-end read, PE for Paired-end read');
#parser.add_argument('-len', '--readLength', dest='readLength', required=True, type=int, choices=range(1,10000), help='Read length. Positive integer');
parser.add_argument('--aligner', type=str, dest='aligner', required=True, help='Aligner to use')
parser.add_argument('--genomeIndex', type=str, dest='genomeIndex', required=True, help='Genome indexes for the aligner')
parser.add_argument('-fa', '--fasta', type=str, dest='fasta', required=True, help='Genome sequence. e.g., hg38.fa')

parser.add_argument('-ref', '--refseq', type=str, dest='refseq', default=None, help='Transcriptome in refseq format. e.g., hg38.ref.txt')
parser.add_argument('--threadNumber', type=int, dest='threadN', default=4, choices=range(1,100), help='Number of threads for multi-threading feature. Positive integer')
parser.add_argument('--deltaPSI', type=float, dest='deltaPSI', default=0.05, help='Delta PSI cutoff. i.e., significant event must show bigger deltaPSI than this cutoff')
parser.add_argument('--highConfidence', type=int, dest='highConfidence', default=1, help='Minimum number of circular junction counts required')
parser.add_argument('--libType', type=str, dest='libType', default='fr-unstranded', choices=['fr-unstranded','fr-firststrand', 'fr-secondstrand'], help='Minimum number of circular junction counts required')
parser.add_argument('--keepTemp', type=str, dest='keepTemp', default='Y', choices=['Y','N'], help='Keep temp files or not')

args = parser.parse_args()

s1=args.s1; s2=args.s2;
gtf = args.gtf;
outDir = args.outDir;
readType = args.readType;
#readLength= args.readLength;
aligner = args.aligner;
genomeIndex = args.genomeIndex;
threadN = args.threadN;
deltaPSI = args.deltaPSI;
highConfidence = args.highConfidence;
libType = args.libType;
keepTemp = args.keepTemp;
fasta = args.fasta;
refseq = args.refseq;


def listToString(x):
  rVal = '';
  for a in x:
    rVal += a+' ';
  return rVal;

def uniq(inlist):
    # order preserving
    uniques = []
    for item in inlist:
        if item not in uniques:
            uniques.append(item)
    return uniques

#
### fastq files or bam file
#
sample_1=s1.split(','); ## each end of a pair is separated by :
sample_2=s2.split(','); ## each end of a pair is separated by :
#
SEPE = readType; ## single-end or paired
#

#
##### checking GTF format ##### it does the minimal checking for now.
#
tempGTF_file = open(gtf); ## open gtf file
for line in tempGTF_file:
  if line.strip()[0]=='#': ## comments, skip this line
    continue;  
  gtfEle = line.strip().split('\t');
  if len(gtfEle)<9: ## may be incorrect gtf format    
    print ("Incorrect GTF file format. Non-comment lines in GTF file must have 9 tab-delimited columns.");
    sys.exit();
  break;  ## just check the first non-comment column
###

os.system('mkdir -p '+ outDir);
oFile = open(outDir+'/commands.txt', 'a'); ## file that will contain list of commands excuted here

### setting up the logging format
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=outDir+'/log.seekCRIT.'+  str(datetime.datetime.now()) + '.txt',
                    filemode='w')

##### Getting Start Time ######
logging.debug('seekCRIT version: %s' % seekCRIT_ver);
logging.debug('Start the program with [%s]\n', listToString(sys.argv));
startTime = time.time();

pythonPath = os.environ['_'];  ## pythonPath that will be used in the remainder of the program

scriptPath = os.path.abspath(os.path.dirname(__file__));  ## absolute script path
binPath = scriptPath + '/bin';  ## absolute bin path
outPath = os.path.abspath(outDir); ## absolute output path

s1Path = outPath + '/SAMPLE_1';
os.system('mkdir -p '+ s1Path);
s2Path = outPath + '/SAMPLE_2';
os.system('mkdir -p '+ s2Path);

## making folders for replicates ##
s1rPath = s1Path+'/REP_';
s2rPath = s2Path+'/REP_';
for rr in range(0,len(sample_1)): ## sample_1
  os.system('mkdir -p '+ s1rPath+str(rr+1));
for rr in range(0,len(sample_2)): ## sample_2
  os.system('mkdir -p '+ s2rPath+str(rr+1));

finalPath = outPath+'/seekCRIT_output';  ## absolute seekCRIT result path
os.system('mkdir -p '+ finalPath);

tempPath = outPath + '/temp';
os.system('mkdir -p '+ tempPath);  ## absolute path for temp results

#
### putting keys in log file
#
logging.debug("################### folder names and associated input files #############");
for fki in range(0,len(sample_1)): ## for each replicate of sample_1
  repTempFolder = "SAMPLE_1\REP_"+str(fki+1);
  associatedFile = sample_1[fki];
  logging.debug(repTempFolder+"\t"+associatedFile);

for fki in range(0,len(sample_2)): ## for each replicate of sample_2
  repTempFolder = "SAMPLE_2\REP_"+str(fki+1);
  associatedFile = sample_2[fki];
  logging.debug(repTempFolder+"\t"+associatedFile);

logging.debug("#########################################################################\n");

if refseq==None:
    logging.debug("converting gtf file to RefSeq file ");
    cmd = pythonPath + ' bin/GTFtoREFSEQ.py '+ gtf +' refseq.txt'
    oFile.write('######    converting gtf file to RefSeq file  #####\n'+cmd+'\n#\n');
    oFile.flush();
    status,output=subprocess.getstatusoutput(cmd);
    logging.debug("converting gtf file to RefSeq file is done with status %s" % status);
    if (int(status)!=0): ## it did not go well
      logging.debug("error in converting gtf file to RefSeq file: %s" % (status));
      logging.debug("error detail: %s" % output);
      raise Exception();
    logging.debug(output);
    refseq='refseq.txt'



########## functions here... ############

def doSTARMapping(): ## do STAR mapping
  logging.debug("mapping the first sample");

  for rr in range(0,len(sample_1)): ## for each replicate of sample_1
    rTempFolder = s1rPath+str(rr+1);
    cmd = 'STAR --chimSegmentMin 2 --outFilterMismatchNmax 3 --runThreadN '+str(threadN)+' --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate ';
    cmd += ' --genomeDir '+genomeIndex+ ' --sjdbGTFfile ' + gtf;
    cmd += ' --outFileNamePrefix ' + rTempFolder + '/ --readFilesIn ';
    if SEPE=='PE': ## paired-end
      cmd += sample_1[rr].split(':')[0]+' '+sample_1[rr].split(':')[1];
    else: ## single-end
      cmd += sample_1[rr];
    oFile.write('######  running STAR for sample_1, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n');
    oFile.flush();
    status,output=subprocess.getstatusoutput(cmd);
    logging.debug("mapping sample_1, rep_"+str(rr+1)+" is done with status %s" % status);
    if (int(status)!=0): ## it did not go well
      logging.debug("error in mapping sample_1, rep_%d: %s" % ((rr+1),status));
      logging.debug("error detail: %s" % output);
      raise Exception();
    logging.debug(output);


  logging.debug("mapping the second sample");
  for rr in range(0,len(sample_2)): ## for each replicate of sample_2
    rTempFolder = s2rPath+str(rr+1);
    cmd = 'STAR --chimSegmentMin 2 --outFilterMismatchNmax 3 --runThreadN '+str(threadN)+' --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate ';
    cmd += ' --genomeDir '+genomeIndex+ ' --sjdbGTFfile ' + gtf;
    cmd += ' --outFileNamePrefix ' + rTempFolder + '/ --readFilesIn ';
    if SEPE=='PE': ## paired-end
      cmd += sample_2[rr].split(':')[0]+' '+sample_2[rr].split(':')[1];
    else: ## single-end
      cmd += sample_2[rr];
    oFile.write('######  running STAR for sample_2, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n');
    oFile.flush();
    status,output=subprocess.getstatusoutput(cmd);
    logging.debug("mapping sample_2, rep_"+str(rr+1)+" is done with status %s" % status);
    if (int(status)!=0): ## it did not go well
      logging.debug("error in mapping sample_2, rep_%d: %s" % ((rr+1),status));
      logging.debug("error detail: %s" % output);
      raise Exception();
    logging.debug(output);

  return;
##### end of doSTARMapping ####


#def doTophatMapping(): ## do tophat mapping, NOT USED FOR NOW
#  logging.debug("mapping the first sample");
#
#  for rr in range(0,len(sample_1)): ## for each replicate of sample_1
#    rTempFolder = s1rPath+str(rr+1);
#    cmd = 'tophat -a '+str(tophatAnchor)+' -m 0 -I 300000 -p 4 -g 20 --library-type ' + libType + ' --no-novel-indels ';
#    cmd += ' -N 3 --segment-mismatches 2 -G '+gtf+' -o '+rTempFolder;
#    if SEPE=='PE': ## paired-end
#      cmd += ' -r '+str(insertLength[rr]) + ' --mate-std-dev ' + str(sigma[rr])+' ' + bIndex +' '+sample_1[rr].split(':')[0]+' '+sample_1[rr].split(':')[1];
#    else: ## single-end
#      cmd += ' ' + bIndex +' '+sample_1[rr];
#    oFile.write('######  running tophat for sample_1, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n');
#    oFile.flush();
#    status,output=subprocess.getstatusoutput(cmd);
#    logging.debug("mapping sample_1, rep_"+str(rr+1)+" is done with status %s" % status);
#    if (int(status)!=0): ## it did not go well
#      logging.debug("error in mapping sample_1, rep_%d: %s" % ((rr+1),status));
#      logging.debug("error detail: %s" % output);
#      raise Exception();
#    logging.debug(output);
#
#  if lite==1:
#    return;
#
#  logging.debug("mapping the second sample");
#  for rr in range(0,len(sample_2)): ## for each replicate of sample_2
#    rTempFolder = s2rPath+str(rr+1);
#    cmd = 'tophat -a '+str(tophatAnchor)+' -m 0 -I 300000 -p 4 -g 20 --library-type '+ libType +' --no-novel-indels ';
#    cmd += ' -N 3 --segment-mismatches 2 -G '+gtf+' -o '+rTempFolder;
#    if SEPE=='PE': ## paired-end
#      cmd += ' -r '+str(insertLength2[rr]) + ' --mate-std-dev ' + str(sigma2[rr])+' ' + bIndex +' '+sample_2[rr].split(':')[0]+' '+sample_2[rr].split(':')[1];
#    else: ## single-end
#      cmd += ' ' + bIndex +' '+sample_2[rr];
#    oFile.write('######  running tophat for sample_2, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n');
#    oFile.flush();
#    status,output=subprocess.getstatusoutput(cmd);
#    logging.debug("mapping sample_2, rep_"+str(rr+1)+" is done with status %s" % status);
#    if (int(status)!=0): ## it did not go well
#      logging.debug("error in mapping sample_2, rep_%d: %s" % ((rr+1),status));
#      logging.debug("error detail: %s" % output);
#      raise Exception();
#    logging.debug(output);
#
#  return;
##### end of doTophatMapping ####

def indexBamFile(): ## indexing bam files to use pysam
  logging.debug("indexing BAM File function..");
  bamFile=0; ## currently not supporting bam file input
  for rr in range(0,len(sample_1)): ## for each replicate of sample_1
    rTempFolder = s1rPath+str(rr+1);
    bam_fn='';
    if bamFile==0:  ## we know the location of the bam file
      bam_fn = rTempFolder+'/Aligned.sortedByCoord.out.bam';
    else: ## bam file is provided
      bam_fn = sample_1[rr];

    if LooseVersion(pysam.version.__samtools_version__) < LooseVersion('1.3'):
      pysam.sort(bam_fn, rTempFolder+'/aligned.sorted'); ## it will make aligned.sorted.bam file
      pysam.index(rTempFolder+'/aligned.sorted.bam'); ## it will make aligned.sorted.bam.bai file
    else:
      pysam.sort(bam_fn, '-o', rTempFolder+'/aligned.sorted.bam'); ## it will make aligned.sorted.bam file
      pysam.index(rTempFolder+'/aligned.sorted.bam'); ## it will make aligned.sorted.bam.bai file

  for rr in range(0,len(sample_2)): ## for each replicate of sample_2
    rTempFolder = s2rPath+str(rr+1);
    bam_fn='';
    if bamFile==0:  ## we know the location of the bam file
      bam_fn = rTempFolder+'/Aligned.sortedByCoord.out.bam';
    else: ## bam file is provided
      bam_fn = sample_2[rr];

    if LooseVersion(pysam.version.__samtools_version__) < LooseVersion('1.3'):
      pysam.sort(bam_fn, rTempFolder+'/aligned.sorted'); ## it will make aligned.sorted.bam file
      pysam.index(rTempFolder+'/aligned.sorted.bam'); ## it will make aligned.sorted.bam.bai file
    else:
      pysam.sort(bam_fn, '-o', rTempFolder+'/aligned.sorted.bam'); ## it will make aligned.sorted.bam file
      pysam.index(rTempFolder+'/aligned.sorted.bam'); ## it will make aligned.sorted.bam.bai file
### end of indexBamFile() function ###


def detectCircRNAs():  ## detecting circular RNAs from chimeric junctions

  
  logging.debug("detecting circular RNAS for sample_1");
  for rr in range(0,len(sample_1)): ## for each replicate of sample_1
    rTempFolder = s1rPath+str(rr+1);
    cjn = rTempFolder+'/Chimeric.out.junction'; ## chimeric junction name from STAR aligner
   
    cmd = pythonPath +' bin/circ_detection.py -j '+cjn+' -g '+fasta+' -r '+refseq+' -o '+ rTempFolder+'/circ.output.txt'; 

    oFile.write('######  detecting circular RNAs for sample_1, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n');
    oFile.flush();
    status,output=subprocess.getstatusoutput(cmd);
    logging.debug("detecting circular RNAs for sample_1, rep_"+str(rr+1)+" is done with status %s" % status);
    if (int(status)!=0): ## it did not go well
      logging.debug("error in detecting circular RNAs for sample_1, rep_%d: %s" % ((rr+1),status));
      logging.debug("error detail: %s" % output);
      raise Exception();
    logging.debug(output);

  logging.debug("detecting circular RNAs for sample_2");
  for rr in range(0,len(sample_2)): ## for each replicate of sample_2
    rTempFolder = s2rPath+str(rr+1);
    cjn = rTempFolder+'/Chimeric.out.junction'; ## chimeric junction name from STAR aligner


    cmd = pythonPath +' bin/circ_detection.py -j '+cjn+' -g '+fasta+' -r '+refseq+' -o '+ rTempFolder+'/circ.output.txt'; 

    oFile.write('######  detecting circular RNAs for sample_2, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n');
    oFile.flush();
    status,output=subprocess.getstatusoutput(cmd);
    logging.debug("detecting circular RNAs for sample_2, rep_"+str(rr+1)+" is done with status %s" % status);
    if (int(status)!=0): ## it did not go well
      logging.debug("error in detecting circular RNAs for sample_2, rep_%d: %s" % ((rr+1),status));
      logging.debug("error detail: %s" % output);
      raise Exception();
    logging.debug(output);

### end of detectCircRNAs() function ###
	
	
	
	
	
def processCircRNAs():  ## processing circular RNAs from two sample groups

  logging.debug("processing circular RNAs from two sample groups");
  
  rTempFolder = s1rPath+str(1)+'/circ.output.txt';
  circ_1=rTempFolder;
  for rr in range(1,len(sample_1)):
    rTempFolder = s1rPath+str(rr+1);
    cir = rTempFolder+'/circ.output.txt'
    circ_1=circ_1+','+cir;
	
  rTempFolder = s2rPath+str(1)+'/circ.output.txt';
  circ_2=rTempFolder;
  for rr in range(1,len(sample_2)):
    rTempFolder = s2rPath+str(rr+1);
    cir = rTempFolder+'/circ.output.txt'
    circ_2=circ_2+','+cir;

	
  rTempFolder = s1rPath+str(1)+'/Aligned.sortedByCoord.out.bam';
  bam_1=rTempFolder;
  for rr in range(1,len(sample_1)):
    rTempFolder = s1rPath+str(rr+1);
    bam = rTempFolder+'/Aligned.sortedByCoord.out.bam'
    bam_1=bam_1+','+bam;
	
  rTempFolder = s2rPath+str(1)+'/Aligned.sortedByCoord.out.bam';
  bam_2=rTempFolder;
  for rr in range(1,len(sample_2)):
    rTempFolder = s2rPath+str(rr+1);
    bam = rTempFolder+'/Aligned.sortedByCoord.out.bam'
    bam_2=bam_2+','+bam;
	
  cmd= pythonPath + ' bin/processCIRC.BAM.py '+ circ_1 +' '+circ_2+' '+bam_1+' '+bam_2+' '+str(highConfidence)+' '+finalPath+' '+SEPE;
  oFile.write('######  processing circular RNAs   #####\n'+cmd+'\n#\n');
  oFile.flush();
  status,output=subprocess.getstatusoutput(cmd);
  logging.debug(" processing circular RNAs is done with status %s" % status);
  if (int(status)!=0): ## it did not go well
    logging.debug("error in processing circular RNAs: %s" % (status));
    logging.debug("error detail: %s" % output);
    raise Exception();
    logging.debug(output);  
  

### end of processCircRNAs() function ###





############################################ main process ###############################################################
def main():
  ####
  #### 1. STAR mapping
  ####
  logging.debug("start mapping..")
  try:
    doSTARMapping();
    pass;
  except:
    logging.debug("There is an exception in mapping");
    logging.debug("Exception: %s" % sys.exc_info()[0]);
    logging.debug("Detail: %s" % sys.exc_info()[1]);
    sys.exit(-1);
  logging.debug("done mapping..");

  ####
  #### 2. index bam files
  ####
  logging.debug("indexing bam files to use pysam");

  try:
    indexBamFile();
    pass;
  except:
    logging.debug("There is an exception in indexing bam files");
    logging.debug("Exception: %s" % sys.exc_info()[0]);
    logging.debug("Detail: %s" % sys.exc_info()[1]);
    sys.exit(-2);
  logging.debug("done indexing bam files..");

  ####
  #### 3. detect circular RNAs from Chimeric output junctions
  ####
  logging.debug("detecitng circRNAs from chimeric junctions");

  try:

    detectCircRNAs();
    pass;
  except:
    logging.debug("There is an exception in detecting circRNAs");
    logging.debug("Exception: %s" % sys.exc_info()[0]);
    logging.debug("Detail: %s" % sys.exc_info()[1]);
    sys.exit(-2);
  logging.debug("done detecting circRNAs..");


  ###
  ### 4. processing circRNAs
  ###
  logging.debug("processing circRNAs from two sample groups");

  try:
    processCircRNAs();
    pass;
  except:
    logging.debug("There is an exception in processing circRNAs");
    logging.debug("Exception: %s" % sys.exc_info()[0]);
    logging.debug("Detail: %s" % sys.exc_info()[1]);
    sys.exit(-2);
  logging.debug("done processing circRNAs..");



  #############
  ## calculate total running time
  #############
  logging.debug("Program ended");
  currentTime = time.time();
  runningTime = currentTime-startTime; ## in seconds
  logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));

  sys.exit(0);


if __name__ == '__main__':
    main()
