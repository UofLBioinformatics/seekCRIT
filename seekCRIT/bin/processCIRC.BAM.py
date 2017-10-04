#
## this program processes CIRCexplorer output and BAM files
## to detect circular counts and linear counts
#

### import necessary libraries
import re,os,sys,logging,time,datetime,scipy,numpy;
import fisher,mne;  ## for p-value and FDR calculation
import pysam; ## use pysam package to access bam/sam files

myVer = "0.0.1";

### checking out the number of arguments
if (len(sys.argv)<8): 
  print('Not enough arguments!!');
  print ('It takes at least 7 arguments.');
  print ('Usage:\n\tProgramName CIRCexplorer_1_1[,CIRCexplorer_1_2]* CIRCexplorer_2_1[,CIRCexplorer_2_2]* BAM_2_1[,BAM_1_2]* BAM_2_1[,BAM_2_2]* circCountCutoff outFolder SEPE');
  print ('Example\n\tProgramName SAMPLE_1_circ.txt SAMPLE_2_circ.txt SAMPLE_1.bam SAMPLE_2.bam 3 CRIT_output SE');
  print ('NOTE:\n\tIt is highly recommended to put bam index files (bai) with bam files in the same folder for a significant speed-up.');
  sys.exit();

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

###
cfs_1 = sys.argv[1].split(','); ## CIRCexplorer output files for sample 1
cfs_2 = sys.argv[2].split(','); ## CIRCexplorer output files for sample 2
bams_1 = sys.argv[3].split(','); ## sam files for sample 1
bams_2 = sys.argv[4].split(','); ## sam files for sample 2
highConf = int(sys.argv[5]); ## cutoff for the high confidence circular RNAs
oDir = sys.argv[6]; ## output file
SEPE=sys.argv[7]; ## single-end or paired-end

#print kf_1,kf_2, sf_1, sf_2, cutoff, oDir

if not os.path.isdir(oDir):
  os.mkdir(oDir);



### setting up the logging format
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=oDir+'/log.processCIRC.SAM.'+  str(datetime.datetime.now()) + '.txt',
                    filemode='w')

##### Getting Start Time ######
logging.debug('Start the program with [%s]\n', listToString(sys.argv));
startTime = time.time();

#
####### input verification
#
## check # of replicates
#
numC1 = len(cfs_1);
numC2 = len(cfs_2);
numB1 = len(bams_1);
numB2 = len(bams_2); 
#
if (numC1 != numB1) : ## different number of replicates for the sample 1
  logging.debug("The number of replicates in the sample 1 does not match. There are %d replicates in CIRCexplorer output while there are %d replicates in bam files" % (numC1,numB1));
  logging.debug("Please check your input parameters.");
  sys.exit(-101);
if (numC2 != numB2) : ## different number of replicates for the sample 2
  logging.debug("The number of replicates in the sample 2 does not match. There are %d replicates in CIRCexplorer output while there are %d replicates in bam files" % (numC2,numB2));
  logging.debug("Please check your input parameters.");
  sys.exit(-102);
#
logging.debug("Number of replicates, Sample_1: %d" % numC1);
logging.debug("Number of replicates, Sample_2: %d" % numC2);
#
### process CIRCexplorer output files to come up with a combined list
#
S1=0; ## constant for sample_1
S2=1; ## constant for sample_2
#
c_dic=[{},{}]; ## dictiionary for circ explorer files. c_dic[0] for sample_1, c_dic[1] for sample_2
##
for i in range(numC1): ## init circDictionary for sample_1
  c_dic[S1][i]={};
for i in range(numC2): ## init circDictionary for sample_2
  c_dic[S2][i]={};
##
combined_circ_d = {};
wholeLine={};
rind=0;
for c1 in cfs_1: ### each CIRCexplorer output file 
  iFile = open(c1); ## open circ file for each replicate
  header = iFile.readline(); ## header, skip it for now
  for line in iFile: ## process each line
    ele = line.strip().split('\t');
    key = ele[0]+'::'+ele[1]+'::'+ele[2];
    jCount=int(ele[12]); ## circular junction count
    value = ele[17].replace(':','_').replace('-','_').split('|'); ### ['chrX_10529328_10532075', 'chrX_10571517_10582168']
    if key in combined_circ_d:
      if combined_circ_d[key] != value: ### same circRNA but different intron. Can it be possible? If so, how can we detect this?
        logging.debug("Same circRNA but different flanking introns: %s, %s, %s" % (key,value,combined_circ_d[key]));
    else:
      combined_circ_d[key]=value; ## put flanking exon information
      wholeLine[key]=line.strip(); ## put whole line here
    if key in c_dic[S1][rind]: ## it should not happend
      logging.debug("Duplicated key in sample_1, replicate_%d: %s" % (rind,key));
    else: ## unique circRNA for this replicate
      c_dic[S1][rind][key]=jCount;
  iFile.close();
  rind+=1; ## increase replicate number
#
rind=0;
for c2 in cfs_2: ### each CIRCexplorer output file 
  iFile = open(c2); ## open circ file for each replicate
  header = iFile.readline(); ## header, skip it for now
  for line in iFile: ## process each line
    ele = line.strip().split('\t');
    key = ele[0]+'::'+ele[1]+'::'+ele[2];
    jCount=int(ele[12]); ## circular junction count
    value = ele[17].replace(':','_').replace('-','_').split('|'); ### ['chrX_10529328_10532075', 'chrX_10571517_10582168']
    if key in combined_circ_d:
      if combined_circ_d[key] != value: ### same circRNA but different intron. Can it be possible? How can we detect this?
        logging.debug("Same circRNA but different flanking introns: %s, %s, %s" % (key,value,combined_circ_d[key]));
    else:
      combined_circ_d[key]=value; ## put flanking exon information
      wholeLine[key]=line.strip(); ## put whole line here
    if key in c_dic[S2][rind]: ## it should not happend
      logging.debug("Duplicated key in sample_2, replicate_%d: %s" % (rind,key));
    else: ## unique circRNA for this replicate
      c_dic[S2][rind][key]=jCount;
  iFile.close();
  rind+=1;
#
for i in range(numC1):
  logging.debug("Number of unique circRNAs in sample_1, repplicate_%d: %d" % (i+1,len(c_dic[S1][i])));
for i in range(numC2):
  logging.debug("Number of unique circRNAs in sample_2, repplicate_%d: %d" % (i+1,len(c_dic[S2][i])));
logging.debug("Done combining circRNAs from all smaples and all replicates. Number of combined circRNAs: %d" % len(combined_circ_d));
##
#
#### now process BAM files to assign linear counts to all of the circRNAs ####
#
def is_unique(r):
  # check if read is uniquely mapped
  for tag in r.tags:
    if tag[0]=='NH': # NH is tag for number of hits
      if int(tag[1])==1: # uniquely mapped if NH=1
        if SEPE=='SE': ## single end, sufficient
          return True;
        elif r.is_proper_pair:
          return True;
  return False
#

### dictionary for linear count per replicate ###
lcs1={}; lcs2={}; ## linear count for sample_1 and count for sample_2

## init lcs dictionaries
for i in range(numC1):
  lcs1[i]={}; ## for each replicate
for i in range(numC2):
  lcs2[i]={}; ## for each replicate

tJunctions={};
junctions={};
#chunk = 10000; ## to speed up the search
bamIndex=0;


def getLinearCounts(bams, lcs):
  i=0; ## replicate number
  for s1 in bams: ## for each bam file
    if len(s1.strip())<1: ## incorrect split. a user might accidently put comma at the end of input bam file list
      continue; ### just skip this entry, probably the last one though
    logging.debug("processing %s" % s1.strip());
    sFile = pysam.AlignmentFile(s1.strip(),'rb'); ## open bam file
    for read in sFile.fetch(until_eof=True):
    #for read in sFile.fetch("chr1"):
      if not is_unique(read): ## not unique, go to the next one
        continue; ## process next read
  
      chr = sFile.getrname(read.tid);
      mc = read.pos+1;
      mString = read.cigarstring; ## mapping string, 50M or aMbNcMdNeMfNgM format, CIGAR string
      if 'D' in mString or 'I' in mString or 'S' in mString or 'H' in mString or 'P' in mString or 'X' in mString or '=' in mString: ## skip
        continue; ## go to next line
  
      ### check to see if the line is either exonic read or junction read
      split_mString = mString.split('M');
      if len(split_mString)==2:  ## exonic read
        continue; ## go to next line
      elif len(split_mString)>=3:    ###### junction read ###########
        jS=mc; jE=mc-1;
        for ec in range(0,len(split_mString)-2): ## for each coordinate
          secondNumber = int(split_mString[ec].split('N')[-1]);
          jS = jE+secondNumber; ## 1-base
          jE = jS+int(split_mString[ec+1].split('N')[0]); ## 0-base
          key = chr+'_'+str(jS)+'_'+str(jE);
  
          ## edge count
          if key in lcs[i]: ## already have it
            lcs[i][key]+=1;
          else: ## new edge
            lcs[i][key]=1;
  
    logging.debug("Done populating edeg count for %s with %d junctions" % (s1.strip(), len(lcs[i])));
    sFile.close();
    edgeFile = open(s1.strip()+'.edgeCount', 'w');
    for k in lcs[i]:
      edgeFile.write(k+'\t'+str(lcs[i][k])+'\n');
    edgeFile.close();
  
    i+=1; ## increase replicate number
  logging.debug("Done running getLinearCounts function");
################## end of getLinearCunts function #####################

fp = {}; ## fisher p-value
sortedKey=[];
counts_d = {};
##### p-value and FDR ####
def computeFisherExact(): ## compute fisher's exact test, use sum for replicates
  for k in combined_circ_d: ## for all unique circRNAs
    cc_1=[]; cc_2=[]; lc_1=[]; lc_2=[];
    ufi = combined_circ_d[k][0]; dfi = "noKey";
    if len(combined_circ_d[k])==2: ## it has two introns listed
      dfi=combined_circ_d[k][1];  ## upstream and downstream flanking introns
 
    ## get cc_1 and lc_1
    for i in range(numC1): ## for each replicate
      cVal=0; lVal=0; 
      if k in c_dic[S1][i]:
        cVal = c_dic[S1][i][k];
        if ufi in lcs1[i]: ## upstream junction count exists
          lVal = lcs1[i][ufi];
        if dfi in lcs1[i]: ## downstream junction count exists
          lVal += lcs1[i][dfi];
      cc_1.append(cVal); lc_1.append(lVal); 
   
    ## get cc_2 and lc_2
    for i in range(numC2): ## for each replicate
      cVal=0; lVal=0;
      if k in c_dic[S2][i]:
        cVal = c_dic[S2][i][k];
        if ufi in lcs2[i]: ## upstream junction count exists
          lVal = lcs2[i][ufi];
        if dfi in lcs2[i]: ## downstream junction count exists
          lVal += lcs2[i][dfi];

      cc_2.append(cVal); lc_2.append(lVal); 

    counts_d[k]=[cc_1,lc_1,cc_2,lc_2]; 

    ##print (cc_1, cc_2, lc_1, lc_2);
    n1=2*sum(cc_1);n2=sum(lc_1);n3=2*sum(cc_2);n4=sum(lc_2);
     
    p = fisher.pvalue(n1,n2,n3,n4);
    fp[k] = p.two_tail; ## saving p-value. it has p.left_tail, p.right_tail, and p.two_tail values

  for k in sorted(fp):
    sortedKey.append(k); ## sort keys for fdr calculation
  logging.debug("Done computing two-tail fisher exact test");


def computeFisherMethod(): ### combined p-value for replicates. not done yet.
  pass;
### end of computeFisherMethod function 


fdr={};
def computeFDR(): ### compute FDR (Benjamini-Hochberg)
  pValList=[];
  for k in sortedKey: ## going through the sorted key
    pValList.append(fp[k]); ## appending p value
  myFDR=list(mne.stats.fdr_correction(pValList)[1]);

  for i in range(len(sortedKey)): ## assign fdr values
    fdr[sortedKey[i]]=myFDR[i];
  logging.debug("Done computing B-H fdr");
    
### end of computeFDR function


def writeResult(): ## write out result
  logging.debug("Writing output..");
  ## sort fp by p-value
  rValue = sorted(fp.items(),key=lambda x:x[1]); ## tuple, sorted by p-value
  outFile = open(oDir + '/circRNAs.pVal.FDR.txt','w');
  header = 'chrom\tcircRNA_start\tcircRNA_end\tstrand\texonCount\texonSizes\texonOffsets\tcircType\tgeneName\tisoformName\texonIndexOrIntronIndex\tFlankingIntrons';
  header +='\tCircularJunctionCount_Sample_1\tLinearJunctionCount_Sample_1\tCircularJunctionCount_Sample_2\tLinearJunctionCount_Sample_2\tPBI_Sample_1\tPBI_Sample_2';
  header +='\tdeltaPBI(PBI_1-PBI_2)\tpValue\tFDR';
  outFile.write(header+'\n'); 

  for vv in rValue:
    outKey=vv[0]; ## chr::start::end
    ele = wholeLine[outKey].strip().split('\t'); ## whole 18 columns here
    c1Str = ','.join([str(xx) for xx in counts_d[outKey][0]]); l1Str = ','.join([str(xx) for xx in counts_d[outKey][1]]);
    c2Str = ','.join([str(xx) for xx in counts_d[outKey][2]]); l2Str = ','.join([str(xx) for xx in counts_d[outKey][3]]);
    c1 = scipy.sum(counts_d[outKey][0]); l1 = scipy.sum(counts_d[outKey][1]); 
    c2 = scipy.sum(counts_d[outKey][2]); l2 = scipy.sum(counts_d[outKey][3]);
    PBI_1=0.0; PBI_2=0.0; ## percent backspliced in
    if (c1+l1)>0: ## sample_1 has count
      PBI_1 = float(2*c1)/float(2*c1+l1);
    if (c2+l2)>0: ## sample_2 has count
      PBI_2 = float(2*c2)/float(2*c2+l2);
    deltaPBI = PBI_1 - PBI_2; 
    myP = fp[outKey]; myFDR=fdr[outKey]; 

    outStr = '\t'.join([ele[0],ele[1],ele[2],ele[5],ele[9],ele[10],ele[11],ele[13],ele[14],ele[15],ele[16],ele[17]]);
    outStr+= '\t'+c1Str+'\t'+l1Str+'\t'+c2Str+'\t'+l2Str+'\t'+str(PBI_1)+'\t'+str(PBI_2)+'\t'+str(deltaPBI)+'\t'+str(myP)+'\t'+str(myFDR);
    outFile.write(outStr+'\n');

  outFile.close();
  logging.debug("Done writing output..");

#
##### processes #######
#
getLinearCounts(bams_1, lcs1);
getLinearCounts(bams_2, lcs2);
#
computeFisherExact();
#
computeFDR();
#
writeResult();
#
#############
## calculate total running time
#############
logging.debug("Program ended");
currentTime = time.time();
runningTime = currentTime-startTime; ## in seconds
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));

sys.exit(0);
