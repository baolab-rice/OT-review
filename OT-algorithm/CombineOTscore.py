#python2.7
#ref: hg38

import sys,math  
import os
import numpy as np
import csv
import pandas as pd
import pickle
import argparse
import re



# for the ROC curve, we only analyze off-targets with certain PAM sites 
# assuming that no software can find the relatively rare PAM sites
# that are not GG/GA/AG
# referred to CRISPOR

#assume 20bp

#Code for Cropit, MIT websitescore, and Hsu aggregation are from CRISPOR.
#Inplement of CCTOP, COSMID by YP
#CFD from original author
#fix COSMID for indels


def find(s, ch):
  return [i for i, ltr in enumerate(s) if ltr == ch]
#reference https://stackoverflow.com/questions/11122291/how-to-find-char-in-string-and-get-all-the-indexes





def findRuns(lst):
  """ yield (start, end) tuples for all runs of ident. numbers in lst 
  >>> list(findRuns([1,1,1,0,0,1,0,1,1,1]))
  [(0, 3), (5, 6), (7, 10)]
  """
  start,end=False,False

  for i,x in enumerate(lst):
    if x and start is False:
      start=i
    if x==0 and start is not False and end is False:
      end=i-1
    if start is not False and end is not False:
       yield start,end+1       #and len is (end-start)
       start,end=False,False
    
  if start is not False:
    yield start,i+1       #and len is (end-start)

def Cropitoff(guideSeq, otSeq):  #from CRISPOR
  """
  see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4605288/ PMID 26032770
  >>> int(Cropitoff("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
  650
  # mismatch in 3' part
  >>> int(Cropitoff("GGGGGGGGGGGGGGGGGGGA","GGGGGGGGGGGGGGGGGGGG"))
  575
  # mismatch in 5' part
  >>> int(Cropitoff("AGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
  642
  # only mismatches -> least likely offtarget
  >>> int(Cropitoff("AAAAAAAAAAAAAAAAAAAA","GGGGGGGGGGGGGGGGGGGG"))
  -27
  """
  if len(guideSeq)==23:
    guideSeq = guideSeq[:20]
    otSeq = otSeq[:20]

  assert(len(guideSeq)==len(otSeq)==20)

  penalties = [5,5,5,5,5,5,5,5,5,5,70,70,70,70,70,50,50,50,50,50]
  score = 0.0

  # do the score only for the non-mism positions
  misList = []
  score = 0.0
  for i in range(0, 20):
    if guideSeq[i]!=otSeq[i]:
      misList.append(1)
    else:
      misList.append(0)
      score += penalties[i]
    
  # get the runs of mismatches and update score for these positions
  consecPos = set()
  singlePos = set()
  for start, end in findRuns(misList):
    if end-start==1:
      score += -penalties[start]/2.0
    else:
      # mean if they happen to fall into different segments
      startScore = penalties[start]
      endScore = penalties[end-1]
      score += -((startScore+endScore)/2.0)
  return score



#########HSU#########
compTable = { "a":"t", "A":"T", "t" :"a", "T":"A", "c":"g", "C":"G", "g":"c", "G":"C", "N":"N", "n":"n", 
        "Y":"R", "R" : "Y", "M" : "K", "K" : "M", "W":"W", "S":"S",
        "H":"D", "B":"V", "V":"B", "D":"H", "y":"r", "r":"y","m":"k",
        "k":"m","w":"w","s":"s","h":"d","b":"v","d":"h","v":"b","y":"r","r":"y" }

nuclFreqs= {('A', 'A'): 0.4819440141370613, ('G', 'G'): 0.6571297543038187, ('U', 'T'): 0.4663759334374003, ('U', 'C'): 0.18755561795635842, ('C', 'T'): 0.3917125484856841, ('G', 'A'): 0.472948896301865, ('G', 'T'): 1.0, ('A', 'G'): 0.2796160896968995, ('U', 'G'): 0.787929779020387, ('C', 'C'): 0.0, ('A', 'C'): 0.6804018833297995, ('C', 'A'): 0.5931243444910334}
posFreqs = [0.294369386, 0.29164666, 0.190210984, 0.306896763, 0.167251773, 0.219909422, 0.169797251, 0.406475982, 0.628680509, 0.588183598, 0.907111342, 0.522909141, 1.256594863, 0.693851359, 0.552888666, 1.158572718, 1.662766602, 1.01548686, 1.428913939]


def complRna(seq):
    " complement the sequence and translate to RNA "
    newseq = []
    for nucl in seq.upper():
        newseq.append( compTable[nucl].replace("T", "U") )
    return "".join(newseq)
#from crispor

def Hsuoff(guideSeq, otSeq):#aggregrate frequencies!!! from crispor
  """
  The Hsu score in a version that only uses the aggregrate frequencies
  >>> calcHsuSuppScore2("GAGTCCGAGCAGAAGAAGAA","GAGTCCGAGCAGAAGAAGAG")
  1.1258838441954209
  >>> calcHsuSuppScore2("GTGTCCGAGCAGAAGAAGAA","GAGTCCGAGCAGAAGAAGAA")
  0.14186956352790206
  >>> calcHsuSuppScore2("GAGTCCGAGCAGAAGAAGAA","GAGTCAGAACAGAAGAACAA")
  0.0
  >>> calcHsuSuppScore2("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG")
  1.0
  >>> calcHsuSuppScore2("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGT")
  0.5597235206124074
  >>> calcHsuSuppScore2("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGC")
  0.0
  >>> calcHsuSuppScore2("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGCGG")
  0.0
  >>> calcHsuSuppScore2("GAGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG")
  0.231942405261347
  """
  guideSeq = guideSeq[1:20]
  otSeq = otSeq[1:20]
  #print guideSeq, len(guideSeq), otSeq, len(otSeq)
  assert(len(guideSeq)==19)
  assert(len(otSeq)==19)# Hsu ignores pos 0

  rnaSeq = complRna(guideSeq)
  # "Predicted cutting frequencies for genome-wide targets were calculated by
  # multiplying, in series: fest = f(1) * g(N1,N1') * f(2) * g(N2,N2') * ... * h
  # with values f(i) and g(Ni, Ni') at position i corresponding,
  # respectively, to the aggregate position- and base-mismatch cutting
  # frequencies for positions and pairings indicated in Fig. 2c"
  mismatchPosList = []
  score = 1.0
  for i in range(0, 19):
    rnaNucl, dnaNucl = rnaSeq[i], otSeq[i]
    # "In case of a match, both were set equal to 1."
    if (rnaNucl, dnaNucl) in [('C', 'G'), ('U', 'A'), ('A', 'T'), ('G', 'C')]:
      f = 1.0
      g = 1.0
    else:
      f = posFreqs[i]
      g = nuclFreqs[(rnaNucl, dnaNucl)]
      mismatchPosList.append(i)
    score *= f * g

  # "the value h meanwhile re-weighted the estimated
  # frequency by the minimum pairwise distance between consecutive mismatches
  # in the target sequence. this distance value, in base-pairs, was divided
  # by 18 to give a maximum value of 1 (in cases where fewer than 2
  # mismatches existed, or where mismatches occurred on opposite ends of
  # the 19 bp target-window"
  if len(mismatchPosList)<2:
    h = 1.0
  else:
    dists = []
    for left, right in zip(mismatchPosList[:-1], mismatchPosList[1:]):
      dists.append(right-left)
    minDist = min(dists)
    h = minDist / 18.0
  score *= h

  return score

#########HSU#########


###############################MITwebsite###############################
hitScoreM = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]
def calcMitScore(string1,string2, startPos=0):
    """
    The MIT off-target score
    see 'Scores of single hits' on http://crispr.mit.edu/about
    startPos can be used to feed sequences longer than 20bp into this function
    the most likely off-targets have a score of 100
    >>> int(calcMitScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
    100
    # mismatches in the first three positions have no effect
    >>> int(calcMitScore("AGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
    100
    # less likely off-targets have lower scores
    >>> int(calcMitScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGA"))
    41
    """
    # The Patrick Hsu weighting scheme
    #print string1, string2
    if len(string1)==len(string2)==23:
        string1 = string1[:20]
        string2 = string2[:20]

    assert(len(string1)==len(string2)==20)

    dists = [] # distances between mismatches, for part 2
    mmCount = 0 # number of mismatches, for part 3
    lastMmPos = None # position of last mismatch, used to calculate distance

    score1 = 1.0
    for pos in range(0, len(string1)):
        if string1[pos]!=string2[pos]:
            mmCount+=1
            if lastMmPos!=None:
                dists.append(pos-lastMmPos)
            score1 *= 1-hitScoreM[pos]
            lastMmPos = pos

    # 2nd part of the score - distribution of mismatches
    if mmCount<2: # special case, not shown in the paper
        score2 = 1.0
    else:
        avgDist = sum(dists)/len(dists)
        score2 = 1.0 / (((19-avgDist)/19.0) * 4 + 1)

    # 3rd part of the score - mismatch penalty
    if mmCount==0: # special case, not shown in the paper
        score3 = 1.0
    else:
        score3 = 1.0 / (mmCount**2)

    score = score1 * score2 * score3 * 100
    return score

###############################MITwebsite###############################


###############################CFD###############################
def get_parser():
    parser = argparse.ArgumentParser(description='Calculates CFD score')
    parser.add_argument('--wt',
        type=str,
        help='WT 23mer sgRNA sequence')
    parser.add_argument('--off',
        type=str,
        help='Off-target 23mer sgRNA sequence')
    return parser

#Reverse complements a given string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

#Unpickle mismatch scores and PAM scores
def get_mm_pam_scores():
    try:
        mm_scores = pickle.load(open('mismatch_score.pkl','rb'))
        pam_scores = pickle.load(open('pam_scores.pkl','rb'))
        return (mm_scores,pam_scores)
    except: 
        raise Exception("Could not find file with mismatch scores or PAM scores")

#Calculates CFD score
def calc_cfd(wt,sg,pam):
    mm_scores,pam_scores = get_mm_pam_scores()
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    print "sg:", sg
    print "DNA", wt
    s_list = list(sg)
    wt_list = list(wt)
    for i,sl in enumerate(s_list):
        if wt_list[i] == sl:
            score*=1
        else:
            key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
            print key
            score*= mm_scores[key]
    score*=pam_scores[pam]
    return (score)

def CFD_main(wt,off):
    mm_scores,pam_scores = get_mm_pam_scores()
    #wt =  "AGTCTGAGTCGGAGCCAGGGGGG"
    #off = "GGTCTGAGTCGGAGCCAGGGCGG"
    m_wt = re.search('[^ATCG]',wt)
    print m_wt
    m_off = re.search('[^ATCG]',off)
    print m_off
    if (m_wt is None) and (m_off is None):
        pam = off[-2:]
        sg = off[:-3]
        cfd_score = calc_cfd(wt,sg,pam)
        return (cfd_score)

###############################CFD###############################

def CCTOPoff(MMpattern):
  guide_len=len(MMpattern)
  MM_index=find(MMpattern,"*")#get a list of MM index
  #print MM_index
  MM_score=[]

  for i in MM_index:
    s=math.pow(1.2, float(i+1))
    #print s
    MM_score[len(MM_score):] = [s]

  CCTOP_score=sum(MM_score)
  #print MM_index
  #print MM_score
  #print CCTOP_score
  return 224.025599547-CCTOP_score

  
#if PAM is different, there will be +20 for the score... 
def COSMIDoff(MMpattern,NRG_):

# if there is indel, calculate with the same way
  scorelist=[0.12,0.13,0.15,0.17,0.19,0.21,0.23,0.27,0.35,0.5,0.7,0.8,1.1,1.3,1.9,2.3,3.0,4.0,5.0,6.0]
  COSMID_score=0


  guide_len=len(MMpattern)#in case len(gRNA)is not 20

# if indel, truncate the left side.
#del:+0.51
#ins:+0.7
  find(MMpattern,"*")
  MM_index=find(MMpattern,"*")#get a list of MM index
  print MM_index


  for i in MM_index:
    COSMID_score+=scorelist[i]

#  else:
#    COSMID_score=COSMID_score+0.51
#
#    div = len(MMpattern)-20
#    MMpattern_left=MMpattern[0:div]
#    MMpattern_right=MMpattern[-20:]
#    print(MMpattern_left,MMpattern_right,MMpattern)

#    MM_index=find(MMpattern,"*")#get a list of MM index


  if (NRG_!=1):
    COSMID_score+=20.0

  #print -COSMID_score
  return 48.4-COSMID_score
  #larger score, higher off-target. 48.4 means perfect match.




####generate matching pattern####

def match_pattern(wt,off):
  if len(wt)==len(off)==23:
    wt = wt[:20]
    off = off[:20]
  assert(len(wt)==len(off)==20)

  pattern=''
  for i in range(0, 20):
    if wt[i]!=off[i]:
      pattern=pattern+'*'
    else:
      pattern=pattern+'.'
  return pattern
    
def match_indel(wt,off):
  assert(len(wt)==len(off))
  pattern=''
  for i in range(0, len(wt)):
    if wt[i]!=off[i]:
      pattern=pattern+'*'
    else:
      pattern=pattern+'.'
  return pattern
  

####generate matching pattern####
 
def main():  


# consider chr index



  #filepath = sys.argv[1]

  #MMpattern="**.......*.....*...."
  #print list(findRuns([1,1,1,0,0,1,0,1,1,1]))
  #print Cropitoff("AGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG")
  #print CCTOPoff("********************")
  #print Hsuoff("GAGTCCGAGCAGAAGAAGAA","GAGTCCGAGCAGAAGAAGAG")


  folder_path = '/home/yp11/Desktop/OT_review_0105.Data/0204_regenerate_everything'

  # generate merged csv
  # use Oseq output for CNN, process seperately
 
  #CRISTA_path = folder_path+"/CRISTA.csv"
  oseq_path= folder_path+"/small_casoffinder/Truecasoff_CNN_Elevation_pred.csv"

#CRISTA use the same index as Tagscan
  #pd_CRISTA = pd.read_csv(CRISTA_path)
  pd_oseq= pd.read_csv(oseq_path)

#only use OT sequence here
  print("here")
  #new_df_1=pd.merge(pd_oseq, pd_CRISTA, how='left', left_on=['pos','start'], right_on = ['chromosome','start_position']) #for single gRNA design
  
  #pd_oseq=pd.merge(new_df_1, pd_predCRISPR,  how='left', left_on=['pos','start'], right_on = ['chromosome','position']) #for single gRNA design
  # generate merged csv



  for row in pd_oseq.itertuples():
   # print(row)
    print(row.crRNA_noind,row.DNA_noind)
    wt=str(row.crRNA_noind).upper()
    wt=wt.replace('R','G')
    ot=str(row.DNA_noind).upper()
    print(wt,ot)


    if ot[-2:]== "GG" or ot[-2:]== "AG" or ot[-2:]== "GA":
      NRG_=1
    else:
      NRG_=0
    
    #print Cropitoff(row.wildtype_seq,row.offtarget_seq)
    #print Hsuoff(row.wildtype_seq,row.offtarget_seq)
    #print CCTOPoff(row.mismatchPos)
    #print COSMIDoff(row.mismatchPos,NGG_)

    if len(wt)==23 and len(ot)==23:
        mismatchPos=match_pattern(wt,ot)
        pd_oseq.at[row.Index,'Cropit']=Cropitoff(wt,ot)
        pd_oseq.at[row.Index,'Hsu']=Hsuoff(wt,ot)
        pd_oseq.at[row.Index,'CCTOP']=CCTOPoff(mismatchPos)
        pd_oseq.at[row.Index,'COSMID']=COSMIDoff(mismatchPos,NRG_)
    #pd_oseq.at[row.Index,'CFD']=CFD_main(wt,ot)
    #print (CFD_main(wt,ot))
        pd_oseq.at[row.Index,'MIT']=calcMitScore(wt,ot)
    elif len(wt)==23 and len(ot)==22: # 1 base del in ot
        ot=str(row.DNA).upper()
        mismatchPos=match_pattern(wt,ot)
        pd_oseq.at[row.Index,'COSMID']=COSMIDoff(mismatchPos,NRG_)-0.51
    elif len(ot)==24: # 1 base extra in ot
        wt=str(row.crRNA).upper()
        wt=wt.replace('R','G')
        ot=str(row.DNA).upper()

#get the info of the first base
        if wt[0]==ot[0]:
            cosmid_score=0.0
        else:
            cosmid_score=-0.1        
        wt=wt[:23]
        ot=ot[:23]
        mismatchPos=match_pattern(wt,ot)
        pd_oseq.at[row.Index,'COSMID']=COSMIDoff(mismatchPos,NRG_)-0.7
  
 # new_df_3=pd.merge(pd_oseq, pd_oseq,  how='left', left_on=['offtargetSeq','chrom','start'], right_on = ['DNA','Chromosome','Position']) #for single 
  pd_oseq.to_csv('/home/yp11/Desktop/OT_review_0105.Data/0204_regenerate_everything/small_casoffinder/Truecasoff_CNN_Elevation_pred_others.csv',index=False)
   


  



#   if not os.path.isfile(filepath):
#       print("File path {} does not exist. Exiting...".format(filepath))
#       sys.exit()



#   fo = open(sys.argv[2], "w") 
#
#   with open(filepath) as fp:  
#      line = fp.readline()
#      for line in fp:
#         if line[0]=="@":
#            fo.write(line)
#         else:
#            linelist=line.split("	")
#            if linelist[2] == "chr" + chrID and chrstart <= int(linelist[3]) and int(linelist[3]) <= chrend :
#               #print ("{}".format(sys.argv[1],sys.argv[2],chrID, chrstart, chrend))
#               fo.write(line)
#   
#   fo.close() 


 # print(calcMitScore("GAGTCCGAGCAGAAGAAGAA","GAGTCCGAGCAGAAGAAGAG"))
 # test=CFD_main("GAGTCCGAGCAGAAGAAGAA","GAGTCCGAGCAGAAGATGGG")
 # match_pattern("GAGTCCGAGCAGAAGAAGAA","GAGTCCGAGCAGAAGATGGG")
  #print(test)
if __name__ == '__main__':  
   main()

