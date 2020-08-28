from __future__ import division
import numpy as np
from numpy import *
import scipy.io as sio
from scipy import stats
import math
import pickle
from sklearn.metrics import roc_auc_score
from svmutil import *
import csv
import sys
import os

def position_binary(seq_on, seq_off):
    L=len(seq_on)
    dic=['AT', 'AC', 'AG', 'TC', 'TG', 'TA', 'GA', 'GT', 'GC', 'CA', 'CT', 'CG']
    bina = np.zeros((1, 12))
    if seq_on[0] == seq_off[0]:
        bina = bina
    else:
        onoff = seq_on[0] + seq_off[0]
        for j in range(0, len(dic)):
            if onoff == dic[j]:
                bina[0, j] = 1
                continue

    for i in range(1, L):
        bina_ini=np.zeros((1,12))
        if seq_on[i]==seq_off[i]:
            bina=np.column_stack((bina,bina_ini))
        else:
            onoff=seq_on[i]+seq_off[i]
            for j in range(0, len(dic)):
                if onoff==dic[j]:
                    bina_ini[0,j]=1
                    bina=np.column_stack((bina,bina_ini))
                    continue
    
    binary=bina
    return binary

def countGC(s, start,end):
#'compute the GC counts and GC contents of a sequences from (start,end)'
    GCcounts=len(s[start:end].replace('A', '').replace('T', ''))
    GCcontent=round(GCcounts/len(s),15)
    S=s[start:end]
    rat=0
    if len(S.replace('A', '').replace('T', ''))!=0:
        GCskew=round((len(S.replace('A', '').replace('T', '').replace('C', ''))-len(S.replace('A', '').replace('T', '').replace('G', '')))/(len(S.replace('A', '').replace('T', ''))),15)
    else:
        GCskew=0
    
    if len(S.replace('C', '').replace('G', ''))!=0:
        ATskew=round((len(S.replace('C', '').replace('G', '').replace('T', ''))-len(S.replace('C', '').replace('G', '').replace('A', '')))/(len(S.replace('C', '').replace('G', ''))),15)
    else:
        ATskew=0
        
    if ATskew!=0:
        rat=round(GCskew/ATskew,15)
        
    return GCcounts, GCcontent, GCskew, ATskew, rat

def deltaGCACration(seq_on, seq_off, start,end):
    GCcounts1, GCcontent1, GCskew1, ATskew1, rat1=countGC(seq_on, start,end)
    GCcounts2, GCcontent2, GCskew2, ATskew2, rat2=countGC(seq_off, start,end)
    GC=np.zeros((1,5))
    GC[0][0]=GCcounts2-GCcounts1
    GC[0][1]=GCcontent2-GCcontent1
    GC[0][2]=GCskew2-GCskew1
    GC[0][3]=ATskew2-ATskew1
    GC[0][4]=rat2-rat1
    return GC
   
def feaAll(seq_ons, seq_offs):
    N=len(seq_ons)
    seq_on0=seq_ons[0]
    seq_off0=seq_offs[0]
    GC0=deltaGCACration(seq_on0, seq_off0, 0,len(seq_on0))
    binary0=position_binary(seq_on0, seq_off0)
    fea_all=np.column_stack((GC0,binary0))
    for i in xrange(1,N):
        seq_on=seq_ons[i]
        seq_off=seq_offs[i]
        GC=deltaGCACration(seq_on, seq_off, 0, len(seq_on))
        binary=position_binary(seq_on, seq_off)
        fea=np.column_stack((GC,binary))
        fea_all=np.row_stack((fea_all,fea))
    return fea_all

def mapminmax(feature, ymax, ymin):
    minmax = np.zeros((2, len(feature[0, :])))
    minmax[0][0] = max(feature[:, 0])
    minmax[1][0] = min(feature[:, 0])

    norm_fea = (ymax - ymin) * (feature[:, 0] - min(feature[:, 0])) / (max(feature[:, 0]) - min(feature[:, 0])) + ymin

    for i in xrange(1, len(feature[0, :])):
        minmax[0][i] = max(feature[:, i])
        minmax[1][i] = min(feature[:, i])
        norm = (ymax - ymin) * (feature[:, i] - min(feature[:, i])) / (max(feature[:, i]) - min(feature[:, i])) + ymin
        norm_fea = np.column_stack((norm_fea, norm))
    return minmax, norm_fea

def train_models(fea_train_All, num_pos, negative_ind):
    fea_pos=fea_train_All[0:num_pos,:]
    fea_neg=fea_train_All[num_pos:len(fea_train_All[:,0]),:]
    num_mod=len(negative_ind[0,:])
    minmaxs={}
    C=1
    G=-4
    para = '-t 2 -b 1' + ' -c ' + str(2 ** C) + ' -g ' + str(2 ** G)
    for i in xrange(0,num_mod):
        ind=negative_ind[:,i]
        fea_tr_neg=fea_neg[ind,:]
        lab_n=np.zeros((len(fea_pos[:,0]),1))
        lab_p=np.ones((len(fea_tr_neg[:,0]),1))
        labels=np.row_stack((lab_p,lab_n))
        fea_inte=np.vstack((fea_pos,fea_tr_neg))
        minmax, norm_fea=mapminmax(fea_inte, 1, 0)
        train_file_path=libsvm_dataformat(norm_fea, labels, [], 'train')
        #os.remove('train_file_path')
        Lab, Fea = svm_read_problem(train_file_path)        
        m = svm_train(Lab, Fea, para)
        minmaxs.update({'minmax_'+str(i): minmax})
        svm_save_model('SVM_Model/model_'+str(i),m)
    #pickle.dump(models, open('Models', 'wb'))
    pickle.dump(minmaxs, open('minmaxs', 'wb'))

def pred_off_target(test_fea):
    #models=pickle.load(open('Models', 'rb'))
    minmaxs=pickle.load(open('minmaxs', 'rb'))
    labels=np.zeros((len(test_fea[:,0]),1))
    n=len(minmaxs)
    scores=np.zeros((len(test_fea[:,0]),1))
    for i in xrange(0, n):
        minmax=minmaxs['minmax_'+str(i)]
        model=svm_load_model('SVM_Model/model_'+str(i))
        test_file_path=libsvm_dataformat(test_fea, labels, minmax, 'test')
        #os.remove('test_file_path')
        Lab, Fea = svm_read_problem(test_file_path)
        #outputfile(Fea,'Fea'+str(i))
        p_label_a, p_acc, p_val = svm_predict(Lab, Fea, model, '-b 1')
        k=0
        if p_label_a[0]==1:
            if p_val[0][0]>p_val[0][1]:
                k=0
            else:
                k=1
        else:
            if p_val[0][0]>p_val[0][1]:
                k=1
            else:
                k=0
                    
        for j in xrange(0, len(test_fea[:,0])):
            scores[j]=scores[j]+p_val[j][k]
    scores=scores/n
    return scores
            

        
def libsvm_dataformat(feature, labels, minmax, flag):
    if flag == 'train':
        file_path = 'train_fea'
        norm_fea=feature
    elif flag == 'test':
        ymax = 1
        ymin = 0
        norm_fea = (ymax - ymin) * (feature[:, 0] - minmax[1, 0]) / (minmax[0, 0] - minmax[1, 0]) + ymin
        for i in xrange(1, len(feature[0, :])):
            norm = (ymax - ymin) * (feature[:, i] - minmax[1, i]) / (minmax[0, i] - minmax[1, i]) + ymin
            norm_fea = np.column_stack((norm_fea, norm))

        file_path = 'test_fea'
    
    f = open(file_path, 'w')
    rows = len(labels)
    columns = len(norm_fea[0, :])
    where_are_NaNs = isnan(norm_fea)
    norm_fea[where_are_NaNs] = 0
    for i in xrange(0, rows):
        f.write(str(labels[i, 0]) + '\t')
        for j in xrange(0, columns):
            f.write(str(j + 1) + ':' + str(norm_fea[i, j]) + '\t')
        f.write('\n')
    f.close()
    return file_path

def output_matrix(matrix,name):
    file_path = name  + '.csv'
    csv_file=open(file_path, 'wb')
    writer = csv.writer(csv_file, delimiter=',')
    for line in matrix:
        writer.writerow(line)

def readCasOffinder(sgRNA, filepath, filetype):
    ## filetype=1: online tool output
    ##             the 3rd column is the potential off-target sequence
    ##             4th: chromosome id
    ##             5th: position
    ##             6th: the direction
    ##             7th: mismatch number
    ## filetype=2: offline tool output
    ##             2nd: chromosome id;
    ##             3rd: position;
    ##             4th: potential off-target sequence;
    ##             5th: strand; 
    ##             6th: mismatch number

    on=sgRNA
    offs=[]
    chros=[]
    poss=[]
    strands=[]
    misnums=[]
    zero_ind=[]
    feas=np.zeros((2,281))
    file=open(filepath)
    i=0
    for line in file:
        s=line.split('\t')
        off=[]
        chro=[]
        pos=[]
        st=[]
        mn=[]
        if filetype==1:
            if i==0:
                i=i+1
                continue
            off=s[2].upper()
            chro=s[3]
            pos=s[4]
            st=s[5]
            mn=s[6]
        elif filetype==2:
            off=s[3].upper()
            chro=s[1]
            pos=s[2]
            st=s[4]
            mn=s[5]
        #print s
        offs.append(off)
        chros.append(chro)
        poss.append(pos)
        strands.append(st)
        misnums.append(int(mn[0]))
        if int(mn[0])==0:
            zero_ind.append(i)
        GC=deltaGCACration(on, off, 0,len(on))
        binary=position_binary(on, off)
        fea=np.column_stack((GC, binary))
        feas=np.vstack((feas,fea))
        i=i+1
    return feas[2:len(feas[:,0]),:], offs, chros, poss, strands, misnums, zero_ind
            

if __name__ == '__main__':
#    train_models(fea_train_All, num_pos, negative_ind)
    sgRNA = sys.argv[1]
    filepath = sys.argv[2]
    filepath = filepath.replace('\\\\', '/')
    filepath = filepath.replace('\\', '/')
    filepath = filepath.replace('//', '/')
    filetype = int(sys.argv[3])
    feas, offs, chros, poss, strands, misnums, zero_ind=readCasOffinder(sgRNA, filepath, filetype)
    scores=pred_off_target(feas)
    ind=np.where(scores>0.5)
    out_offs=list()
    out_chros=list()
    out_poss=list()
    out_strands=list()
    out_misnums=list()
    out_scores=list()
    ind_re=list(set(ind[0]).difference(set(zero_ind)))

    for i in xrange(0, len(ind_re)):
        out_offs.append(offs[ind_re[i]])
        out_chros.append(chros[ind_re[i]])
        out_poss.append(poss[ind_re[i]])
        out_strands.append(strands[ind_re[i]])
        out_misnums.append(misnums[ind_re[i]])
        out_scores.append(scores[ind_re[i]])

    if len(zero_ind)>1:
        out_zeros_offs=list()
        out_zeros_chros=list()
        out_zeros_poss=list()
        out_zeros_strands=list()
        out_zeros_misnums=list()
        out_zeros_scores=list()
        for j in xrange(0, len(zero_ind)):
            out_zeros_offs=offs[zero_ind[j]]
            out_zeros_chros=chros[zero_ind[j]]
            out_zeros_poss=poss[zero_ind[j]]
            out_zeros_strands=strands[zero_ind[j]]
            out_zeros_misnums=misnums[zero_ind[j]]
            out_zeros_scores=scores[zero_ind[j]]
    score_index=np.zeros((len(out_scores),1))

    for i in xrange(0, len(score_index)):
        score_index[i,0]=i

    Out=np.hstack((out_scores, score_index))
    out = sorted(Out, key=lambda row: row[0], reverse=True)
    file_path = 'predict_results/predict_results.csv'
    with open(file_path, 'w') as csv_file:
        fieldnames = ['off-target','chromosome', 'position', 'direction', 'mismatch number', 'score']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        if len(zero_ind)>1:
            for i in xrange(0, len(out_zeros_scores)):
                ind = int(zero_ind[i])
                writer.writerow(
                    {'off-target': out_zeros_offs[ind], 'chromosome':out_zeros_chros[ind],'position': out_zeros_poss[ind], 'direction': out_zeros_strands[ind],
                     'mismatch number': out_zeros_misnums[ind], 'score': '*'})

        for i in xrange(0, len(out_scores)):
            ind = int(out[i][1])
            writer.writerow(
                {'off-target': out_offs[ind], 'chromosome':out_chros[ind],'position': out_poss[ind], 'direction': out_strands[ind],
                 'mismatch number': out_misnums[ind], 'score': out_scores[ind]})

