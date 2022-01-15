import os
import pandas as pd
from collections import Counter
import regex as re
import math
import numpy as np
import concurrent.futures
#----------------------------------------------------------------------------

def pdbs_seq_info(path):
    file=open(path,'r')
    pdb=[]
    seq=[]
    for each in file:
        if each.startswith('>'):
         pdb.append(each.strip())
        else:
            seq.append(each.strip())

    file.close()
    return pdb,seq

#----------------------------------------------------------------------------------
def pi_abs(path_to_all_pi_files):

    seqscore={}
    for i in range(5,16):
        numb=int(i)
        file = os.listdir(path_to_all_pi_files +str('_')+ str(numb))
        backstring = file[0][6:]
        temp_seqscore=[]
        for j in range(len(pdb)):
            temp_seqscore.append([])
            file_2=open(path_to_all_pi_files + str('_')+str(numb)+str('\\')+pdb[j][1:]+backstring,'r')
            l = []
            for x in file_2:
                l.append(x)

            for x in l[7:]:
                temp_seqscore[j].append(float(x[38:50].strip()))
            l = []
            file_2.close()
        seqscore[str(i)]=temp_seqscore


    ##absolute
    a=[]
    b=[]
    c=[]
    d=[]
    e=[]
    f=[]
    g=[]
    h=[]
    i=[]
    j=[]
    k=[]

    for each in seqscore['5']:
        for eacheach in each:
            a.append(eacheach)

    for each in seqscore['6']:
        for eacheach in each:
            b.append(eacheach)


    for each in seqscore['7']:
        for eacheach in each:
            c.append(eacheach)

    for each in seqscore['8']:
        for eacheach in each:
            d.append(eacheach)


    for each in seqscore['9']:
        for eacheach in each:
            e.append(eacheach)


    for each in seqscore['10']:
        for eacheach in each:
            f.append(eacheach)


    for each in seqscore['11']:
        for eacheach in each:
            g.append(eacheach)


    for each in seqscore['12']:
        for eacheach in each:
            h.append(eacheach)


    for each in seqscore['13']:
        for eacheach in each:
            i.append(eacheach)



    for each in seqscore['14']:
        for eacheach in each:
            j.append(eacheach)


    for each in seqscore['15']:
        for eacheach in each:
            k.append(eacheach)

    return a,b,c,d,e,f,g,h,i,j,k

if __name__=="__main__":
    pdb, seq = pdbs_seq_info('F:\\sab_new_bound_features\\2.0\\2.0_latest\\latest_bench_test_own\\2.0_own_test_latest.fasta')
    a,b,c,d,e,f,g,h,i,j,k=pi_abs('F:\\sab_new_bound_features\\2.0\\2.0_latest\\new_test\\pi\\pi')
    df = pd.DataFrame(list(zip(*[a,b,c,d,e,f,g,h,i,j,k])))
    df.columns = ['pi_abs_2.0_bound_5','pi_abs_2.0_bound_6','pi_abs_2.0_bound_7','pi_abs_2.0_bound_8','pi_abs_2.0_bound_9','pi_abs_2.0_bound_10','pi_abs_2.0_bound_11','pi_abs_2.0_bound_12','pi_abs_2.0_bound_13','pi_abs_2.0_bound_14','pi_abs_2.0_bound_15']
    df.to_csv('F:\\sab_new_bound_features\\2.0\\2.0_latest\\new_test\\features\\pi_abs_test_bench.csv', index=False)

