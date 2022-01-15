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

def di_abs(path_to_di_values,common_backstring_on_di_values):
    seqscore={}
    count=0
    seqscore=[]
    for i in range(len(seq)):
        print(i)
        seqscore.append([])
        file=open(path_to_di_values+pdb[i][1:].strip().lower()+str(common_backstring_on_di_values),'r')
        l = []
        for x in file:
            l.append(x)

        for x in l[7:]:
            seqscore[i].append(float(x[38:50].strip()))
        l = []
        file.close()
    diabsolutescore=[]
    for i in range(len(seqscore)):
        for j in range(len(seqscore[i])):
            diabsolutescore.append(seqscore[i][j])
    return diabsolutescore,seqscore


#------------------------------------------------------------------------------------

if __name__=="__main__":
    pdb, seq = pdbs_seq_info('F:\\sab_new_bound_features\\2.0\\2.0_latest\\latest_bench_test_own\\2.0_own_test_latest.fasta')
    diabsolutescore,seqscore=di_abs('F:\\sab_new_bound_features\\2.0\\2.0_latest\\new_test\\di\\','_202107312216_unbound.tbl')
    df = pd.DataFrame(diabsolutescore)
    df.columns=['di_abs_2.0_bound']
    #df.to_csv('F:\\sab_new_bound_features\\2.0\\2.0_latest\\new_test\\features\\di_abs_test_bench.csv', index=False)

