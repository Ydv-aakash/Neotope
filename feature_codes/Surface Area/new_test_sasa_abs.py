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
def sasa_abs(path_to_sasa_files_from_psaia,common_tag_from_psaia):

    l=[]
    episasa=[]
    testing=[]
    for i in range(len(seq)):
        print(i)
        episasa.append([])
        f_in = open(path_to_sasa_files_from_psaia + pdb[i][1:].strip()+common_tag_from_psaia,'r')
        ch=[]
        for x in f_in:
            ch.append(x.strip())
        for z in ch[12:]:
            episasa[i].append(z[110:122].strip())
        testing=ch.copy()
        ch=[]
    seqscore=[]
    for i in range(len(episasa)):
        seqscore.append([])
        for j in range(len(episasa[i])):
            seqscore[i].append(float(episasa[i][j]))


    sasaabsolutescore=[]

    for i in range(len(seqscore)):
        for j in range(len(seqscore[i])):
            sasaabsolutescore.append(seqscore[i][j])
    return sasaabsolutescore
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if __name__=='__main__':
    pdb, seq = pdbs_seq_info('F:\\sab_new_bound_features\\2.0\\2.0_latest\\latest_bench_test_own\\2.0_own_test_latest.fasta')
    sasaabsolutescore=sasa_abs('F:\\sab_new_bound_features\\2.0\\2.0_latest\\new_test\\sasa\\','_202107312213_unbound.tbl')
    df = pd.DataFrame(sasaabsolutescore)
    df.columns=['sasa_abs_2.0_bound']

    df.to_csv('F:\\sab_new_bound_features\\2.0\\2.0_latest\\new_test\\sasa_abs_test_bench.csv', index=False)
