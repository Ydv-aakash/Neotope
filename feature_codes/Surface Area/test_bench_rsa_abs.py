
import os
import pandas as pd
from collections import Counter
import regex as re
import math
import numpy as np
from pdbs_seq_coordinates import pdbs_seq_info,ca_coordinates

pdb,seq=pdbs_seq_info('F:\\sab_new_bound_features\\2.0\\bench_test\\new_2.0_benchmark_test.fasta')

res_id,res_name,ca_xcoordinate,ca_ycoordinate,ca_zcoordinate=ca_coordinates(seq,'F:\\sab_new_bound_features\\2.0\\bench_test\\split_2.0_bench\\',pdb)
#---------------------------------------------------------------------------------------------------------------------------------
df=pd.read_csv('F:\\sab_new_bound_features\\2.0\\bench_test\\test_bench_pdbs.csv')
pdbs_sasa=[]
for i in range(len(df)):
    pdbs_sasa.append(df['pdbs'][i][7:8])


df2=pd.read_csv('F:\\sab_new_bound_features\\2.0\\bench_test\\features\\sasa\\sasa_abs_test_bench.csv')
score=[]
for i in range(len(df2)):
    score.append(df2['sasa_abs_2.0_bound'][i])


rsa=[]
for i in range(len(pdbs_sasa)):
    if (pdbs_sasa[i]=='ALA' or pdbs_sasa[i]=='a'or pdbs_sasa[i]=='A'):
        a=score[i]
        b=float(a)
        c=b/121.0
        rsa.append(c)
    elif pdbs_sasa[i] == 'ARG' or pdbs_sasa[i]=='r'or pdbs_sasa[i]=='R':
        a = score[i]
        b = float(a)
        c = b / 265.0
        rsa.append(c)
    elif pdbs_sasa[i]=='ASN'or pdbs_sasa[i]=='n'or pdbs_sasa[i]=='N':
        a=score[i]
        b=float(a)
        c=b/187.0
        rsa.append(c)


    elif pdbs_sasa[i]=='ASP'or pdbs_sasa[i]=='d'or pdbs_sasa[i]=='D':
        a=score[i]
        b=float(a)
        c=b/187.0
        rsa.append(c)
    elif pdbs_sasa[i]=='CYS'or pdbs_sasa[i]=='c'or pdbs_sasa[i]=='C':
        a=score[i]
        b=float(a)
        c=b/148.0
        rsa.append(c)
    elif pdbs_sasa[i]=='GLU'or pdbs_sasa[i]=='e'or pdbs_sasa[i]=='E':
        a=score[i]
        b=float(a)
        c=b/214.0
        rsa.append(c)
    elif pdbs_sasa[i]=='GLN'or pdbs_sasa[i]=='q'or pdbs_sasa[i]=='Q':
        a=score[i]
        b=float(a)
        c=b/214.0
        rsa.append(c)
    elif pdbs_sasa[i]=='GLY'or pdbs_sasa[i]=='g'or pdbs_sasa[i]=='G':
        a=score[i]
        b=float(a)
        c=b/97.0
        rsa.append(c)

    elif pdbs_sasa[i]=='HIS'or pdbs_sasa[i]=='h'or pdbs_sasa[i]=='H':
        a=score[i]
        b=float(a)
        c=b/216.0
        rsa.append(c)
    elif pdbs_sasa[i]=='ILE'or pdbs_sasa[i]=='i'or pdbs_sasa[i]=='I':
        a=score[i]
        b=float(a)
        c=b/195.0
        rsa.append(c)
    elif pdbs_sasa[i]=='LEU'or pdbs_sasa[i]=='l'or pdbs_sasa[i]=='L':
        a=score[i]
        b=float(a)
        c=b/191.0
        rsa.append(c)
    elif pdbs_sasa[i]=='LYS'or pdbs_sasa[i]=='k'or pdbs_sasa[i]=='K':
        a=score[i]
        b=float(a)
        c=b/230.0
        rsa.append(c)
    elif pdbs_sasa[i]=='MET'or pdbs_sasa[i]=='m'or pdbs_sasa[i]=='M':
        a=score[i]
        b=float(a)
        c=b/203.0
        rsa.append(c)
    elif pdbs_sasa[i]=='PHE'or pdbs_sasa[i]=='f'or pdbs_sasa[i]=='F':
        a=score[i]
        b=float(a)
        c=b/228.0
        rsa.append(c)
    elif pdbs_sasa[i]=='PRO'or pdbs_sasa[i]=='p'or pdbs_sasa[i]=='P':
        a=score[i]
        b=float(a)
        c=b/154.0
        rsa.append(c)
    elif pdbs_sasa[i]=='SER'or pdbs_sasa[i]=='s'or pdbs_sasa[i]=='S':
        a=score[i]
        b=float(a)
        c=b/143.0
        rsa.append(c)
    elif pdbs_sasa[i]=='THR'or pdbs_sasa[i]=='t'or pdbs_sasa[i]=='T':
        a=score[i]
        b=float(a)
        c=b/163.0
        rsa.append(c)
    elif pdbs_sasa[i]=='TRP'or pdbs_sasa[i]=='w'or pdbs_sasa[i]=='W':
        a=score[i]
        b=float(a)
        c=b/264.0
        rsa.append(c)
    elif pdbs_sasa[i]=='TYR'or pdbs_sasa[i]=='y'or pdbs_sasa[i]=='Y':
        a=score[i]
        b=float(a)
        c=b/255.0
        rsa.append(c)
    elif pdbs_sasa[i]=='VAL'or pdbs_sasa[i]=='v'or pdbs_sasa[i]=='V':
        a=score[i]
        b=float(a)
        c=b/165.0
        rsa.append(c)

    a=0
    b=0
    c=0

df_final=pd.DataFrame(rsa)
df_final.columns=['rsa_2.0_bound_abs']
df_final.to_csv('F:\\test_bench_rsa_abs.csv',index=False)
