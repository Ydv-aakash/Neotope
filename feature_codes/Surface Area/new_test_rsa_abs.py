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
def rsa_abs(path_to_sasa_abs_csv_file,residue):
    df2 = pd.read_csv(path_to_sasa_abs_csv_file)
    score = []
    for i in range(len(df2)):
        score.append(df2['sasa_abs_2.0_bound'][i])

    rsa = []
    for i in range(len(residue)):
        if (residue[i] == 'ALA' or residue[i] == 'a' or residue[i] == 'A'):
            a = score[i]
            b = float(a)
            c = b / 121.0
            rsa.append(c)
        elif residue[i] == 'ARG' or residue[i] == 'r' or residue[i] == 'R':
            a = score[i]
            b = float(a)
            c = b / 265.0
            rsa.append(c)
        elif residue[i] == 'ASN' or residue[i] == 'n' or residue[i] == 'N':
            a = score[i]
            b = float(a)
            c = b / 187.0
            rsa.append(c)


        elif residue[i] == 'ASP' or residue[i] == 'd' or residue[i] == 'D':
            a = score[i]
            b = float(a)
            c = b / 187.0
            rsa.append(c)
        elif residue[i] == 'CYS' or residue[i] == 'c' or residue[i] == 'C':
            a = score[i]
            b = float(a)
            c = b / 148.0
            rsa.append(c)
        elif residue[i] == 'GLU' or residue[i] == 'e' or residue[i] == 'E':
            a = score[i]
            b = float(a)
            c = b / 214.0
            rsa.append(c)
        elif residue[i] == 'GLN' or residue[i] == 'q' or residue[i] == 'Q':
            a = score[i]
            b = float(a)
            c = b / 214.0
            rsa.append(c)
        elif residue[i] == 'GLY' or residue[i] == 'g' or residue[i] == 'G':
            a = score[i]
            b = float(a)
            c = b / 97.0
            rsa.append(c)

        elif residue[i] == 'HIS' or residue[i] == 'h' or residue[i] == 'H':
            a = score[i]
            b = float(a)
            c = b / 216.0
            rsa.append(c)
        elif residue[i] == 'ILE' or residue[i] == 'i' or residue[i] == 'I':
            a = score[i]
            b = float(a)
            c = b / 195.0
            rsa.append(c)
        elif residue[i] == 'LEU' or residue[i] == 'l' or residue[i] == 'L':
            a = score[i]
            b = float(a)
            c = b / 191.0
            rsa.append(c)
        elif residue[i] == 'LYS' or residue[i] == 'k' or residue[i] == 'K':
            a = score[i]
            b = float(a)
            c = b / 230.0
            rsa.append(c)
        elif residue[i] == 'MET' or residue[i] == 'm' or residue[i] == 'M':
            a = score[i]
            b = float(a)
            c = b / 203.0
            rsa.append(c)
        elif residue[i] == 'PHE' or residue[i] == 'f' or residue[i] == 'F':
            a = score[i]
            b = float(a)
            c = b / 228.0
            rsa.append(c)
        elif residue[i] == 'PRO' or residue[i] == 'p' or residue[i] == 'P':
            a = score[i]
            b = float(a)
            c = b / 154.0
            rsa.append(c)
        elif residue[i] == 'SER' or residue[i] == 's' or residue[i] == 'S':
            a = score[i]
            b = float(a)
            c = b / 143.0
            rsa.append(c)
        elif residue[i] == 'THR' or residue[i] == 't' or residue[i] == 'T':
            a = score[i]
            b = float(a)
            c = b / 163.0
            rsa.append(c)
        elif residue[i] == 'TRP' or residue[i] == 'w' or residue[i] == 'W':
            a = score[i]
            b = float(a)
            c = b / 264.0
            rsa.append(c)
        elif residue[i] == 'TYR' or residue[i] == 'y' or residue[i] == 'Y':
            a = score[i]
            b = float(a)
            c = b / 255.0
            rsa.append(c)
        elif residue[i] == 'VAL' or residue[i] == 'v' or residue[i] == 'V':
            a = score[i]
            b = float(a)
            c = b / 165.0
            rsa.append(c)

        a = 0
        b = 0
        c = 0

    return rsa


if __name__=='__main__':
    pdb, seq = pdbs_seq_info('F:\\sab_new_bound_features\\2.0\\2.0_latest\\latest_bench_test_own\\2.0_own_test_latest.fasta')
    residue=[]
    for each in seq:
        for eachres in each:
            residue.append(eachres)

    rsa=rsa_abs('F:\\sab_new_bound_features\\2.0\\2.0_latest\\new_test\\features\\sasa_abs_test_bench.csv',residue)
    df_final = pd.DataFrame(rsa)
    df_final.columns = ['rsa_2.0_bound_abs']
    df_final.to_csv('F:\\sab_new_bound_features\\2.0\\2.0_latest\\new_test\\features\\test_bench_rsa_abs.csv', index=False)
