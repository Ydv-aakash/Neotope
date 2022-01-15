import os
import pandas as pd

file=open('F:\\sab_new_bound_features\\2.0\\2.0_latest\\latest_bench_test_own\\2.0_own_test_latest.fasta','r')

pdb=[]
seq=[]
for each in file:
    if each.startswith('>'):
     pdb.append(each.strip())
    else:
        seq.append(each.strip())
file.close()


for i in range(len(pdb)):
    file_3=open('F:\\sab_new_bound_features\\2.0\\2.0_latest\\new_test\\separate_fasta_files\\'+pdb[i][1:]+'.fasta','w')
    file_3.write(pdb[i])
    file_3.write('\n')
    file_3.write(seq[i])
    file_3.close()

for i in range(len(pdb)):
    print(i)
    command='psiblast -query F:\\sab_new_bound_features\\2.0\\2.0_latest\\new_test\\separate_fasta_files\\'+pdb[i][1:]+'.fasta -db C:\WINDOWS\system32\swissprot -num_iterations 3 -out F:\\sab_new_bound_features\\2.0\\2.0_latest\\new_test\\separate_fasta_files\\'+pdb[i][1:] + ' -out_ascii_pssm F:\\sab_new_bound_features\\2.0\\2.0_latest\\new_test\\separate_fasta_files\\'+pdb[i][1:]+'_pssm -evalue 0.001'
    os.system(command)


#############################################################################################################################################
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
def ca_coordinates(seq,path_to_complexes,pdb):

    res_id=[]
    res_name=[]
    ca_xcoordinate=[]
    ca_ycoordinate=[]
    ca_zcoordinate=[]
    l = []
    for i in range(len(seq)):
        ca_xcoordinate.append([])
        ca_ycoordinate.append([])
        ca_zcoordinate.append([])
        res_id.append([])
        res_name.append([])
        f_in = open(path_to_complexes + pdb[i][1:]+ '.pdb', 'r')
        for x in f_in:
            if x.startswith('ATOM') and x[11:15].strip() == 'CA':
                l.append(x)
        positions = []
        for p in l:
            temp=p[22:27].strip()
            if any(c.isalpha() for c in temp)==True:
                res = [re.findall(r'(\w+?)(\d+)', temp)[0]]
                res_2=res[0]
                res_3="".join(res_2)
                res_4=int(res_3)
                positions.append(int(res_4))
            else:
                positions.append(int(temp))
        count_positions = Counter(positions)
        x = 0
        while x <= (len(l) - 1):
            res_type_check=l[x][22:27].strip()
            if x == (len(l) - 1):
                ca_xcoordinate[i].append(l[x][30:38].strip())
                ca_ycoordinate[i].append(l[x][38:46].strip())
                ca_zcoordinate[i].append(l[x][47:55].strip())
                res_id[i].append(l[x][22:27].strip())
                res_name[i].append(l[x][15:20].strip())
                x += 1

            elif x != (len(l) - 1) and x < (len(l) - 1):
                if (count_positions[positions[x]]>= 1 and len(l[x][15:20].strip()) == 3 ):
                    ca_xcoordinate[i].append(l[x][30:38].strip())
                    ca_ycoordinate[i].append(l[x][38:46].strip())
                    ca_zcoordinate[i].append(l[x][47:55].strip())
                    res_id[i].append(l[x][22:27].strip())
                    res_name[i].append(l[x][15:20].strip())
                    x += 1
                #elif any(c.isalpha() for c in res_type_check)==True:
                 #   ca_xcoordinate[i].append(l[x][30:39].strip())
                 #   ca_ycoordinate[i].append(l[x][38:46].strip())
                  #  ca_zcoordinate[i].append(l[x][47:55].strip())
                  #  res_id[i].append(l[x][22:27].strip())
                  #  res_name[i].append(l[x][15:20].strip())
                  #  x += 1

                elif (count_positions[positions[x]] >= 1 and len(l[x][15:20].strip()) > 3):
                    if count_positions[positions[x]] != 3:
                        if (len(l[x][15:20].strip()) > 3 and len(l[x + 1][15:20].strip()) > 3):
                            if (float(l[x][55:60].strip()) >= float(l[x + 1][55:60].strip())):
                                ca_xcoordinate[i].append(l[x][30:38].strip())
                                ca_ycoordinate[i].append(l[x][38:46].strip())
                                ca_zcoordinate[i].append(l[x][47:55].strip())
                                res_id[i].append(l[x][22:27].strip())
                                res_name[i].append(l[x][15:20].strip())
                                x += 2
                            elif (float(l[x][55:60].strip()) < float(l[x + 1][55:60].strip())):
                                ca_xcoordinate[i].append(l[x + 1][30:38].strip())
                                ca_ycoordinate[i].append(l[x + 1][38:46].strip())
                                ca_zcoordinate[i].append(l[x + 1][47:55].strip())
                                res_id[i].append(l[x+ 1][22:27].strip())
                                res_name[i].append(l[x+1][15:20].strip())
                                x += 2
                    elif count_positions[positions[x]] == 3:
                        if (len(l[x][15:20].strip()) > 3 and len(l[x + 1][15:20].strip()) > 3 and len(l[x + 2][15:20].strip()) > 3 ):
                            if (float(l[x][55:60].strip()) >= (float(l[x + 1][55:60].strip()) and float(l[x + 2][55:60].strip()))):
                                ca_xcoordinate[i].append(l[x][30:38].strip())
                                ca_ycoordinate[i].append(l[x][38:46].strip())
                                ca_zcoordinate[i].append(l[x][47:55].strip())
                                res_id[i].append(l[x][22:27].strip())
                                res_name[i].append(l[x][15:20].strip())
                                x += 3
                            elif ((float(l[x][55:60].strip()) and float(l[x + 2][55:60].strip())) < float(l[x + 1][55:60].strip())):
                                ca_xcoordinate[i].append(l[x + 1][30:38].strip())
                                ca_ycoordinate[i].append(l[x + 1][38:46].strip())
                                ca_zcoordinate[i].append(l[x + 1][47:55].strip())
                                res_id[i].append(l[x+ 1][22:27].strip())
                                res_name[i].append(l[x+1][15:20].strip())
                                x += 3
                            elif ((float(l[x][55:60].strip()) and float(l[x + 1][55:60].strip())) < float(l[x + 2][55:60].strip())):
                                ca_xcoordinate[i].append(l[x + 2][30:38].strip())
                                ca_ycoordinate[i].append(l[x + 2][38:46].strip())
                                ca_zcoordinate[i].append(l[x + 2][47:55].strip())
                                res_id[i].append(l[x+ 2][22:27].strip())
                                res_name[i].append(l[x+2][15:20].strip())
                                x += 3
                        elif (len(l[x][15:20].strip()) > 3 and len(l[x + 1][15:20].strip()) > 3 and len(l[x + 2][15:20].strip()) == 3):
                            if (float(l[x][55:60].strip()) >= float(l[x + 1][55:60].strip())):
                                ca_xcoordinate[i].append(l[x][30:38].strip())
                                ca_ycoordinate[i].append(l[x][38:46].strip())
                                ca_zcoordinate[i].append(l[x][47:55].strip())
                                res_id[i].append(l[x][22:27].strip())
                                res_name[i].append(l[x][15:20].strip())
                                x += 2
                            elif (float(l[x][55:60].strip()) < float(l[x + 1][55:60].strip())):
                                ca_xcoordinate[i].append(l[x + 1][30:38].strip())
                                ca_ycoordinate[i].append(l[x + 1][38:46].strip())
                                ca_zcoordinate[i].append(l[x + 1][47:55].strip())
                                res_id[i].append(l[x + 1][22:27].strip())
                                res_name[i].append(l[x + 1][15:20].strip())
                                x += 2

        positions=[]
        l=[]
        f_in.close()
    return res_id,res_name,ca_xcoordinate,ca_ycoordinate,ca_zcoordinate


#---------------------------------------------------------------------------------------------------------------------------------------------------------
def pssm_val(path_to_separate_fasta_files):
    pssm_values=[]
    for i in range(len(pdb)):
        pssm_values.append([])
        l = []
        file = open(path_to_separate_fasta_files+pdb[i][1:] +'_pssm', 'r')
        for each in file:
            l.append(each.strip())
        file.close()
        temp_1=[]
        for j in range(3,(len(seq[i])+3)):
            temp_2=l[j].split(' ')
            while ("" in temp_2):
                temp_2.remove("")
            temp_1.append(temp_2[2:22])

        pssm_values[i].append(temp_1)
        l=[]
    return pssm_values

#-----------------------------------------------------------------------------------------------------------------------

def separate_amino_pssm(type_amino,pssm_values):

    pssm_sep_values=[]
    for i in range(len(pssm_values)):
        pssm_sep_values.append([])
        for j in range(len(seq[i])):
            if type_amino=='ala':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][0]))
            elif type_amino == 'gly':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][7]))
            elif type_amino == 'ile':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][9]))
            elif type_amino == 'leu':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][10]))
            elif type_amino == 'val':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][19]))
            elif type_amino == 'met':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][12]))
            elif type_amino == 'phe':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][13]))
            elif type_amino == 'trp':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][17]))
            elif type_amino == 'pro':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][14]))
            elif type_amino == 'cys':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][4]))
            elif type_amino == 'ser':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][15]))
            elif type_amino == 'thr':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][16]))
            elif type_amino == 'tyr':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][18]))
            elif type_amino == 'asn':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][2]))
            elif type_amino == 'gln':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][5]))
            elif type_amino == 'his':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][8]))
            elif type_amino == 'lys':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][11]))
            elif type_amino == 'arg':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][1]))
            elif type_amino == 'asp':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][3]))
            elif type_amino == 'glu':
                pssm_sep_values[i].append(float(pssm_values[i][0][j][6]))

    return pssm_sep_values


def pssm_flattened(not_flatt_pssm_values):
    pssm_flatt=[]
    for i in range(len(seq)):
        for j in range(len(seq[i])):
            pssm_flatt.append(not_flatt_pssm_values[i][j])
    return pssm_flatt


#------------------------------------------------------------------------------------------------------
def pssm_cal(path_to_save_pssm_abs_file):
    acids = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    acids_keys = list(acids.keys())
    acids_vals = list(acids.values())

    sep_amino_pssm={}
    for each in acids_vals:
        sep_amino_pssm['pssm_'+each.lower()+'values']=separate_amino_pssm(acids_keys[acids_vals.index(each)].lower(),pssm_values)
    sep_amino_pssm_flatt={}
    for each in sep_amino_pssm:
        sep_amino_pssm_flatt[each]=pssm_flattened(sep_amino_pssm[each])
    for each in sep_amino_pssm_flatt:
        df00= pd.DataFrame(sep_amino_pssm_flatt[each])
        df00.columns = ['2.0_bound_pssm_'+each[5]+'_abs']
        df00.to_csv(path_to_save_pssm_abs_file+'new_test_pssm_'+each[5]+'_abs.csv', index=False)




if __name__=='__main__':
    pdb, seq = pdbs_seq_info('F:\\sab_new_bound_features\\2.0\\2.0_latest\\latest_bench_test_own\\2.0_own_test_latest.fasta')
    res_id, res_name, ca_xcoordinate, ca_ycoordinate, ca_zcoordinate = ca_coordinates(seq,
                                                                                      'F:\\sab_new_bound_features\\2.0\\entire_bound_benchtest_2.0_split_pdbs\\',
                                                                                      pdb)
    pssm_values=pssm_val('F:\\sab_new_bound_features\\2.0\\2.0_latest\\new_test\\separate_fasta_files\\')
    pssm_cal('F:\\sab_new_bound_features\\2.0\\2.0_latest\\new_test\\pssm\\')

