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
def bfactor_calc(pdb,path_to_complexes):
    bfactor = []

    for i in range(len(pdb)):
        bfactor.append([])
        file_in = open(path_to_complexes + pdb[i][1:] + '.pdb', 'r')
        l = []
        for x in file_in:
            if x.startswith('ATOM'):
                l.append(x)
        file_in.close()
        for j in range(len(seq[i])):
            temp = []
            for k in range(len(l)):
                if (res_id[i][j] == l[k][22:27].strip() and res_name[i][j] == l[k][16:20].strip()):
                    temp.append(float(l[k][60:68].strip()))
            bfactor[i].append(temp)

    final_bfactor = []

    for i in range(len(bfactor)):
        final_bfactor.append([])
        for j in range(len(bfactor[i])):
            temp = (sum(bfactor[i][j]) / len(bfactor[i][j]))
            final_bfactor[i].append(temp)
    return final_bfactor

#--------------------------------------------------------------------------------------------------

def featurecalculation_bfac(separation):
    scores = []

    eucdistance=[]
    for i in range(len(seq)):
        eucdistance.append([])

        for j in range(len(ca_xcoordinate[i])):
            eucdistance[i].append([])
            for k in range(len(ca_xcoordinate[i])):
                d=math.sqrt(((float(ca_xcoordinate[i][j])-float(ca_xcoordinate[i][k]))**2)+ ((float(ca_ycoordinate[i][j])-float(ca_ycoordinate[i][k]))**2)+ ((float(ca_zcoordinate[i][j])-float(ca_zcoordinate[i][k]))**2))
                eucdistance[i][j].append(d)
    atomindexeucdistance=[]
    for i in range(len(seq)):
        atomindexeucdistance.append([])
        for j in range(len(eucdistance[i])):
            atomindexeucdistance[i].append([])
            for k in range(len(eucdistance[i])):
                if eucdistance[i][j][k]<=separation:
                    atomindexeucdistance[i][j].append(k)

    finalscoreofeach=[]
    for i in range(len(seq)):
        finalscoreofeach.append([])
        try:
            for j in range(len(atomindexeucdistance[i])):
                addscore = 0
                for k in range(len(atomindexeucdistance[i][j])):
                    addscore+=bfac[i][atomindexeucdistance[i][j][k]]
                c = addscore / len(atomindexeucdistance[i][j])
                finalscoreofeach[i].append(c)

        except: print('i is'+ str(i) + 'j is' + str(j) + 'k is' + str(k) + ' has problem')
    scores.extend(finalscoreofeach)

    final_feature=[]
    for i in range(len(scores)):
        for j in range(len(scores[i])):
            final_feature.append(scores[i][j])

    return final_feature

#-------------------------------------------------------------------------------------------------------------

def bfactor_avg(path_to_save_bfactor_file):
    final_scores = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        separation = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        results = executor.map(featurecalculation_bfac, separation)
        for result in results:
            final_scores.append(result)

    df = pd.DataFrame(final_scores[0])
    for i in range(len(final_scores)):
        if i != 0:
            df1 = pd.DataFrame(final_scores[i])
            df = pd.concat([df, df1], axis=1)

    column_names = []
    for each in separation:
        column_names.append('2.0_bfac_avg_' + str(each))

    df.columns = column_names
    df.to_csv(path_to_save_bfactor_file, index=False)

def bfac_abs(final_bfactor,path_to_save_bfactor_absolute):
    bfeature = []
    for i in range(len(final_bfactor)):
        for j in range(len(final_bfactor[i])):
            bfeature.append(final_bfactor[i][j])

    df = pd.DataFrame(bfeature)
    df.columns = ['2.0_bound_bfac']
    df.to_csv(path_to_save_bfactor_absolute, index=False)


#------------------------------------------------------------------------------
if __name__=="__main__":
    pdb, seq = pdbs_seq_info('/home/19bt60r19/latest_test/2.0_own_test_latest.fasta')
    res_id, res_name, ca_xcoordinate, ca_ycoordinate, ca_zcoordinate = ca_coordinates(seq,'/home/19bt60r19/entire_bound_benchtest_2.0_split_pdbs/',pdb)
    bfac = bfactor_calc(pdb, '/home/19bt60r19/entire_bound_benchtest_2.0_split_pdbs/')

    bfactor_avg('//home//19bt60r19//new_2.0_test//test_bench_bfac_avg.csv')
    bfac_abs(bfac,'//home//19bt60r19//new_2.0_test//test_bench_bfac_abs.csv')
