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

#----------------------------------------------------------------------------------
def rsa_avg(path_to_sasa_files_from_psaia,common_tag_from_psaia):

    l=[]
    episasa=[]
    testing=[]
    residues = []
    for i in range(len(seq)):

        episasa.append([])
        residues.append([])
        f_in = open(path_to_sasa_files_from_psaia + pdb[i][1:].strip()+common_tag_from_psaia,'r')
        ch=[]
        for x in f_in:
            ch.append(x.strip())
        for z in ch[12:]:
            episasa[i].append(z[110:122].strip())
            residues[i].append(z[95:105].strip())
        testing=ch.copy()
        ch=[]
    rsa=[]
    for i in range(len(episasa)):
        rsa.append([])
        for r in range(len(episasa[i])):
            if residues[i][r]=='ALA':
                a=episasa[i][r]
                b=float(a)
                c=b/121.0
                rsa[i].append(c)
            elif residues[i][r] == 'ARG':
                a = episasa[i][r]
                b = float(a)
                c = b / 265.0
                rsa[i].append(c)
            elif residues[i][r]=='ASN':
                a=episasa[i][r]
                b=float(a)
                c=b/187.0
                rsa[i].append(c)


            elif residues[i][r]=='ASP':
                a=episasa[i][r]
                b=float(a)
                c=b/187.0
                rsa[i].append(c)
            elif residues[i][r]=='CYS':
                a=episasa[i][r]
                b=float(a)
                c=b/148.0
                rsa[i].append(c)
            elif residues[i][r]=='GLU':
                a=episasa[i][r]
                b=float(a)
                c=b/214.0
                rsa[i].append(c)
            elif residues[i][r]=='GLN':
                a=episasa[i][r]
                b=float(a)
                c=b/214.0
                rsa[i].append(c)
            elif residues[i][r]=='GLY':
                a=episasa[i][r]
                b=float(a)
                c=b/97.0
                rsa[i].append(c)

            elif residues[i][r]=='HIS':
                a=episasa[i][r]
                b=float(a)
                c=b/216.0
                rsa[i].append(c)
            elif residues[i][r]=='ILE':
                a=episasa[i][r]
                b=float(a)
                c=b/195.0
                rsa[i].append(c)
            elif residues[i][r]=='LEU':
                a=episasa[i][r]
                b=float(a)
                c=b/191.0
                rsa[i].append(c)
            elif residues[i][r]=='LYS':
                a=episasa[i][r]
                b=float(a)
                c=b/230.0
                rsa[i].append(c)
            elif residues[i][r]=='MET':
                a=episasa[i][r]
                b=float(a)
                c=b/203.0
                rsa[i].append(c)
            elif residues[i][r]=='PHE':
                a=episasa[i][r]
                b=float(a)
                c=b/228.0
                rsa[i].append(c)
            elif residues[i][r]=='PRO':
                a=episasa[i][r]
                b=float(a)
                c=b/154.0
                rsa[i].append(c)
            elif residues[i][r]=='SER':
                a=episasa[i][r]
                b=float(a)
                c=b/143.0
                rsa[i].append(c)
            elif residues[i][r]=='THR':
                a=episasa[i][r]
                b=float(a)
                c=b/163.0
                rsa[i].append(c)
            elif residues[i][r]=='TRP':
                a=episasa[i][r]
                b=float(a)
                c=b/264.0
                rsa[i].append(c)
            elif residues[i][r]=='TYR':
                a=episasa[i][r]
                b=float(a)
                c=b/255.0
                rsa[i].append(c)
            elif residues[i][r]=='VAL':
                a=episasa[i][r]
                b=float(a)
                c=b/165.0
                rsa[i].append(c)

            a=0
            b=0
            c=0
    seqscore=[]
    for i in range(len(episasa)):
        seqscore.append([])
        for j in range(len(rsa[i])):
            seqscore[i].append(float(rsa[i][j]))

    return seqscore
#------------------------------------------------------------------------------------------------------

def featurecalculation_rsa(separation):
    scorez = []
    eucdistance = []
    for i in range(len(seq)):
        eucdistance.append([])

        for j in range(len(ca_xcoordinate[i])):
            eucdistance[i].append([])
            for k in range(len(ca_xcoordinate[i])):
                d = math.sqrt(((float(ca_xcoordinate[i][j]) - float(ca_xcoordinate[i][k])) ** 2) + (
                            (float(ca_ycoordinate[i][j]) - float(ca_ycoordinate[i][k])) ** 2) + (
                                          (float(ca_zcoordinate[i][j]) - float(ca_zcoordinate[i][k])) ** 2))
                eucdistance[i][j].append(d)
    atomindexeucdistance = []
    for i in range(len(seq)):
        atomindexeucdistance.append([])
        for j in range(len(eucdistance[i])):
            atomindexeucdistance[i].append([])
            for k in range(len(eucdistance[i])):
                if eucdistance[i][j][k] <= separation:
                    atomindexeucdistance[i][j].append(k)

    finalchamscoreofeach = []
    for i in range(len(seq)):
        finalchamscoreofeach.append([])
        try:
            for j in range(len(atomindexeucdistance[i])):
                addscore = 0
                for k in range(len(atomindexeucdistance[i][j])):
                    addscore += seqscore[i][atomindexeucdistance[i][j][k]]
                c = addscore / len(atomindexeucdistance[i][j])
                finalchamscoreofeach[i].append(c)

        except:
            print('i is' + str(i) + 'j is' + str(j) + 'k is' + str(k) + ' has problem')
    scorez.extend(finalchamscoreofeach)

    final_feature = []
    for i in range(len(scorez)):
        for j in range(len(scorez[i])):
            final_feature.append(scorez[i][j])

    return final_feature


#------------------------------------------------------------------------------------------------------

if __name__=='__main__':
    pdb, seq = pdbs_seq_info('/home/19bt60r19/latest_test/2.0_own_test_latest.fasta')
    res_id, res_name, ca_xcoordinate, ca_ycoordinate, ca_zcoordinate = ca_coordinates(seq,
                                                                                      '/home/19bt60r19/entire_bound_benchtest_2.0_split_pdbs/',
                                                                                      pdb)
    seqscore = rsa_avg('/home/19bt60r19/new_2.0_test/sasa_files/','_202107312213_unbound.tbl')
    final_scores = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        separation = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        results = executor.map(featurecalculation_rsa, separation)
        for result in results:
            final_scores.append(result)

    df = pd.DataFrame(final_scores[0])
    for i in range(len(final_scores)):
        if i != 0:
            df1 = pd.DataFrame(final_scores[i])
            df = pd.concat([df, df1], axis=1)

    column_names = []
    for each in separation:
        column_names.append('2.0_rasa_avg_' + str(each))
    df.columns = column_names
    df.to_csv('//home//19bt60r19//new_2.0_test//new_test_rasa_avg.csv', index=False)

