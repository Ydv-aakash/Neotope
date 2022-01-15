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

def featurecalculation(separation):
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
                if 0<eucdistance[i][j][k]<=separation:
                    atomindexeucdistance[i][j].append(k)

    lengthofcontacts = []
    for i in range(len(seq)):
        lengthofcontacts.append([])
        for j in range(len(seq[i])):
            lengthofcontacts[i].append(len(atomindexeucdistance[i][j]))

    meanoflengthofcontacts=[]
    for i in range(len(lengthofcontacts)):
        meanoflengthofcontacts.append((sum(lengthofcontacts[i]))/len(lengthofcontacts[i]))

    relativecontactnumber=[]
    for i in range(len(lengthofcontacts)):
        relativecontactnumber.append([])
        for j in range(len(lengthofcontacts[i])):
            relativecontactnumber[i].append(lengthofcontacts[i][j]-meanoflengthofcontacts[i])

    final_feature=[]
    for i in range(len(seq)):
        for j in range(len(relativecontactnumber[i])):
            final_feature.append(relativecontactnumber[i][j])

    return final_feature
#----------------------------------------------------------------------------------------------------------------------

if __name__=="__main__":
    pdb, seq = pdbs_seq_info('F:\\sab_new_bound_features\\2.0\\2.0_latest\\latest_bench_test_own\\2.0_own_test_latest.fasta')
    res_id, res_name, ca_xcoordinate, ca_ycoordinate, ca_zcoordinate = ca_coordinates(seq,
                                                                                      'F:\\sab_new_bound_features\\2.0\\2.0_latest\\new_test\\pdbs\\',
                                                                                      pdb)

    names = ['zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten', 'eleven', 'twelve','thirteen', 'fourteen', 'fifteen', 'sixteen','seventeen','eighteen','nineteen','twenty']
    relativecn = {}
    for sep in range(4, 21):
        relativecn['2.0_rcn_'+names[sep]+'_abs_bound'] = featurecalculation(separation=sep)
    df=pd.DataFrame(relativecn)
    df.to_csv('F:\\sab_new_bound_features\\2.0\\2.0_latest\\new_test\\features\\new_test_rcn_abs_bound.csv',index=False)




