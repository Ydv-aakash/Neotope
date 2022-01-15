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

def featurecalculation_hse(seqscore,separation):
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
#-------------------------------------------------------------------------------------------------------------------------------

def differentradius(radii):
    residuenames = []
    hsevalues = []
    aminoacids = ['ALA', 'CYS', 'MET', 'GLU', 'GLN', 'PRO', 'ILE', 'LEU', 'ARG', 'ASP', 'ASN', 'LYS', 'VAL', 'TYR',
                  'THR', 'TRP', 'PHE', 'GLY', 'SER', 'HIS']
    a = []
    b = []
    from Bio import PDB
    for i in range(len(seq)):

        a.append([])
        b.append([])
        pdbfile = open(path_to_complexes + pdb[i][1:] + '.pdb')
        p = PDB.PDBParser()
        s = p.get_structure('X', pdbfile)
        m = s[0]
        RADIUS = radii
        hse = PDB.HSExposureCB(m, RADIUS)
        residue_list = PDB.Selection.unfold_entities(m, 'R')
        hse = PDB.HSExposureCA(m, RADIUS)
        residue_list = PDB.Selection.unfold_entities(m, 'R')
        for r in residue_list:
            a[i].append(r.get_resname())
            b[i].append(r.xtra)

    for i in range(len(a)):
        residuenames.append([])
        hsevalues.append([])
        for j in range(len(a[i])):
            if a[i][j] in aminoacids:
                residuenames[i].append(a[i][j])
                hsevalues[i].append(b[i][j])

    seqscore = []
    for i in range(len(seq)):
        seqscore.append([])
        for j in range(len(seq[i])):
            try:
                if type_of_hse in hsevalues[i][j]:
                    seqscore[i].append(hsevalues[i][j][type_of_hse])
                else:
                    seqscore[i].append(int(0))
            except:
                print(str('i is ') + str(i) + 'and' + str('j is ') + str(j) + str(hsevalues[i][j]))
    scoring_input_to_fc = {}
    for each in range(4, 21):
        scoring_input_to_fc['hse_'+type_hse + '_' + str(radii) + '_' + 'avg' + '_' + str(each)] = featurecalculation_hse(seqscore,
                                                                                                              separation=each)
    return scoring_input_to_fc


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def hse_avg(type_hse,path_to_save_file):
    final_scores = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        radius_info = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        results = executor.map(differentradius, radius_info)
        for result in results:
            final_scores.append(result)

    names = ['four', 'five', 'six', 'seven', 'eight', 'nine', 'ten', 'eleven', 'twelve', 'thirteen', 'fourteen',
             'fifteen', 'sixteen', 'seventeen', 'eighteen', 'nineteen', 'twenty']

    for i in range(len(final_scores)):
        df = pd.DataFrame(final_scores[i])
        df.to_csv(path_to_save_file+type_hse+'_avg_' + names[i] + '.csv', index=False)


if __name__=='__main__':
    pdb, seq = pdbs_seq_info('/home/19bt60r19/latest_test/2.0_own_test_latest.fasta')
    res_id, res_name, ca_xcoordinate, ca_ycoordinate, ca_zcoordinate = ca_coordinates(seq,
                                                                                      '/home/19bt60r19/entire_bound_benchtest_2.0_split_pdbs/',
                                                                                      pdb)
    path_to_complexes = '/home/19bt60r19/entire_bound_benchtest_2.0_split_pdbs/'
    hse_types=['EXP_HSE_A_U','EXP_HSE_A_D','EXP_HSE_B_U','EXP_HSE_B_D']
    for each in hse_types:
        if each=='EXP_HSE_A_U':
            type_of_hse=each
            type_hse='au'
            hse_avg(type_hse, '//home//19bt60r19//new_2.0_test//new_test_bench_hse_')
        elif each=='EXP_HSE_A_D':
            type_of_hse=each
            type_hse = 'ad'
            hse_avg(type_hse, '//home//19bt60r19//new_2.0_test//new_test_bench_hse_')
        elif each=='EXP_HSE_B_U':
            type_of_hse=each
            type_hse = 'bu'
            hse_avg(type_hse, '//home//19bt60r19//new_2.0_test//new_test_bench_hse_')
        elif each=='EXP_HSE_B_D':
            type_of_hse=each
            type_hse = 'bd'
            hse_avg(type_hse, '//home//19bt60r19//new_2.0_test//new_test_bench_hse_')
        
