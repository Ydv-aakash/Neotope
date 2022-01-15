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
def compo(seq):
    epico = []
    nonepico = []
    for i in range(len(seq)):
        count_epi = 0
        count_nonepi = 0
        epico.append([])
        nonepico.append([])
        for j in range(len(seq[i])):
            if seq[i][j].isupper():
                count_epi += 1
            else:
                count_nonepi += 1
        epico[i].append(count_epi)
        nonepico[i].append(count_nonepi)

    ala_epi = []
    arg_epi = []
    asn_epi = []
    asp_epi = []
    cys_epi = []
    glu_epi = []
    gln_epi = []
    gly_epi = []
    his_epi = []
    ile_epi = []
    leu_epi = []
    lys_epi = []
    met_epi = []
    phe_epi = []
    pro_epi = []
    ser_epi = []
    thr_epi = []
    trp_epi = []
    tyr_epi = []
    val_epi = []

    for i in range(len(seq)):
        count_epi_ala = 0
        count_epi_arg = 0
        count_epi_asn = 0
        count_epi_asp = 0
        count_epi_cys = 0
        count_epi_glu = 0
        count_epi_gln = 0
        count_epi_gly = 0
        count_epi_his = 0
        count_epi_ile = 0
        count_epi_leu = 0
        count_epi_lys = 0
        count_epi_met = 0
        count_epi_phe = 0
        count_epi_pro = 0
        count_epi_ser = 0
        count_epi_thr = 0
        count_epi_trp = 0
        count_epi_tyr = 0
        count_epi_val = 0

        ala_epi.append([])
        arg_epi.append([])
        asn_epi.append([])
        asp_epi.append([])
        cys_epi.append([])
        glu_epi.append([])
        gln_epi.append([])
        gly_epi.append([])
        his_epi.append([])
        ile_epi.append([])
        leu_epi.append([])
        lys_epi.append([])
        met_epi.append([])
        phe_epi.append([])
        pro_epi.append([])
        ser_epi.append([])
        thr_epi.append([])
        trp_epi.append([])
        tyr_epi.append([])
        val_epi.append([])

        for j in range(len(seq[i])):
            if seq[i][j] == 'A':
                count_epi_ala += 1
            if seq[i][j] == 'L':
                count_epi_leu += 1
            if seq[i][j] == 'R':
                count_epi_arg += 1
            if seq[i][j] == 'K':
                count_epi_lys += 1
            if seq[i][j] == 'N':
                count_epi_asn += 1
            if seq[i][j] == 'M':
                count_epi_met += 1
            if seq[i][j] == 'D':
                count_epi_asp += 1
            if seq[i][j] == 'F':
                count_epi_phe += 1
            if seq[i][j] == 'C':
                count_epi_cys += 1
            if seq[i][j] == 'P':
                count_epi_pro += 1
            if seq[i][j] == 'Q':
                count_epi_gln += 1
            if seq[i][j] == 'S':
                count_epi_ser += 1
            if seq[i][j] == 'E':
                count_epi_glu += 1
            if seq[i][j] == 'T':
                count_epi_thr += 1
            if seq[i][j] == 'G':
                count_epi_gly += 1
            if seq[i][j] == 'W':
                count_epi_trp += 1
            if seq[i][j] == 'H':
                count_epi_his += 1
            if seq[i][j] == 'Y':
                count_epi_tyr += 1
            if seq[i][j] == 'I':
                count_epi_ile += 1
            if seq[i][j] == 'V':
                count_epi_val += 1

        ala_epi[i].append(count_epi_ala)
        arg_epi[i].append(count_epi_arg)
        asn_epi[i].append(count_epi_asn)
        asp_epi[i].append(count_epi_asp)
        cys_epi[i].append(count_epi_cys)
        glu_epi[i].append(count_epi_glu)
        gln_epi[i].append(count_epi_gln)
        gly_epi[i].append(count_epi_gly)
        his_epi[i].append(count_epi_his)
        ile_epi[i].append(count_epi_ile)
        leu_epi[i].append(count_epi_leu)
        lys_epi[i].append(count_epi_lys)
        met_epi[i].append(count_epi_met)
        phe_epi[i].append(count_epi_phe)
        pro_epi[i].append(count_epi_pro)
        ser_epi[i].append(count_epi_ser)
        thr_epi[i].append(count_epi_thr)
        trp_epi[i].append(count_epi_trp)
        tyr_epi[i].append(count_epi_tyr)
        val_epi[i].append(count_epi_val)

    ala_nonepi = []
    arg_nonepi = []
    asn_nonepi = []
    asp_nonepi = []
    cys_nonepi = []
    glu_nonepi = []
    gln_nonepi = []
    gly_nonepi = []
    his_nonepi = []
    ile_nonepi = []
    leu_nonepi = []
    lys_nonepi = []
    met_nonepi = []
    phe_nonepi = []
    pro_nonepi = []
    ser_nonepi = []
    thr_nonepi = []
    trp_nonepi = []
    tyr_nonepi = []
    val_nonepi = []
    for i in range(len(seq)):
        ala_nonepi.append([])
        arg_nonepi.append([])
        asn_nonepi.append([])
        asp_nonepi.append([])
        cys_nonepi.append([])
        glu_nonepi.append([])
        gln_nonepi.append([])
        gly_nonepi.append([])
        his_nonepi.append([])
        ile_nonepi.append([])
        leu_nonepi.append([])
        lys_nonepi.append([])
        met_nonepi.append([])
        phe_nonepi.append([])
        pro_nonepi.append([])
        ser_nonepi.append([])
        thr_nonepi.append([])
        trp_nonepi.append([])
        tyr_nonepi.append([])
        val_nonepi.append([])

        count_nonepi_ala = 0
        count_nonepi_arg = 0
        count_nonepi_asn = 0
        count_nonepi_asp = 0
        count_nonepi_cys = 0
        count_nonepi_glu = 0
        count_nonepi_gln = 0
        count_nonepi_gly = 0
        count_nonepi_his = 0
        count_nonepi_ile = 0
        count_nonepi_leu = 0
        count_nonepi_lys = 0
        count_nonepi_met = 0
        count_nonepi_phe = 0
        count_nonepi_pro = 0
        count_nonepi_ser = 0
        count_nonepi_thr = 0
        count_nonepi_trp = 0
        count_nonepi_tyr = 0
        count_nonepi_val = 0
        for j in range(len(seq[i])):
            if seq[i][j] == 'a':
                count_nonepi_ala += 1
            if seq[i][j] == 'l':
                count_nonepi_leu += 1
            if seq[i][j] == 'r':
                count_nonepi_arg += 1
            if seq[i][j] == 'k':
                count_nonepi_lys += 1
            if seq[i][j] == 'n':
                count_nonepi_asn += 1
            if seq[i][j] == 'm':
                count_nonepi_met += 1
            if seq[i][j] == 'd':
                count_nonepi_asp += 1
            if seq[i][j] == 'f':
                count_nonepi_phe += 1
            if seq[i][j] == 'c':
                count_nonepi_cys += 1
            if seq[i][j] == 'p':
                count_nonepi_pro += 1
            if seq[i][j] == 'q':
                count_nonepi_gln += 1
            if seq[i][j] == 's':
                count_nonepi_ser += 1
            if seq[i][j] == 'e':
                count_nonepi_glu += 1
            if seq[i][j] == 't':
                count_nonepi_thr += 1
            if seq[i][j] == 'g':
                count_nonepi_gly += 1
            if seq[i][j] == 'w':
                count_nonepi_trp += 1
            if seq[i][j] == 'h':
                count_nonepi_his += 1
            if seq[i][j] == 'y':
                count_nonepi_tyr += 1
            if seq[i][j] == 'i':
                count_nonepi_ile += 1
            if seq[i][j] == 'v':
                count_nonepi_val += 1

        ala_nonepi[i].append(count_nonepi_ala)
        arg_nonepi[i].append(count_nonepi_arg)
        asn_nonepi[i].append(count_nonepi_asn)
        asp_nonepi[i].append(count_nonepi_asp)
        cys_nonepi[i].append(count_nonepi_cys)
        glu_nonepi[i].append(count_nonepi_glu)
        gln_nonepi[i].append(count_nonepi_gln)
        gly_nonepi[i].append(count_nonepi_gly)
        his_nonepi[i].append(count_nonepi_his)
        ile_nonepi[i].append(count_nonepi_ile)
        leu_nonepi[i].append(count_nonepi_leu)
        lys_nonepi[i].append(count_nonepi_lys)
        met_nonepi[i].append(count_nonepi_met)
        phe_nonepi[i].append(count_nonepi_phe)
        pro_nonepi[i].append(count_nonepi_pro)
        ser_nonepi[i].append(count_nonepi_ser)
        thr_nonepi[i].append(count_nonepi_thr)
        trp_nonepi[i].append(count_nonepi_trp)
        tyr_nonepi[i].append(count_nonepi_tyr)
        val_nonepi[i].append(count_nonepi_val)

    seqscore=[]
    for i in range(len(seq)):
        seqscore.append([])
        for j in range(len(seq[i])):
            if seq[i][j]=='A' or seq[i][j]=='a':
                try:
                    ala = ((ala_epi[i][0] / epico[i][0]) / (ala_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(ala)

                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'L' or seq[i][j] == 'l':
                try:
                    leu = ((leu_epi[i][0] / epico[i][0]) / (leu_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(leu)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'R' or seq[i][j] == 'r':
                try:
                    arg = ((arg_epi[i][0] / epico[i][0]) / (arg_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(arg)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'K' or seq[i][j] == 'k':
                try:

                    lys = ((lys_epi[i][0] / epico[i][0]) / (lys_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(lys)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'N' or seq[i][j] == 'n':
                try:
                    asn = ((asn_epi[i][0] / epico[i][0]) / (asn_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(asn)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'M' or seq[i][j] == 'm':
                try:
                    met = ((met_epi[i][0]/ epico[i][0]) / (met_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(met)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'D' or seq[i][j] == 'd':
                try:
                    asp = ((asp_epi[i][0] / epico[i][0]) / (asp_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(asp)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'F' or seq[i][j] == 'f':
                try:
                    phe = ((phe_epi[i][0] / epico[i][0]) / (phe_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(phe)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'C' or seq[i][j] == 'c':
                try:
                    cys = ((cys_epi[i][0] / epico[i][0]) / (cys_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(cys)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'P' or seq[i][j] == 'p':
                try:

                    pro = ((pro_epi[i][0] / epico[i][0]) / (pro_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(pro)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'Q' or seq[i][j] == 'q':
                try:
                    gln = ((gln_epi[i][0] / epico[i][0]) / (gln_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(gln)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'S' or seq[i][j] == 's':
                try:
                    ser = ((ser_epi[i][0] / epico[i][0]) / (ser_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(ser)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'E' or seq[i][j] == 'e':
                try:
                    glu = ((glu_epi[i][0] / epico[i][0]) / (glu_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(glu)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'T' or seq[i][j] == 't':
                try:
                    thr = ((thr_epi[i][0] / epico[i][0]) / (thr_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(thr)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'G' or seq[i][j] == 'g':
                try:
                    gly = ((gly_epi[i][0] / epico[i][0]) / (gly_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(gly)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'W' or seq[i][j] == 'w':
                try:
                    trp = ((trp_epi[i][0] / epico[i][0]) / (trp_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(trp)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'H' or seq[i][j] == 'h':
                try:
                    his = ((his_epi[i][0] / epico[i][0]) / (his_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(his)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'Y' or seq[i][j] == 'y':
                try:
                    tyr = ((tyr_epi[i][0] / epico[i][0]) / (tyr_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(tyr)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'I' or seq[i][j] == 'i':
                try:
                    ile = ((ile_epi[i][0] / epico[i][0]) / (ile_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(ile)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

            if seq[i][j] == 'V' or seq[i][j] == 'v':
                try:

                    val = ((val_epi[i][0] / epico[i][0]) / (val_nonepi[i][0] / nonepico[i][0]))
                    seqscore[i].append(val)
                except ZeroDivisionError:
                    seqscore[i].append(float(0))

    return seqscore

#----------------------------------------------------------------------------------------------------



def featurecalculation_compo(separation):
    scorez = []
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

    finalchamscoreofeach=[]
    for i in range(len(seq)):
        finalchamscoreofeach.append([])
        try:
            for j in range(len(atomindexeucdistance[i])):
                addscore = 0
                for k in range(len(atomindexeucdistance[i][j])):
                    addscore+=seqscore[i][atomindexeucdistance[i][j][k]]
                c = addscore / len(atomindexeucdistance[i][j])
                finalchamscoreofeach[i].append(c)

        except: print('i is'+ str(i) + 'j is' + str(j) + 'k is' + str(k) + ' has problem')
    scorez.extend(finalchamscoreofeach)

    final_feature=[]
    for i in range(len(scorez)):
        for j in range(len(scorez[i])):
            final_feature.append(scorez[i][j])
    
    return final_feature
#---------------------------------------------------------------------------------------------------------


def compo_avg(path_to_save_bfactor_file):
    final_scores = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        separation = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        results = executor.map(featurecalculation_compo, separation)
        for result in results:
            final_scores.append(result)

    df = pd.DataFrame(final_scores[0])
    for i in range(len(final_scores)):
        if i != 0:
            df1 = pd.DataFrame(final_scores[i])
            df = pd.concat([df, df1], axis=1)

    column_names = []
    for each in separation:
        column_names.append('2.0_new_aacompo_avg_' + str(each))

    df.columns = column_names

    df.to_csv(path_to_save_bfactor_file,index=False)


def compo_abs(seqscore,path_tosave_abs_file):
    scores = []
    for i in range(len(seqscore)):
        for j in range(len(seqscore[i])):
            scores.append(seqscore[i][j])

    df2 = pd.DataFrame(scores)
    df2.columns = ['2.0_new_aacompo_abs']
    df2.to_csv(path_tosave_abs_file,index=False)



if __name__=="__main__":
    pdb, seq = pdbs_seq_info('/home/19bt60r19/latest_test/2.0_own_test_latest.fasta')
    res_id, res_name, ca_xcoordinate, ca_ycoordinate, ca_zcoordinate = ca_coordinates(seq,'/home/19bt60r19/entire_bound_benchtest_2.0_split_pdbs/',pdb)
    seqscore = compo(seq)

    compo_avg('//home//19bt60r19//new_2.0_test//new_2.0_aacompo_test_avg.csv')
    compo_abs(seqscore, '//home//19bt60r19//new_2.0_test//new_2.0_test_aacompo_abs.csv')





