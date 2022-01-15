import os
import pandas as pd
from collections import Counter
import regex as re
import math
import concurrent.futures
import numpy as np

file = open('//home//19bt60r19//2.0_bound//bench_test//new_2.0_benchmark_test.fasta', 'r')

pdb = []
seq = []
for each in file:
    if each.startswith('>'):
        pdb.append(each.strip())
    else:
        seq.append(each.strip())

file.close()

res_id = []
res_name = []
xcoordinate = []
ycoordinate = []
zcoordinate = []
l = []
for i in range(len(seq)):
    print(i)
    xcoordinate.append([])
    ycoordinate.append([])
    zcoordinate.append([])
    res_id.append([])
    res_name.append([])
    f_in = open('//home//19bt60r19//2.0_bound//bench_test//complexes//' + pdb[i][1:] + '.pdb', 'r')
    for x in f_in:
        if x.startswith('ATOM') and x[11:15].strip() == 'CA':
            l.append(x)
    positions = []
    for p in l:
        temp = p[22:27].strip()
        if any(c.isalpha() for c in temp) == True:
            res = [re.findall(r'(\w+?)(\d+)', temp)[0]]
            res_2 = res[0]
            res_3 = "".join(res_2)
            res_4 = int(res_3)
            positions.append(int(res_4))
        else:
            positions.append(int(temp))
    count_positions = Counter(positions)
    x = 0
    while x <= (len(l) - 1):
        res_type_check = l[x][22:27].strip()
        if x == (len(l) - 1):
            xcoordinate[i].append(l[x][30:38].strip())
            ycoordinate[i].append(l[x][38:46].strip())
            zcoordinate[i].append(l[x][47:55].strip())
            res_id[i].append(l[x][22:27].strip())
            res_name[i].append(l[x][15:20].strip())
            x += 1

        elif x != (len(l) - 1) and x < (len(l) - 1):
            if (count_positions[positions[x]] >= 1 and len(l[x][15:20].strip()) == 3):
                xcoordinate[i].append(l[x][30:38].strip())
                ycoordinate[i].append(l[x][38:46].strip())
                zcoordinate[i].append(l[x][47:55].strip())
                res_id[i].append(l[x][22:27].strip())
                res_name[i].append(l[x][15:20].strip())
                x += 1
            # elif any(c.isalpha() for c in res_type_check)==True:
            #   xcoordinate[i].append(l[x][30:39].strip())
            #   ycoordinate[i].append(l[x][38:46].strip())
            #  zcoordinate[i].append(l[x][47:55].strip())
            #  res_id[i].append(l[x][22:27].strip())
            #  res_name[i].append(l[x][15:20].strip())
            #  x += 1

            elif (count_positions[positions[x]] >= 1 and len(l[x][15:20].strip()) > 3):
                if count_positions[positions[x]] != 3:
                    if (len(l[x][15:20].strip()) > 3 and len(l[x + 1][15:20].strip()) > 3):
                        if (float(l[x][55:60].strip()) >= float(l[x + 1][55:60].strip())):
                            xcoordinate[i].append(l[x][30:38].strip())
                            ycoordinate[i].append(l[x][38:46].strip())
                            zcoordinate[i].append(l[x][47:55].strip())
                            res_id[i].append(l[x][22:27].strip())
                            res_name[i].append(l[x][15:20].strip())
                            x += 2
                        elif (float(l[x][55:60].strip()) < float(l[x + 1][55:60].strip())):
                            xcoordinate[i].append(l[x + 1][30:38].strip())
                            ycoordinate[i].append(l[x + 1][38:46].strip())
                            zcoordinate[i].append(l[x + 1][47:55].strip())
                            res_id[i].append(l[x + 1][22:27].strip())
                            res_name[i].append(l[x + 1][15:20].strip())
                            x += 2
                elif count_positions[positions[x]] == 3:
                    if (len(l[x][15:20].strip()) > 3 and len(l[x + 1][15:20].strip()) > 3 and len(
                            l[x + 2][15:20].strip()) > 3):
                        if (float(l[x][55:60].strip()) >= (
                                float(l[x + 1][55:60].strip()) and float(l[x + 2][55:60].strip()))):
                            xcoordinate[i].append(l[x][30:38].strip())
                            ycoordinate[i].append(l[x][38:46].strip())
                            zcoordinate[i].append(l[x][47:55].strip())
                            res_id[i].append(l[x][22:27].strip())
                            res_name[i].append(l[x][15:20].strip())
                            x += 3
                        elif ((float(l[x][55:60].strip()) and float(l[x + 2][55:60].strip())) < float(
                                l[x + 1][55:60].strip())):
                            xcoordinate[i].append(l[x + 1][30:38].strip())
                            ycoordinate[i].append(l[x + 1][38:46].strip())
                            zcoordinate[i].append(l[x + 1][47:55].strip())
                            res_id[i].append(l[x + 1][22:27].strip())
                            res_name[i].append(l[x + 1][15:20].strip())
                            x += 3
                        elif ((float(l[x][55:60].strip()) and float(l[x + 1][55:60].strip())) < float(
                                l[x + 2][55:60].strip())):
                            xcoordinate[i].append(l[x + 2][30:38].strip())
                            ycoordinate[i].append(l[x + 2][38:46].strip())
                            zcoordinate[i].append(l[x + 2][47:55].strip())
                            res_id[i].append(l[x + 2][22:27].strip())
                            res_name[i].append(l[x + 2][15:20].strip())
                            x += 3
                    elif (len(l[x][15:20].strip()) > 3 and len(l[x + 1][15:20].strip()) > 3 and len(
                            l[x + 2][15:20].strip()) == 3):
                        if (float(l[x][55:60].strip()) >= float(l[x + 1][55:60].strip())):
                            xcoordinate[i].append(l[x][30:38].strip())
                            ycoordinate[i].append(l[x][38:46].strip())
                            zcoordinate[i].append(l[x][47:55].strip())
                            res_id[i].append(l[x][22:27].strip())
                            res_name[i].append(l[x][15:20].strip())
                            x += 2
                        elif (float(l[x][55:60].strip()) < float(l[x + 1][55:60].strip())):
                            xcoordinate[i].append(l[x + 1][30:38].strip())
                            ycoordinate[i].append(l[x + 1][38:46].strip())
                            zcoordinate[i].append(l[x + 1][47:55].strip())
                            res_id[i].append(l[x + 1][22:27].strip())
                            res_name[i].append(l[x + 1][15:20].strip())
                            x += 2

    positions = []
    l = []
    f_in.close()




#------------------------------------------------------------------------------------------------------------------------------------------------------

def featurecalculation(seqscore, separation):
    cham830score = []
    eucdistance = []
    for i in range(len(seq)):
        eucdistance.append([])

        for j in range(len(xcoordinate[i])):
            eucdistance[i].append([])
            for k in range(len(xcoordinate[i])):
                d = math.sqrt(((float(xcoordinate[i][j]) - float(xcoordinate[i][k])) ** 2) + ((float(ycoordinate[i][j]) - float(ycoordinate[i][k])) ** 2) + ((float(zcoordinate[i][j]) - float(zcoordinate[i][k])) ** 2))
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
    cham830score.extend(finalchamscoreofeach)

    simplelistcham = []
    for i in range(len(cham830score)):
        for j in range(len(cham830score[i])):
            simplelistcham.append(cham830score[i][j])
    cham830feature = simplelistcham
    return cham830feature

def differentradii(radii):
    eucdistance = []
    for i in range(len(seq)):
        eucdistance.append([])
        for j in range(len(xcoordinate[i])):
            eucdistance[i].append([])
            for k in range(len(xcoordinate[i])):
                d = math.sqrt(((float(xcoordinate[i][j]) - float(xcoordinate[i][k])) ** 2) + ((float(ycoordinate[i][j]) - float(ycoordinate[i][k])) ** 2) + ((float(zcoordinate[i][j]) - float(zcoordinate[i][k])) ** 2))
                eucdistance[i][j].append(d)

    atomindexeucdistance = []
    for i in range(len(seq)):
        atomindexeucdistance.append([])
        for j in range(len(eucdistance[i])):
            atomindexeucdistance[i].append([])
            for k in range(len(eucdistance[i])):
                if 0 < eucdistance[i][j][k] <= radii:
                    atomindexeucdistance[i][j].append(k)

    lengthofcontacts = []
    for i in range(len(seq)):
        lengthofcontacts.append([])
        for j in range(len(seq[i])):
            lengthofcontacts[i].append(len(atomindexeucdistance[i][j]))

    meanoflengthofcontacts = []
    for i in range(len(lengthofcontacts)):
        meanoflengthofcontacts.append((sum(lengthofcontacts[i])) / len(lengthofcontacts[i]))

    relativecontactnumber = []
    for i in range(len(lengthofcontacts)):
        relativecontactnumber.append([])
        for j in range(len(lengthofcontacts[i])):
            relativecontactnumber[i].append(lengthofcontacts[i][j] - meanoflengthofcontacts[i])
    scoring_input_to_fc = {}
    for each in range(4, 21):
        scoring_input_to_fc['rcn' + '_' + str(radii) + '_' + 'avg' + '_' + str(each)] = featurecalculation(relativecontactnumber,separation=each)
    return scoring_input_to_fc



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

final_scores=[]
with concurrent.futures.ProcessPoolExecutor() as executor:
    radius_info=[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    results=executor.map(differentradii,radius_info)
    for result in results:
        final_scores.append(result)

names = ['four', 'five', 'six', 'seven', 'eight', 'nine', 'ten', 'eleven', 'twelve','thirteen', 'fourteen', 'fifteen', 'sixteen','seventeen','eighteen','nineteen','twenty']


for i in range(len(final_scores)):

  df=pd.DataFrame(final_scores[i])
  df.to_csv('//home//19bt60r19//test_bench_rcn_avg_'+names[i]+'.csv',index=False)

