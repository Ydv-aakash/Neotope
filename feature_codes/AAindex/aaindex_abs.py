import os
import pandas as pd
from collections import Counter
import regex as re
import math
import numpy as np
import concurrent.futures



file=open('//home//19bt60r19//2.0_bound//dataset//2.0_ultimate_bound_train.fasta','r')

pdb=[]
seq=[]
for each in file:
    if each.startswith('>'):
     pdb.append(each.strip())
    else:
        seq.append(each.strip())

file.close()


res_id=[]
res_name=[]
xcoordinate=[]
ycoordinate=[]
zcoordinate=[]
l = []
for i in range(len(seq)):
    print(i)
    xcoordinate.append([])
    ycoordinate.append([])
    zcoordinate.append([])
    res_id.append([])
    res_name.append([])
    f_in = open('//home//19bt60r19//2.0_bound//complexes//' + pdb[i][1:]+ '.pdb', 'r')
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
            xcoordinate[i].append(l[x][30:38].strip())
            ycoordinate[i].append(l[x][38:46].strip())
            zcoordinate[i].append(l[x][47:55].strip())
            res_id[i].append(l[x][22:27].strip())
            res_name[i].append(l[x][15:20].strip())
            x += 1

        elif x != (len(l) - 1) and x < (len(l) - 1):
            if (count_positions[positions[x]]>= 1 and len(l[x][15:20].strip()) == 3 ):
                xcoordinate[i].append(l[x][30:38].strip())
                ycoordinate[i].append(l[x][38:46].strip())
                zcoordinate[i].append(l[x][47:55].strip())
                res_id[i].append(l[x][22:27].strip())
                res_name[i].append(l[x][15:20].strip())
                x += 1
            #elif any(c.isalpha() for c in res_type_check)==True:
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
                            res_id[i].append(l[x+ 1][22:27].strip())
                            res_name[i].append(l[x+1][15:20].strip())
                            x += 2
                elif count_positions[positions[x]] == 3:
                    if (len(l[x][15:20].strip()) > 3 and len(l[x + 1][15:20].strip()) > 3 and len(l[x + 2][15:20].strip()) > 3 ):
                        if (float(l[x][55:60].strip()) >= (float(l[x + 1][55:60].strip()) and float(l[x + 2][55:60].strip()))):
                            xcoordinate[i].append(l[x][30:38].strip())
                            ycoordinate[i].append(l[x][38:46].strip())
                            zcoordinate[i].append(l[x][47:55].strip())
                            res_id[i].append(l[x][22:27].strip())
                            res_name[i].append(l[x][15:20].strip())
                            x += 3
                        elif ((float(l[x][55:60].strip()) and float(l[x + 2][55:60].strip())) < float(l[x + 1][55:60].strip())):
                            xcoordinate[i].append(l[x + 1][30:38].strip())
                            ycoordinate[i].append(l[x + 1][38:46].strip())
                            zcoordinate[i].append(l[x + 1][47:55].strip())
                            res_id[i].append(l[x+ 1][22:27].strip())
                            res_name[i].append(l[x+1][15:20].strip())
                            x += 3
                        elif ((float(l[x][55:60].strip()) and float(l[x + 1][55:60].strip())) < float(l[x + 2][55:60].strip())):
                            xcoordinate[i].append(l[x + 2][30:38].strip())
                            ycoordinate[i].append(l[x + 2][38:46].strip())
                            zcoordinate[i].append(l[x + 2][47:55].strip())
                            res_id[i].append(l[x+ 2][22:27].strip())
                            res_name[i].append(l[x+2][15:20].strip())
                            x += 3
                    elif (len(l[x][15:20].strip()) > 3 and len(l[x + 1][15:20].strip()) > 3 and len(l[x + 2][15:20].strip()) == 3):
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

    positions=[]
    l=[]
    f_in.close()


#------------------------------------------------------------------------------------------------------------------------------------------------------
file2=open('//home//19bt60r19//2.0_bound//aaindex_na.txt','r')

l=[]
for x in file2:
    l.append(x.strip())
file2.close()
featurenames=[]
fvalues=[]
count=0
for x in range(len(l)):
    if x == 0:
        if (l[x][0]=='H' and l[x+1][0]=='D'):
            featurenames.append(l[x][1:].strip())
        if (l[x][0] == 'I' and l[x + 3] == '//'):
            fvalues.append([])
            fvalues[count].append(l[x + 1].strip())
            fvalues[count].append(l[x + 2].strip())
            count += 1
    if x>1:
        if (l[x][0]=='H' and l[x+1][0]=='D'and l[x-1]=='//'):
            featurenames.append(l[x][1:].strip())
        if (l[x][0]=='I' and l[x+3]=='//'):
            fvalues.append([])
            fvalues[count].append(l[x+  1].strip())
            fvalues[count].append(l[x + 2].strip())
            count += 1
#------------------------------------------------------------------------------------------------------------------------------------------------------------
final=[]
featurehaving_na=[]

for i in range(len(fvalues)):
    try:
        temp2 = []
        final.append([])
        a = "  ".join(fvalues[i])
        for j in range(len(a)):
            temp=list(a)
            if (a[j]==' ' and a[j-1] != ' '):
                v="".join(temp2)
                final[i].append(float(v.strip()))
                temp2 = []
            elif (type(a[j])==str or a[j]=='.'):
                temp2.append((a[j]))
            if j==(len(a)-1):
                temp3="".join(temp2)
                final[i].append(float(temp3.strip()))
    except:
        featurehaving_na.append(i)
        print('featurehaving_na is = '+str(i))


feature_indices_required=[]
for i in range(len(fvalues)):
    if i not in featurehaving_na:
        feature_indices_required.append(i)

featuresvaluesafter_removed_na= {}
test=[]
for i in feature_indices_required:
    try:
        temp2 = []
        finaltemp=[]
        a = "  ".join(fvalues[i])
        for j in range(len(a)):
            temp=list(a)
            if (a[j]==' ' and a[j-1] != ' '):
                v="".join(temp2)
                finaltemp.append(float(v.strip()))
                featuresvaluesafter_removed_na[featurenames[i]]=finaltemp
                temp2 = []
            elif (type(a[j])==str or a[j]=='.'):
                temp2.append((a[j]))
            if j==(len(a)-1):
                temp3="".join(temp2)
                finaltemp.append(float(temp3.strip()))
                featuresvaluesafter_removed_na[featurenames[i]]=finaltemp

    except:
        test.append(i)
        print('error'+ str(i))


###################################################################################

def sequencescore(score):
    seqscore=[]
    for i in range(len(seq)):
        seqscore.append([])
        for j in range(len(seq[i])):
            if seq[i][j]=='A' or seq[i][j]=='a':
                seqscore[i].append(score[0])
            if seq[i][j] == 'L' or seq[i][j] == 'l':
                seqscore[i].append(score[10])
            if seq[i][j] == 'R' or seq[i][j] == 'r':
                seqscore[i].append(score[1])
            if seq[i][j] == 'K' or seq[i][j] == 'k':
                seqscore[i].append(score[11])
            if seq[i][j] == 'N' or seq[i][j] == 'n':
                seqscore[i].append(score[2])
            if seq[i][j] == 'M' or seq[i][j] == 'm':
                seqscore[i].append(score[12])
            if seq[i][j] == 'D' or seq[i][j] == 'd':
                seqscore[i].append(score[3])
            if seq[i][j] == 'F' or seq[i][j] == 'f':
                seqscore[i].append(score[13])
            if seq[i][j] == 'C' or seq[i][j] == 'c':
                seqscore[i].append(score[4])
            if seq[i][j] == 'P' or seq[i][j] == 'p':
                seqscore[i].append(score[14])
            if seq[i][j] == 'Q' or seq[i][j] == 'q':
                seqscore[i].append(score[5])
            if seq[i][j] == 'S' or seq[i][j] == 's':
                seqscore[i].append(score[15])
            if seq[i][j] == 'E' or seq[i][j] == 'e':
                seqscore[i].append(score[6])
            if seq[i][j] == 'T' or seq[i][j] == 't':
                seqscore[i].append(score[16])
            if seq[i][j] == 'G' or seq[i][j] == 'g':
                seqscore[i].append(score[7])
            if seq[i][j] == 'W' or seq[i][j] == 'w':
                seqscore[i].append(score[17])
            if seq[i][j] == 'H' or seq[i][j] == 'h':
                seqscore[i].append(score[8])
            if seq[i][j] == 'Y' or seq[i][j] == 'y':
                seqscore[i].append(score[18])
            if seq[i][j] == 'I' or seq[i][j] == 'i':
                seqscore[i].append(score[9])
            if seq[i][j] == 'V' or seq[i][j] == 'v':
                seqscore[i].append(score[19])
    return seqscore


def featurecalculation(seqscore,separation):
    cham830score = []

    eucdistance=[]
    for i in range(len(seq)):
        eucdistance.append([])

        for j in range(len(xcoordinate[i])):
            eucdistance[i].append([])
            for k in range(len(xcoordinate[i])):
                d=math.sqrt(((float(xcoordinate[i][j])-float(xcoordinate[i][k]))**2)+ ((float(ycoordinate[i][j])-float(ycoordinate[i][k]))**2)+ ((float(zcoordinate[i][j])-float(zcoordinate[i][k]))**2))
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
    cham830score.extend(finalchamscoreofeach)

    simplelistcham=[]
    for i in range(len(cham830score)):
        for j in range(len(cham830score[i])):
            simplelistcham.append(cham830score[i][j])
    cham830feature=simplelistcham
    return cham830feature



a4={}
a5={}
a6={}
a7={}
a8={}
a9={}
a10={}
a11={}
a12={}
a13={}
a14={}
a15={}
a16={}

a17={}
a18={}
a19={}
a20={}
countofloop=0
for i in feature_indices_required:
    a4[featurenames[i]+str('_avg_four')+'2.0_bound']=featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]),separation=4)
df = pd.DataFrame(a4)
df.to_csv('//home//19bt60r19//2.0_aaindi(4)_bound_avg.csv', index=False)

for i in feature_indices_required:
    a5[featurenames[i] + str('_avg_five')+'2.0_bound'] = featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]), separation=5)
df1 = pd.DataFrame(a5)
df1.to_csv('//home//19bt60r19//2.0_aaindi(5)_bound_avg.csv', index=False)

for i in feature_indices_required:
    a6[featurenames[i]+str('_avg_six')+'2.0_bound']=featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]),separation=6)
df2 = pd.DataFrame(a6)
df2.to_csv('//home//19bt60r19//2.0_aaindi(6)_bound_avg.csv', index=False)

for i in feature_indices_required:
    a7[featurenames[i] + str('_avg_seven')+'2.0_bound'] = featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]), separation=7)
df3 = pd.DataFrame(a7)
df3.to_csv('//home//19bt60r19//2.0_aaindi(7)_bound_avg.csv', index=False)

for i in feature_indices_required:
    a8[featurenames[i]+str('_avg_eight')+'2.0_bound']=featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]),separation=8)
df4 = pd.DataFrame(a8)
df4.to_csv('//home//19bt60r19//2.0_aaindi(8)_bound_avg.csv', index=False)

for i in feature_indices_required:
    a9[featurenames[i] + str('_avg_nine')+'2.0_bound'] = featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]), separation=9)
df5 = pd.DataFrame(a9)
df5.to_csv('//home//19bt60r19//2.0_aaindi(9)_bound_avg.csv', index=False)

for i in feature_indices_required:
    a10[featurenames[i]+str('_avg_ten')+'2.0_bound']=featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]),separation=10)
df6 = pd.DataFrame(a10)
df6.to_csv('//home//19bt60r19//2.0_aaindi(10)_bound_avg.csv', index=False)

for i in feature_indices_required:
    a11[featurenames[i] + str('_avg_eleven')+'2.0_bound'] = featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]), separation=11)
df7 = pd.DataFrame(a11)
df7.to_csv('//home//19bt60r19//2.0_aaindi(11)_bound_avg.csv', index=False)

for i in feature_indices_required:
    a12[featurenames[i]+str('_avg_twelve')+'2.0_bound']=featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]),separation=12)
df8 = pd.DataFrame(a12)
df8.to_csv('//home//19bt60r19//2.0_aaindi(12)_bound_avg.csv', index=False)

for i in feature_indices_required:
    a13[featurenames[i] + str('_avg_thirteen')+'2.0_bound'] = featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]), separation=13)
df9 = pd.DataFrame(a13)
df9.to_csv('//home//19bt60r19//2.0_aaindi(13)_bound_avg.csv', index=False)

for i in feature_indices_required:
    a14[featurenames[i]+str('_avg_fourteen')+'2.0_bound']=featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]),separation=14)
df10 = pd.DataFrame(a14)
df10.to_csv('//home//19bt60r19//2.0_aaindi(14)_bound_avg.csv', index=False)

for i in feature_indices_required:
    a15[featurenames[i] + str('_avg_fifteen')+'2.0_bound'] = featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]), separation=15)

df11 = pd.DataFrame(a15)
df11.to_csv('//home//19bt60r19//2.0_aaindi(15)_bound_avg.csv', index=False)

for i in feature_indices_required:
    a16[featurenames[i]+str('_avg_sixteen')+'2.0_bound']=featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]),separation=16)

df12 = pd.DataFrame(a16)
df12.to_csv('//home//19bt60r19//2.0_aaindi(16)_bound_avg.csv', index=False)

for i in feature_indices_required:
    a17[featurenames[i]+str('_avg_seventeen')+'2.0_bound']=featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]),separation=17)
df13 = pd.DataFrame(a17)
df13.to_csv('//home//19bt60r19//2.0_aaindi(17)_bound_avg.csv', index=False)

for i in feature_indices_required:
    a18[featurenames[i]+str('_avg_eighteen')+'2.0_bound']=featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]),separation=18)
df14 = pd.DataFrame(a18)
df14.to_csv('//home//19bt60r19//2.0_aaindi(18)_bound_avg.csv', index=False)

for i in feature_indices_required:
    a19[featurenames[i]+str('_avg_nineteen')+'2.0_bound']=featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]),separation=19)
df15 = pd.DataFrame(a19)
df15.to_csv('//home//19bt60r19//2.0_aaindi(19)_bound_avg.csv', index=False)

for i in feature_indices_required:
    a20[featurenames[i]+str('_avg_twenty')+'2.0_bound']=featurecalculation(sequencescore(featuresvaluesafter_removed_na[featurenames[i]]),separation=20)
df16 = pd.DataFrame(a20)
df16.to_csv('//home//19bt60r19//2.0_aaindi(20)_bound_avg.csv', index=False)
