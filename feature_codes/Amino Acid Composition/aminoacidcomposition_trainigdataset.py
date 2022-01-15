import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
file=open('F:\\sab_new_bound_features\\publish\\datasets\\final_datasets\\final_2.0_latest_train.fasta')

seq=[]
pdbs=[]
for eachinfo in file:
    if eachinfo.startswith('>'):
        pdbs.append(eachinfo.strip())
    else:seq.append(eachinfo.strip())

file.close()
#---------------------------------------------------------------
epitopes=[]
count=0
for i in range(len(seq)):
    for j in range(len(seq[i])):
        if seq[i][j].isupper():
            epitopes.append(seq[i][j])
            count+=1
            epitopes[i:i+count] = [''.join(epitopes[i:i+count])]
    count=0

nonepitopes=[]
count=0
for i in range(len(seq)):
    for j in range(len(seq[i])):
        if seq[i][j].islower():
            nonepitopes.append(seq[i][j])
            count+=1
            nonepitopes[i:i+count] = [''.join(nonepitopes[i:i+count])]
    count=0

count=0
for i in range(len(epitopes)):
    count=count+len(epitopes[i])

print("total epitopes are- " + str(count))
count = 0
for i in range(len(nonepitopes)):
    count = count + len(nonepitopes[i])

print("total nonepitopes are- " + str(count))


epico=''.join(epitopes)
nonepico=''.join(nonepitopes)

countepiala=0
countepiarg=0
countepiasn=0
countepiasp=0
countepicys=0
countepiglu=0
countepigln=0
countepigly=0
countepihis=0
countepiile=0
countepileu=0
countepilys=0
countepimet=0
countepiphe=0
countepipro=0
countepiser=0
countepithr=0
countepitrp=0
countepityr=0
countepival=0

for i in range(len(epico)):
    if epico[i]=='A':
        countepiala+=1
    elif epico[i]=='R':
        countepiarg+=1
    elif epico[i]=='N':
        countepiasn+=1
    elif epico[i]=='D':
        countepiasp+=1
    elif epico[i]=='C':
        countepicys+=1
    elif epico[i]=='E':
        countepiglu+=1
    elif epico[i]=='Q':
        countepigln+=1
    elif epico[i]=='G':
        countepigly+=1
    elif epico[i]=='H':
        countepihis+=1
    elif epico[i]=='I':
        countepiile+=1
    elif epico[i]=='L':
        countepileu+=1
    elif epico[i]=='K':
        countepilys+=1
    elif epico[i]=='M':
        countepimet+=1
    elif epico[i]=='F':
        countepiphe+=1
    elif epico[i]=='P':
        countepipro+=1
    elif epico[i]=='S':
        countepiser+=1
    elif epico[i]=='T':
        countepithr+=1
    elif epico[i]=='W':
        countepitrp+=1
    elif epico[i]=='Y':
        countepityr+=1
    elif epico[i]=='V':
        countepival+=1
print(countepiala+countepiarg+countepiasn+countepiasp+countepicys+countepiglu+countepigln+countepigly+countepihis+countepiile+countepileu+countepilys+countepimet+countepiphe+countepipro+
countepiser+countepithr+countepitrp+countepityr+countepival)

countnonepiala=0
countnonepiarg=0
countnonepiasn=0
countnonepiasp=0
countnonepicys=0
countnonepiglu=0
countnonepigln=0
countnonepigly=0
countnonepihis=0
countnonepiile=0
countnonepileu=0
countnonepilys=0
countnonepimet=0
countnonepiphe=0
countnonepipro=0
countnonepiser=0
countnonepithr=0
countnonepitrp=0
countnonepityr=0
countnonepival=0

for i in range(len(nonepico)):
    if nonepico[i]=='a':
        countnonepiala+=1
    elif nonepico[i]=='r':
        countnonepiarg+=1
    elif nonepico[i]=='n':
        countnonepiasn+=1
    elif nonepico[i]=='d':
        countnonepiasp+=1
    elif nonepico[i]=='c':
        countnonepicys+=1
    elif nonepico[i]=='e':
        countnonepiglu+=1
    elif nonepico[i]=='q':
        countnonepigln+=1
    elif nonepico[i]=='g':
        countnonepigly+=1
    elif nonepico[i]=='h':
        countnonepihis+=1
    elif nonepico[i]=='i':
        countnonepiile+=1
    elif nonepico[i]=='l':
        countnonepileu+=1
    elif nonepico[i]=='k':
        countnonepilys+=1
    elif nonepico[i]=='m':
        countnonepimet+=1
    elif nonepico[i]=='f':
        countnonepiphe+=1
    elif nonepico[i]=='p':
        countnonepipro+=1
    elif nonepico[i]=='s':
        countnonepiser+=1
    elif nonepico[i]=='t':
        countnonepithr+=1
    elif nonepico[i]=='w':
        countnonepitrp+=1
    elif nonepico[i]=='y':
        countnonepityr+=1
    elif nonepico[i]=='v':
        countnonepival+=1
print(countnonepiala+countnonepiarg+countnonepiasn+countnonepiasp+countnonepicys+countnonepiglu+countnonepigln+countnonepigly+countnonepihis+countnonepiile+countnonepileu+countnonepilys+countnonepimet+countnonepiphe+countnonepipro+
countnonepiser+countnonepithr+countnonepitrp+countnonepityr+countnonepival)



labels=['GLU','LYS','SER','ASP','ASN','THR','LEU','ARG','GLN', 'PRO','GLY','ALA','TYR','VAL','ILE','PHE', 'HIS', 'TRP', 'CYS','MET']

e=[countepiglu,countepilys,countepiser,countepiasp,countepiasn,countepithr,countepileu,countepiarg,countepigln,countepipro,countepigly,countepiala,countepityr,countepival,
   countepiile,countepiphe,countepihis,countepitrp,countepicys,countepimet]
n=[countnonepiglu,countnonepilys,countnonepiser,countnonepiasp,countnonepiasn,countnonepithr,countnonepileu,countnonepiarg,countnonepigln,
   countnonepipro,countnonepigly,countnonepiala,countnonepityr,countnonepival,countnonepiile,countnonepiphe,countnonepihis,countnonepitrp,countnonepicys,countnonepimet]
e2=[]
for i in e:
    i=((i/len(epico))*100)
    e2.append(i)
n2=[]
for i in n:
    i=((i/len(nonepico))*100)
    n2.append(i)

x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, e2, width, label='Epitopes')
rects2 = ax.bar(x + width/2, n2, width, label='Non-Epitopes')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Frequency(%)')
ax.set_title('Amino Acids')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()


fig.tight_layout()

plt.show()


#-----------------------------
#g-test for charge and polar aa
import numpy as np
from scipy.stats import chi2_contingency
obs=np.array([[countepiglu,countepilys,countepiser,countepithr,countepiasp,countepiasn,countepihis,countepiarg,countepigln,countepityr],
              [countnonepiglu,countnonepilys,countnonepiser,countnonepithr,countnonepiasp,countnonepiasn,countnonepihis,countnonepiarg,countnonepigln,countnonepityr]])

g, p, dof, expctd = chi2_contingency(obs, lambda_="log-likelihood")

print(g,p)
#------------------------------------------------------------