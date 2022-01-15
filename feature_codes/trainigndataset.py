import pandas as pd
file_chains=pd.read_csv('F:\\sab_new_bound_features\\publish\\epitope3d\\chainsepitope3d.csv')
file_pdbs=pd.read_csv('F:\\sab_new_bound_features\\publish\\epitope3d\\pdbsepitope3d.csv')
chains=[]
pdbs=[]
for eachchain in file_chains['chains']:
    chains.append(eachchain)
for eachpdb in file_pdbs['pdbs']:
    pdbs.append(eachpdb)
######################################################################################
#fetching fasta seq
#pdb_fasta=[]
#
#for eachprotein in pdbs:
#    file_pro=open('F:\\sab_new_bound_features\\publish\\epitope3d\\all\\'+eachprotein+'.fasta')
all_info={}

for prot in range(len(pdbs)):
    fasta_seq = []
    fasta_chain = []
    seq_dict={}
    file_pro = open('F:\\sab_new_bound_features\\publish\\epitope3d\\chains_fasta\\all\\'+pdbs[prot]+'.fasta')

    file_pro_parts=[]
    for each in file_pro:
        file_pro_parts.append(each)
    file_pro.close()

    temp_pdb=[]
    temp_seq=[]
    count=0
    tempseqq=[]
    for each in range(len(file_pro_parts)):
        if file_pro_parts[each].startswith('>'):
            temp_pdb.append(file_pro_parts[each].strip())
            count+=1
        else:tempseqq.append(file_pro_parts[each].strip())
        if (count==1 and len(temp_pdb)>1 and each!=(len(file_pro_parts)-1)):
            newtemp="".join(tempseqq)
            tempseqq=[]
            temp_seq.append(newtemp)
        count=0
        if each==(len(file_pro_parts)-1):
            tempseqq.append(file_pro_parts[each].strip())
            newtemp = "".join(tempseqq)
            temp_seq.append(newtemp)
            tempseqq = []

    for eachchain in range(len(temp_pdb)):
        if temp_pdb[eachchain][-1]==chains[prot]:
            fasta_chain.append(temp_pdb[eachchain][-1])
            fasta_seq.append(temp_seq[eachchain])

    seq_dict[chains[prot]] = fasta_seq
    all_info[pdbs[prot]] = seq_dict
#############################################################################################
lengthofall=[]
for eachprot in range(len(pdbs)):
    lengthofall.append(len(all_info[pdbs[eachprot]][chains[eachprot]]))

######################################################################################
final_info={}
for i in range(len(pdbs)):
    final_temp_seq=[]
    seq_again_dict={}
    for eachseq in all_info[pdbs[i]][chains[i]]:
        temp_seq=[]
        for amino in eachseq:
            if amino !='X':
                temp_seq.append(amino)
        final_temp_seq.append("".join(temp_seq))
    while "" in final_temp_seq:
        final_temp_seq.remove("")
    seq_again_dict[chains[i]]=final_temp_seq
    final_info[pdbs[i]]=seq_again_dict



lengthofall_new=[]
for eachprot in range(len(pdbs)):
    lengthofall_new.append(len(final_info[pdbs[eachprot]][chains[eachprot]]))
####################################################################################

file=open('F:\\sab_new_bound_features\\publish\\epitope3d\\epitope3d_training_dataset.fasta','w')

for i in range(len(final_info)):
    file.write('>'+pdbs[i]+'_'+chains[i])
    file.write('\n')
    file.write(final_info[pdbs[i]][chains[i]][0].lower())
    file.write('\n')

file.close()





