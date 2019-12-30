''' Shannon's Entropy as a measure of Diversity in DNA and peptides '''
''' author @rathsidd'''

'''---------------------------------------------------------------Import Block-----------------------------------------------------------------'''

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

'''-------------------------------------------------------------Define Functions---------------------------------------------------------------'''

def SeqProb(Seq_set,column_name,length):
    import numpy as np
    import pandas as pd
    print('extracting sequences and converting to dict')
    temp=Seq_set[str(column_name)].apply(lambda x:list(x)).to_dict()
    print('converting to dataframe')
    temp=pd.DataFrame.from_dict(temp,orient='index',columns=np.linspace(1,length,length).astype(int))
    freq=temp.apply(pd.value_counts).fillna(value=0)
    prob=freq.apply(lambda x:x/len(temp))
    ln_prob=prob.apply(lambda x:np.log(x))
    return(prob,ln_prob)

def SeqEntropy(Prob, logProb):
    import numpy as np
    import pandas as pd
    entropy=list((-1*Prob*logProb).fillna(value=0).sum())
    return(entropy)
        

'''-----------------------------------------------------------------Read Data------------------------------------------------------------------'''

dna_data = {}
pep_data = {}
#dna_N=pd.DataFrame()
sets = np.linspace(1,3,num=3).astype(int)
d=os.getcwd()
for i in sets:
    temp=pd.read_csv(d+'DNA Sequences/set_'+str(i)+'all_DNA_counts.csv')
#    for j in range(len(temp)):
#        if 'N' in list(temp['Unnamed: 0'].loc[j]):
#            dna_N.append(temp.iloc[j,:])
#            temp.drop([j])
    dna_data['set_'+str(i)]=temp[temp['Total']>2]
    tempd=pd.read_csv(d+'Peptide Sequences/NGS/All_peptides_Set'+str(i)+'.csv')
    tempd['Total']=tempd['CP1']+tempd['CP2']+tempd['CP3']+tempd['CE']
    pep_data['set_'+str(i)]=tempd

'''-----Nucleotide/Amino-Acid Position Specific probabilities (and logarithms) among unique sequences (not weighted by population or pans)-----'''

cols=list(dna_data['set_1'].columns)
cols[0]='Seq'
for i in sets:
    dna_data['set_'+str(i)].columns=cols
    

dna_prob={}
log_dna_prob={}
pep_prob={}
log_pep_prob={}
print('---------Calculating probabilities & Entropies of all DNA & peptides---------')
for i in sets:
    key = 'set_'+str(i)
    dna_prob[key],log_dna_prob[key] = SeqProb(dna_data[key],'Seq',length=36)
    pep_prob[key],log_pep_prob[key] = SeqProb(pep_data[key],'AA_seq',length=12)
    
# Locationwise Entropy for DNA and Pep unique sequences (not weighted by population or pans)
dna_ent={}
pep_ent={}
average_ent=[]
for i in sets:
    key = 'set_'+str(i)
    dna_ent[key] = SeqEntropy(dna_prob[key],log_dna_prob[key])
    pep_ent[key]= SeqEntropy(pep_prob[key],log_pep_prob[key])
    average_ent.append([np.mean(dna_ent[key]),np.mean(pep_ent[key])])    
dna_ent=pd.DataFrame.from_dict(dna_ent)
pep_ent=pd.DataFrame.from_dict(pep_ent)

#Locationwise Entropy based on pans (unique Seqs only, disregarding populations)
pan1_dna= {}
pan2_dna = {}
pan3_dna = {}
Elu_dna = {}
pan1_pep= {}
pan2_pep = {}
pan3_pep = {}
Elu_pep = {}

for i in sets:
    key = 'set_'+str(i)
    pan1_dna[key]=dna_data[key][dna_data[key]['pan1']>0]
    pan1_pep[key]=pep_data[key][pep_data[key]['CP1']>0]
    pan2_dna[key]=dna_data[key][dna_data[key]['pan2']>0]
    pan2_pep[key]=pep_data[key][pep_data[key]['CP2']>0]
    pan3_dna[key]=dna_data[key][dna_data[key]['pan3']>0]
    pan3_pep[key]=pep_data[key][pep_data[key]['CP3']>0]
    Elu_dna[key]=dna_data[key][dna_data[key]['Eluate']>0]
    Elu_pep[key]=pep_data[key][pep_data[key]['CE']>0]

''' Pan 1 '''

print('---------Calculating probabilities & Entropies of all in Pan1---------')
dna_pan1_prob={}
log_dna_pan1_prob={}
pep_pan1_prob={}
log_pep_pan1_prob={}
for i in sets:
    key = 'set_'+str(i)
    dna_pan1_prob[key],log_dna_pan1_prob[key] = SeqProb(pan1_dna[key],'Seq',length=36)
    pep_pan1_prob[key],log_pep_pan1_prob[key] = SeqProb(pan1_pep[key],'AA_seq',length=12)
    
# Locationwise Entropy for DNA and Pep unique sequences (not weighted by population or pans)
dna_pan1_ent={}
pep_pan1_ent={}
average_pan1_ent=[]
for i in sets:
    key = 'set_'+str(i)
    dna_pan1_ent[key] = SeqEntropy(dna_pan1_prob[key],log_dna_pan1_prob[key])
    pep_pan1_ent[key] = SeqEntropy(pep_pan1_prob[key],log_pep_pan1_prob[key])
    average_pan1_ent.append([np.mean(dna_pan1_ent[key]),np.mean(pep_pan1_ent[key])])    
dna_pan1_ent=pd.DataFrame.from_dict(dna_pan1_ent)
pep_pan1_ent=pd.DataFrame.from_dict(pep_pan1_ent)

''' pan 2 '''

print('---------Calculating probabilities & Entropies of all in Pan2---------')
dna_pan2_prob={}
log_dna_pan2_prob={}
pep_pan2_prob={}
log_pep_pan2_prob={}
for i in sets:
    key = 'set_'+str(i)
    dna_pan2_prob[key],log_dna_pan2_prob[key] = SeqProb(pan2_dna[key],'Seq',length=36)
    pep_pan2_prob[key],log_pep_pan2_prob[key] = SeqProb(pan2_pep[key],'AA_seq',length=12)
    
# Locationwise Entropy for DNA and Pep unique sequences (not weighted by population or pans)
dna_pan2_ent={}
pep_pan2_ent={}
average_pan2_ent=[]
for i in sets:
    key = 'set_'+str(i)
    dna_pan2_ent[key] = SeqEntropy(dna_pan2_prob[key],log_dna_pan2_prob[key])
    pep_pan2_ent[key] = SeqEntropy(pep_pan2_prob[key],log_pep_pan2_prob[key])
    average_pan2_ent.append([np.mean(dna_pan2_ent[key]),np.mean(pep_pan2_ent[key])])    
dna_pan2_ent=pd.DataFrame.from_dict(dna_pan2_ent)
pep_pan2_ent=pd.DataFrame.from_dict(pep_pan2_ent)

''' pan 3'''

print('---------Calculating probabilities & Entropies of all in Pan3---------')
dna_pan3_prob={}
log_dna_pan3_prob={}
pep_pan3_prob={}
log_pep_pan3_prob={}
for i in sets:
    key = 'set_'+str(i)
    dna_pan3_prob[key],log_dna_pan3_prob[key] = SeqProb(pan3_dna[key],'Seq',length=36)
    pep_pan3_prob[key],log_pep_pan3_prob[key] = SeqProb(pan3_pep[key],'AA_seq',length=12)
    
# Locationwise Entropy for DNA and Pep unique sequences (not weighted by population or pans)
dna_pan3_ent={}
pep_pan3_ent={}

average_pan3_ent=[]
for i in sets:
    key = 'set_'+str(i)
    dna_pan3_ent[key] = SeqEntropy(dna_pan3_prob[key],log_dna_pan3_prob[key])
    pep_pan3_ent[key] = SeqEntropy(pep_pan3_prob[key],log_pep_pan3_prob[key])
    average_pan3_ent.append([np.mean(dna_pan3_ent[key]),np.mean(pep_pan3_ent[key])])    
dna_pan3_ent=pd.DataFrame.from_dict(dna_pan3_ent)
pep_pan3_ent=pd.DataFrame.from_dict(pep_pan3_ent)

''' Eluate '''
print('---------Calculating probabilities & Entropies of all in Eluate---------')
dna_elu_prob={}
log_dna_elu_prob={}
pep_elu_prob={}
log_pep_elu_prob={}
for i in sets:
    key = 'set_'+str(i)
    dna_elu_prob[key],log_dna_elu_prob[key] = SeqProb(Elu_dna[key],'Seq',length=36)
    pep_elu_prob[key],log_pep_elu_prob[key] = SeqProb(Elu_pep[key],'AA_seq',length=12)
    
# Locationwise Entropy for DNA and Pep unique sequences (not weighted by population or pans)
dna_elu_ent={}
pep_elu_ent={}
average_elu_ent=[]
for i in sets:
    key = 'set_'+str(i)
    dna_elu_ent[key] = SeqEntropy(dna_elu_prob[key],log_dna_elu_prob[key])
    pep_elu_ent[key] = SeqEntropy(pep_elu_prob[key],log_pep_elu_prob[key])
    average_elu_ent.append([np.mean(dna_elu_ent[key]),np.mean(pep_elu_ent[key])])    
dna_elu_ent=pd.DataFrame.from_dict(dna_elu_ent)
pep_elu_ent=pd.DataFrame.from_dict(pep_elu_ent)
