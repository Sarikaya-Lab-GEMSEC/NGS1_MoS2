''' Correlation between DNA sequences from Biological and Technical Replicates '''
''' author @rathsidd '''

'''---------------------------------------------------------------Import Block-----------------------------------------------------------------'''

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

'''-----------------------------------------------------------------Read Data------------------------------------------------------------------'''

dna_data = {}
sets = ['set_1','set_2','set_3']
pans = ['pan1','pan2','pan3','Eluate']
reps = ['a','b']

d = os.getcwd()

for i in sets:
    dna_data[i]={}
    for j in pans:
        dna_data[i][j]={}
        for k in reps:
            dna_data[i][j][k]=pd.read_csv(d+'/DNA_NGS_1/'+i+'_'+j+'_'+k+'.tsv.tsv',sep='\t')
            col = list(dna_data[i][j][k].columns)
            col[0]='Seq'
            dna_data[i][j][k].columns=col
            dna_data[i][j][k].set_index('Seq',inplace=True)

'''--------------------------------------------------------------merge tech reps---------------------------------------------------------------'''

dna_ab={}
corr_techreps={}
for i in dna_data.keys():
    dna_ab[i]={}
    corr_techreps[i]={}
    for j in dna_data[i].keys():
        dna_ab[i][j]=pd.concat([dna_data[i][j]['a'],dna_data[i][j]['b']],axis=1,sort=False).fillna(value=0)
        dna_ab[i][j].columns=['count_a','count_b']
        dna_ab[i][j]['total_ab']=dna_ab[i][j].sum(axis=1)
        corr_techreps[i][j]=dna_ab[i][j].corr()
#        sns.scatterplot(dna_ab[i][j]['count_a'],dna_ab[i][j]['count_b'],alpha=0.5)
#        plt.show()
#        plt.close()

'''--------------------------------------------------------------merge bio reps---------------------------------------------------------------'''

dna_ss ={}
corr_bioreps={}
for j in pans:
    dna_ss[j]=pd.concat([dna_ab['set_1'][j]['total_ab'],dna_ab['set_2'][j]['total_ab'],dna_ab['set_3'][j]['total_ab']],axis=1,sort=False).fillna(value=0)
    dna_ss[j].columns=['total_ab_1','total_ab_2','total_ab_3']
    corr_bioreps[j]=dna_ss[j].corr()
    sns.scatterplot(dna_ss[j]['total_ab_1'],dna_ss[j]['total_ab_2'],alpha=0.5)
    plt.show()
    plt.close()
    sns.scatterplot(dna_ss[j]['total_ab_2'],dna_ss[j]['total_ab_3'],alpha=0.5)
    plt.show()
    plt.close()
    sns.scatterplot(dna_ss[j]['total_ab_1'],dna_ss[j]['total_ab_3'],alpha=0.5)
    plt.show()
    plt.close()