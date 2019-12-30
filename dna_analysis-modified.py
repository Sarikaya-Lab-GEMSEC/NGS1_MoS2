''' analysis of DNA reads from NextSeq NGS '''
''' @author: Siddharth Rath'''

import numpy as np
import pandas as pd
import seaborn as sns
#import math
import matplotlib.pyplot as plt
#import scipy as scp
''' read sets '''
sets={}
for i in ['set_1','set_2','set_3']:
    sets[i]={}
    for j in ['pan1','pan2','pan3','Eluate']:
        sets[i][j]={}
        for k in ['a','b']:
            sets[i][j][k]=pd.read_csv(str(i)+'_'+str(j)+'_'+str(k)+'.tsv.tsv',delimiter='\t')
            sets[i][j][k].columns=['DNA_seq','counts']
            sets[i][j][k].set_index('DNA_seq',inplace=True)

''' consolidate technical replicates'''
techreps={}
for i in sets.keys():
    techreps[i]={}
    for j in sets[i].keys():
        techreps[i][j]=pd.concat([sets[i][j]['a'],sets[i][j]['b']],axis=1,sort=False).fillna(value=0)
        techreps[i][j].columns=['count_a','count_b']
        techreps[i][j][str(j)]=techreps[i][j]['count_a']+techreps[i][j]['count_b']
        techreps[i][j].drop(columns=['count_a','count_b'],axis=1,inplace=True)

''' consolidate biological replicates'''
bioreps={}
total_sets={}
for i in techreps.keys():
    bioreps[i]=pd.concat([techreps[i]['pan1'],techreps[i]['pan2'],techreps[i]['pan3'],techreps[i]['Eluate']],axis=1,sort=False).fillna(value=0)
    bioreps[i]['Total']=bioreps[i].sum(axis=1)
    #bioreps[i].to_csv(str(i)+'_all_DNA_counts.csv')
    total_sets[i]=bioreps[i].Total

log_total={}
for i in total_sets.keys():
    log_total[i]=np.log2(total_sets[i])
    
sns.kdeplot(log_total['set_1'],bw=0.09,shade=True)
plt.title('log_count_distribution_set_1')
#plt.xlim((0,1e4))
plt.xlabel('log2(Counts)')
plt.ylabel('kernel_density')
plt.show()
plt.close()
sns.kdeplot(log_total['set_2'],bw=0.09,shade=True)
plt.title('log_count_distribution_set_2')
#plt.xlim((0,1e4))
plt.xlabel('log2(Counts)')
plt.ylabel('kernel_density')
plt.show()
plt.close()
sns.kdeplot(log_total['set_3'],bw=0.09,shade=True)
plt.title('log_count_distribution_set_3')
#plt.xlim((0,1e4))
plt.xlabel('log2(Counts)')
plt.ylabel('kernel_density')