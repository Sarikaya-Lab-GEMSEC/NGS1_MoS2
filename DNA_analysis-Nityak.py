# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 14:41:55 2019

@author: knity
"""

import pandas as pd
import matplotlib.pyplot as plt

#class DNA_Analysis:
#    'read in tsv files'
#    'takes in a file location (does not include name of file)'
#    def loadTSV(location):
#        for i in ['set_1','set_2','set_3']:
#            for j in ['pan1','pan2','pan3','Eluate']:
#                for k in ['a','b']:
#                    name = str(i) + '_' + str(j) + '_' + str(k)
#                    data = pd.read_csv(location + str(i)+'_'+str(j)+'_'+str(k) + '.tsv.tsv',delimiter='\t')
#                    data.columns=['DNA_seq','counts']
#                    data.set_index('DNA_seq',inplace=True)
#                    sets[name] = data
# 
#    def reps():
#        for i in ['set_1','set_2','set_3']:
#            for j in ['pan1','pan2','pan3','Eluate']:
#                k = ['a', 'b']
#                name = str(i) + '_' + str(j)
#                name_a = name + '_' + str(k[0])
#                name_b = name + '_' + str(k[1])
#                reps = sets[name_a].merge(sets[name_b], left_on='DNA_seq', right_on='DNA_seq')
#                #reps.columns=['sequence', 'a', 'b']
#                tech_reps[name] = reps['counts_x'] + tech_reps['counts_y']
#                tech_reps.columns=['sequence', 'a', 'b', name]
#                #tech_reps.drop(columns=[name_a,name_b],axis=1,inplace=True)
#    
#    loadTSV('./DNA_NGS_1/')
#    reps()


# Load in data
sets = {}
for i in ['set_1','set_2','set_3']:
    sets[i]={}
    for j in ['pan1','pan2','pan3','Eluate']:
        sets[i][j]={}
        for k in ['a','b']:
            name = str(i) + '_' + str(j) + '_' + str(k)
            data = pd.read_csv('./DNA_NGS_1/' + str(i)+'_'+str(j)+'_'+str(k) + '.tsv.tsv',delimiter='\t')
            data.columns=['DNA_seq','counts']
            data.set_index('DNA_seq',inplace=True)
            sets[i][j][k] = data

# Counts after consolidating technical replicates
merged_sets = {}
for i in ['set_1','set_2','set_3']:
    merged_sets[i] = {}
    for j in ['pan1','pan2','pan3','Eluate']:
        merged_sets[i][j] = {}
        m_data = sets[i][j]['a'].merge(sets[i][j]['b'], left_on='DNA_seq', right_on='DNA_seq')
        merged_sets[i][j] = m_data

# Technical replicates count distribution
for i in ['set_1','set_2','set_3']:
    for j in ['pan1','pan2','pan3','Eluate']:
        plt.figure()
        merged_sets[i][j].head(500).plot()
        plt.xticks([])
        plt.xlabel('DNA Sequences')
        plt.ylabel('Count')
        plt.title(f'{i} {j} A and B count distributions')
        plt.savefig(f'./DNA_NGS_1/plots/TechRepsDistZoom/{i}_{j}ABCounts_Zoom.png')

# Similarity in technical replicates for all pans
for i in ['set_1','set_2','set_3']:
    fig, ax = plt.subplots()
    ax.scatter(sets[i]['pan1']['a'].head(500), sets[i]['pan1']['b'].head(500), label = 'Pan 1', alpha=0.6, marker = 'o')
    ax.scatter(sets[i]['pan2']['a'].head(500), sets[i]['pan2']['b'].head(500), label = 'Pan 2', alpha=0.6, marker = '^')
    ax.scatter(sets[i]['pan3']['a'].head(500), sets[i]['pan3']['b'].head(500), label = 'Pan 3', alpha=0.6, marker = '*')
    t = ax.scatter(sets[i]['Eluate']['a'].head(500), sets[i]['Eluate']['b'].head(500), label = 'Eluate', alpha=0.6, marker = 's')
    plt.xlabel('a count')
    plt.ylabel('b count')
    plt.title(i)
    plt.legend(bbox_to_anchor=(0.9, .05), loc='lower right')
    ax.grid(True)
    plt.savefig(f'./DNA_NGS_1/plots/AvsB_allPans/{i}_AvsBCounts.png')

 
# Counts after consolidating biological replicates  
for i in ['pan1','pan2','pan3', 'Eluate']:
    fig, ax = plt.subplots()
    ax.scatter(sets['set_1'][i]['a'].head(500), sets['set_1'][i]['b'].head(500), label = 'Set 1', alpha=0.6, marker = 'o')
    ax.scatter(sets['set_2'][i]['a'].head(500), sets['set_3'][i]['b'].head(500), label = 'Set 2', alpha=0.6, marker = '^')
    t = ax.scatter(sets['set_3'][i]['a'].head(500), sets['set_3'][i]['b'].head(500), label = 'Set 3', alpha=0.6, marker = 's')
    plt.xlabel('a count')
    plt.ylabel('b count')
    plt.xlim()
    plt.title(i)
    plt.legend(bbox_to_anchor=(0.9, .05), loc='lower right')
    ax.grid(True)
    plt.savefig(f'./DNA_NGS_1/plots/CompareSets/{i}_Counts_per_set.png')    
    

    
    
    
    
    
    
    
    
    
    
    
    
    