# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 15:15:50 2019

@author: jtfl2
"""

import numpy as np
import os
import matplotlib.pyplot as plt
#import seaborn as sns
import math
#import random
#
full_sets = {}
d = os.getcwd()
for i in os.listdir(d + '/sets/'):
    print(i)
    set_name = i.replace('all_DNA_counts.csv','').replace('_', ' ')
    full_sets[set_name] = {}
    with open(d + '/sets/' + i, 'r') as g:
        g.readline()
        string = g.readline()
        while string:
            l = string.split(',')
            full_sets[set_name][l[0]] = {}
            full_sets[set_name][l[0]]['Pan 1'] = float(l[1])
            full_sets[set_name][l[0]]['Pan 2'] = float(l[2])
            full_sets[set_name][l[0]]['Pan 3'] = float(l[3])
            full_sets[set_name][l[0]]['Eluate'] = float(l[4])
            full_sets[set_name][l[0]]['Total'] = float(l[5])
            string = g.readline()

pans = ['Pan 1', 'Pan 2', 'Pan 3', 'Eluate', 'Total']
DNA = ['A', 'T', 'C', 'G']
y_pos = np.arange(len(pans))

for set_name in list(full_sets.keys()):
    print(set_name)
    mat = {}
    set_entropy = []
    for pan in pans:
        mat[pan] = np.zeros((4,36))
        print(pan)
        sequences = []
        try:
            for i in full_sets[set_name]:
                for dna in range(len(i)):
                    mat[pan][DNA.index(i[dna])][dna] += full_sets[set_name][i][pan]
        except:
            continue
    entropy_values = []
    for pan in pans:
        mat[pan] = mat[pan]/mat['Total']
        plt.matshow(mat[pan])
        plt.savefig(d + '/5/Probability of Nucleotides Per Location for' + set_name + ', ' + pan + '.png')
        plt.close()
        entropy = 0
        for i in range(len(mat[pan])):
            for j in range(len(mat[pan][0])):
                entropy += -mat[pan][i][j]*math.log(mat[pan][i][j])
        entropy_values.append(entropy)
    if pan != 'Total':
        set_entropy.append(sum(entropy_values))
    plt.bar(y_pos, entropy_values, align='center', alpha=0.5)
    plt.xticks(y_pos, pans)
    plt.ylabel('Entropy')
    plt.title('Entropy per Pan for '+ set_name)
    plt.savefig(d + '/5/Entropy per pan for '+ set_name + '.png')
    plt.close()
    
    
setnames = list(full_sets.keys())
y_pos = np.arange(len(setnames))
plt.bar(y_pos, set_entropy, align='center', alpha=0.5)
plt.xticks(y_pos, setnames)
plt.ylabel('Entropy')
plt.title('Entropy per Set')
plt.savefig(d + '/5/Entropy per Set.png')
plt.close()

    
        
        
    
                