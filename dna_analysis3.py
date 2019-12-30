# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 14:43:03 2019

@author: jtfl2
"""

import numpy as np
import os
import matplotlib.pyplot as plt
#import seaborn as sns
import math
import random
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
print('saving figures')
for set_name in list(full_sets.keys()):
    print(set_name)
    for pan in pans:
        print(pan)
        sequences = []
        for i in full_sets[set_name]:
            for j in range(int(full_sets[set_name][i][pan])):
                sequences.append(i)
        print('choosing at random')
        value = 2
        step = 10
        unique = []
        values = []
        for exp in range(7):
            temp = {}
            new_value = value * (step)**exp
            for _ in range(new_value):
                temp[random.choice(sequences)] = 1
            entropy = 0
            for seq in temp:
                prob = full_sets[set_name][seq][pan]/full_sets[set_name][seq]['Total']
                entropy += -prob*math.log(prob)
            unique.append(entropy)
            values.append(new_value)
            
        y_pos = np.arange(len(values))
        plt.bar(y_pos, unique, align='center', alpha=0.5)
        plt.xticks(y_pos, values)
        plt.title(set_name + ', ' + pan)
        plt.xlabel('Selection Size')
        plt.ylabel('Amount Unique')
        plt.savefig(d + '/6/Entropy of Unique Sequence Count from Random Selection of ' + set_name + ', ' + pan + '.png')
        plt.close()