# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 12:24:18 2019

@author: jtfl2
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import math
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


#set_name = list(full_sets.keys())[0]
#pan = 'pan1'

pans = ['Pan 1', 'Pan 2', 'Pan 3', 'Eluate', 'Total']
print('saving figures')
for set_name in list(full_sets.keys()):

    print(set_name)
    for pan in pans:
        print(pan)
        counts = []
        for i in full_sets[set_name]:
            counts.append(full_sets[set_name][i][pan])
        
        value = int(max(counts)/500)
        if value < 50:
            value = 50
        print('full ' + str(value))
        plt.hist(counts, bins = value)
        #plt.hist(counts, bins = int(max(counts)))
        plt.title('Count Distribution of ' + set_name + ', ' + pan)
        #plt.xlim((0, 50))
        plt.xlabel('Counts')
        plt.ylabel('Density')
#        plt.show()
        plt.savefig(d + '/1/Count Distribution of ' + set_name + ', ' + pan + ' full.png')
        plt.close()
        
        print('zoom ' + str(max(counts)))
        plt.hist(counts, bins = int(max(counts)))
        plt.title('Count Distribution of ' + set_name + ', ' + pan)
        plt.xlim((0, 50))
        plt.xlabel('Counts')
        plt.ylabel('Density')
        plt.savefig(d + '/1/Count Distribution of ' + set_name + ', ' + pan + ' zoom.png')
        plt.close()
        
        print('log')
        counts = []
        for i in full_sets[set_name]:
            if full_sets[set_name][i][pan] != 0:
                counts.append(math.log2(full_sets[set_name][i][pan]))
            
        print('full')
        plt.hist(counts, bins = int(max(counts)))
        #plt.hist(counts, bins = int(max(counts)))
        plt.title('Log2 Count Distribution of ' + set_name + ', ' + pan)
        #plt.xlim((0, 50))
        plt.xlabel('Counts')
        plt.ylabel('Density')
#        plt.show()
        plt.savefig(d + '/2/Log2 Count Distribution of ' + set_name + ', ' + pan + ' full.png')
        plt.close()
        
