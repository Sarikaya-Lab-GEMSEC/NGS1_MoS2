""" Translation of DNA sequences into proteins from NGS"""

import numpy as np
import pandas as pd
import math
import itertools as it
import os
from time import time
#from time import time



def logb2(a):
    b = math.log(a,2)
    return(b)



def translate(seq):
    """
    Translate DNA to protein
    """
    table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'Q',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',}

    ncontain = [''.join(x) for x in it.product('ACTGN',repeat=3)]
    correct = [''.join(x) for x in it.product('ACTG',repeat=3)]
    rejects=list(set(ncontain)-set(correct))
    peptide = ""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i : i+3]
            if codon in rejects:
                break
            else:
                peptide += table[codon]
    return peptide

vlog = np.vectorize(logb2)
vtranslate = np.vectorize(translate)
dir1 = os.getcwd()
for file in os.listdir():
    now = time()
    print(file)
    DNAfile = pd.read_csv(dir1 + '\\' + file + '\\tsv\\' + file + '_lib' + '\\main_barcodes_counts.tsv',delimiter='\t')
    columns=['DNA_seq','Counts']
    DNAfile.columns=columns
    
    pepcount = vlog(DNAfile.Counts.values)
    print(time() - now)
    pepseq = vtranslate(DNAfile.DNA_seq.values)
    print(time() -now)
    PEPseq = pd.DataFrame(pepseq,index=DNAfile.index,columns=['AA_seq'])
    Copynum = pd.DataFrame(pepcount,index=DNAfile.index,columns=['counts'])
    PEPfile = pd.concat([PEPseq,Copynum],axis=1)
    
    pepfile = PEPfile[PEPfile['AA_seq'].apply(lambda x: len(str(x))==12)&
                       PEPfile['AA_seq'].apply(lambda x: "_" not in str(x))&
                       PEPfile['counts'].apply(lambda x: x>=1.0)]
    pepfile.to_csv(dir1 + '\\' + file + '\\'+file+'.csv')
        