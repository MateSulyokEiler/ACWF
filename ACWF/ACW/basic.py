#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 11:10:59 2022

@author: raxis
"""

#to import degree, mutate, compare, difference, code, textin, textin2, predict, floatify, floatify_lowcut, mass, letters, letterz, lettersorg, apr, nt, gk, bases, masstable, predictor, median, lettercodes, ch_n, compare_shift, ranger

import math
from collections import defaultdict

def median(lst):
    n = len(lst)
    s = sorted(lst)
    return (s[n//2-1]/2.0+s[n//2]/2.0, s[n//2])[n % 2] if n else None


def degree(rad):

    if rad == None:
        deg = rad
    else:
        deg = math.degrees(rad)

    return (deg)


def mutate(peptide, newaa, location):
    peptide = peptide[0:location-1]+newaa+peptide[location:]
    return peptide

def compare(peptide1, peptide2):
    i=0
    for n in range(len(peptide1)):
        if peptide1[n]==peptide2[n]:
            i+=1
    return i

def compare_shift(peptide1, peptide2, shift):
    i=0
    for n in range(len(peptide1)-shift):
        if peptide1[n]==peptide2[n+shift]:
            i+=1
    if i == len(peptide1)-shift:
        return 1
    i=0
    for n in range(len(peptide2)-shift):
        if peptide2[n]==peptide1[n+shift]:
            i+=1
    if i == len(peptide1)-shift:
        return -1
    return 0

def difference(peptide1, peptide2):
    for n in range(len(peptide1)):
        if peptide1[n]!=peptide2[n]:
            break
    return peptide1[n], peptide2[n], n

def code(peptide1, peptide2):
    code=''
    for n in range(len(peptide1)):
        if peptide1[n] == peptide2[n]:
            code=code+"X"
        else:
            code=code+"0"
    return code


def textin(textin, n):
    with open(textin, 'r+') as f:
        textout = []
        textinlist = f.read().splitlines()
        for i in textinlist:
            textout.append(i[n:n+6])
    return textout

def textin2(textin):
    #meglévő útvonalhoz beolvasó
    with open(textin, 'r+') as f:
        textout=[]
        textinlist = f.read().splitlines()
        for i in textinlist:
            textoutline=[]
            for n in range(int((len(i)+2)/17)):
                textoutline.append(i[n*17:n*17+6])
            textout.append(textoutline)
    return textout

def predict(peptide):

    i = 1.08

    for n,aa in enumerate(peptide):
        i += predictor[aa][n]
    i = round(i,2)
    return i


def functionify(lista, function):
    listb = []
    for n in lista:
        if isinstance(n, list):
            listb.append(functionify(n))
        else:
            try:
                listb.append(function(n))
            except (ValueError, TypeError):
                listb.append(n)
    return listb


def floatify(lista):
    listb = []
    for n in lista:
        if isinstance(n, list):
            listb.append(floatify(n))
        else:
            try:
                listb.append(float(n))
            except (ValueError, TypeError):
                listb.append(n)
    return listb


def floatify_lowcut(lista, cut):
    listb = []
    for n in lista:
        if isinstance(n, list):
            listb.append(floatify_lowcut(n, cut))
        else:
            try:
                if float(n) < cut: 
                    listb.append(0)
                else:
                    listb.append(float(n))
            except (ValueError, TypeError):
                listb.append(n)
    return listb


def mass(peptide):
    mass = 17.999
    for aa in peptide:
        mass += masstable[aa]
    return mass 


def ch_n(count):
    # chain letters from numbers for convenience
    
    ch_let=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
            'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    
    return ch_let[count]


def ranger(start, end, step):
    output = [start]
    while start < end:
        start += step
        output.append(start)
    return output
    
    


letters = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X']
letterz = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
lettersorg = ['F', 'L', 'M', 'I','V', 'C', 'Y', 'W', 'A', 'G', 'T', 'S', 'N', 'H', 'Q', 'P', 'D', 'E', 'K', 'R', 'X']

apr = ['F', 'L', 'M', 'I','V', 'C', 'Y', 'W']
nt = ['A', 'G', 'T', 'S', 'N', 'H', 'Q']
gk = ['P', 'D', 'E', 'K', 'R']

bases = ['A','T','G','C']


lettercodes = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 
     'TYC': 'X', 'PTR': 'X'}

lettercodes = defaultdict(lambda: 'X', lettercodes)


masstable = {
    #average mass
    'A':  71.0788, #Ala
    'C': 103.1388, #Cys
    'D': 115.0886, #Asp
    'E': 129.1155, #Glu
    'F': 147.1766, #Phe
    'G':  57.0519, #Gly
    'H': 137.1411, #His
    'I': 113.1594, #Ile
    'K': 128.1741, #Lys
    'L': 113.1594, #Leu
    'M': 131.1926, #Met
    'N': 114.1038, #Asn
    'P':  97.1167, #Pro
    'Q': 128.1307, #Gln
    'R': 156.1875, #Arg
    'S':  87.0782, #Ser
    'T': 101.1051, #Thr
    'V':  99.1326, #Val
    'W': 186.2132, #Trp
    'Y': 163.1760, #Tyr
    }

predictor = {
    'A': {0:-0.26, 1:-0.32, 2:-0.27, 3:-0.14, 4:-0.43, 5:-0.22},
    'C': {0:-0.09, 1:-0.21, 2: 0.03, 3:-0.05, 4:-0.17, 5:-0.05},
    'D': {0:-0.49, 1:-0.43, 2:-0.56, 3:-0.41, 4:-0.56, 5:-0.36},
    'E': {0:-0.51, 1:-0.41, 2:-0.43, 3:-0.30, 4:-0.61, 5:-0.39},
    'F': {0:-0.13, 1:-0.11, 2: 0.05, 3:-0.05, 4:-0.13, 5:-0.11},
    'G': {0:-0.23, 1:-0.37, 2:-0.46, 3:-0.37, 4:-0.30, 5:-0.33},
    'H': {0:-0.32, 1:-0.26, 2:-0.26, 3:-0.30, 4:-0.35, 5:-0.25},
    'I': {0:-0.06, 1:-0.08, 2: 0.26, 3: 0.09, 4:-0.06, 5:-0.07},
    'K': {0:-0.39, 1:-0.45, 2:-0.51, 3:-0.35, 4:-0.59, 5:-0.32},
    'L': {0:-0.10, 1:-0.18, 2: 0.02, 3: 0.04, 4:-0.22, 5:-0.13},
    'M': {0:-0.17, 1:-0.25, 2:-0.02, 3:-0.10, 4:-0.19, 5:-0.18},
    'N': {0:-0.40, 1:-0.34, 2:-0.49, 3:-0.27, 4:-0.46, 5:-0.30},
    'P': {0:-0.56, 1:-0.38, 2:-0.56, 3:-0.51, 4:-0.42, 5:-0.45},
    'Q': {0:-0.37, 1:-0.30, 2:-0.36, 3:-0.34, 4:-0.48, 5:-0.32},
    'R': {0:-0.45, 1:-0.41, 2:-0.46, 3:-0.33, 4:-0.52, 5:-0.35},
    'S': {0:-0.37, 1:-0.35, 2:-0.41, 3:-0.30, 4:-0.48, 5:-0.23},
    'T': {0:-0.34, 1:-0.33, 2:-0.28, 3:-0.23, 4:-0.40, 5:-0.23},
    'V': {0:-0.05, 1:-0.14, 2: 0.19, 3: 0.14, 4:-0.19, 5: 0.01},
    'W': {0:-0.17, 1:-0.17, 2:-0.09, 3:-0.06, 4:-0.12, 5:-0.16},
    'Y': {0:-0.23, 1:-0.11, 2:-0.13, 3:-0.06, 4:-0.18, 5:-0.15}}


GRAVY_data = {
    'A': 1.80,
    'C': 2.50,
    'D': -3.50,
    'E': -3.50,
    'F': 2.80,
    'G': -0.40,
    'H': -3.20,
    'I': 4.50,
    'K': -3.90,
    'L': 3.80,
    'M': 1.90,
    'N': -3.50,
    'P': -1.60,
    'Q': -3.50,
    'R': -4.50,
    'S': -0.80,
    'T': -0.70,
    'V': 4.20,
    'W': -0.90,
    'Y': -1.30,
    'X': 0
}

def gravy_score(seq):
    
    gravy_sum = sum([GRAVY_data[aa] for aa in seq])
    
    return gravy_sum/len(seq)


"(2N4Q6N)+(2N4Q6N)"
"(1M3G,4G6V)+(1M4G5V,2V3G6V)"


def edit_contribution(cont):
    
    cont = cont.split('+')
    A = cont[0][1:-1].replace(',','')
    B = cont[1][1:-1].replace(',','')
    A = ''.join(i for i in A if not i.isdigit())
    B = ''.join(i for i in B if not i.isdigit())
    seq = A + B
    
    return seq, A, B

