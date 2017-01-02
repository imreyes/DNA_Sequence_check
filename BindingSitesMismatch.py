# -*- coding: utf-8 -*-
"""
Spyder Editor

Title: Function Gallery for Finding DNA Replicating Origin
Author: Guang Yang
Date: Dec. 23, 2016.
"""

# Print out lists.
def prtout(text):
    for i in text:
        print(i, end = ' ')

# DNA pattern to number
def ptn(text):
    def gtn(char):
        if char=='A' or char=='a':
            return 0
        elif char=='C' or char=='c':
            return 1
        elif char=='G' or char=='g':
            return 2
        elif char=='T' or char=='t':
            return 3
    num=0
    for char in text:
        num=num*4+gtn(char)
    return num


# DNA code number to pattern.
def ntp(num,k):
    def ntg(num):
        if num==0:
            return 'A'
        elif num==1:
            return 'C'
        elif num==2:
            return 'G'
        elif num==3:
            return 'T'
    seq=''
    count=0
    while count<k:
        seq+=ntg(num%4)
        num=num//4
        count+=1
    if num!=0:
        print('Warning: residual number detected; code not fully transfered!')
    return seq[::-1]


# Get compliment nucleotide.
def compnuc(char):
    if char=='A' or char=='a':
        return 'T'
    elif char=='T' or char=='t':
        return 'A'
    elif char=='C' or char=='c':
        return 'G'
    elif char=='G' or char=='g':
        return 'C'
    else:
        print('Error: wrong neucliotide letter!')
        return ''

    
# Return compliment strand of a given DNA sequence.
def compstr(text):
    comp=''
    for i in text:
        comp=compnuc(i)+comp
    return comp


# Generate a skew-list from DNA sequence.
# The list gets smaller when encountering more C, and larger with more T.
def skew(genome):
    counter = 0
    skewlist = list()
    skewlist.append(counter)
    for i in genome:
        if i == 'C' or i == 'c':
            counter -= 1
        elif i == 'G' or i == 'g':
            counter += 1
        skewlist.append(counter)
    return skewlist

# Find the minimal position of the skew list of a DNA sequence.
# The turning point shows possible origin location of replication.
def minskew(genome):
    counter = 0
    mins = 0
    idx = list()
    for i, j in enumerate(genome):
        if j == 'C' or j == 'c':
            counter -= 1
        elif j == 'G' or j == 'g':
            counter += 1
        if counter < mins:
            mins = counter
            idx = list()
        if counter == mins:
            idx.append(i+1)
    return idx


# Finding Hamming distance between the two strings of equal lengths.
# Specifically, finding distance between two genomes of same lengths.
def hamd(genome1, genome2, d = -1):
    if len(genome1) != len(genome2):
        print("Error: strings differ in length!")
        exit()
    dist = 0
    length = len(genome1)
    for i in range(length):
        if genome1[i] != genome2[i]:
            dist += 1
            if d >=0 and dist > d:
                return dist
    return dist


# Approximated matching of DNA sequence.
def appmatch(pattern, genome, d):
    k = len(pattern)
    n = len(genome)
    i = 0
    idx = list()
    while i <= n-k:
        if hamd(pattern, genome[i:(i+k)]) <= d:
            idx.append(i)
        i += 1
    return idx
    

# Find neighbors close in Hamming Distance d, of a given pattern.
def neighbor(pattern, d):
    nbhood = list()
    cache = list()
    seed = ['A', 'T', 'C', 'G']
    for elem in seed:
        nbhood.append(elem)
    idx = 1
    k = len(pattern)
    while idx < k:
        sub = pattern[0:idx]
        for elem in nbhood:
            if hamd(sub, elem, d) < d:
                for nuc in seed:
                    cache.append(elem + nuc)
            else:
                cache.append(elem + pattern[idx])
        nbhood = cache
        cache = list()
        idx += 1
    return nbhood
    
def countfreqapp(genome, pattern, d):
    k = len(pattern)
    i = 0
    n = len(genome)
    count = 0
    while i <= n-k:
        sub = genome[i:(i+k)]
        if hamd(sub, pattern, d) <= d:
            count += 1
        i += 1
    return count

# Find the most frequent approximated matches of a DNA sequence
# with given pattern and hamming distance.
def freqpatapp(genome, k, d):
    flist = [0]*(4**k)
    i = 0
    n = len(genome)
    while i <= n-k:
        pattern = genome[i:(i+k)]
        nbhood = neighbor(pattern, d)
        for elem in nbhood:
            flist[ptn(elem)] += 1
        i += 1
    m = max(flist)
    maxes = [i for i, j in enumerate(flist) if j == m]
    fplist = list()
    for idx in maxes:
        fplist.append(ntp(idx, k))
    fplist.sort()
    return fplist


def freqpatapp_comp(genome, k, d):
    flist = [0]*(4**k)
    i = 0
    n = len(genome)
    while i <= n-k:
        pattern = genome[i:(i+k)]
        nbhood = neighbor(pattern, d)
        nbhood.extend(neighbor(compstr(pattern), d))
        for elem in nbhood:
            flist[ptn(elem)] += 1
        i += 1
    m = max(flist)
    maxes = [i for i, j in enumerate(flist) if j == m]
    fplist = list()
    for idx in maxes:
        fplist.append(ntp(idx, k))
    fplist.sort()
    return fplist
