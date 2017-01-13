# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 13:50:31 2017

Function Gallery for Real-World DNA Sequence Problems.
@author: Guang Yang
"""
import numpy as np
import math
import random

# Print out lists.
def prtout(text):
    for i in text:
        print(i)


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
   

# Find k-mers that appear in all strings of DNA colletion.
def checkpat(DNA, k, d):
    def freqlist(genome, k, d):
        flist = [0]*(4**k)
        i = 0
        n = len(genome)
        while i <= n-k:
            pattern = genome[i:(i+k)]
            nbhood = neighbor(pattern, d)
            for elem in nbhood:
                flist[ptn(elem)] += 1
            i += 1
        return flist
    fmat = np.zeros([len(DNA), 4**k])
    for i, genome in enumerate(DNA):
        fmat[i] = freqlist(genome, k, d)
    checklist = np.ones(4**k)
    for lines in fmat:
        checklist *= lines
    patterns = list()
    for i, j in enumerate(checklist):
        if j:
            patterns.append(ntp(i, k))
    patterns.sort()
    return patterns


# Calculate entropy.
def entropy(DNA):
    denom = len(DNA)
    entropy = 0
    for i in range(len(DNA[0])):
        chars = {'A':0, 'T':0, 'C':0, 'G':0}
        for genome in DNA:
            chars[genome[i]] += 1/denom
        for val in chars.values():
            val = round(val, 1)
            if val==0:
                continue
            entropy -= round(val * math.log(val, 2),5)      
    return entropy


# Find distance between designated sequence pattern
# and strings of DNA.
# By: summing minima of all distances between pattern and Dna strings.
def dpatdna(pattern, DNA):
    k = len(pattern)
    distance = 0
    for Dna in DNA:
        n = len(Dna)
        dist = k
        i = 0
        while i <= n - k:
            d = hamd(pattern, Dna[i:(i+k)])
            if dist > d:
                dist = d
            if dist == 0:
                break
            i += 1
        distance += dist
    return distance


# Find median string that minimizes distance score.
# The string has the smallest summed distances with all strings of DNA.
# Input: DNA as list of strings containing genome (ATCG), k as designated k-mer to search,
# greedy to find all (True) or only one (False).
def medstr(DNA, k, greedy = False):
    distance = math.inf
    if greedy:
        flist = [0] * (4**k)
    if not greedy:
        median = ''
    for i in range(4**k):
        pattern = ntp(i, k)
        dist = dpatdna(pattern, DNA)
        if greedy:
            flist[i] = dist
        if (not greedy) and distance > dist:
            distance = dist
            median = pattern
    if not greedy:
        return median
    else:
        m = min(flist)
        return [ntp(i, k) for i, j in enumerate(flist) if j == m]


# Find the most probable k-mer pattern in a DNA string,
# given the profile probability of the k-mer.
# Input: DNA as single string, integer k as designated k-mer to search,
# profile as array, indicating probability for genome(ACGT, row-wise)
# at each location of the k-mer (column).
def profilemostk(DNA, k, profile):
    n = len(DNA)
    maxprob = 0
    maxidx = 0
    i = 0
    def getidx(char):
        if char == 'A':
            return 0
        elif char == 'C':
            return 1
        elif char == 'G':
            return 2
        elif char == 'T':
            return 3
    while i <= n - k:
        pattern = DNA[i:(i+k)]
        prob = 1
        for j in range(k):
            idx = getidx(pattern[j])
            p = float(profile[idx][j])
            if p == 0:
                prob = 0
                break
            prob *= p
        if maxprob < prob:
            maxprob = prob
            maxidx = i
        i += 1
    return DNA[maxidx:(maxidx+k)]


# Return profile matrix of motif collection.
# Columns represent ith genome position in the k-mer;
# rows represent probability of 'A', 'C', 'G', 'T'.
def profile(motifs):
    k = len(motifs[0])
    t = len(motifs)
    mat = [[], [], [], []]
    gdict = {'A':0, 'C':1, 'G':2, 'T':3}
    for i in range(k):
        chars = {'A':0, 'T':0, 'C':0, 'G':0}
        for j in range(t):
            char = motifs[j][i]
            chars[char] += round(1/t, 3)
        for key, val in chars.items():
            mat[gdict[key]].append(val)
    return mat


# Find concensus string (the most probable k-mer in the DNA string, judged by profile)
def consensus(motifs):
    k = len(motifs[0])
    t = len(motifs)
    cons = list()
    for i in range(k):
        chars = {'A':0, 'T':0, 'C':0, 'G':0}
        for j in range(t):
            char = motifs[j][i]
            chars[char] += 1
        v = list(chars.values())
        k = list(chars.keys())
        cons.append(k[v.index(max(v))])
    cons = ''.join(cons)
    return cons


# Return score of motif collection.
def score(motifs):
    k = len(motifs[0])
    t = len(motifs)
    scores = 0
    for i in range(k):
        chars = {'A':0, 'T':0, 'C':0, 'G':0}
        for j in range(t):
            char = motifs[j][i]
            chars[char] += 1
        scores += (t - max(chars.values()))
    return scores


# Upgrade the profile() function, applying Laplace's Rule.
# When counting genome appearances in motifs, add the entire matrix by one.
def laprofile(motifs):
    k = len(motifs[0])
    t = len(motifs)
    mat = [[], [], [], []]
    gdict = {'A':0, 'C':1, 'G':2, 'T':3}
    for i in range(k):
        chars = {'A':0, 'T':0, 'C':0, 'G':0}
        for j in range(t):
            char = motifs[j][i]
            chars[char] += round(1/t, 3)
        for key, val in chars.items():
            mat[gdict[key]].append(val + 1/t)
    return mat
        

# Search for best motif collection for the DNA collections.
# Not an accurate solution.
# Input: list of DNA strings, integer k for size of k-mer.
def greedymotif(DNA, k, laplace = True):
    t = len(DNA)
    bestmotifs = []
    for i in DNA:
        bestmotifs.append(i[:k])
    string1 = DNA[0]
    for i in range(len(string1) - k + 1):
        motifs = [string1[i:(i+k)]]
        for j in range(1, t):
            if laplace:
                pro = laprofile(motifs)
            else:
                pro = profile(motifs)
            motifs.append(profilemostk(DNA[j], k, pro))
        if score(motifs) < score(bestmotifs):
            bestmotifs = motifs
    return bestmotifs


# Search motif gallery using randomized algorithm.
def randmotifs(DNA, k):
    n = len(DNA[0])
    motifs = []                     # Initialize motifs collection
    for string in DNA:
        i = random.randint(0, n-k)
        motifs.append(string[i:(i+k)])
    bestmotifs = motifs.copy()
    minscore = score(bestmotifs)
    #minscore = entropy(bestmotifs)
    counter = 0
    while True:
        pro = laprofile(motifs)
        motifs = []
        for string in DNA:
            motif = profilemostk(string, k, pro)
            motifs.append(motif)
        newscore = score(motifs)
        #newscore = entropy(motifs)
        if minscore > newscore:     # For one optimization, try to get min-scored motifs.
            bestmotifs = motifs.copy()
            minscore = newscore
        else:
            return bestmotifs


# Refining functions for randmotifs - get better chance to reach 'global minimum'.
def randsearch(DNA, k, rep = 1000):
    i = 0
    t = len(DNA)
    bestmotifs = randmotifs(DNA, k)
    minscore = score(bestmotifs)
    while i < rep:
        motifs = randmotifs(DNA, k)
        newscore = score(motifs)
        if minscore > newscore:
            minscore = newscore
            bestmotifs = motifs.copy()
            print(minscore, i)
        i += 1
    return bestmotifs


# Motif generator based on profile probabilities.
# Not useful here; just to keep a reference.
def RandProGen(profile):
    k = len(profile[0])
    t = len(profile)
    string = []
    probs = []
    glist = ['A', 'C', 'G', 'T']
    for j in range(k):
        for i in range(t):
            probs.append(profile[i][j])
        string.append(np.random.choice(glist, p = probs))
    return ''.join(string)


# Gibbs sampler which optimizes only 1 motif, instead of entire motif colelction.
def GibbsSampler(DNA, k, N):
    n = len(DNA[0])
    t = len(DNA)
    motifs = []
    for string in DNA:
        i = random.randint(0, n-k)
        motifs.append(string[i:(i+k)])
    bestmotifs = motifs.copy()
    minscore = score(bestmotifs)
    for rep in range(N):
        i = random.randint(0, t-1)
        cache = motifs.copy()
        cache.pop(i)
        pro = profile(cache)
        rand = random.randint(0, n-k)
        motif = DNA[i][rand:(rand+k)]
        motifs[i] = motif
        newscore = score(motifs)
        if minscore > newscore:
            bestmotifs = motifs.copy()
            minscore = newscore
            print(minscore, i)
    return bestmotifs


# Refining functions for GibbsSampler - get better chance to reach 'global minimum'.
def randGibbs(DNA, k, N, rep = 20):
    i = 0
    t = len(DNA)
    bestmotifs = GibbsSampler(DNA, k, N)
    minscore = score(bestmotifs)
    while i < rep:
        motifs = GibbsSampler(DNA, k, N)
        newscore = score(motifs)
        if minscore > newscore:
            minscore = newscore
            bestmotifs = motifs.copy()
            print(minscore, i)
        i += 1
    return bestmotifs
