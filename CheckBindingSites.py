# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
# Count replicate numbers of a given DNA pattern.
def pcount(text,pattern):
    count=0
    tlength=len(text)
    plength=len(pattern)
    pos=0
    while pos <= tlength-plength:
        if text[pos:(pos+plength)]==pattern:
            count=count+1
        pos=pos+1
    return count

# Find the most frequent k--mers in a DNA sequence.
def freqpat(text, k):
    maxcount=0
    fplist=list()
    pos=0
    tlength=len(text)
    while pos<= tlength-k:
        pattern=text[pos:(pos+k)]
        count=pcount(text,pattern)
        if count>=maxcount and not pattern in fplist:
            if count>maxcount:
                fplist=list()
                maxcount=count
            fplist.append(pattern)
        pos+=1
    fplist.sort()
    return fplist


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


# Use list to store times of appereance of patterns.
# This returns the entire list.
def freqlist(text,k):
    flist=[0]*(4**k)
    pos=0
    tlength=len(text)
    while pos<=tlength-k:
        pattern=text[pos:(pos+k)]
        idx=ptn(pattern)
        flist[idx]+=1
        pos+=1
    return flist


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
        exit


        
# Return compliment strand of a given DNA sequence.
def compstr(text):
    comp=''
    for i in text:
        comp=compnuc(i)+comp
    return comp



# Find positions of a string where designated pattern is found.
def posstr(text, pattern):
    tlength=len(text)
    plength=len(pattern)
    pos=0
    plist=list()
    while pos <= tlength-plength:
        if text[pos:(pos+plength)]==pattern:
            plist.append(pos)
        pos=pos+1
    return plist



# Find clump - in a section of genome (text) of length L, search for all k-mer patterns
# that appears at least t times in the section.
def findclump(text, k, L, t):
    clump = [0] * (4**k)
    pos = 0
    tlength = len(text)
    plist = list()
    subtext = text[pos:(pos+L)]
    flist = freqlist(subtext, k)
    for i in range(len(flist)):
        if flist[i] >= t and clump[i] == 0:
            clump[i] = 1
    pos += 1
    while pos <= tlength - L:
        pattern1 = text[(pos-1):(pos+k-1)]
        idx1 = ptn(pattern1)
        pattern2 = text[(pos+L-k):(pos+L)]
        idx2 = ptn(pattern2)
        flist[idx1] -= 1
        flist[idx2] += 1
        if flist[idx2] >= t and clump[idx2] == 0:
            clump[idx2] = 1
        pos += 1
    for i in range(len(clump)):
        if clump[i] == 1:
            pattern = ntp(i, k)
            plist.append(pattern)
    plist.sort()
    return plist


# Print answers.
def prtout(text):
    for i in text:
        print(i,end=' ')
