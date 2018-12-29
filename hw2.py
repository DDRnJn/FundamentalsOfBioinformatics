#Bioinformatics HW 2
#Dhruv Ranjan

import sys
import os
import random
import copy

def motifEnumeration(Dna, k, d):

    patterns = []
    for seq in Dna:
        patterns += patternGen(seq, k)
    finalPatterns = []
    for pattern in patterns:
        neighborhood = neighbors(pattern, d)
        for neighbor in neighborhood:
            inSeq = False
            inSeqs = True
            for seq in Dna:
                dnaPatterns = patternGen(seq, k)
                for pattern2 in dnaPatterns:
                    if hammingDistance(neighbor,pattern2)<=d:
                        inSeq = True
                if inSeq != True:
                    inSeqs = False
                inSeq = False
            if inSeqs == True:
                finalPatterns += [neighbor]
            inSeqs = True
    finalPatterns = list(set(finalPatterns))
    #print finalPatterns
    return " ".join(finalPatterns)

def patternGen(text, k):
    
    patterns = [] 
    for i in xrange(0, len(text)-k+1):
        patterns += [text[i:(k+i)]]
    return patterns 

def hammingDistance(s1,s2):

    if len(s1)==len(s2)==0:
        return 0
    else:
        if s1[0]==s2[0]:
            return hammingDistance(s1[1:],s2[1:])
        else:
            return hammingDistance(s1[1:],s2[1:]) + 1

#returns all neighbors of pattern at most d hamming distance away 
def neighbors(pattern, d):

    if d == 0:
        return pattern
    if len(pattern) == 1:
        return ["A", "C", "G", "T"]
    neighborhood = []
    suffixNeighbors = neighbors(pattern[1:],d)
    for text in suffixNeighbors:
        if hammingDistance(pattern[1:], text) < d:
            neighborhood += ["A" + text]
            neighborhood += ["T" + text]
            neighborhood += ["C" + text]
            neighborhood += ["G" + text]
        else:
            neighborhood += [pattern[0] + text]
    return neighborhood

def medianString(k, Dna):

    distance = sys.maxint
    for i in xrange(0,(4**k)):
        pattern = numberToPattern(i, k)
        if distance > distanceBetweenPatternAndString(pattern, Dna):
            distance = distanceBetweenPatternAndString(pattern, Dna)
            median = pattern
    return median


def numberToPattern(index, k):
    if k==1:
        return numberToSymbol(index)
    else:
        prefixIndex = index/4
        remainder = index % 4
        symbol = numberToSymbol(remainder)
        return numberToPattern(prefixIndex,k-1) + symbol

def numberToSymbol(n):
    if n == 0:
        return "A"
    elif n==1:
        return "C"
    elif n==2:
        return "G"
    else:
        return "T"

def distanceBetweenPatternAndString(pattern, Dna):
    
    k = len(pattern)
    distance = 0
    for text in Dna:
        hd = sys.maxint
        patterns = patternGen(text, k)
        for kmer in patterns:
            if hd > hammingDistance(pattern,kmer):
                hd = hammingDistance(pattern,kmer)
        distance += hd
    return distance

def profileMostProbableKmer(text,k,matrix):

    patterns = patternGen(text,k)
    bestPattern = patterns[0]
    bestScore = 0
    for pattern in patterns:
        score = profileScore(pattern,matrix)
        if score > bestScore:
            bestScore = score
            bestPattern = pattern
    return bestPattern

def profileScore(pattern, matrix):

    score = 1
    for i in xrange(0,len(pattern)):
        if pattern[i]=="A":
            score *= matrix[0][i]
        elif pattern[i]=="C":
            score *= matrix[1][i]
        elif pattern[i]=="G":
            score *= matrix[2][i]
        elif pattern[i]=="T":
            score *= matrix[3][i]
    return score

def readFile(filename, mode="rb"): 
    
    fin = contents = None
    try:
        fin = open(filename, mode)
        contents = fin.read()
    finally:
        if (fin != None): fin.close()
    return contents

def firstKmerMatrix(Dna, k):

    motifs = []
    for seq in Dna:
        currentKmer = seq[0:k]
        motifs += [currentKmer]
    return motifs

def greedyMotifSearch(Dna, k, t):

    bestMotifs = firstKmerMatrix(Dna, k)
    firstStringMotifs = patternGen(Dna[0],k)
    for motif in firstStringMotifs:
        motifs = []
        motifs += [motif]
        for i in xrange(1,t):
            currentProfile = formProfile(motifs)
            currentSeq = Dna[i]
            mostProbable = profileMostProbableKmer(currentSeq, k, currentProfile)
            motifs += [mostProbable]
        if motifScore(motifs) < motifScore(bestMotifs):
            bestMotifs = motifs
    return bestMotifs

def greedyMotifSearchFileWrapper(fileName):

    contents=open(fileName).readlines()
    (k, t) = contents[0].split()
    contents.pop(0)
    Dna = contents
    bestMotifs = greedyMotifSearch(Dna, int(k), int(t))
    print " ".join(bestMotifs)

def motifScore(motifs):

    totalScore = 0
    t = len(motifs)
    k = len(motifs[0])
    for i in xrange(0,k):
        countA = 0
        countC = 0
        countG = 0
        countT = 0
        for j in xrange(0,t):
            currentBase = motifs[j][i]
            if currentBase=="A":
                countA += 1
            elif currentBase=="C":
                countC += 1
            elif currentBase=="G":
                countG += 1
            else:
                countT += 1
        total = countA+countC+countG+countT
        maxCount = max([countA,countC,countG,countT])
        totalScore += (total - maxCount)
    return totalScore

#working as intended
def formCountMatrix(motifs):

    t = len(motifs)
    k = len(motifs[0])
    countMatrix = [[],[],[],[]]
    for i in xrange(0,k):
        countA = 0
        countC = 0
        countG = 0
        countT = 0
        for j in xrange(0,t):
            currentBase = motifs[j][i]
            if currentBase=="A":
                countA += 1
            elif currentBase=="C":
                countC += 1
            elif currentBase=="G":
                countG += 1
            else:
                countT += 1
        countMatrix[0]+=[countA]
        countMatrix[1]+=[countC]
        countMatrix[2]+=[countG]
        countMatrix[3]+=[countT]
    return countMatrix

def formCountMatrixPseudo(motifs):

    t = len(motifs)
    k = len(motifs[0])
    countMatrix = [[],[],[],[]]
    for i in xrange(0,k):
        countA = 1
        countC = 1
        countG = 1
        countT = 1
        for j in xrange(0,t):
            currentBase = motifs[j][i]
            if currentBase=="A":
                countA += 1
            elif currentBase=="C":
                countC += 1
            elif currentBase=="G":
                countG += 1
            else:
                countT += 1
        countMatrix[0]+=[countA]
        countMatrix[1]+=[countC]
        countMatrix[2]+=[countG]
        countMatrix[3]+=[countT]
    return countMatrix


def formProfile(motifs):

    countMatrix = formCountMatrix(motifs)
    k = len(countMatrix[0])
    profileMatrix = [[],[],[],[]]
    for i in xrange(0,k):
        countA = countMatrix[0][i]
        countC = countMatrix[1][i]
        countG = countMatrix[2][i]
        countT = countMatrix[3][i]
        total = countA + countC + countG + countT
        percentA = float(countA)/total
        percentC = float(countC)/total
        percentG = float(countG)/total
        percentT = float(countT)/total
        profileMatrix[0]+=[percentA]
        profileMatrix[1]+=[percentC]
        profileMatrix[2]+=[percentG]
        profileMatrix[3]+=[percentT]
    return profileMatrix

def formProfilePseudo(motifs):

    countMatrix = formCountMatrixPseudo(motifs)
    k = len(countMatrix[0])
    profileMatrix = [[],[],[],[]]
    for i in xrange(0,k):
        countA = countMatrix[0][i]
        countC = countMatrix[1][i]
        countG = countMatrix[2][i]
        countT = countMatrix[3][i]
        total = countA + countC + countG + countT
        percentA = float(countA)/total
        percentC = float(countC)/total
        percentG = float(countG)/total
        percentT = float(countT)/total
        profileMatrix[0]+=[percentA]
        profileMatrix[1]+=[percentC]
        profileMatrix[2]+=[percentG]
        profileMatrix[3]+=[percentT]
    return profileMatrix

def greedyMotifSearchPseudo(Dna, k, t):

    bestMotifs = firstKmerMatrix(Dna, k)
    firstStringMotifs = patternGen(Dna[0],k)
    for motif in firstStringMotifs:
        motifs = []
        motifs += [motif]
        for i in xrange(1,t):
            currentProfile = formProfilePseudo(motifs)
            currentSeq = Dna[i]
            mostProbable = profileMostProbableKmer(currentSeq, k, currentProfile)
            motifs += [mostProbable]
        if motifScore(motifs) < motifScore(bestMotifs):
            bestMotifs = motifs
    return bestMotifs

def greedyMotifSearchFileWrapperPseudo(fileName):

    contents=open(fileName).readlines()
    (k, t) = contents[0].split()
    contents.pop(0)
    Dna = contents
    bestMotifs = greedyMotifSearchPseudo(Dna, int(k), int(t))
    print " ".join(bestMotifs)

def randomizedMotifSearch(Dna, k, t):
    
    motifs = []
    for seq in Dna:
        kmers = patternGen(seq, k)
        randomIndex = random.randint(0, len(kmers)-1)
        motifs += [kmers[randomIndex]]
    bestMotifs = motifs
    while True:
        profile = formProfilePseudo(motifs)
        motifs = []
        for i in xrange(0,t):
            currentSeq = Dna[i]
            mostProbable = profileMostProbableKmer(currentSeq, k, profile)
            motifs += [mostProbable]
        if motifScore(motifs) < motifScore(bestMotifs):
            bestMotifs = motifs
        else:
            return bestMotifs
        
#not working yet. Will fix later.
def gibbsSampler(Dna, k, t, N):
    motifs = []
    for seq in Dna:
        kmers = patternGen(seq, k)
        randomIndex = random.randint(0,len(kmers)-1)
        motifs += [kmers[randomIndex]]
    bestMotifs = motifs
    for j in xrange(0, N):
        i = random.randint(0,t)
        profile = formProfilePseudo(motifs.pop(i))
        motifi = profileRandomlyGeneratedKmer(profile,Dna[i],k)
        motifs += [motifi]
        if motifScore(motifs) < motifScore(bestMotifs):
            bestMotifs = motifs
    return bestMotifs

def profileRandomlyGeneratedKmer(profile, seq, k):

    kmers = patternGen(seq, k)
    probabilities = []
    probabilitiesCorrected = []
    for pattern in kmers:
        probabilities += [profileScore(pattern,profile)]
    totalProb = sum(probabilities)
    for prob in probabilities:
        probabilitiesCorrected += [float(prob)/totalProb]
    randomNum = random.randint(0,totalProb)
    total = 0
    for i in xrange(len(probabilitiesCorrected)):
        total += probabilitiesCorrected[i]
        if randomNum < total:
            return seq[i]
        
def randomizedMotifSearchWrapper(fileName):

    contents=open(fileName).readlines()
    (k, t) = contents[0].split()
    contents.pop(0)
    Dna = contents
    allMotifs = []
    bestMotifs = []
    for i in xrange(1000):
        bestMotif = randomizedMotifSearch(Dna, int(k), int(t))
        allMotifs += [bestMotif]
    bestMotifs = allMotifs[0]
    for motifs in allMotifs:
        if motifScore(motifs) < motifScore(bestMotifs):
            bestMotifs = motifs
    for m in bestMotifs:
        print m
    #print " ".join(bestMotifs)

def randomizedMotifSearchWrapper2(Dna, k, t):

    allMotifs = []
    for i in xrange(1000):
        bestMotif = randomizedMotifSearch(Dna, k, t)
        allMotifs += [bestMotif]
    bestMotifs = allMotifs[0]
    for motifs in allMotifs:
        if motifScore(motifs) < motifScore(bestMotifs):
            bestMotifs = motifs
    print " ".join(bestMotifs)


    



def profileMostProbableTest():

    text = "CCCAACTAGGCAAGTCTTCCAACCGCATCTATAGTGATCTAAACGACTCTTCACGTCCGTCTGGCGCATTTTCCTATCAACTAGAGTTGGATCACATGCAACCGGTCCACTCTGGGGAGTGGATGCAGGTCCGAGTTGCGCAAAATCATGACTGCCGGCCTTTGCTAACTAGCCCCATAGGTGGGAGGATGTAGATTCGTCTCTTGTGTTTGGGACGTCGACGTCATCGATATGGCTGATGCAAGGACTTAGAGGGAACATCGTGGAATTGGCGGCAGCTTTACGGGAGTGAGCTTAGACTGAACCGCCGCAGTATCGAACTTTCCACTGGTCACAGTGCGCCACGGGGATGCAGCTTTTATCTGACGCGACAATGTAAGGTATACGCGAATTCCATGAAAACCTCAGCTTGAGTAGTGGGTGGGCACACAGGACTCATCCGCTGGTCAGTGGCATCTTTCATGGTCAACAGACTCTTTATACAATCTTTTAAATCGTTGCGCAGGTACGGCCTCAATTATATTCAAGCGGCGTAGAACGCAAGTGTGACAGACCGTAGGTCACGCAAGGTCAGCGCAACGTTACGAACTAAGCATCATTTAATGGGTAGTCACATGGCTGAGATGACCCGTCACCAGTGGATGCCACTTAGCACACTTGCTCGGCCGAGGCAACGTCCACGGAGACAGTGAACTCCTGCTCGTATGCTTTCTAATACCCTCTATCTCCGCCCTTTGCACCTACTAGGCGTCACAATCTACAATACTAAATACTTAGAGTAGTGCTGATACGTTTTGAGCACTTCTTTATGGCTTTACCGATGATGATTAAACGGGCTCCTCTACCTATCGCTACCCAAAGCTCTCACATGAGCGTCAGCTGTGCGGACCCATTAGCTCGCGTACGCGCAGAGGAGTCACGGACGCACGCACAGACCATAACCCATGGGGGAAAGGTCTTGTCCGGGGGACCGCTAAATCGGGTCAGCCACGATAGACTT"
    k = 14
    matrix = [[0.268,0.268,0.155,0.282,0.225,0.211, 0.169, 0.31, 0.183, 0.282, 0.31, 0.254, 0.268, 0.268],
              [0.268,0.282,0.225,0.338,0.324,0.239,0.268,0.239,0.324,0.239,0.254,0.211,0.324,0.268],
              [0.225,0.183,0.352,0.225,0.239,0.268,0.254,0.197,0.254,0.225,0.183,0.324,0.211,0.269],
              [0.239,0.268,0.268,0.155,0.211,0.282,0.31,0.254,0.239,0.254,0.254,0.211,0.197,0.255]]
    return profileMostProbableKmer(text,k,matrix)

sys.setrecursionlimit(10000000)
