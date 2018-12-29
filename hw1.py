#Bioinformatics HW1
#Dhruv Ranjan

import copy
import os 
import sys

#count instances of str text in str pattern
def patternCount(text, pattern):
    count = 0
    k = len(pattern)
    for i in range(0, len(text) - len(pattern)+1):
        if text[i:(k+i)]==pattern:
            count += 1
    return count

#find most frequent words of length k in str text
def frequentWords(text, k):
    frequentPatterns = []
    counts = [] 
    for i in range(len(text)-k+1):
        pattern = text[i:(k+i)]
        counts += [patternCount(text, pattern)]
    maxCount = max(counts)
    for i in range(0, len(text)-k+1):
        if counts[i] == maxCount:
            frequentPatterns += [text[i:(k+i)]]
    frequentPatterns = list(set(frequentPatterns))
    return frequentPatterns

def fasterFrequentWords(text, k):
    frequentPatterns = []
    frequencyArray = computeFrequencies(text, k)
    maxCount = max(frequencyArray)
    for i in range(0,4**k):
        if frequencyArray[i] == maxCount:
            pattern = numberToPattern(i, k)
            frequentPatterns += [pattern]
    return list(set(frequentPatterns))
    

#find the reverse complement of the str pattern
def reverseComplement(pattern):

    newPattern = []
    for i in range(len(pattern)):
        if pattern[i] == "A":
            newPattern += "T"
        elif pattern[i] == "T":
            newPattern += "A"
        elif pattern[i] == "C":
            newPattern += "G"
        else:
            newPattern += "C"
    newPattern = newPattern[::-1]
    return "".join(newPattern)

#find the starting positions(0 indexed) of str pattern in str genome

def patternMatch(pattern, genome):

    positions = []
    for i in range(0,len(genome)-len(pattern)+1):
        if genome[i:len(pattern)+i] == pattern:
            positions += [i]
    return " ".join(str(e) for e in positions)

def patternMatchFromFile(pattern):

    positions = []
    genome = readFile("texts/vibrioCholeraeGenome.txt")
    for i in range(0,len(genome)-len(pattern)+1):
        if genome[i:len(pattern)+i] == pattern: #or genome[i:len(pattern)+i] == reverseComplement(pattern):
            positions += [i]
    return " ".join(str(e) for e in positions)


def readFile(filename, mode="rb"): 
    
    fin = contents = None
    try:
        fin = open(filename, mode)
        contents = fin.read()
    finally:
        if (fin != None): fin.close()
    return contents


#returns all k-mers forming (L,t)-clumps in genome
def findClumps(genome, k, L, t):

    frequentPatterns = []
    clumps = []
    frequencyArray = [] 
    for i in range(0, (4**k)):
        clumps += [0]
    for i in range(0,len(genome)-L+1):
        text = genome[i:(L+i)]
        frequencyArray = computeFrequencies(text, k)
        for j in range(0,(4**k)):
            if frequencyArray[j]>=t:
                clumps[j]=1
    for i in range(0,(4**k)):
        if clumps[i]==1:
            pattern = numberToPattern(i, k)
            frequentPatterns += [pattern]
    return " ".join(frequentPatterns) 

def betterClumpFinding(genome, k, L, t):
    frequentPatterns = []
    clumps = []
    frequencyArray = [] 
    for i in range(0, 4**k):
        clumps += [0]
    text = genome[0:L]
    frequencyArray = computeFrequencies(text, k)
    for i in range(0, 4**k):
        if frequencyArray[i]>= t:
            clumps[i] = 1
    for i in range(1,len(genome)-L+1):
        firstPattern = genome[i-1:(k+i-1)]
        index = patternToNumber(firstPattern)
        frequencyArray[index] -= 1
        lastPattern = genome[i+L-k:(i+L-k+k)]
        index = patternToNumber(lastPattern)
        frequencyArray[index] += 1
        if frequencyArray[index] >= t:
            clumps[index] = 1
    for i in range(0,4**k):
        if clumps[i]==1:
            pattern = numberToPattern(i, k)
            frequentPatterns += [pattern]
    return " ".join(frequentPatterns)

def betterClumpFindingFromFile(k, L, t):
    frequentPatterns = []
    clumps = []
    frequencyArray = []
    genome = readFile("texts/E-coli.txt")
    for i in range(0, 4**k):
        clumps += [0]
    text = genome[0:L]
    frequencyArray = computeFrequencies(text, k)
    for i in range(0, 4**k):
        if frequencyArray[i]>= t:
            clumps[i] = 1
    for i in range(1,len(genome)-L+1):
        firstPattern = genome[i-1:(k+i-1)]
        index = patternToNumber(firstPattern)
        frequencyArray[index] -= 1
        lastPattern = genome[i+L-k:(i+L-k+k)]
        index = patternToNumber(lastPattern)
        frequencyArray[index] += 1
        if frequencyArray[index] >= t:
            clumps[index] = 1
    for i in range(0,4**k):
        if clumps[i]==1:
            pattern = numberToPattern(i, k)
            frequentPatterns += [pattern]
    return " ".join(frequentPatterns)

def computeFrequencies(text,k):

    frequencyArray = []
    for i in range(0,(4**k)):
        frequencyArray += [0]
    for i in range(0, len(text)-k+1):
        pattern = text[i:(i+k)]
        j = patternToNumber(pattern)
        frequencyArray[j] += 1
    return frequencyArray #" ".join(str(e) for e in frequencyArray) 
    
def patternToNumber(pattern):
    if len(pattern)==0:
        return 0
    else:
        symbol = pattern[len(pattern)-1]
        prefix = pattern[:-1]
        return 4 * patternToNumber(prefix) + symbolToNumber(symbol)

def numberToPattern(index, k):
    if k==1:
        return numberToSymbol(index)
    else:
        prefixIndex = index/4
        remainder = index % 4
        symbol = numberToSymbol(remainder)
        #prefixPattern = numberToPattern(prefixIndex,k-1)
        #return prefixPattern + symbol
        return numberToPattern(prefixIndex,k-1) + symbol

def symbolToNumber(symbol):
    if symbol == "A":
        return 0
    elif symbol == "C":
        return 1
    elif symbol == "G":
        return 2
    else:
        return 3

def numberToSymbol(n):
    if n == 0:
        return "A"
    elif n==1:
        return "C"
    elif n==2:
        return "G"
    else:
        return "T"

#list of the skews as you run down the genome
def skew(genome):

    skews = []
    for i in range(0, len(genome)+1):
        skews += [0]
    for i in range(0,len(genome)):
        if genome[i] == "C":
            skews[i+1] = skews[i] - 1
        elif genome[i] == "G":
            skews[i+1] = skews[i] + 1
        else:
            skews[i+1]=skews[i]
    return skews    
    #return " ".join(str(e) for e in skews)

#indices of the minimum skews 
def minSkews(genome):

    skews = skew(genome)
    minSkewList = [] 
    currentMin = sys.maxsize
    for i in range(0,len(genome)):
        if skews[i]<currentMin:
            currentMin = skews[i]
            minSkewList = []
            minSkewList += [i]
        elif skews[i] == currentMin:
            minSkewList += [i]
    return " ".join(str(e) for e in minSkewList)
    

def tests():
    print(frequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4))
    print(reverseComplement("AAAACCCGGT"))
    print(patternMatch("TTTTTATTT", """TGGTTTTTATTTTTTATTACGTTTTTATTGTTTTTATGCTTTTTTATCCATTTTTTTATCGGCTTTTTATTTATGAGCCCCCCTTTTTATTTTTTTATGGTTTTTTATTTTTTTATTTTTTATCCTTTTTATTAATTTTTTATATTTTTATTTTTTATCATTTTTTATACTGTTTTTATGTTTTTATACAGAATCCGTTTTTATGTTTTTTATCTTTTTATACTTTTTTTATTTTTTATACACTGTTTTTATTTTTTATTTTTTATTTTTTATCCTTTTTATTTTTTATGGTTTTTATAGATGAACTTTTTTATATTTTTATATTTTTATTTTTTTATTTTTTATCTTTTTATTTTTTATCCAGTTTTTATTACGTTTTTTTATTCAATTTTTTATATTTTTATCATTTTTATGCTTTTTTATCTTTAAATTTTTATTCTTTTTATCGTTTTTATTTTTTATAAGGGTTTTTATAATTTTTATGAATTTTTTATTGGATTTTTATGAAATTTTTATCAAATTTTTATTTTTTATTTTTTATCGTCTTTTTTTTATGGTTTTTATGCTTTTTTATTTTTTATTTTTTTTTATGTTTTTATTAAATTTTTATCGACATTTTTATTGTTTTTATCATTAAGTGTAAGGAGACGGTTTTTATTGCTTTTTATGCGCCGTTTTTTATCGCGGTAGTCGCGGATTTTTATTTTTTTTTATCGGTTTTTATGAATATTTTTATCCGAGATTTTTATTCTTTTTATCCAAACATTTTTATGCTGAGTTTTTATAATTTTTTATAATCTCCTTTTTATACTCTCGGCTTTTTATCCTATAGCTTTTTATATTTTTATGTTTTTATAAGGTTTTTTATATTTTTATTTTTTATTTTTTATCTTTTTATGATATTTTTTATTCAGTTTTTATCTTTTTATTTTTTATCTTTTTATCTATTTTTATTTTTTATTTTTTATATTTTTATCAGGGCCATTTTTATATTTTTATAGCCTTTTTATCGCAGTCTTTTTATGTTTTTTATCCTTTTTTATAACCTTAGTCGGTTTTTATTTTTTATAAGTTTTTATAATTTTTATTCGCCACTTTTTATCATTTTTATTTTTTATTGTTTTTATTTTTTATATTTTTATGATTTTTATATTTTTATAGATTTTTATCTTTTTATGTTGATTTTTTATGTTTTTATTTTTTATGTTTTTATTAATTTTTATTTTTTATTTTTTATTTTTTATTTTTTATTTTTTATAACTTTTTATTTTTTATCTTTTTATGAAGTTTTTATGACCTTTTTATTATTTTTATTTTTTATTTTTTATAGCATTTTTATATTTTTATGGTGTTTTTATGTTTTTATAAATGATTTTTATCTTTTTATAGCCATTTTTTATTTTTTATTTTTTATGATTTTTATTTTTTATTTTTTATTGATTTTTATTAACTTTTTATCCCAGTTTTTATACTTTTTATGTCGGCGTATTTTTATTTTTTATTTTTTATATTTTTATTTTTTATCGATTTTTATCTTTTTATTTTTTATCCAGGTTTTTTATTGGGATTTTTATTTTTTATTGCTTTTTATATTTTTATAGTTTTTATGTGTTTTTATGTTTTTATCTTTTTATACCCATTTTTATCCGTTTTTATTTTTTATGTCATTTTTATTTTTTATTTTTTATATCAAGTACCGGGCTCTGTTTTTATTCGTTTTTTATCTTTTTATGTTTTTATATTTTTATTTTTTATTTTTTATGTTAATTTTTATTTTTTATCTTTTTTATCGGTTTTTATCGACATTTTTTATGTTTTTTATCTTTTTATCTTTTTATACAGTTGGTTTTTTATTGTTTTTATATTTTTATCTTTTTATTTTTTATTTTTTATGTGGAAAGCTTTTTATTTTTTATTCTACTTTTTATATTTTTATTGCACATGACCAGCTTTTTATTTTTTATTTTTTTATTTTTTATTTTTTTATTTTTTTATTTTTTATTTTTTATGTTTTTTATCATTTTTATTTTTTATGGGACCATTTTTATCTTTTTATCCTTTTTATGCTTTTTATACCTTTTTATGTTTTTATACCATTTTTATGGTTCTTTTTATGGTGGGTTTTTATAGATAACGGCTTTTTATACGACAATTTTTATTTTTTATTATTTTTATCTGGTGGTTTTTTATTAGTTCATTTTTTATTCAACTTTTTATTTTTTATTTTTTATCAGCCGTTTTTTTATTTTTTATCCTTTTTATACTTTTTATTTTTTATAGTTTTTATTTTTTTATTTTTTATTTTTTATTTTTTATTTTTTATGCTAGAGGTTTTTATTACGCTAATTTTTATTATTTTTATTTTTTATTTTTTATATTTTTATGTCCTCAATTTTTATCGATTTTTATTAAATTTTTATTTTTTATCCGTTTTTTATTTTTTTATTTTTTATATTTTTATTTTTTATGATCTGATTTTTTTATGGAATTTTTATTTTTTATTCTTTTTTATTTTTTATTCTTTTTATTTTTTTATTTTTTATGGCTACCGATGTTTTTATATTTTTTATTTTTTATTCTTTTTATTTTTTATCTTTTTATAATTTTTTATTTTTTATTCTTTTTATTTTTTATTTTTTATGTTTTTTTATTTTTTATATTTTTATTTTTTATTGGTTTTTATGTTTTTATGCCGATTATTTTTATTTTTTATGAGTATTTTTATACTTTTTATATTTTTATTTTTTATTTTTTTATCTTTTTATGCACAGTTTTTTATAATTTTTATTGCAGCTTAAAGGACTCAAGCATTTTTATATTTTTATGTTTTTATGTTTTTATTTTTTATAGCGTCACTACCCGGATGGATTTTTATCGTTTTTATGTTTTTTATTTTTTATAGTTTTTATTTTTTATTTGGTTTTTATCCCTTTTTATTTTTTTATATATGTCAGTTTTTATTGCTTTTTATACAACTTTTTTATGCTCTTTTTATCCGTTTTTATTTTTTATAATCTTTTTATTTCACTTTTTATTTTTTATCGTTTTTATTTTTTATTGGTTTTTATCTTTTTATAATTTTTATGGTTTTTATGCTAATGCCGCTTTTTATATTTTTATGTTTTTATTTTTTTATTTTTTATCATTCTTTTTATCGGATTTTTATCCTTCTGCCTTTTTTATTTGTTTTTTATTTTTTATTTTTTATTTTTTATTTTTTTTATCAGTTTTTATGTAATTTTTATGTTTTTATACTTTTTATTAATTTTTATTTTTTATTTTTTATTTTTTATCATTTTTTATGTTTTTATTTTTTTATGTTTTTATGTGTACTTAGTCGTGGTTTTTATCGTTTTTATGATTGTTTTTATGGCTAGTCAAATATTTTTATACCCCCAGGTTTTTATAACGTTTTTATCAAGACATTTTTATTTTTTATCAGCACAGATTTTTTATAGTTTTTTTATCCTGTGAGACACCACTTTTTATTCTTTTTATTTTTTATGCTTTTTATTTCTTTTTATCTTTTTATTTTTTATTTTTTATGATTGGTTTTTATCTCACTGCTGTTTTTATCATTTTTATCGGATGATTGTTTTTATTTTTTATCGTGTTAGGTTTTTATATTTTTATTCTTTTTATTATATTTTTATGGTTTTTATAGTTTTTATGGTTTTTATTTTTTATTTTTTTATTTTTTATACTTTTTATATTTTTATAGTTTTTATTCTTTTTTATACAGACTTTTTTATATTTTTTATTTTTTATCAGTTTTTATTTTTTATTTTGGTTTTTATGGCTTTTTATGATTTTTATTTTTTATGCTTTTTATTCCACGTTTTTATTTTTTATTTTTTATGAATAGCTTTTTATGTTTTTATGGGCTTTTTATTCTTTTTATTTTTTATTTTTTATGCTACAGTTTTTATTAACCCTTCGGCGCTTTTTATTTTTTTTATCTTTTTATTTTTTATTTCAGGGTTTTTTTATCTTTTTATTTTTTATGATTTTTATAGATTTTTATTTTTTATTATTTATTTTTATTTTTTATTTTTTTATGGAATTTTTATTATTTTTATCCGTTTTTATCATTTTTATCTTTTTTATATACGCGATTTTTATGTTTTTTATTGTTTTTATTTTTTATGCTTTTTATTTTTTATCTTTTTATTTTTTATAGAGCTTTTTATGACTTTTTATGTTTTTATTGCACGTTTTTATCGAAATTTTTATTTTTTATTCTTGCGGATTTTTTATATTTTTTATTTTTTTATATTGGAGTTTTTATTTTTTTATGTTTTTTATAATATGAGTTTTTATTTTTTATAAGTTTTTTATTTACTGACGCATTTTTATAGGGTTTTTTATGCTTATTTTTATCTGTCTGTTTTTTATTACTTTTTATATTTTTTATATAATTTTTTATATTTTTATTTTTTTATCTCTTTTTATGTTTTTATTTTTTATTTTTTTATATTTTTATCTTTTTATTTTTTATGTGTTTTTATGCGTTTTTATCTTTTTATTTTTTATTTTTTATTGTAGATATTTTTATTTTTTATTTTTTATTTTTTATGTTTTTATATTTTTATCTTTTTATCCTGTTTTTATCTTTTTATTTTTTATTTTTTATAATGTTTTTATTTTTTATCCGTTTTTTATGGAGGCTACATTTTTATATTTTTTATCATTTTTATTCCCTTTTTATTGCTTTTTTATTTTTTTATGGGTTTTTATCGTTTCGTGTTTTTATTGTTTTTATTTTTTATCGCTATTTTTATTTTTTTATCTTTTTATTTTTTATGCGCTTTTTATGGCTTTTTATAACGCACTTTTTTATACATTCTGGTTTTTTATATTTTTATGACGATTACGCCTTTTTATGTTTTTATTTTTTATGACTATTTATTTTTATGGTTTTTATTTTTTATGAACTTTTTATCCCTCTTTTTTATTTTTTATTAAGACTGTTTTTATCCCCTGGTTTTTATTGCACTCAGATGGATCTTTTTATATCTTTTTATCTTTTTATTTTTTTATTTTTTATTTTTTATTTTTTATTTTTTATTGATTTTTATTCTGTGTTTTTTATTTTTTATGATCTTTTTATGTTTTTTATGAGTTTTTATGATTTTTATTACTTTTTATCGTACTCTTTTTATTATTTTTTATATTTTTATTTTTTATCTTTTTATTGCCTTTTTATTTTTTTTATGAGCGCGACCGGTTTTTATCAGTTTTTATGTTTTTATTTTTTATGTTTTTATTTTTTTTATGAGTTTTTATTTTTTATTTTTTATATTTTTATTTTTTATAGGTGGAATCTGTCTTTTTATCCGGCTTTTTATTTTTTATTTTTTATAGGGCTCTTTTTTATCCGGAATCTTGTTTTTTATGTTTTTATTTTTTATTGATAGTATTTTTATATTTTTATCCTCTTTTTATCCTTCTTTTTATTTTTTATTTTTTATTTGGTTTTTATTTTTTATAGGTCACTTTTTATTTTTTATCTTTTTATGCTAAGCTTTTTATTTTTTATTTTTTATATCATTTTTATACTTTTTTATAGATTTTTTATGCGTTTTTTTATGTTTTTATTTTTTATGGCGACCTGAGGTGTTTTTTATCAACCAAGCTTTTTATCATTTTTATTTTTTATTTTTTATTTTTTATATTTTTTATATTTGTTTTTATTTTTTTTATCTTTTTATTTTTTATTCCTTTATTTTTATGCTTTTTATTTTTTATCTTTTTTTATAGCGTTTTTATCTTTTTATCTCTAAGGCAAGCCTTTTTTATATTTTTATTTTTTATAACCCATTTTTTTATCCCCCTTTTTTATCTTTTTATTATTTTTTATTGGTTTTTATTCCTTTTTATTTTTTATCGCTAGTTTTTATTTTTTATTTTTTATGTACCATTTTTATAGTTTTTATTGCCTTCCTTTTTATAGTTTTTTATAGGCTTTTTATGTGTTTTTATTTTTTATAGATTTTTATGTATTTTTATCTTTTTTTATTTTTTTATGGCTTTTTATTGTTTTTTATTTTTTATTATTTTTATAAGGTTTTTATATTTTTATATTTTTTATAGGTTCAAGTTTTTATTTTTTATTTTTTTATTTTTTATATTTTTATAACGGTTTTTATGTGGTGTTTTTTATATTTTTTATTTTTTATACTTTTTATAATTTTTATAATCTTTTTTATTTTTTATTTTTTTATTTGGTTTTTATTTACTATTTTTATTTTTTTATTTTTTATTGGTTTTTATTTTTTATTTTTTATCCTTTTTATCATAACTCCAGTTTTTTATATTGGTTCATTTTTATAACCTTTTTATTAGGGTTTTTATTTTTTATAGAGTTTTTATAGATCGGGGCGTTAATTTTTTATTTTTTATCTAGCATTTTTATTTTTTATTCGGACATTTTTTTATATTTTTATGCTTTTTATCTTTTTATTGGTTTTTATACAGTTTTTTATGTTTTTATGCGCATGTGGGGGTTTTTATCGTACTTTTTATTGTTTTTATTTTTTATTTTTTATATTTTTATTTGTTTTTATATTTTTATTTTTTATTCATATTTTTATCTTTTTTATGCGTTTTTTATAGTTTTTATGGGACGCTTTTTATCTTTTTATGCATTTTTATTTTTTATTTTTTATTTTTTATTTTTTATCCTTTTTATGGTGCTATTTTTATTCTTTTTATATCCTTTTTATTTTTTATGTTTTTATTCTTTTTATTGGGCAGTACCCTTTTTATCAGTTGCTTTTTATTTTTTATCCCATTTTTATTTTTTATTGCTTTTTATGGACTTTTTATGCGTTTTTATCCTTTTTTATAACGCCCTTTTTATAAATGGTTTTTATACCCCTTTTTATTTTTTATTGTTTTTATTGACTTTTTATTTTTTATTTTTTATATTTTTATGTAATGCCGGTTTTTATTTTTTATTTTTTATATTTTTATGTTTTTATGGGTCGTTTTTTATGGTATTTTTATTTTTTATGCTGATTGACTTTTTATTTTTTATGACTAATTTTTATCCAGATTTTTATTTTTTATTAATTTTTATCAAGTGTGTTTTTTATTTTTTATTTTTTATACTTCTTTTTATTTTTTATATTAACCTTTTTATTGTTTTTATCTTTTTATGTTTTTATAGACTTTTTATTTTTTATTTCAGAGTTTTTATTCATAGGGTTTTTATGTTTTTATTGAAGATTCTTTTTTATTAGCTTTTTATTTTTTATAGAGGCGTTTTTATTATAGGAGTAATTTTTATTGCATTTTTATCTTTTTTATATTTTTATGTTTTTATATTCTTTTTTTATTGGTCCTTTTTATCGTTTTTATGGCAACCGTTTTTATGTTCCAATGACTTTTTATTTAACTCGTTTTTATTTTTTATTAGGGATTTTTTATCTTTTTATATTTTTATGCATTACCCTTTTTATTTTTTATTGATTTTTATAAGAACGGGATTTTTTATCCCTCTTTTTATGACGAGTTTTTATCAATTTTTATCATTTTTATGGCTTTTTATTTAGGTTTTTATGATTATTTTTATGTTTTTATTTTTTATGTTTTTATATGTTTTTATAGTTTTTATTTTTTATGGTTTTTTTATGTTTTTATGCTTTTTTATTTTTTTTATGGTTTTTTATGTGTTTTTATAGCCAATTTTTTATTGGTTGTTTTTATCTTTTTATCTTTTTTATGGCGCCTTTTTATTGTTTTTATCGTTTTTATGCGTGAGTTTTTATTTTTTATATTTTTATCATTTTTATTTTTTTATGTTTTTATTTTTTATATTTTTATTTTTTATAATTTTTATTGTTTTTATCTTTTTATCTTTTTATATACTTTTTATTTTTTATCGAAATTTTTATTTTTTTATTTTTTATTAATTTTTTATTTTTTATACTTTATTTTTATCTTTATATGGGCTTTTTTATATTTTTTTATTTCTAGCCGGCGATTTTTATTAGCTTTTTATCAATTTTTATTTTTTATCAGTTTTTTATCTTTTTATTTTTTATCCTTTTTTTATCTTTTTATGATCGTTTTTATCGAGTTTTTATCTCACATTTTTATACTCACTTTTTTATTTTTTATATTTTTATGTATTTTTTATCTTTTTATATTTTTATTCGGTTTTTATGCGTTTTTATTATTTTTATACCCTACTTTTTATTTTTTATGTTTTTATGATTTTTATGATTTTTATTTTTTATCCCTTTTTATGGCCCCTTTTTATGTTTTTATTTTTTTATTTTTTTATATTTTTATCTGCTTTTTATACCTGAATTTTTATTAGTAATGCCGGCTTTTTTTATGTTTTTATGTCGTTTTTATCATATTTTTATTTTTTATCTTTTTTATTTTTTTATGATTTTTATAGGATTTTTTATTTTTTATTTTTTATTCTTTTTTATCTCGGTTGATTGCATTCACTAGATTTTTATAGTAGATGCGAGCTTCTAGCCCCCACGATTTTTATACACGCCCTTTTTATTCTTTTTATTAAACGCCCTTTTTATCATTTTTATATTTTTATAATTTTTATCAGTTTTTATTTTTTATTCGCCTTTTTATATTTTTATTATTTTTATTTTTTATAGCCAAGTGTTTTTATCATGAATTTTTTATATTTTTATCCGTTTTTATTGCGCCTATTTTTATTACTGAACTATTTTTTATTTTTAATTTTTATTTTTTTATTTTTTATCCTTTTTATTTTTTATTTTTTATATTTTTATATTTTTATAATTTTTTATCATTTTTATTTTTTTATATTTTTTATTTTTTTATCGTATTTTTATCCATTTTTATGTTTTTTTTATGTGTTTTTATGTTTTTATTTTTTATTGGTCAATTTTTATTTTTTATGGTAGCGACCTTTTTATACTTTTTATCACGATCTTTTTATTGATTTTTATCGCACTTTTTATTATTTTTATCGGGCATTTTTATCTTTTTATTTTTTATTGTGCAAATCTTTTTATGGTTTTTATCGCCATTTTTATGTTTTTATGGTATCTTTTTATGTTTTTATGCAGTCCAGTTTTTATACTGAATTTTTATCTTTTTATTGTTTTTATCGCGATATTTTTATCACCTTTTTATTGACTTTTTTTATTTTTTATTTTGGTTTTTTATCAGTTTTTATGTTATTTTTATATTTTTATTTTTTTATGTTTTTATTACTTTTTATTTTTTATACATTTTTATTTTTTATTTTTTATTTTTTATGCGTTTTTATGACCATTTTTATGCGGTTTTTATGAATTTTTATCTTTTTATATTTTTATGTTTTTATGGACTCTTTTTATATGTGTTTTTATATTTTTATGGTTTTTATGCCCCGACAAGTTTTTATCTTTTTAT"""))

def hammingDistance(s1,s2):

    if len(s1)==len(s2)==0:
        return 0
    else:
        if s1[0]==s2[0]:
            return hammingDistance(s1[1:],s2[1:])
        else:
            return hammingDistance(s1[1:],s2[1:]) + 1


#d is the maximum allowed mismatch (inclusive) 
def approximatePatternMatching(pattern, genome, d):
    
    positions = []
    for i in range(0,len(genome)-len(pattern)+1):
        if hammingDistance(genome[i:len(pattern)+i],pattern)<=d:
            positions += [i]
    return " ".join(str(e) for e in positions)
    return len(positions)

def approximatePatternCount(pattern, genome, d):
    
    positions = []
    for i in range(0,len(genome)-len(pattern)+1):
        if hammingDistance(genome[i:len(pattern)+i],pattern)<=d:
            positions += [i]
    return len(positions)

def frequentWordsWithMismatches(text, k, d):

    frequentPatterns = []
    close = []
    frequencyArray = [] 
    for i in range(0, 4**k):
        close += [0]
        frequencyArray += [0]
    for i in range(0,len(text)-k+1):
        neighborhood = neighbors(text[i:(k+i)],d)
        for pattern in neighborhood:
            index = patternToNumber(pattern)
            close[index] = 1
    for i in range(0,4**k):
        if close[i]==1:
            pattern = numberToPattern(i,k)
            frequencyArray[i] = approximatePatternCount(pattern, text, d)
    maxCount = max(frequencyArray)
    for i in range(0, 4**k):
        if frequencyArray[i] == maxCount:
            pattern = numberToPattern(i,k)
            frequentPatterns += [pattern]
    #return frequentPatterns
    return " ".join(frequentPatterns)

#Doesn't quite work. I'm not sure why. Will look at it again later to try and figure it out. 
def frequentWordsMismatchesReverse(text, k, d):

    frequentPatterns = []
    close = []
    frequencyArray = [] 
    for i in range(0, 4**k):
        close += [0]
        frequencyArray += [0]
    for i in range(0,len(text)-k+1):
        neighborhood = neighbors(text[i:(k+i)],d)
        neighborhood2 = neighbors(reverseComplement(text[i:(k+i)]),d)
        for pattern in neighborhood:
            index = patternToNumber(pattern)
            close[index] = 1
        for pattern in neighborhood2:
            index = patternToNumber(pattern)
            close[index] = 1 
    for i in range(0,4**k):
        if close[i]==1:
            pattern = numberToPattern(i,k)
            frequencyArray[i] = approximatePatternCount(pattern, text, d)
    patternsToAdd = []
    maxCount = 0
    for i in range(0, 4**k):
        reverseCount = frequencyArray[patternToNumber(reverseComplement(numberToPattern(frequencyArray[i],k)))]
        if frequencyArray[i] + reverseCount > maxCount:
            patternsToAdd = [] 
            pattern = numberToPattern(i,k)
            patternsToAdd += [pattern]
            if reverseCount != 0:
                pattern = reverseComplement(numberToPattern(frequencyArray[i],k))
                patternsToAdd += [pattern]
            maxCount = frequencyArray[i] + frequencyArray[patternToNumber(reverseComplement(numberToPattern(frequencyArray[i],k)))]
        elif frequencyArray[i] + reverseCount == maxCount:
            pattern = numberToPattern(i,k)
            patternsToAdd += [pattern]
            if reverseCount != 0:
                pattern = reverseComplement(numberToPattern(frequencyArray[i],k))
                patternsToAdd += [pattern]
            maxCount = frequencyArray[i] + frequencyArray[patternToNumber(reverseComplement(numberToPattern(frequencyArray[i],k)))]
    #return frequentPatterns
    return " ".join(list(set(patternsToAdd)))
        
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

#tests()
sys.setrecursionlimit(10000000)
