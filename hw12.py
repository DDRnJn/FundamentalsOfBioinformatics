#Bioinformatics hw10
#Dhruv Ranjan

from collections import defaultdict
import sys
import math
import copy

def BWTWrapper(fileName):

    contents = open(fileName).readlines()
    text = contents[0]
    b = BWT(text)
    output = "".join(str(i) for i in b)
    print(output)

def BWT(text):

    rotations = []
    bMatrix = []
    for i in range(len(text)):
        current = text[-i:]+text[:-i]
        rotations.append(current)
    bMatrix = sorted(rotations)
    bwt = map(lambda x: x[-1], bMatrix)
    return bwt

def BWTToFirst(bwt):

    first = sorted(bwt)
    return first

def BWTInverseWrapper(fileName):

    contents = open(fileName).readlines()
    text = contents[0]
    bInverse = BWTInverse(text)
    print(bInverse)

def BWTInverse(text):

    currentIndex = 0
    inverse = ""
    fCol = "".join(BWTToFirst(text))
    currentChar = fCol[currentIndex]
    inverse+=currentChar
    currentIndex = text.find(currentChar)
    for i in range(len(text)-1):
        currentChar = fCol[currentIndex]
        inverse+=currentChar
        nIndex = findCharN(fCol,currentChar,currentIndex)
        currentIndex = findNChar(text, currentChar, nIndex)
    final = inverse[1:]+inverse[0]
    return final

def findCharN(text, char, n):

    return text[:n+1].count(char)-1

def findNChar(text, char, n):

    split = text.split(char,n+1)
    return len(text)-len(split[-1])-1

def BWMatchingWrapper(fileName):

    contents = open(fileName).readlines()
    text = contents[0]
    p = contents[1]
    patterns = p.split()
    counts = BWMatching(text,patterns)
    fout = open("texts/hw12/BWMatchingAnswer.txt","wt")
    output = " ".join(str(i) for i in counts)
    #print output
    fout.write(output)

def BWMatching(BW, patterns):

    fCol = "".join(BWTToFirst(BW))
    lastToFirst = LTF(BW,fCol)
    counts = []
    for i in patterns:
        counts.append(patternMatch(BW,i,lastToFirst))
    return counts

def patternMatch(lCol,pattern,lastToFirst):

    top = 0
    bottom = len(lCol)-1
    while top <= bottom:
        if len(pattern)!=0:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            currentPositions = lCol[top:bottom+1]
            if symbol in currentPositions:
                topIndex = top + currentPositions.index(symbol)
                bottomIndex = bottom - currentPositions[::-1].index(symbol)
                top = lastToFirst[topIndex]
                bottom = lastToFirst[bottomIndex]
            else:
                return 0
        else:
            return bottom-top+1

def LTF(BW,fCol):

    lastToFirst = []
    for i in range(len(BW)-1):
        currentIndex = findCharN(BW,BW[i],i)
        lastToFirst.append(findNChar(fCol,BW[i],currentIndex)-1)
    return lastToFirst
                        
def betterBWMatchingWrapper(fileName):

    contents = open(fileName).readlines()
    text = contents[0].strip()
    p = contents[1]
    patterns = p.split()
    counts = betterBWMatching(text,patterns)
    fout = open("texts/hw12/betterBWMatchingAnswer.txt","wt")
    output = " ".join(str(i) for i in counts)
    fout.write(output)

def betterBWMatching(BW, patterns):

    fCol = "".join(BWTToFirst(BW))
    counts = []
    countDict = getCD(BW)
    FO = firstOccurence(fCol)
    for i in patterns:
        counts.append(betterPatternMatch(FO,BW,i,countDict))
    return counts

def getCD(BW):

    alphabet = list(set(BW))
    countDict = defaultdict(lambda: defaultdict(int))
    for i in alphabet:
        countDict[0][i]=0
    for i in range(len(BW)):
        current = copy.copy(countDict[i])
        current[BW[i]]+=1
        countDict[i+1]=copy.copy(current)
    return countDict

def firstOccurence(fCol):

    FO = defaultdict(int)
    alphabet = list(set(fCol))
    for i in alphabet:
        FO[i]=fCol.index(i)
    return FO

def betterPatternMatch(FO,lCol,pattern,count):

    top = 0
    bottom = len(lCol)-1
    while top <= bottom:
        if len(pattern)!=0:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            currentPositions = lCol[top:bottom+1]
            if symbol in currentPositions:
                top = FO[symbol]+count[top][symbol]
                bottom = FO[symbol]+count[bottom+1][symbol]-1
            else:
                return 0
        else:
            return bottom-top+1

#multiplePatternMatchingWrapper("texts/hw12/multiplePatternMatching.txt")
def multiplePatternMatchingWrapper(fileName):

    contents = open(fileName).readlines()
    text = contents.pop(0).strip()+"$"
    patterns = [i.strip() for i in contents]
    counts = multiplePatternMatching(text,patterns)
    counts = sorted(map(str,counts))
    fout = open("texts/hw12/multipleMatchingAnswer.txt","wt")
    output = " ".join(str(i) for i in counts)
    print(output)
    print(len(output))
    fout.write(output)

def multiplePatternMatching(text, patterns):

    BW = BWT(text)
    SA = suffixArray(text)
    fCol = "".join(BWTToFirst(BW))
    counts = []
    countDict = getCD(BW)
    FO = firstOccurence(fCol)
    for i in patterns:
        counts+=mPatternMatch(FO,BW,i,countDict,SA)
    return counts

def mPatternMatch(FO,lCol,pattern,count,SA):

    top = 0
    bottom = len(lCol)-1
    res = []
    while top <= bottom:
        if len(pattern)!=0:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            currentPositions = lCol[top:bottom+1]
            if symbol in currentPositions:
                top = FO[symbol]+count[top][symbol]
                bottom = FO[symbol]+count[bottom+1][symbol]-1
            else:
                return res
        else:
            for i in range(top,bottom+1):
                res.append(SA[i])
            return res

def suffixArray(text):

    suffixes = []
    positions = []
    for i in range(len(text)):
        suffixes.append(text[-i:])
    positions+=[0]
    for i in range(1,len(text)):
        positions.append(len(text)-i)
    positions2 = sorted(zip(suffixes,positions), key=lambda x: x[0])
    output = [i[1] for i in positions2]
    return output

def approxMatchingWrapper(fileName):

    contents = open(fileName).readlines()
    text = contents.pop(0)
    d = contents.pop(len(contents)-1)
    patterns = contents.split().strip()
    positions = approxMatching(text,patterns,d)

def approxMatching(text,patterns,d):

    BW = BWT(text)
    SA = suffixArray(text)
    fCol = "".join(BWTToFirst(BW))
    counts = []
    countDict = getCD(BW)
    FO = firstOccurence(fCol)
    n = len(text)
    k = n/(d+1)
    for i in patterns:
        n = len(text)
        k = n/(d+1)
        
























        
    
