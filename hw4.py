#Bioinformatics hw4
#Dhruv Ranjan

import copy
from operator import itemgetter
from collections import Counter

def rnaCodonTable(fileName):

    contents = open(fileName).readlines()
    rnaCodonDict = {}
    for line in contents:
        (codon, AA) = line.split()
        rnaCodonDict[codon] = AA
    return rnaCodonDict

def aminoAcidWeightTable(fileName):

    contents = open(fileName).readlines()
    aminoAcidWeightDict = {}
    for line in contents:
        (AA, weight) = line.split()
        aminoAcidWeightDict[AA] = weight
    return aminoAcidWeightDict

#where geneticCode is a codon -> amino acid dict        
def rnaTranslate(pattern, geneticCode):

    peptide = ""
    codons = getCodons(pattern)
    for s in codons:
        AA = geneticCode.get(s)
        if AA != "STOP":
            peptide += AA
    return peptide

def getCodons(pattern):

    codons = []
    codonNumber = len(pattern)-2
    for s in xrange(0,codonNumber,3):
        currentCodon = pattern[s] + pattern[s+1] + pattern[s+2]
        codons += [currentCodon]
    return codons

def rnaTranslateWrapper(fileName):

    geneticCode = rnaCodonTable("texts/hw4/rnaCodonTable.txt")
    contents = open(fileName).read()
    pattern = contents
    peptide = rnaTranslate(pattern, geneticCode)
    fout = open("texts\hw4\hw4q1Answer.txt", "wt")
    fout.write(peptide)

#find the reverse complement of the str pattern
def reverseComplement(pattern):

    newPattern = []
    for i in xrange(len(pattern)):
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

def reverseComplementRna(pattern):

    newPattern = []
    for i in xrange(len(pattern)):
        if pattern[i] == "A":
            newPattern += "U"
        elif pattern[i] == "U":
            newPattern += "A"
        elif pattern[i] == "C":
            newPattern += "G"
        elif pattern[i] == "G":
            newPattern += "C"
    newPattern = newPattern[::-1]
    return "".join(newPattern)

def peptideEncoding(text, peptide, geneticCode):
    rna = dnaToRna(text)
    rnaSubstrings = []
    dnaSubstrings = []
    codonLength = 3
    substringLength = len(peptide)*codonLength
    for s in xrange(0,len(rna)):
        substring = rna[s:s+substringLength]
        if len(substring) == substringLength:
            currentPeptide = rnaTranslate(substring, geneticCode)
            currentReversePeptide = rnaTranslate(reverseComplementRna(substring), geneticCode)
            if currentPeptide == peptide:
                rnaSubstrings += [substring]
            if currentReversePeptide == peptide:
                rnaSubstrings += [substring]
    for s in rnaSubstrings:
        dnaSubstrings += [rnaToDna(s)]
    return dnaSubstrings

def peptideEncodingWrapper(fileName):

    geneticCode = rnaCodonTable("texts/hw4/rnaCodonTable.txt")
    contents = open(fileName).readlines()
    text = contents[0].rstrip()
    contents.pop(0)
    peptide = contents[0]
    dnaSubstrings = peptideEncoding(text, peptide, geneticCode)
    print str(dnaSubstrings)
    fout = open("texts\hw4\hw4q2Answer.txt", "wt")
    for s in dnaSubstrings:
        fout.write(s)
        fout.write("\n")

def dnaToRna(text):

    rna = ""
    for s in text:
        if s == "T":
            rna +=  "U"
        else:
            rna += s
    return rna

def rnaToDna(text):

    dna = ""
    for s in text:
        if s == "U":
            dna += "T"
        else:
            dna += s
    return dna

def numCyclicSubseqs(n):

    if n==1:
        return 1
    elif n==2:
        return 2
    else:
        return (n-1) * (n)

def cyclicSpectrum(peptide, AAMassDict):

    prefixMasses = []
    prefixMasses += [0]
    for i in xrange(len(peptide)):
        prefixMasses += [prefixMasses[i] + int(AAMassDict.get(peptide[i]))]
    peptideMass = prefixMasses[len(peptide)]
    cyclicSpectrum = []
    cyclicSpectrum += [0]
    for i in xrange(0,len(peptide)-1):
        for j in xrange(i+1, len(peptide)):
            cyclicSpectrum += [prefixMasses[j]-prefixMasses[i]]
            if i >= 0 and j < len(peptide):
                cyclicSpectrum += [peptideMass - (prefixMasses[j] - prefixMasses[i])]
    return sorted(cyclicSpectrum) + [peptideMass]

def linearSpectrum(peptide, AAMassDict):

    prefixMasses = []
    prefixMasses += [0]
    for i in xrange(len(peptide)):
        prefixMasses += [prefixMasses[i] + int(AAMassDict.get(peptide[i]))]
    linearSpectrum = []
    linearSpectrum += [0]
    for i in xrange(0,len(peptide)-1):
        for j in xrange(i+1, len(peptide)):
            linearSpectrum += [prefixMasses[j]-prefixMasses[i]]
    return sorted(linearSpectrum) + [peptideMass]

def linearSpectrumNum(peptide, AAMassDict):

    prefixMasses = []
    prefixMasses += [0]
    for i in xrange(len(peptide)):
        prefixMasses += [prefixMasses[i] + int(peptide[i])]
    peptideMass = prefixMasses[len(peptide)]
    linearSpectrum = []
    linearSpectrum += [0]
    for i in xrange(0,len(peptide)-1):
        for j in xrange(i+1, len(peptide)):
            linearSpectrum += [prefixMasses[j]-prefixMasses[i]]
    return sorted(linearSpectrum) + [peptideMass]

def cyclicSpectrumNum(peptide, AAMassDict):

    prefixMasses = []
    prefixMasses += [0]
    for i in xrange(len(peptide)):
        prefixMasses += [prefixMasses[i] + int((peptide[i]))]
    peptideMass = prefixMasses[len(peptide)]
    cyclicSpectrum = []
    cyclicSpectrum += [0]
    for i in xrange(0,len(peptide)-1):
        for j in xrange(i+1, len(peptide)):
            cyclicSpectrum += [prefixMasses[j]-prefixMasses[i]]
            if i >= 0 and j < len(peptide):
                cyclicSpectrum += [peptideMass - (prefixMasses[j] - prefixMasses[i])]
    return sorted(cyclicSpectrum) + [peptideMass]

def cyclicSpectrumWrapper(fileName):

    contents = open(fileName).read()
    peptide = contents
    AAMassDict = aminoAcidWeightTable("texts/hw4/aminoAcidWeightTable.txt")
    cyclicSpec = cyclicSpectrum(peptide, AAMassDict)
    fout = open("texts\hw4\hw4q3Answer.txt", "wt")
    print " ".join(str(i) for i in cyclicSpec)
    fout.write(" ".join(str(i) for i in cyclicSpec))

def numLinearSubseqs(n):

    if n==1:
        return 1
    elif n ==2:
        return 4
    else:
        return (((n-1)*((n-1)+3))/2) + 2

def cyclopeptideSequencing(spectrum, AAWeightList, AAMassDict):

    peptides = [[]]
    peptidesToKeep = []
    finalPeptides = []
    parentMass = int(spectrum[len(spectrum)-1])
    while len(peptides) != 0:
        peptides = expand(peptides, AAWeightList)
        for peptide in peptides:
            keep = 0
            if (sum(map(int, peptide))) == parentMass:
                if map(str, cyclicSpectrumNum(peptide, AAMassDict)) == spectrum:
                    finalPeptides += [peptide]
                    keep = 1
            elif checkConsistency(peptide, spectrum, AAMassDict)==False:
                keep = 1
            if keep == 0:
                peptidesToKeep += [peptide]
        peptides = peptidesToKeep
        peptidesToKeep = []
    return finalPeptides

def cyclopeptideSequencingWrapper(fileName):

    fin = open(fileName, "r")
    contents = fin.readlines()
    contents = contents[0]
    spectrum = []
    currentNum = ""
    for s in contents:
        if s == " ":
            spectrum += [currentNum]
            currentNum = ""
        else:
            currentNum += s
    spectrum += [currentNum]
    spectrum[len(spectrum)-1] = spectrum[len(spectrum)-1].rstrip()
    spectrum = map(int, spectrum)
    spectrum.sort
    spectrum = map(str, spectrum)
    AAWeights = open("texts/hw4/aminoAcidWeightTable.txt").readlines()
    AAWeightList = []
    for line in AAWeights:
        (AA, weight) = line.split()
        AAWeightList += [weight]
    AAWeightList = sorted(list(set(AAWeightList)))
    AAMassDict = aminoAcidWeightTable("texts/hw4/aminoAcidWeightTable.txt")
    peptides = cyclopeptideSequencing(spectrum, AAWeightList, AAMassDict)
    fout = open("texts\hw4\hw4q4Answer.txt", "wt")
    peptides = formatPeptides(peptides)
    print " ".join(peptides)
    fout.write (" ".join(peptides))
    fin.close()
    fout.close() 

def formatPeptides(peptides):

    newPeptides = []
    for peptide in peptides:
        if len(peptide) == 1:
            newPeptides += [peptide[0]]
        else:
            first = peptide.pop(0)
            last = peptide.pop(len(peptide)-1)
            newPeptides += [first + "-" + "-".join(peptide) + "-" + last]
    return newPeptides

def checkConsistency(peptide, spectrum, AAMassDict):

    linSpec = linearSpectrumNum(peptide, AAMassDict)
    linSpec = map(str, linSpec)
    for AA in linSpec:
        if linSpec.count(AA) > spectrum.count(AA):
            return False
    return True

def mass(peptide):

    sum = 0
    for AA in peptide:
        sum += int(AA)
    return str(sum)

def parentMass(spectrum):

    return spectrum[len(spectrum)-1]

def expand(peptides, AAWeightList):

    currentLength = len(peptides[0])
    finalLength = currentLength + 1
    newPeptides = []
    for peptide in peptides:
        for j in xrange(0,len(AAWeightList)):
            newPeptides += [peptide + [AAWeightList[j]]]
    for i in xrange(len(newPeptides)):
        if newPeptides[i] == []:
            newPeptides.pop(i)
    return newPeptides

def cyclopeptideScore(peptide, spectrum, AAMassDict):

    score = 0
    newSpec = copy.copy(spectrum)
    theoreticalSpec = cyclicSpectrum(peptide, AAMassDict)
    for AA in theoreticalSpec:
        if str(AA) in newSpec:
            score += 1
            newSpec.remove(str(AA))
    return score

def cyclopeptideScoreNum(peptide, spectrum, AAMassDict):

    score = 0
    newSpec = copy.copy(spectrum)
    theoreticalSpec = cyclicSpectrumNum(peptide, AAMassDict)
    for AA in theoreticalSpec:
        if str(AA) in newSpec:
            score += 1
            newSpec.remove(str(AA))
    return score

def linearScore(peptide, spectrum, AAMassDict):
    
    score = 0
    newSpec = copy.copy(spectrum)
    theoreticalSpec = linearSpectrumNum(peptide, AAMassDict)
    for AA in theoreticalSpec:
        if str(AA) in newSpec:
            score += 1
            newSpec.remove(str(AA))
    return score
            
def cyclopeptideScoreWrapper(fileName):

    AAMassDict = aminoAcidWeightTable("texts/hw4/aminoAcidWeightTable.txt")
    contents = open(fileName).readlines()
    peptide = (contents.pop(0)).rstrip()
    contents = contents[0]
    spectrum = []
    currentNum = ""
    for s in contents:
        if s == " ":
            spectrum += [currentNum]
            currentNum = ""
        else:
            currentNum += s
    spectrum += [currentNum]
    spectrum[len(spectrum)-1] = spectrum[len(spectrum)-1].rstrip()
    score = cyclopeptideScore(peptide, spectrum, AAMassDict)
    print score

def leaderboardCyclopeptideSequencing(spectrum, N, AAMassDict, AAWeightList):

    leaderboard = [[]]
    leaderPeptide = []
    peptidesToKeep = []
    finalPeptides = []
    parentMass = int(spectrum[len(spectrum)-1])
    while len(leaderboard) != 0:
        leaderboard = expand(leaderboard, AAWeightList)
        for peptide in leaderboard:
            keep = 0
            if (sum(map(int, peptide))) == parentMass:
                if cyclopeptideScoreNum(peptide, spectrum, AAMassDict) > cyclopeptideScoreNum(leaderPeptide, spectrum, AAMassDict):
                    leaderPeptide = peptide
            elif (sum(map(int,peptide))) > parentMass:
                keep = 1
            if keep == 0:
                peptidesToKeep += [peptide]
        leaderboard = peptidesToKeep
        leaderboard = trim(leaderboard, spectrum, N, AAMassDict)
        peptidesToKeep = []
    return leaderPeptide

def trim(leaderboard, spectrum, N, AAMassDict):

    linearScores = []
    for j in xrange(len(leaderboard)):
        peptide = leaderboard[j]
        linearScores += [(peptide, linearScore(peptide, spectrum, AAMassDict))]
    linearScores = sorted(linearScores,key=itemgetter(1))
    linearScores.reverse()
    newLinearScores = copy.copy(linearScores)
    for j in xrange(N+1, len(linearScores)):
        if linearScores[j] < linearScores[N]:
            newLinearScores[j:len(newLinearScores)] = []
    newLeaderboard = []
    for j in newLinearScores:
        newLeaderboard += [j[0]]
    return newLeaderboard
            
def leaderboardCyclopeptideSequencingWrapper(fileName):

    AAMassDict = aminoAcidWeightTable("texts/hw4/aminoAcidWeightTable.txt")
    contents = open(fileName).readlines()
    N = contents[0].rstrip()
    contents.pop(0)
    N = int(N)
    spectrum = []
    currentNum = ""
    contents = contents[0]
    for s in contents:
        if s == " ":
            spectrum += [currentNum]
            currentNum = ""
        else:
            currentNum += s
    spectrum += [currentNum]
    spectrum[len(spectrum)-1] = spectrum[len(spectrum)-1].rstrip()
    AAWeights = open("texts/hw4/aminoAcidWeightTable.txt").readlines()
    AAWeightList = []
    for line in AAWeights:
        (AA, weight) = line.split()
        AAWeightList += [weight]
    AAWeightList = sorted(list(set(AAWeightList)))
    finalPeptide = leaderboardCyclopeptideSequencing(spectrum, N, AAMassDict, AAWeightList)
    finalPeptide = formatPeptides([finalPeptide])
    print finalPeptide[0]
    
def spectralConvolution(spectrum):

    convolution = []
    for i in xrange(len(spectrum)):
        for j in xrange(len(spectrum)):
            if spectrum[i] - spectrum[j] > 0:
                convolution += [spectrum[i]-spectrum[j]]
    return convolution

def spectralConvolutionWrapper(fileName):

    contents = open(fileName).read()
    spectrum = []
    currentNum = ""
    for s in contents:
        if s == " ":
            spectrum += [currentNum]
            currentNum = ""
        else:
            currentNum += s
    spectrum += [currentNum]
    spectrum[len(spectrum)-1] = spectrum[len(spectrum)-1].rstrip()
    spectrum = (map(int, spectrum))
    spectrum.sort()
    print spectrum
    convolution = spectralConvolution(spectrum)
    print " ".join(map(str, convolution))
    fout = open("texts\hw4\spectralConvolutionAnswer.txt", "wt")
    fout.write(" ".join(map(str,convolution)))

def convolutionCyclopeptideSequencing(M, N, spectrum):

    convolution = spectralConvolution(map(int,spectrum))
    alphabet = mMostFrequent(M, convolution)
    finalPeptide = leaderboardCyclopeptideSequencing(spectrum, N, alphabet, alphabet)
    return finalPeptide

def mMostFrequent(M, convolution):

    alphabet = []
    counts = []
    seen = []
    for AA in convolution:
        if AA >= 57 and AA <=200 and AA not in seen:
            count = convolution.count(AA)
            seen += [AA]
            counts += [(AA, count)]
    counts = sorted(counts,key=itemgetter(1))
    counts.reverse()
    newCounts = copy.copy(counts)
    for j in xrange(M+1, len(counts)):
        if counts[j] < counts[M]:
            newCounts[j:len(newCounts)] = []
    for j in newCounts:
        alphabet += [j[0]]
    print map(str, sorted(alphabet))
    return map(str, sorted(alphabet))

def convolutionCyclopeptideSequencingWrapper(fileName):

    contents = open(fileName).readlines()
    M = contents[0].rstrip()
    N = contents[1].rstrip()
    M = int(M)
    N = int(N)
    contents.pop(0)
    contents.pop(0)
    contents = contents[0]
    spectrum = []
    currentNum = ""
    for s in contents:
        if s == " ":
            spectrum += [currentNum]
            currentNum = ""
        else:
            currentNum += s
    spectrum += [currentNum]
    print spectrum
    spectrum[len(spectrum)-1] = spectrum[len(spectrum)-1].rstrip()
    spectrum = (map(int, spectrum))
    spectrum.sort()
    spectrum = (map(str, spectrum))
    finalPeptide = convolutionCyclopeptideSequencing(M, N, spectrum)
    finalPeptide = formatPeptides([map(str,finalPeptide)])
    print finalPeptide[0]

def tester():

    AAMassDict = aminoAcidWeightTable("texts/hw4/aminoAcidWeightTable.txt")
    s = [0,71,113,129,147,200,218,260,313,331,347,389,460]
    s = map(str, s)
    print linearSpectrumNum([71,147], AAMassDict)
    print linearScore([71,147],s,AAMassDict)
