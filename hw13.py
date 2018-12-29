#Bioinformatics hw10
#Dhruv Ranjan

from collections import defaultdict
import sys
import math
import copy


def hiddenPathProbWrapper(fileName):

    contents = open(fileName).readlines()
    path = contents.pop(0).strip()
    contents.pop(0) #pop dashed line
    states = contents.pop(0).split()
    contents.pop(0) #pop dashed line
    matrix = defaultdict(lambda: defaultdict(int))
    keys = contents.pop(0).split()
    for i in xrange(len(contents)):
        currentRow = contents[i].split()
        currentSym = currentRow.pop(0)
        for j in xrange(len(currentRow)):
            matrix[keys[j]][currentSym] = float(currentRow[j])
    print str(matrix)
    prob = hiddenPathProb(path,states,matrix)
    print prob

def hiddenPathProb(path,states,matrix):

    initialProb = float(1)/len(states)
    for i in xrange(len(path)-1):
        initialProb *= matrix[path[i]][path[i+1]]
    return initialProb

def HMMOutcomeWrapper(fileName):

    contents = open(fileName).readlines()
    x = contents.pop(0).strip()
    contents.pop(0) #pop dashed line
    alphabet = contents.pop(0).split()
    contents.pop(0) #pop dashed line
    path = contents.pop(0).strip()
    contents.pop(0) #pop dashed line
    states = contents.pop(0)
    contents.pop(0) #pop dashed line
    matrix = defaultdict(lambda: defaultdict(int))
    keys = contents.pop(0).split()
    for i in xrange(len(contents)):
        currentRow = contents[i].split()
        currentSym = currentRow.pop(0)
        for j in xrange(len(currentRow)):
            matrix[keys[j]][currentSym] = float(currentRow[j])
    print str(matrix)
    prob = HMMOutcome(x, path, matrix)
    print prob
    
def HMMOutcome(x, path, matrix):

    initialProb = 1.0
    for i in xrange(len(x)):
        initialProb *= matrix[x[i]][path[i]]
    return initialProb

def viterbiDecodingWrapper(fileName):

    contents = open(fileName).readlines()
    x = contents.pop(0).strip()
    contents.pop(0)
    alphabet = contents.pop(0).split()
    contents.pop(0)
    states = contents.pop(0).split()
    contents.pop(0)
    transitionMatrix = defaultdict(lambda: defaultdict(int))
    keys = contents.pop(0).split()
    for i in xrange(len(contents)):
        currentRow = contents.pop(0).split()
        if len(currentRow)==1:
            break
        currentSym = currentRow.pop(0)
        for j in xrange(len(currentRow)):
            transitionMatrix[keys[j]][currentSym] = float(currentRow[j])
    print str(transitionMatrix)
    emissionMatrix = defaultdict(lambda: defaultdict(int))
    keys2 = contents.pop(0).split()
    for i in xrange(len(contents)):
        currentRow = contents[i].split()
        currentSym = currentRow.pop(0)
        for j in xrange(len(currentRow)):
            emissionMatrix[currentSym][keys2[j]] = float(currentRow[j])
    print str(emissionMatrix)
    path = viterbiDecoding(x, states, transitionMatrix, emissionMatrix)
    print path

def viterbiDecoding(x, states, transitionMatrix, emissionMatrix):

    n = len(x)
    initialProb = 1.0
    pathDict = defaultdict(lambda: defaultdict(float))
    prevDict = defaultdict(lambda: defaultdict(str))
    source = 0
    for state in states:
        pathDict[source][state] = initialProb*emissionMatrix[state][x[0]]*len(states)
    for i in xrange(1,len(x)):
        for state in states:
            m = -sys.maxint-1
            current = ""
            for state2 in states:
                currentProb = pathDict[i-1][state2]*transitionMatrix[state][state2]*emissionMatrix[state][x[i]]
                if currentProb > m:
                    m = currentProb
                    current = state2
            pathDict[i][state]=m
            prevDict[i][state]=current
    m = -sys.maxint-1
    current = ""
    for state in states:
        if pathDict[len(x)-1][state]>m:
            m=pathDict[len(x)-1][state]
            current=state
    path = current
    for i in reversed(range(1,len(x))):
        current = prevDict[i][current]
        path += current
    return path[::-1]

def outcomeLikelihoodWrapper(fileName):

    contents = open(fileName).readlines()
    contents = open(fileName).readlines()
    x = contents.pop(0).strip()
    contents.pop(0)
    alphabet = contents.pop(0).split()
    contents.pop(0)
    states = contents.pop(0).split()
    contents.pop(0)
    transitionMatrix = defaultdict(lambda: defaultdict(int))
    keys = contents.pop(0).split()
    for i in xrange(len(contents)):
        currentRow = contents.pop(0).split()
        if len(currentRow)==1:
            break
        currentSym = currentRow.pop(0)
        for j in xrange(len(currentRow)):
            transitionMatrix[keys[j]][currentSym] = float(currentRow[j])
    print str(transitionMatrix)
    emissionMatrix = defaultdict(lambda: defaultdict(int))
    keys2 = contents.pop(0).split()
    for i in xrange(len(contents)):
        currentRow = contents[i].split()
        currentSym = currentRow.pop(0)
        for j in xrange(len(currentRow)):
            emissionMatrix[currentSym][keys2[j]] = float(currentRow[j])
    print str(emissionMatrix)
    prob = outcomeLikelihood(x,states,transitionMatrix,emissionMatrix)
    print prob

def outcomeLikelihood(x,states,transitionMatrix,emissionMatrix):

    initialProb = 1.0
    source = 0
    pathDict = defaultdict(lambda: defaultdict(float))
    for state in states:
        pathDict[state][source] = (initialProb/len(states)) * emissionMatrix[state][x[0]]
    for i in xrange(1,len(x)):
        for state in states:
            for state2 in states:
                    pathDict[state][i] += pathDict[state2][i-1]*transitionMatrix[state][state2]*emissionMatrix[state][x[i]]
    prob = 0.0
    for i in range(len(states)):
        prob += pathDict[states[i]][len(x)-1]
    return prob

def profileHMMWrapper(fileName):

    contents = open(fileName).readlines()
    threshold = float(contents.pop(0))
    contents.pop(0)
    alphabet = contents.pop(0).split()
    contents.pop(0)
    alignment = [str(line.strip()) for line in contents]
    print threshold
    print alphabet
    print alignment
    (transitionMatrix, emissionMatrix) = profileHMM(threshold,alphabet,alignment)

def profileHMM(threshold, alphabet, alignment):

    seed = seedAlignment(alignment,threshold)
    profile = formProfile(seed)
    transitionMatrix = defaultdict(lambda: defaultdict(float))
    emissionMatrix = defaultdict(lambda: defaultdict(float))
    numMatchStates = len(seed[0])
    numInsertStates = numMatchStates+2
    matchStates = []
    insertStates = []
    deleteStates = []
    #Set possible states
    for i in xrange(1,numMatchStates+1):
        matchStates.append("M"+str(i))
        deleteStates.append("D"+str(i))
        insertStates.append("I"+str(i))
    '''for i in xrange(numInsertStates):
        insertStates.append("I"+str(i))'''
    #Set initial emission values
    for i in alphabet:
        emissionMatrix["S"][i]=0
        emissionMatrix["E"][i]=0
        emissionMatrix["I0"][i]=0
    for i in xrange(len(matchStates)):
        for j in alphabet:
            emissionMatrix[matchStates[i]][j]=0
            emissionMatrix[insertStates[i]][j]=0
            emissionMatrix[deleteStates[i]][j]=0
    #Set initial transition values
    transitionMatrix["S"]["I0"]=0
    transitionMatrix["S"]["D1"]=0
    transitionMatrix["S"]["M1"]=0
    transitionMatrix["I0"]["I0"]=0
    transitionMatrix["I0"]["D1"]=0
    transitionMatrix["I0"]["M1"]=0
        
def seedAlignment(alignment,threshold):

    seed = copy.copy(alignment)
    n = len(alignment[0])
    for i in xrange(n):
        for j in alignment:
            count = 0 
            if j[i]=="-":
                count+=1
        if count/float(n) < threshold:
            removeI(seed,i)
    return seed

def removeI(seed,i):
    
    for j in xrange(len(seed)):
        seed[j] = seed[j][:i]+seed[j][(i+1):]
    return seed

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
            elif currentBase=="T":
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

def HMMParameterEstimationWrapper(fileName):

    contents = open(fileName).readlines()
    x = contents.pop(0).strip()
    contents.pop(0)
    alphabet = contents.pop(0).split()
    contents.pop(0)
    path = contents.pop(0).strip()
    contents.pop(0)
    states = contents.pop(0).split()
    (transitionMatrix, emissionMatrix) = HMMParameterEstimation(x,alphabet,path,states)
    print str(transitionMatrix)
    print str(emissionMatrix)
    transitionOutput = ""
    emissionOutput = ""
    #print output in matrix format
    transitionOutput += " " + " ".join(str(i) for i in states)
    for state in states:
        currentLine = "\n"
        currentLine += str(state)
        for state2 in states:
            currentLine += " " + str(transitionMatrix[state][state2])
        transitionOutput+=currentLine
    emissionOutput += " " + " ".join(str(i) for i in alphabet)
    for state in states:
        currentLine = "\n"
        currentLine += str(state)
        for letter in alphabet:
            currentLine += " " + str(emissionMatrix[state][letter])
        emissionOutput+=currentLine
    output = transitionOutput + "\n" + "--------" + "\n" + emissionOutput
    print output

def HMMParameterEstimation(x,alphabet,path,states):

    transitionMatrix = defaultdict(lambda: defaultdict(float))
    emissionMatrix = defaultdict(lambda: defaultdict(float))
    emissionMatrix = calcEmission(x,path,emissionMatrix,alphabet,states)
    transitionMatrix = calcTransition(x,path,transitionMatrix,states)
    return (transitionMatrix, emissionMatrix)

def calcEmission(x,path,emissionMatrix,alphabet,states):

    for state in states:
        for letter in alphabet:
            emissionMatrix[state][letter]=0.0
    for i in xrange(len(path)):
        emissionMatrix[path[i]][x[i]]+=1.0
    for state in emissionMatrix:
        currentSum = getSum(emissionMatrix[state])
        for state2 in emissionMatrix[state]:
            if currentSum!=0.0:
                emissionMatrix[state][state2]=float(emissionMatrix[state][state2]/currentSum)
            else:
                emissionMatrix[state][state2]=float(1.0/len(alphabet))
    return emissionMatrix

def calcTransition(x,path,transitionMatrix,states):

    for state in states:
        for state2 in states:
            transitionMatrix[state][state2]=0.0
    for i in xrange(len(path)-1):
        transitionMatrix[path[i]][path[i+1]]+=1.0
    for state in transitionMatrix:
        currentSum = getSum(transitionMatrix[state])
        for state2 in transitionMatrix[state]:
            if currentSum!=0.0:
                transitionMatrix[state][state2]=float(transitionMatrix[state][state2]/currentSum)
            else:
                transitionMatrix[state][state2]=float(1.0/len(states))
    return transitionMatrix

def getSum(d):

    return sum(d.values())

def viterbiLearningWrapper(fileName):

    contents = open(fileName).readlines()
    iterations = int(contents.pop(0).strip())
    contents.pop(0)
    x = contents.pop(0).strip()
    contents.pop(0)
    alphabet = contents.pop(0).split()
    contents.pop(0)
    states = contents.pop(0).split()
    contents.pop(0)
    transitionMatrix = defaultdict(lambda: defaultdict(int))
    keys = contents.pop(0).split()
    for i in xrange(len(contents)):
        currentRow = contents.pop(0).split()
        if len(currentRow)==1:
            break
        currentSym = currentRow.pop(0)
        for j in xrange(len(currentRow)):
            transitionMatrix[keys[j]][currentSym] = float(currentRow[j])
    emissionMatrix = defaultdict(lambda: defaultdict(int))
    keys2 = contents.pop(0).split()
    for i in xrange(len(contents)):
        currentRow = contents[i].split()
        currentSym = currentRow.pop(0)
        for j in xrange(len(currentRow)):
            emissionMatrix[currentSym][keys2[j]] = float(currentRow[j])
    (transitionMatrix, emissionMatrix) = viterbiLearning(iterations, x, alphabet, states, transitionMatrix, emissionMatrix)
    transitionOutput = ""
    emissionOutput = ""
    #print output in matrix format
    transitionOutput += " " + " ".join(str(i) for i in states)
    for state in states:
        currentLine = "\n"
        currentLine += str(state)
        for state2 in states:
            currentLine += " " + str(transitionMatrix[state][state2])
        transitionOutput+=currentLine
    emissionOutput += " " + " ".join(str(i) for i in alphabet)
    for state in states:
        currentLine = "\n"
        currentLine += str(state)
        for letter in alphabet:
            currentLine += " " + str(emissionMatrix[state][letter])
        emissionOutput+=currentLine
    output = transitionOutput + "\n" + "--------" + "\n" + emissionOutput
    print output
    

def viterbiLearning(iterations, x, alphabet, states, transitionMatrix, emissionMatrix):

    for i in xrange(iterations):
        optimalPath = viterbiDecoding(x, states, transitionMatrix, emissionMatrix)
        (transitionMatrix,emissionMatrix) = HMMParameterEstimation(x,alphabet,optimalPath,states)
    return (transitionMatrix, emissionMatrix)

def softDecodingWrapper(fileName):

    contents = open(fileName).readlines()
    x = contents.pop(0).strip()
    contents.pop(0)
    alphabet = contents.pop(0).split()
    contents.pop(0)
    states = contents.pop(0).split()
    contents.pop(0)
    transitionMatrix = defaultdict(lambda: defaultdict(int))
    keys = contents.pop(0).split()
    for i in xrange(len(contents)):
        currentRow = contents.pop(0).split()
        if len(currentRow)==1:
            break
        currentSym = currentRow.pop(0)
        for j in xrange(len(currentRow)):
            transitionMatrix[keys[j]][currentSym] = float(currentRow[j])
    #print str(transitionMatrix)
    emissionMatrix = defaultdict(lambda: defaultdict(int))
    keys2 = contents.pop(0).split()
    for i in xrange(len(contents)):
        currentRow = contents[i].split()
        currentSym = currentRow.pop(0)
        for j in xrange(len(currentRow)):
            emissionMatrix[currentSym][keys2[j]] = float(currentRow[j])
    #print str(emissionMatrix)
    finalDict = softDecoding(x,alphabet,states,transitionMatrix,emissionMatrix)
    print finalDict

def softDecoding(x,alphabet,states,transitionMatrix,emissionMatrix):

    frontDict = defaultdict(lambda: defaultdict(float))
    backDict = defaultdict(lambda: defaultdict(float))
    frontDict = calcFrontDict(states, x, emissionMatrix, transitionMatrix, frontDict)
    backDict = calcBackDict(states, x, emissionMatrix, transitionMatrix, backDict)
    print str(frontDict)
    print "\n"
    print str(backDict)
    print "\n"
    finalDict = defaultdict(lambda: defaultdict(float))
    for i in xrange(len(x)):
        currentSum = 0.0
        for state in states:
            current = frontDict[i][state]*backDict[i][state]
            finalDict[i][state] = current
            currentSum += current
        for state in states:
            finalDict[i][state] = finalDict[i][state]/currentSum
    return finalDict

def calcFrontDict(states, x, emissionMatrix, transitionMatrix, frontDict):

    initialProb = 1.0
    for state in states:
        frontDict[0][state] = emissionMatrix[state][x[0]]*(initialProb/len(states))
    for i in xrange(1,len(x)):
        for state in states:
            for state2 in states:
                frontDict[i][state]+=frontDict[i-1][state2]*transitionMatrix[state2][state]*emissionMatrix[state][x[i]]
    return frontDict

def calcBackDict(states, x, emissionMatrix, transitionMatrix, backDict):

    initialProb = 1.0
    #back = x[::-1]
    for state in states:
        backDict[len(x)-1][state] = 1.0
    for i in xrange(len(x)-2,-1,-1):
        for state in states:
            for state2 in states:
                backDict[i][state] += backDict[i+1][state2]*transitionMatrix[state2][state]*emissionMatrix[state][x[i+1]]

    return backDict

















