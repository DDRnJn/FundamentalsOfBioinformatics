#Bioinformatics hw9
#Dhruv Ranjan

from collections import defaultdict
import sys
import copy

def neighborJoiningWrapper(fileName):

    contents = open(fileName).readlines()
    n = int(contents.pop(0).strip())
    distanceDict = defaultdict(lambda: defaultdict(int))
    for j in range(len(contents)):
        current = map(int, contents[j].split())
        for i in range(len(current)):
            distanceDict[i][j]=current[i]
    finalDict = {}
    treeDict = neighborJoining(n, distanceDict)
    fout = open("texts/hw9/neighborJoiningAnswer.txt","wt")
    output = ""
    for key in treeDict:
        for node in treeDict[key]:
            if key != node: 
                output += str(key) + "->" + str(node) + ":" + str(treeDict[key][node]) + "\n"
    fout.write(output)

def totaldistance(distanceDict, node):

    current = distanceDict[node]
    sum = 0
    for key in current:
        sum += current[key]
    return sum

def makeNJM(distanceDict, n):

    NJM = defaultdict(lambda: defaultdict(int))
    for i in distanceDict:
        for j in distanceDict:
            if i != j:
                NJM[i][j] = (n-2) * distanceDict[i][j] - totaldistance(distanceDict, i) - totaldistance(distanceDict, j)
            else:
                NJM[i][j] = 0
    return NJM

def minNonDiagonal(NJM):

    minElem = sys.maxint
    (minI, minJ) = (-1,-1)
    for i in NJM:
        for j in NJM[i]:
            if i != j:
                current = NJM[i][j]
                if current < minElem:
                    minElem = current
                    minI = i
                    minJ = j
    return (minI, minJ)

def neighborJoining(n, distanceDict):

    if n == 2:
        return distanceDict
    NJM = makeNJM(distanceDict, n)
    (minI, minJ) = minNonDiagonal(NJM)
    delta = (totaldistance(distanceDict, minI) - totaldistance(distanceDict, minJ))/(float(n-2))
    limbLengthI = (0.5)*(distanceDict[minI][minJ]+delta)
    limbLengthJ = (0.5)*(distanceDict[minI][minJ]-delta)
    newNode = max(distanceDict.keys(), key=int) + 1
    distanceDict[newNode][newNode] = 0
    for k in distanceDict:
        distanceDict[k][newNode] = (0.5) * (distanceDict[k][minI] + distanceDict[k][minJ] - distanceDict[minI][minJ])
        distanceDict[newNode][k] = (0.5) * (distanceDict[k][minI] + distanceDict[k][minJ] - distanceDict[minI][minJ])
    distanceDict[newNode][newNode] = 0
    distanceDict.pop(minI)
    distanceDict.pop(minJ)
    for key in distanceDict:
        distanceDict[key].pop(minI)
        distanceDict[key].pop(minJ)
    T = neighborJoining(n-1, distanceDict)
    T[newNode][minI] = limbLengthI
    T[newNode][minJ] = limbLengthJ
    T[minI][newNode] = limbLengthI
    T[minJ][newNode] = limbLengthJ
    return T
    
def hammingDistance(s1,s2):

    if len(s1)==len(s2)==0:
        return 0
    else:
        if s1[0]==s2[0]:
            return hammingDistance(s1[1:],s2[1:])
        else:
            return hammingDistance(s1[1:],s2[1:]) + 1

def getI(leaves,i):

    chars = []
    for k in leaves: 
        chars.append(k[i])
    return chars
        
def smallParsimonyWrapper(fileName):

    contents = open(fileName).readlines()
    n = int(contents.pop(0))
    treeDict = defaultdict(list)
    leaves = []
    for line in contents:
        current = line.split("->")
        key = current[0].strip()
        to = current[1].strip()
        if not to.isdigit():
            leaves.append(to)
        treeDict[key].append(to)
    outputDict = defaultdict(str)
    resultDict = defaultdict(str)
    for i in range(n-1):
        outputDict[i]=""
        outputDict[n+i]=""
    outputDict[i+1]=""
    resultDict = copy.copy(outputDict)
    leafList = range(0,n)
    for i in range(((len(treeDict))/2)+1):
        if len(treeDict[str(n+i)][0])>1:
            treeDict[str(n+i)][0] = leafList.pop(0)
            treeDict[str(n+i)][1] = leafList.pop(0)
    parsimonyScore = 0
    nLength = len(leaves[0])
    for i in range(nLength):
        ithChars = getI(leaves,i)
        for j in range(len(leaves)):
            outputDict[j]=ithChars[j]
        (score, labelling, scoreDict) = smallParsimony(outputDict,treeDict,n)
        parsimonyScore += score
        resultDict = populateResults(resultDict, labelling)
    mappingDict = computeMappings(resultDict, treeDict, n)
    fout = open("texts/hw9/smallParsimonyAnswer.txt","wt")
    output = ""
    for key in mappingDict:
        for node in mappingDict[key]: 
            output += str(key) + "->" + str(node) + ":" + str(mappingDict[key][node]) + "\n"
    fout.write(output)
    return parsimonyScore

def computeMappings(resultDict, treeDict, n):
    mappingDict = defaultdict(lambda: defaultdict(int))
    for key in treeDict:
        if int(key) >= n:
            children = treeDict[key]
            daughter = int(children[0])
            son = int(children[1])
            keySeq = resultDict[int(key)]
            daughterSeq = resultDict[daughter]
            sonSeq = resultDict[son]
            mappingDict[keySeq][daughterSeq] = hammingDistance(keySeq, daughterSeq)
            mappingDict[keySeq][sonSeq] = hammingDistance(keySeq, sonSeq)
            mappingDict[daughterSeq][keySeq] = hammingDistance(keySeq, daughterSeq)
            mappingDict[sonSeq][keySeq] = hammingDistance(keySeq, sonSeq)
    return mappingDict

def populateResults(resultDict, labelling):

    for key in resultDict:
        resultDict[key] += labelling[key]
    return resultDict

def smallParsimony(T,treeDict,n):

    alphabet = ["A", "T", "C", "G"]
    tagDict = defaultdict(int)
    scoreDict = defaultdict(lambda: defaultdict(int))
    for v in T:
        tagDict[v] = 0
        if v < n: #isLeaf(v,n):
            tagDict[v] = 1
            for k in alphabet:
                if T[v]==k:
                    scoreDict[v][k] = 0
                else:
                    scoreDict[v][k] = sys.maxint
    ripeNodes = findRipeNodes(T, tagDict, treeDict, n)
    while len(ripeNodes) != 0:
        v = ripeNodes.pop(0)
        tagDict[v] = 1
        daughter = int(treeDict[str(v)][0])
        son = int(treeDict[str(v)][1])
        dMin = findMins(daughter, scoreDict) 
        sMin = findMins(son, scoreDict)
        dMinScore = scoreDict[daughter].get(dMin[0])
        sMinScore = scoreDict[son].get(sMin[0])
        for k in alphabet:
            deltaik = 1
            deltajk = 1
            if k in dMin:
                deltaik = 0
            if k in sMin:
                deltajk = 0
            scoreDict[v][k] = dMinScore + deltaik + sMinScore + deltajk
        ripeNodes = findRipeNodes(T, tagDict, treeDict, n)
    rootMin = min(scoreDict[v],key=scoreDict[v].get)
    rootMinScore = scoreDict[v][rootMin]
    T[v] = rootMin
    T = findInternals(treeDict, v, scoreDict, T, n)
    return (rootMinScore, T, scoreDict)

def findInternals(treeDict, root, scoreDict, T, n):

    if not(int(root) < n):
        daughter = int(treeDict[str(root)][0])
        son = int(treeDict[str(root)][1])
        dMin = findMins(daughter, scoreDict) 
        sMin = findMins(son, scoreDict)
        dMinScore = scoreDict[daughter].get(dMin[0])
        sMinScore = scoreDict[son].get(sMin[0])
        dLabel = ""
        sLabel = ""
        if len(dMin)==1:
            dLabel=dMin[0]
        else:
            dLabel = minHamming(dMin, T[root])
        if len(sMin)==1:
            sLabel=sMin[0]
        else:
            sLabel = minHamming(sMin, T[root])
        T[daughter]=dLabel
        T[son]=sLabel
        findInternals(treeDict, daughter, scoreDict, T, n)
        findInternals(treeDict, son, scoreDict, T, n)
    return T
        
def minHamming(dsMin, root):

    res = ""
    for node in dsMin:
        if node==root:
            res=node
    if res == "":
        res = dsMin[0]
    return res

def findMins(node, scoreDict):

    mins = []
    currentMin = sys.maxint
    for key in scoreDict[node]:
        if scoreDict[node][key] <= currentMin:
            currentMin = scoreDict[node][key]
    for key in scoreDict[node]:
        if scoreDict[node][key]==currentMin:
            mins.append(key)
    return mins

def findRipeNodes(T, tagDict, treeDict, n):

    ripeNodes = []
    for key in T:
        if tagDict[int(key)]==0:
            children = treeDict[str(key)]
            daughter = int(children[0])
            son = int(children[1])
            if tagDict[daughter] == tagDict[son] == 1:
                ripeNodes.append(int(key))
    return ripeNodes
