#Bioinformatics hw7
#Dhruv Ranjan

import copy
import math

def greedySorting(P):

    approxReversalDistance = 0
    permList = []
    for k in range(1,len(P)+1):
        if int(P[k-1][1:]) != int(k):
            newP = kSort(k, P)
            approxReversalDistance += 1
            permList += [copy.copy(newP)]
        if P[k-1][0] == "-":
            P[k-1] = "+" + P[k-1][1:]
            approxReversalDistance += 1
            permList += [copy.copy(newP)]
        P = newP
    return (approxReversalDistance, permList)

def kSort(k, P):

    if "+" + str(k) in P:
        currentIndex = P.index("+" + str(k))
    else:
        currentIndex = P.index("-" + str(k))
    revPart = P[(k-1):currentIndex+1]
    revPart.reverse()
    for i in range(len(revPart)):
        if revPart[i][0]=="+":
           revPart[i] = "-" + revPart[i][1:]
        else:
            revPart[i] = "+" + revPart[i][1:]
    P[(k-1):currentIndex+1] = revPart
    final = P 
    #print "final: " + str(final)
    return final

def greedySortingWrapper(fileName):

    contents = open(fileName).readlines()
    pList = [str(i) for i in contents[0].strip("{").strip("}").split()]
    pList[0] = pList[0][1:]
    #removes the parens from the first and last element
    pList[len(pList)-1] = pList[len(pList)-1][:len(pList[len(pList)-1])-1]
    print(pList)
    (approxReversalDistance, permList) = greedySorting(pList)
    output = []
    for i in permList:
        output += ["(" + " ".join(str(j) for j in i) + ")"]
    fout = open("texts\hw7\greedySortingAnswer.txt", "wt")
    fout.write("\n".join(str(i) for i in output))

def breakpointCount(permutation):

    n = len(permutation)
    count = 0
    for i in range(n-1):
        if permutation[i+1] != permutation[i] + 1:
            count += 1
    return count

def breakpointCountWrapper(fileName):

    contents = open(fileName).readlines()
    pList = [str(i) for i in contents[0].strip("{").strip("}").split()]
    pList[0] = pList[0][1:]
    #removes the parens from the first and last element
    pList[len(pList)-1] = pList[len(pList)-1][:len(pList[len(pList)-1])-1]
    pList2 = []
    for i in pList:
        if i[0]=="+":
            pList2 += [int(i[1:])]
        else:
            pList2 += [int(i[1:])*-1]
    #print pList2
    count = breakpointCount(pList2)
    print(count)

#for the sample data set, works for range(597), one error with range(598)
#found the issue. Things work now :) 
def checkAns():

    contents1 = open("texts/hw7/greedySortingAnswer.txt").readlines()
    contents2 = open("texts/hw7/ans2.txt").readlines()
    wrongList = []
    wrongList2 = []
    print(len(contents1))
    print(len(contents2))
    for i in range(598):
        if contents2[i] != contents1[i]:
            wrongList += [contents1[i]]
            wrongList2 += [contents2[i]]
    print(len(wrongList))
    print(contents1[596])
    print(wrongList)
    print(contents2[596])
    print(wrongList2)

def chromosomeToCycle(chromosome):

    cycle = []
    for i in range(2*len(chromosome)):
        cycle += [0]
    for j in range(0,len(chromosome)):
        current = chromosome[j]
        if current>0:
            cycle[(2*j)] = 2*(current)-1
            cycle[(2*j)+1] = 2*current
        else:
            cycle[(2*j)] = -2*current
            cycle[(2*j)+1] = -2*current -1
    return cycle

def chromosomeToCycleWrapper(fileName):

    contents = open(fileName).readlines()
    pList = [str(i) for i in contents[0].strip("{").strip("}").split()]
    pList[0] = pList[0][1:]
    #removes the parens from the first and last element
    pList[len(pList)-1] = pList[len(pList)-1][:len(pList[len(pList)-1])-1]
    pList2 = []
    for i in pList:
        if i[0]=="+":
            pList2 += [int(i[1:])]
        else:
            pList2 += [int(i[1:])*-1]
    print(pList2)
    cycle = chromosomeToCycle(pList2)
    output = "(" + " ".join(str(j) for j in cycle) + ")"
    print(output)

def cycleToChromosome(cycle):

    chromosome = []
    for i in range(len(cycle)/2):
        chromosome += [0]
    for j in range(len(cycle)/2):
        if cycle[2*j] < cycle[(2*j)+1]:
            chromosome[j] = math.ceil(cycle[2*j]/2.0)
        else:
            chromosome[j] = -1*math.ceil(cycle[2*j]/2.0)
    return chromosome

def cycleToChromosomeWrapper(fileName):

    contents = open(fileName).readlines()
    pList = [str(i) for i in contents[0].split()]
    pList[0] = pList[0][1:]
    pList[len(pList)-1] = pList[len(pList)-1][:len(pList[len(pList)-1])-1]
    #contents is in form (1 2 4 3)
    #pList is in form [1, 2, 4, 3]
    pList = map(int, pList)
    print(pList)
    chromosome = cycleToChromosome(pList)
    output = map(int, chromosome)
    output = map(str, output)
    print(output)
    for i in range(len(output)):
        if output[i][0] != "-":
            output[i] = "+" + output[i]
    output = "(" + " ".join(str(j) for j in output) + ")"
    print(output)

def coloredEdges(P):

    edges = []
    for chromosome in P:
        nodes = chromosomeToCycle(chromosome)
        for j in range(len(chromosome)):
            edges += [(nodes[(2*j)-1],nodes[(2*j)])]
    return edges

def coloredEdgesWrapper(fileName):

    contents = open(fileName).read()
    tupleList = contents.split(")")
    tupleList.pop(len(tupleList)-1)
    for i in range(len(tupleList)):
        tupleList[i] = tupleList[i][1:]
    tupleList2 = []
    for i in range(len(tupleList)):
        current = [str(k) for k in tupleList[i].split()]
        toAdd = []
        for j in current:
            if j[0]=="+":
                toAdd += [int(j[1:])]
            else:
                toAdd += [int(j[1:])*-1]
        tupleList2 += [toAdd]
    edges = coloredEdges(tupleList2)
    output = " ".join(str(j) + "," for j in edges)
    print(output[:len(output)-1])
            
def graphToGenome(genomeGraph):

    P = []
    for nodes in genomeGraph:
        print(nodes)
        chromosome = cycleToChromosome(nodes)
        P += chromosome
    return P
        

def graphToGenomeWrapper(fileName):

    contents = open(fileName).read()
    graph1 = contents.split("),")
    graph2 = []
    graph3 = []
    for i in range(len(graph1)-1):
        graph1[i] = graph1[i].strip() + ")"
    graph1[-1] = graph1[-1].strip()
    for i in range(len(graph1)):
        graph1[i] = graph1[i].replace(",", "")
    for i in range(len(graph1)):
        graph2 += [str(j) for j in graph1[i].split()]
        graph2[0] = graph2[0][1:]
        graph2[-1] = graph2[-1][:-1]
        graph2 = map(int, graph2)
        graph3 += [graph2]
        graph2 = []
    print(graph3)
    P = graphToGenome(graph3)
    print(P)

def twoBreakDistance(P,Q):

    breakGraph = {}
    blockNum = 0
    cycleNum = 0
    for i in P:
        blockNum += len(i)
    totalGenome = P+Q
    colEdges = coloredEdges(totalGenome)
    for i in range(len(colEdges)):
        breakGraph[i+1] = []
    for i in range(len(colEdges)):
        currentPair = colEdges[i]
        start = currentPair[0]
        end = currentPair[1]
        breakGraph[start] = breakGraph[start] + [end]
        breakGraph[end] = breakGraph[end] + [start]
    for i in range(len(colEdges)):
        breakGraph[i+1] = list(set(breakGraph[i+1]))
    allNodes = [(i+1) for i in range(len(colEdges))]
    seen = []
    while len(seen) < len(colEdges):
        cycleNum += 1
        currentNodes = [allNodes.pop(0)]
        seen += [currentNodes[0]]
        while len(currentNodes) > 0:
            currentNode = currentNodes.pop(0)
            neighbors = breakGraph.get(currentNode, [])
            for i in neighbors:
                if i in allNodes:
                    currentNodes += [i]
            allNodes2 = [j for j in allNodes if j not in currentNodes]
            allNodes = allNodes2
            seen += currentNodes
            seen = list(set(seen))
    return blockNum - cycleNum
    
def twoBreakDistanceWrapper(fileName):

    contents = open(fileName).readlines()
    #print contents
    P = contents[0]
    Q = contents[1]
    P = P.strip()
    Q = Q.strip()
    #print P
    #print Q
    P = P.strip('(').strip(')')
    Q = Q.strip('(').strip(')')
    #print P
    #print Q
    P = P.split(')(')
    Q = Q.split(')(')
    #print P
    #print Q
    P = [e.split() for e in P]
    Q = [e.split() for e in Q]
    #print P
    #print Q
    P = [map(int, i) for i in P]
    Q = [map(int, i) for i in Q]
    #print "P: " + str(P)
    #print "Q: " + str(Q)
    twoBD = twoBreakDistance(P,Q)
    print(twoBD)
    
#find the reverse complement of the str pattern
def reverseComplement(pattern):

    newPattern = ""
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
    return newPattern
    #return "".join(newPattern)

def patternGen(text, k):
    
    patterns = []
    for i in range(0, len(text)-k+1):
        patterns += [text[i:(k+i)]]
    return patterns

'''def sharedKmers(k, str1, str2):

    kmers = []
    str1Kmers = patternGen(str1, k)
    str2Kmers = patternGen(str2, k)
    for i in range(len(str1Kmers)):
        for j in range(len(str2Kmers)):
            if str1Kmers[i]==str2Kmers[j] or reverseComplement(str1Kmers[i])==str2Kmers[j]:
                kmers += [(i,j)]
    return kmers'''

def sharedKmers(k, str1, str2):

    kmers = []
    str1Kmers = patternGen(str1, k)
    str2Kmers = patternGen(str2, k)
    kmerDict = {}
    for i in range(len(str1Kmers)):
        kmerDict[str1Kmers[i]] = []
    for i in range(len(str1Kmers)):
        kmerDict[str1Kmers[i]] = kmerDict[str1Kmers[i]] + [i]
    for i in range(len(str2Kmers)):
        if str2Kmers[i] in kmerDict:
            newKmers = [(j,i) for j in kmerDict[str2Kmers[i]]]
            kmers += newKmers
        elif reverseComplement(str2Kmers[i]) in kmerDict:
            newKmers = [(j,i) for j in kmerDict[reverseComplement(str2Kmers[i])]]
            kmers += newKmers
    return kmers

def sharedKmersWrapper(fileName):

    contents = open(fileName).readlines()
    k = int(contents[0])
    str1 = contents[1].strip()
    str2 = contents[2].strip()
    kmerList = sharedKmers(k, str1, str2)
    fout = open("texts/hw7/sharedKmersAnswer.txt", "wt")
    #print "\n".join(str(i) for i in kmerList)
    fout.write("\n".join(str(i) for i in kmerList))
