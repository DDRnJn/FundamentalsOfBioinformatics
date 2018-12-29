#Bioinformatics HW3
#Dhruv Ranjan

import random
import copy

def patternGen(text, k):
    
    patterns = []
    for i in xrange(0, len(text)-k+1):
        patterns += [text[i:(k+i)]]
    return patterns

def patternGenWrapper(fileName, fileName2):

    contents = open(fileName).readlines()
    k = contents[0]
    contents.pop(0)
    seq = contents[0]
    motifs = patternGen(str(seq), int(k))
    fout = open(fileName2, "wt")
    for m in motifs:
        fout.write("%s\n" % m)
        
def stringFromGenomePath(kmers):

    finalString = ""
    initialString = kmers[0]
    finalString = initialString
    kmers.pop(0)
    for kmer in kmers:
        finalString += kmer[-1:]
    return finalString

def stringFromGenomePathWrapper(fileName):

    contents = [line.rstrip("\n") for line in open(fileName)]
    seq = stringFromGenomePath(contents)
    fout = open ("texts\hw3q2Answer.txt", "wt")
    fout.write(seq)
    return seq

def makeOverlapGraph(kmers):

    graph = {}
    for kmer in kmers:
        for kmer2 in kmers:
            if kmer != kmer2:
                if kmer[-(len(kmer)-1):] == kmer2[:len(kmer2)-1]:
                    graph[kmer] = graph.get(kmer,[]) + [kmer2]
    #for key in graph:
        #print key + " -> " + ", ".join(graph.get(key))
    return graph
        

def makeOverlapGraphWrapper(fileName):

    contents = [line.rstrip("\n") for line in open(fileName)]
    graph = makeOverlapGraph(contents)
    fout = open("texts\hw3q3Answer.txt", "wt")
    for key in graph:
        fout.write(key + " -> " + ", ".join(graph.get(key)))
        fout.write("\n")
    
def makeDeBruijnGraph(k, text):

    kmers = patternGen(text,k)
    deBruijnGraph = makePrelimGraph(kmers)
    return deBruijnGraph
    #for key in deBruijnGraph:
        #print key + " -> " + ", ".join(deBruijnGraph.get(key))

def makePrelimGraph(kmers):

    graph = {}
    newKmers = []
    for kmer in kmers:
        graph[kmer[:-1]] = graph.get(kmer[:-1],[]) + [kmer[1:]]
    return graph

def makeDeBruijnGraphWrapper(fileName):

    contents = open(fileName).readlines()
    k = int(contents[0])
    contents.pop(0)
    text = contents[0]
    graph = makeDeBruijnGraph(k, text)
    fout = open("texts\hw3q4Answer.txt", "wt")
    for key in sorted(graph):
        fout.write(key + " -> " + ",".join(graph.get(key)))
        fout.write("\n")
                    
def deBruinFromKmers(kmers):

    deBruijnGraph = makePrelimGraph(kmers)
    return deBruijnGraph
        
def deBruinFromKmersWrapper(fileName):

    contents = [line.rstrip("\n") for line in open(fileName)]
    graph = deBruinFromKmers(contents)
    fout = open("texts\hw3q5Answer.txt", "wt")
    for key in graph:
        fout.write(key + " -> " + ",".join(graph.get(key)))
        fout.write("\n")

def eulerianCycle(graph):

    randomNode = random.choice(graph.keys())
    #(randomNode,endNode) = chooseNode(graph)
    (initialCycle, newGraph) = randomWalkNoRepeats(graph, randomNode)
    mergedCycle = initialCycle
    while len(newGraph) != 0:
        while True:
            randomNode = random.choice(initialCycle)
            if randomNode in newGraph:
                currentNode = randomNode
                break
        (newCycle, newGraph) = randomWalkNoRepeats(newGraph, currentNode)
        mergedCycle = mergeCycles(initialCycle, newCycle)
        initialCycle = mergedCycle
    return mergedCycle

def mergeCycles(initialCycle, newCycle):

    startingNode = newCycle[0]
    mergedCycle = []
    newInitial = copy.copy(initialCycle)
    appendToEnd = []
    for i in xrange(len(initialCycle)):
        if initialCycle[i] == startingNode:
            break
        appendToEnd += [initialCycle[i]]
        newInitial.pop(0)
    newInitial.pop(len(newInitial)-1)
    mergedCycle = newInitial + appendToEnd + newCycle
    return mergedCycle

def randomWalkNoRepeats(graph, randomNode):
    
    cycle = []
    currentNode = copy.copy(randomNode)
    while True:
        cycle += [currentNode]
        print currentNode
        print graph.get(currentNode)
        randomIndex = random.randint(0,len(graph.get(currentNode))-1)
        randomNeighbor = graph.get(currentNode).pop(randomIndex)
        if graph.get(currentNode)==[]:
            del graph[currentNode]
        currentNode = randomNeighbor
        if currentNode == randomNode:
            cycle += [currentNode]
            break
    return (cycle, graph)

def randomWalkNoRepeatsWrapper(fileName):
    
    contents = open(fileName).readlines()
    graph = {}
    for line in contents:
        (key,arrow,value) = line.partition(" -> ")
        key = key.rstrip()
        value = value.rstrip()
        currentValues = []
        currentString = ""
        for s in value:
            if s != ",":
                currentString += s
            else:
                currentValues += currentString
                currentString = ""
        currentValues += [currentString]
        graph[key] = graph.get(key,[]) + currentValues
    cycle = randomWalkNoRepeats(graph,random.choice(graph.keys()))
    print "->".join(cycle)

#This nonsense doesn't work.
def eulerianCycle2(graph):

    seen = []
    cycle = []
    randomNode = random.choice(graph.keys())
    currentNode = randomNode
    while True:
        if graph.get(currentNode)==[]:
            cycle += [currentNode]
            if len(seen) != 0:
                currentNode = seen.pop(0)
            else:
                break
        else:
            seen.insert(0,currentNode)
            randomIndex = random.randint(0,len(graph.get(currentNode))-1)
            randomNeighbor = graph.get(currentNode).pop(randomIndex)
            currentNode= randomNeighbor
    print "cycle: " + str(cycle)
    return cycle

def chooseNode(graph):

    degreeGraph = {}
    for key in graph:
        for key2 in graph.get(key):
            degreeGraph[key] = (degreeGraph.get(key,(0,0))[0] + 1, degreeGraph.get(key,(0,0))[1])
            degreeGraph[key2] = (degreeGraph.get(key2,(0,0))[0],degreeGraph.get(key2,(0,0))[1] + 1)
    startNode = ""
    endNode = ""
    for key in degreeGraph:
        if degreeGraph.get(key,0)[0] > degreeGraph.get(key,0)[1]:
            startNode = key
        if degreeGraph.get(key,0)[0] < degreeGraph.get(key,0)[1]:
            endNode = key
    #if startNode != "":
        #node = oddNode
    #else:
        #node = random.choice(graph.keys())
    return (startNode, endNode)

def chooseNodeWrapper(fileName):

    contents = open(fileName).readlines()
    graph = {}
    for line in contents:
        (key,arrow,value) = line.partition(" -> ")
        key = key.rstrip()
        value = value.rstrip()
        currentValues = []
        currentString = ""
        for s in value:
            if s != ",":
                currentString += s
            else:
                currentValues += [currentString]
                currentString = ""
        currentValues += [currentString]
        graph[key] = graph.get(key,[]) + currentValues
    (startNode, endNode) = chooseNode(graph)
    print startNode, endNode

def eulerianCycleWrapper(fileName):

    contents = open(fileName).readlines()
    graph = {}
    for line in contents:
        (key,arrow,value) = line.partition(" -> ")
        key = key.rstrip()
        value = value.rstrip()
        currentValues = []
        currentString = ""
        for s in value:
            if s != ",":
                currentString += s
            else:
                currentValues += [currentString]
                currentString = ""
        currentValues += [currentString]
        graph[key] = graph.get(key,[]) + currentValues
    cycle = eulerianCycle(graph)
    fout = open("texts\hw3q6Answer.txt", "wt")
    fout.write("->".join(cycle))

def eulerianPath(graph):

    (startNode, endNode) = chooseNode(graph)
    print startNode, endNode
    graph[endNode] = graph.get(endNode,[]) + [startNode]
    print str(graph[endNode])
    cycle = eulerianCyclePath(graph, startNode, endNode)
    return cycle

def eulerianCyclePath(graph, startNode, endNode):

    randomNode = startNode
    (initialCycle, newGraph) = randomWalkNoRepeatsPath(graph, randomNode, endNode)
    mergedCycle = initialCycle
    while len(newGraph) != 0:
        while True:
            randomNode = random.choice(initialCycle)
            if randomNode in newGraph:
                currentNode = randomNode
                break
        (newCycle, newGraph) = randomWalkNoRepeatsPath(newGraph, currentNode, endNode)
        mergedCycle = mergeCycles(initialCycle, newCycle)
        initialCycle = mergedCycle
    return mergedCycle

def randomWalkNoRepeatsPath(graph, randomNode, endNode):
    
    cycle = []
    currentNode = copy.copy(randomNode)
    finalTrigger = 0
    while True:
        cycle += [currentNode]
        print currentNode
        print graph.get(currentNode)
        print graph.get(endNode)
        while True: #make sure endNode is only used up in the final case
            randomIndex = random.randint(0,len(graph.get(currentNode))-1)
            if graph.get(currentNode)[randomIndex] != endNode:
                randomNeighbor = graph.get(currentNode).pop(randomIndex)
                break
            elif graph.get(currentNode)[randomIndex] == endNode and len(graph.get(endNode))>1:
                randomNeighbor = graph.get(currentNode).pop(0)
                break
            elif graph.get(graph.get(currentNode)[randomIndex]) == [randomNode] and len(graph)==1:
                finalTrigger = 1
                break
        if finalTrigger != 1:
            if graph.get(currentNode)==[]:
                del graph[currentNode]
            currentNode = randomNeighbor
            if currentNode == randomNode:
                cycle += [currentNode]
                break
        else:
            graph = {}
            cycle += randomNeighbor
            break
    return (cycle, graph)

"""def eulerianCycleWithStart(graph, startNode, endNode):

    randomNode = startNode
    (initialCycle, newGraph) = randomWalkNoRepeatsPath(graph, randomNode, endNode)
    mergedCycle = initialCycle
    while len(newGraph) != 0:
        while True:
            randomNode = random.choice(initialCycle)
            if randomNode in newGraph:
                currentNode = randomNode
                break
        (newCycle, newGraph) = randomWalkNoRepeatsPath(newGraph, currentNode, endNode)
        mergedCycle = mergeCyclesPath(initialCycle, newCycle)
        initialCycle = mergedCycle
    return mergedCycle

def randomWalkNoRepeatsPath(graph, randomNode, endNode):
    
    cycle = []
    currentNode = copy.copy(randomNode)
    breakTrigger = 0
    while True:
        cycle += [currentNode]
        while True:
            randomIndex = random.randint(0,len(graph.get(currentNode))-1)
            if graph.get(currentNode)[randomIndex] != endNode:
                break
            elif len(graph.get(endNode,[])) != 0:
                break
            elif graph.get(endNode,[]) == [] and len(graph.get(currentNode))<=1:
                cycle += [endNode]
                breakTrigger = 1
                break
        if breakTrigger == 1:
            del graph[currentNode]
            break
        randomNeighbor = graph.get(currentNode).pop(randomIndex)
        if graph.get(currentNode)==[]:
            del graph[currentNode]
        currentNode = randomNeighbor
        if currentNode == randomNode:
            cycle += [currentNode]
            break
    return (cycle, graph)

def mergeCyclesPath(initialCycle, newCycle):

    startingNode = newCycle[0]
    mergedCycle = []
    newInitial = copy.copy(initialCycle)
    appendToEnd = []
    for i in xrange(len(initialCycle)):
        if initialCycle[i] == startingNode:
            break
        appendToEnd += [initialCycle[i]]
        newInitial.pop(0)
    newInitial.pop(len(newInitial)-1)
    mergedCycle = newInitial + appendToEnd + newCycle
    return mergedCycle"""

def eulerianPathWrapper(fileName):

    contents = open(fileName).readlines()
    graph = {}
    for line in contents:
        (key,arrow,value) = line.partition(" -> ")
        key = key.rstrip()
        value = value.rstrip()
        currentValues = []
        currentString = ""
        for s in value:
            if s != ",":
                currentString += s
            else:
                currentValues += [currentString]
                currentString = ""
        currentValues += [currentString]
        graph[key] = graph.get(key,[]) + currentValues
    path = eulerianPath(graph)
    fout = open("texts\hw3q7Answer.txt", "wt")
    fout.write("->".join(path))
    print "->".join(path)
    

def makeEdgeList(graph):

    edgeList = []
    for key in graph:
        for kmer in graph.get(key):
            edgeList += (key,kmer)
    return edgeList
    





