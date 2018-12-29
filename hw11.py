#Bioinformatics hw10
#Dhruv Ranjan

from collections import defaultdict
import sys
import math

def trieConstructionWrapper(fileName):

    contents = open(fileName).readlines()
    patterns = []
    for line in contents:
        patterns.append(line.strip())
    trie = trieConstruction(patterns)
    fout = open("texts/hw11/trieConstructionAnswer.txt", "wt")
    output = ""
    for key in trie:
        for node in trie[key]:
            output += str(key) + "->" + str(node) + ":" + str(trie[key][node]) + "\n"
    #print output
    fout.write(output)

def trieConstruction(patterns):

    trie = defaultdict(lambda: defaultdict(str))
    currentNode = 0
    newNode = 0
    for pattern in patterns:
        currentNode = 0
        for i in range(len(pattern)):
            currentSymbol = pattern[i]
            if currentSymbol in trie[currentNode].values():
                currentNode = getKey(trie[currentNode], currentSymbol)
            else:
                newNode += 1
                trie[currentNode][newNode] = currentSymbol
                currentNode = newNode
    return trie 
                

def getKey(trieDict, symbol):

    for key, sym in trieDict.iteritems():
        if sym == symbol:
            return key

def addRoot(trie, patterns):

    currentNode = 0
    symbols = []
    for i in range(len(patterns)):
        currentSymbol = patterns[i][0]
        if currentSymbol not in symbols:
            symbols.append(currentSymbol)
    for i in symbols:
        trie[0][currentNode+1] = i
        currentNode += 1
    return (trie, patterns)

def trieMatchingWrapper(fileName):

    contents = open(fileName).readlines()
    text = contents.pop(0).strip()
    patterns = []
    for line in contents:
        patterns.append(line.strip())
    positions = trieMatching(text, patterns)
    fout = open("texts/hw11/trieMatchingAnswer.txt", "wt")
    output = " ".join(str(i) for i in positions)
    fout.write(output)

def prefixTrieMatching(text, trie):

    symbol = text[0]
    text = text[1:]
    v = 0
    path = symbol
    while True:
        if v not in trie:
            return path[:-1]
        elif symbol in trie[v].values():
            v = getKey(trie[v], symbol)
            if len(text)!=0:
                symbol = text[0]
            else:
                return path
            text = text[1:]
            path += symbol
        else:
            #print "No matches found"
            return ""

def trieMatching(text, patterns):

    trie = trieConstruction(patterns)
    positions = []
    counter = 0
    while len(text) != 0:
        path = prefixTrieMatching(text, trie)
        if path in patterns:
            positions.append(counter)
        text = text[1:]
        counter += 1
    return positions

def suffixTreeWrapper(fileName):

    contents = open(fileName).readlines()
    text = contents[0]
    (tree, positionDict, labelDict, lengthDict) = suffixTreeConstruction(text)
    suffixes = []
    for v in positionDict:
        for w in positionDict[v]:
            pos = positionDict[v][w]
            length = lengthDict[v][w]
            current = text[pos:pos+length]
            suffixes.append(current)
    fout = open("texts/hw11/suffixTree.txt", "wt")
    output = "\n".join(str(i) for i in suffixes)
    #print output
    fout.write(output)
    
def suffixTrieConstruction(text):

    trie = defaultdict(lambda: defaultdict(str))
    positionDict = defaultdict(lambda: defaultdict(int))
    labelDict = defaultdict(int)
    currentNode = 0
    newNode = 0
    for i in range(len(text)):
        currentNode = 0
        for j in range(i,len(text)):
            currentSymbol = text[j]
            if currentSymbol in trie[currentNode].values():
                currentNode = getKey(trie[currentNode], currentSymbol)
            else:
                newNode += 1
                trie[currentNode][newNode] = currentSymbol
                positionDict[currentNode][newNode] = j
                currentNode = newNode
        if currentNode not in trie:
            labelDict[currentNode]=i
    return (trie, positionDict, labelDict)

def suffixTreeConstruction(text):

    (trie, positionDict, labelDict) = suffixTrieConstruction(text)
    paths = maxNonBranchingPaths(trie)
    lengthDict = defaultdict(lambda: defaultdict(int))
    positionDict2 = defaultdict(lambda: defaultdict(int))
    for path in paths:
        (v,u) = (path[0],path[-1])
        positionDict2[v][u] = positionDict[v][path[1]]
        lengthDict[v][u] = len(path)-1
    return (trie, positionDict2, labelDict, lengthDict)

def maxNonBranchingPaths(trie):

    paths = []
    for v in trie:
        if len(trie[v]) > 1:
            for w in trie[v]:
                path = [v,w]
                if w in trie:
                    while w in trie and len(trie[w])==1:
                        current = trie[w].keys()[0]
                        path.append(current)
                        w=current
                paths.append(path)
    return paths
    
def suffixArrayWrapper(fileName):

    contents = open(fileName).readlines()
    text = contents[0]
    sArray = suffixArray(text)
    output = str(sArray).lstrip("[").rstrip("]")
    fout = open("texts/hw11/suffixArrayAns.txt", "wt")
    fout.write(output)

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

def longestRepeatWrapper(fileName):

    contents = open(fileName).readlines()
    text = contents[0]
    longestRep = longestRepeat(text)
    print(longestRep)

def longestRepeat(text):

    (tree, positionDict, labelDict, lengthDict) = suffixTreeConstruction(text+"$")
    longest = ""
    root = 0
    seen = []
    queue = [root]
    lDict = defaultdict(int)
    sDict = defaultdict(str)
    while len(queue) != 0:
        v = queue.pop(0)
        for w in positionDict[v]:
            if len(positionDict[w])>1:
                queue.append(w)
                pos = positionDict[v][w]
                length = lengthDict[v][w]
                current = text[pos:pos+length]
                sDict[w]+=sDict[v]+current
                lDict[w]+=lengthDict[v][w] + lDict[v]            
    maxKey = max(lDict, key=lDict.get)
    return sDict[maxKey]

def longestSharedWrapper(fileName):

    contents = open(fileName).readlines()
    s1 = contents[0]
    s2 = contents[1]
    longestSharedAns = longestShared(s1,s2)
    print(longestSharedAns)

def longestShared(s1, s2):

    text = s1+"#"+s2+"$"
    (tree, positionDict, labelDict, lengthDict) = suffixTreeConstruction(text)
    #(tree2, positionDict2, labelDict2, lengthDict2) = suffixTreeConstruction(s2+"$")
    print(max(labelDict.values()))
    print(len(s1))
    print(len(s2))
    longest = ""
    root = 0
    queue = [root]
    lDict = defaultdict(int)
    sDict = defaultdict(str)
    while len(queue) != 0:
        v = queue.pop(0)
        for w in positionDict[v]:
            if len(positionDict[w])>1:
                queue.append(w)
                pos = positionDict[v][w]
                length = lengthDict[v][w]
                current = text[pos:pos+length]
                sDict[w]+=sDict[v]+current
                lDict[w]+=lengthDict[v][w] + lDict[v]
    maxKey = max(lDict, key=lDict.get)
    return sDict[maxKey]

def shortestNonSharedWrapper(fileName):

    contents = open(fileName).readlines()
    s1 = contents[0]
    s2 = contents[1]
    output = shortestNonShared(s1,s2)
    print(output)

def shortestNonShared(s1,s2):

    text = s1+"#"+s2+"$"
    (tree, positionDict, labelDict, lengthDict) = suffixTreeConstruction(text)
    longest = ""
    root = 0
    queue = [root]
    lDict = defaultdict(int)
    sDict = defaultdict(str)
    while len(queue) != 0:
        v = queue.pop(0)
        for w in positionDict[v]:
            #if len(positionDict[w])>1:
            queue.append(w)
            pos = positionDict[v][w]
            length = lengthDict[v][w]
            current = text[pos:pos+length]
            if (sDict[w]+sDict[v]+current) not in s2:
                #return sDict[w]+sDict[v]+current
                sDict[w]+=sDict[v]+current
                lDict[w]+=lengthDict[v][w] + lDict[v]
    maxKey = min(lDict, key=lDict.get)
    print(lDict)
    return sDict[maxKey]
