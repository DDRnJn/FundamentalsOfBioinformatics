#Bioinformatics hw8
#Dhruv Ranjan

from collections import defaultdict
import copy
import sys

def distanceBetweenLeaves(n, tree):

    resultDict = defaultdict(list)
    leaves = findLeaves(tree)
    #print leaves
    for leaf in leaves:
        resultDict[leaf] = [0]*len(leaves) 
    for leaf in leaves:
        distanceDict = BFS(tree,leaf)
        distanceDict[leaf]=0
        for i in range(len(leaves)):
            resultDict[leaf][i]=distanceDict[i]   
    #print str(resultDict)
    return resultDict

def BFS(tree,leaf):

    queue = []
    visitedDict = defaultdict(int)
    seen = []
    queue += [leaf]
    visitedDict[leaf]=0
    while len(queue) != 0:
        current = queue.pop(0)
        for v in tree.get(current):
            print(v)
            if v[0] not in seen:
                seen += [v[0]]
                queue += [v[0]]
                visitedDict[v[0]] = visitedDict[current] + v[1]
    return visitedDict

def findLeaves(tree):

    leaves = []
    #print str(tree)
    for key in tree:
        if len(tree.get(key)) == 1:
            leaves += [key]
    return leaves 

def distanceBetweenLeavesWrapper(fileName):
    
    contents = open(fileName).readlines()
    n = int(contents.pop(0).strip())
    tree = defaultdict(list) 
    for line in contents:
        current = line.split("->")
        key = int(current[0].strip())
        to = int(current[1].split(":")[0].strip())
        weight = int(current[1].split(":")[1].strip())
        tree[key] += [(to, weight)]
    distanceDict = distanceBetweenLeaves(n, tree)
    fout= open("texts/hw8/distanceBetweenLeavesAnswer.txt","wt")
    output = ""
    for key in distanceDict:
        current = distanceDict.get(key)
        output += " ".join(str(i) for i in current) + "\n"
    fout.write(output)

def limbLengthProblem(n,j,distanceDict):

    limbLength = sys.maxint
    for i in distanceDict:
        for k in distanceDict:
            if i != j and k != j: 
                dij = distanceDict[i][j]
                djk = distanceDict[j][k]
                dik = distanceDict[i][k]
                val = (dij + djk - dik)/2.0
                if val < limbLength:
                    limbLength = val
    return limbLength

def limbLengthProblemWrapper(fileName):

    contents = open(fileName).readlines()
    n = int(contents.pop(0).strip())
    m = int(contents.pop(0).strip())
    distanceDict = defaultdict(int)
    for i in range(n):
        distanceDict[i] = [0]*n
    for j in range(len(contents)):
        current = map(int, contents[j].split())
        for i in range(len(current)):
            distanceDict[i][j]=current[i]
    limbLength = limbLengthProblem(n,m,distanceDict)
    return limbLength

def additivePhylogeny(n, actual_n, distanceDict, finalDict):

    if n == 2:
        #return the tree consisting of a single edge of length D1,2
        #return distanceDict
        T = defaultdict(lambda: defaultdict(int))
        for key in distanceDict:
            for weight in range(len(distanceDict.get(key))):
                T[key][weight] = distanceDict[key][weight]
        print(str(T))
        return T
    print(str(distanceDict))
    limbLength = limbLengthProblem(n-1,n-1,distanceDict)
    for j in range(n):
        distanceDict[j][n-1] = distanceDict[j][n-1] - limbLength
        distanceDict[n-1][j] = distanceDict[j][n-1]
    print(str(distanceDict))
    (i2,n2,k2) = (0,0,0)
    for i in range(n):
        for k in range(n):
            din = distanceDict[i][n-1]
            dnk = distanceDict[n-1][k]
            dik = distanceDict[i][k]
            if dik == (din + dnk) and i != k: 
                (i2,n2,k2) = (i,n,k)
    print(i2,n2,k2)
    x = distanceDict[i2][n-1]
    distanceDict.pop(n-1)
    for key in distanceDict:
        if distanceDict[key] != []:
            distanceDict[key].pop(n-1)
    print(str(distanceDict))
    T = additivePhylogeny(n-1, actual_n, distanceDict, finalDict)
    v = -1
    if i2 in T: 
        distances = BFS2(T, i2)
        for key in distances:
            if distances[key] == x:
                v = key
    if v== -1:
        v = actual_n + 1
        actual_n += 1
    T[i2][v] = x
    T[i2][k2] = limbLength-x
    T[v][n2] = limbLength
    return T

def BFS2(tree,leaf):

    queue = []
    visitedDict = defaultdict(int)
    seen = []
    queue += [leaf]
    visitedDict[leaf]=0
    print(str(tree))
    print(str(queue))
    while len(queue) != 0:
        current = queue.pop(0)
        for v in range(len(tree.get(current))):
            if v not in seen:
                seen += [v]
                queue += [v]
                visitedDict[v] = visitedDict[current] + tree.get(current)[v]
    return visitedDict

def additivePhylogenyWrapper(fileName):

    contents = open(fileName).readlines()
    n = int(contents.pop(0).strip())
    distanceDict = defaultdict(int)
    for i in range(n):
        distanceDict[i] = [0]*n
    for j in range(len(contents)):
        current = map(int, contents[j].split())
        for i in range(len(current)):
            distanceDict[i][j]=current[i]
    finalDict = {}
    treeDict = additivePhylogeny(n,n,distanceDict,finalDict)
    fout = open("texts/hw8/additivePhylogenyAnswer.txt","wt")
    output = ""
    for key in treeDict:
        for node in treeDict[key]: 
            output += str(key) + "->" + str(node) + ":" + str(treeDict[key][node]) + "\n"
    fout.write(output)
    print(str(treeDict))
    
