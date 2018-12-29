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
        for i in xrange(len(leaves)):
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
            print v
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
    for i in xrange(n):
        distanceDict[i] = [0]*n
    for j in xrange(len(contents)):
        current = map(int, contents[j].split())
        for i in xrange(len(current)):
            distanceDict[i][j]=current[i]
    limbLength = limbLengthProblem(n,m,distanceDict)
    return limbLength

def additivePhylogeny(n, distanceDict, finalDict):

    if n == 2:
        #return the tree consisting of a single edge of length D1,2
        '''T = defaultdict(list)
        T[0] = [0]
        T[0][0]=distanceDict[0][1]
        T[1] = [0]
        T[1]= distanceDict[0][1]
        print "T is: " + str(T)
        return T'''
        return distanceDict
    print str(distanceDict)
    limbLength = limbLengthProblem(n-1,n-1,distanceDict)
    for j in xrange(n):
        distanceDict[j][n-1] = distanceDict[j][n-1] - limbLength
        distanceDict[n-1][j] = distanceDict[j][n-1]
    print str(distanceDict)
    (i2,n2,k2) = (0,0,0)
    for i in xrange(n):
        for k in xrange(n):
            din = distanceDict[i][n-1]
            dnk = distanceDict[n-1][k]
            dik = distanceDict[i][k]
            if dik == (din + dnk) and i != k: 
                (i2,n2,k2) = (i,n,k)
    print i2,n2,k2
    x = distanceDict[i2][n-1]
    distanceDict.pop(n-1)
    for key in distanceDict:
        if distanceDict[key] != []:
            distanceDict[key].pop(n-1)
    print str(distanceDict)
    T = additivePhylogeny(n-1, distanceDict, finalDict)
    print str(distanceDict)
    print str(i2)
    D2 = BFS2(distanceDict, i2)
    print "D2 is: " + str(D2)
    v = -1
    for key in D2:
        if D2[key] == x:
            v = key
    print v
    print "x is: " + str(x)
    print T
    print n-1
    if T[n-1] == 0:
        T[n-1] = [0]*n
    for extension in xrange(n-len(T[n-1])+1):
        T[n-1].append(0)
    print T
    if (v==-1):
        v = n
    T[n-1][v] = limbLength
    if T[v] == 0:
        T[v] = [0]*n
    for extension in xrange(n-len(T[v])):
        T[v].append(0)
    T[v][n-1] = limbLength
    print T
    return T

def BFS2(tree,leaf):

    queue = []
    visitedDict = defaultdict(int)
    seen = []
    queue += [leaf]
    visitedDict[leaf]=0
    while len(queue) != 0:
        current = queue.pop(0)
        for v in xrange(len(tree.get(current))):
            if v not in seen:
                seen += [v]
                queue += [v]
                visitedDict[v] = visitedDict[current] + tree.get(current)[v]
    return visitedDict

def additivePhylogenyWrapper(fileName):

    contents = open(fileName).readlines()
    n = int(contents.pop(0).strip())
    distanceDict = defaultdict(int)
    for i in xrange(n):
        distanceDict[i] = [0]*n
    for j in xrange(len(contents)):
        current = map(int, contents[j].split())
        for i in xrange(len(current)):
            distanceDict[i][j]=current[i]
    finalDict = {}
    treeDict = additivePhylogeny(n,distanceDict,finalDict)
    print str(treeDict)
    '''fout = open("texts/hw8/additivePhylogenyAnswer.txt","wt")
    output = ""
    for key in treeDict:
        current = treeDict.get(key)
        output += str(key) + "->" + str(current[0]) + ":" + str(current[1]) + "\n"
    fout.write(output)'''

graph2 = {1: [(2,5)],
          2: [(1,5),(3,21)],
          3: [(2,21)]}

graph = {1: [2, 3],
         2: [1, 4, 5, 6],
         3: [1, 4],
         4: [2, 3, 5],
         5: [2, 4, 6],
         6: [2, 5]}


'''def getDistance33(leaf,leaf2,seen,tree):

    if leaf==leaf2:
        print "DONE"
        return 0
    neighbors = tree.get(leaf)
    print leaf,leaf2
    if leaf in seen:
        print "FAIL"
        return 0 
    seen += [leaf] 
    for i in neighbors:
        current = neighbors
        
        dist = getDistance33(i[0],leaf2,seen,tree)
        return dist + int(i[1])

#distanceBetweenLeavesWrapper("texts/hw8/distanceBetweenLeaves.txt")

def getDistance(leaf,leaf2,seen,tree,total,dists):

    if leaf==leaf2:
        print "DONE"
        return total
    neighbors = tree.get(leaf)
    print leaf,leaf2
    if leaf in seen:
        print "FAIL"
        return 0 
    seen += [leaf] 
    for i in neighbors:
        current = neighbors
        total += int(i[1])
        dists += [total] 
        dist = getDistance(i[0],leaf2,seen,tree,total,dists)
        currentmax = 0
        for d in dists:
            if d > currentmax:
                currentmax = d
        return currentmax

def getDistance234(leaf,leaf2,seen,tree,total,dists):

    if leaf==leaf2:
        print "DONE"
        return total
    neighbors = tree.get(leaf)
    print leaf,leaf2
    if leaf in seen:
        print "FAIL"
        return 0 
    seen += [leaf] 
    for i in neighbors:
        current = neighbors
        total += int(i[1])
        dists += [total] 
        dist = getDistance(i[0],leaf2,seen,tree,total,dists)
        if dist != 0:
            return dist

def BFS2(tree,leaf,leaf2):
    queue = []
    seen = []
    path = [] 
    queue = [leaf] + queue
    total = 0
    while not len(queue)==0:
        
        current = queue.pop(len(queue)-1)
        if current not in path:
            path.append(current)
        if current == leaf2:
           return seen
        elif current not in seen:
            for edge in tree.get(current):
                if edge not in seen:
                    queue = [edge] + queue
            seen.append(current)
            
    return False

def BFS3(tree,leaf,leaf2):
    queue = []
    path = [] 
    queue += [leaf]
    total = 0
    while not len(queue)==0:
        
        current = queue.pop(len(queue)-1)
        if current not in path:
            path.append(current)
        if current == leaf2:
           return path
        for edge in tree.get(current):
            if edge not in path:
                queue += [edge]
    return False

def BFS4v1(tree,leaf,leaf2):

    queue = []
    visitedDict = defaultdict(int)
    seen = []
    queue += [leaf]
    visitedDict[leaf]=0
    while len(queue) != 0:
        current = queue.pop(0)
        for v in tree.get(current):
            if v not in seen:
                seen += [v]
                queue += [v]
                visitedDict[v] = visitedDict[current] + 1
    print str(visitedDict)

def getDistance21(tree,leaf,leaf2):

    total = 0
    seenQueue = []
    seenQueue += [leaf] 
    while len(seenQueue)!=0:
        
        new = [seenQueue.pop(len(seenQueue)-1)]
        other = new[len(new)-1]
        print new
        print other
        if leaf2 == other:
            print new
            return total
        neighbors = tree.get(other)
        for node in neighbors:
            if node[0] not in new:
                new2 = []
                new2 = new + [node[0]]
                seenQueue += new2

def populateDict(leaf, leaves, tree, resultDict):

    tree2 = copy.copy(tree) 
    currentDistance = 0
    parent = tree.get(leaf) #parent in form of (key, weight)
    neighbors = tree.get(parent[0][0]) #neighbors is [(key,weight), (key,weight) etc.] of parent
    currentDistance += int(parent[0][1]) #add weight of parent
    seen= []
    while len(neighbors) > 0: 
        currentNode = neighbors.pop(0)
        neighbors2 = tree.get(currentNode[0]) #neighbor of current neighbor (initially neighbors of parent)
        for i in neighbors2:
            currentDistance2 = currentDistance
            if i[0] in leaves:
                resultDict[leaf] += [currentDistance2 + int(i[1])]
            elif i[0] not in seen:
                neighbors += [i[0]]
                seen += [i[0]]
    return resultDict'''

