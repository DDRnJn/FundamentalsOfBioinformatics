#Bioinformatics hw10
#Dhruv Ranjan

from collections import defaultdict
import sys
import math

def farthestFirstTraversalWrapper(fileName):

    contents = open(fileName).readlines()
    first = contents.pop(0).split()
    k = int(first[0])
    m = int(first[1])
    points = []
    for point in contents:
        points += [map(float, point.split())]
    centers = farthestFirstTraversal(k,m,points)
    fout = open("texts/hw10/farthestFirstTraversalAnswer.txt","wt")
    output = ""
    for point in centers:
        line = " ".join(str(i) for i in point)
        output += line + "\n"
    print output
    fout.write(output)

def farthestFirstTraversal(k,m,points):

    centers = []
    firstCenter = points[0]
    centers.append(firstCenter)
    while len(centers) < k:
        maxDistPoint = maxDist(centers, points, m)
        centers.append(maxDistPoint)
    return centers

def maxDist(centers, points, m):

    maxDist = -sys.maxint-1
    maxDistPoint = None
    for i in xrange(len(points)):
        distance = 0.0
        minDist = sys.maxint
        for center in centers:
            distance = euclideanDistance(points[i],center, m)
            if distance < minDist:
                minDist = distance
        if minDist > maxDist:
            maxDist = minDist
            maxDistPoint = points[i]
    return maxDistPoint

def euclideanDistance(point, center, m):

    distance = 0.0
    for i in xrange(m):
        distance += (point[i]-center[i])**2
    return math.sqrt(distance)

def squaredErrorDistortionWrapper(fileName):

    contents = open(fileName).readlines()
    first = contents.pop(0).split()
    k = int(first[0])
    m = int(first[1])
    centers = []
    points = []
    for center in xrange(k):
        current = contents.pop(0)
        centers += [map(float, current.split())]
    contents.pop(0) #pop the "-------" line thing
    for point in contents:
        points += [map(float, point.split())]
    SED = squaredErrorDistortion(k,m,centers,points)
    return SED

def squaredErrorDistortion(k,m,centers,points):

    SED = 0.0
    for i in xrange(len(points)):
        distance = 0.0
        minDist = sys.maxint
        for center in centers:
            distance = euclideanDistance(points[i],center, m)
            if distance < minDist:
                minDist = distance
        SED += minDist**2
    return SED/float(len(points))
    
def lloydAlgorithmWrapper(fileName):

    contents = open(fileName).readlines()
    first = contents.pop(0).split()
    k = int(first[0])
    m = int(first[1])
    points = []
    for point in contents:
        points += [map(float, point.split())]
    centers = lloydAlgorithm(k,m,points)
    fout = open("texts/hw10/lloydAlgorithmAnswer.txt","wt")
    output = ""
    for point in centers:
        line = " ".join(str(i) for i in point)
        output += line + "\n"
    print output
    fout.write(output)
    
def lloydAlgorithm(k,m,points):

    clusters = defaultdict(list)
    centers = points[:k]
    while True:
        for i in xrange(k):
            clusters[i]=[]
        for i in xrange(len(points)):
            distance = 0.0
            minDist = sys.maxint
            minDistCenter = None
            cluster = -1
            for j in xrange(k):
                distance = euclideanDistance(points[i],centers[j], m)
                if distance < minDist:
                    minDist = distance
                    minDistCenter = centers[j]
                    cluster = j
            clusters[cluster].append(points[i])
        newCenters = []
        for cluster in clusters:
            newCenter = centerOfGravity(clusters[cluster],m)
            newCenters.append(newCenter)
        if centers == newCenters:
            break
        else:
            centers = newCenters
    return centers

def centerOfGravity(cluster,m):

    center = []
    for i in xrange(m):
        coord = 0.0
        for j in xrange(len(cluster)):
            coord += cluster[j][i]
        center += [coord/(float(len(cluster)))]
    return center

def EMaximizationWrapper(fileName):

    contents = open(fileName).readlines()
    first = contents.pop(0).split()
    k = int(first[0])
    m = int(first[1])
    stiffness = float(contents.pop(0))
    points = []
    for point in contents:
        points += [map(float, point.split())]
    centers = EMaximization(k,m,stiffness,points)
    fout = open("texts/hw10/EMaxAnswer.txt","wt")
    output = ""
    for point in centers:
        line = " ".join(str(i) for i in point)
        output += line + "\n"
    print output
    fout.write(output)

def EMaximization(k,m,stiffness,points):

    clusters = defaultdict(list)
    centers = points[:k]
    n = len(points)
    for e in xrange(100):
        newCenters = []
        hiddenMatrix = makeHM(k,m,stiffness,centers,points) #E step
        for i in xrange(k):
            currentCenter = [] 
            for j in xrange(m):
                numerator = getNum(hiddenMatrix, points, i, j)
                denominator = getDenom(hiddenMatrix,i,n)
                currentCenter+=[(numerator/denominator)]
            newCenters.append(currentCenter)
        centers = newCenters
    return centers

def getNum(hiddenM, points, i, j):

    res = 0.0
    for m in xrange(len(points)):
        res+=hiddenM[i][m]*points[m][j]
    return res

def getDenom(hiddenM, i, n):

    res = 0.0
    for m in xrange(n):
        res += hiddenM[i][m]
    return res

def makeHM(k,m,stiffness,centers,points): #Statistical partition way. 

    hm = defaultdict(lambda: defaultdict(float))
    n = len(points)
    for i in xrange(k):
        for j in xrange(n):
            numerator = math.exp(-stiffness*euclideanDistance(points[j],centers[i],m))
            denominator = sum([math.exp(-stiffness*euclideanDistance(points[j],x,m)) for x in centers])
            hm[i][j] = numerator/denominator
    return hm

def makeHM2(k,m,stiffness,centers,points): #The Newtonian way. 

    hm = defaultdict(lambda: defaultdict(float))
    n = len(points)
    for i in xrange(k):
        for j in xrange(n):
            if i!=j:
                numerator = 1/(euclideanDistance(points[j],centers[i],m))**2
                denominator = sum([1/(euclideanDistance(points[j],x,m))**2 for x in centers if x!=points[j]])
                hm[i][j] = numerator/denominator
            else:
                hm[i][j] = 0
    return hm

def hClusteringWrapper(fileName):

    contents = open(fileName).readlines()
    k = int(contents.pop(0))
    dMatrix = []
    for line in contents:
        dMatrix += [map(float, line.split())]
    
