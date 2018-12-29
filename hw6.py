#Bioinformatics hw6
#Dhruv Ranjan

import sys

def importMatrix(fileName):

    contents = open(fileName).readlines()
    keys = [str(i) for i in contents[0].strip().split()]
    contents.pop(0)
    for i in range(len(contents)):
        line = contents[i][1:]
        line.strip()
        contents[i] = line
    blosum62 = [[int(i) for i in line.split()] for line in contents]
    return (keys, blosum62)
    
def globalAlignment(s, t, indelCost, blosum62Dict, lookupDict):

    if not(len(s))>0:
        return (indelCost*len(t), '-'*len(t),t)
    if not(len(t))>0:
        return (indelCost*len(s), s, '-'*len(s))
    currentSuffixS = s[1:]
    currentSuffixT = t[1:]
    (solution1, lookupDict) = getSolution(currentSuffixS, t, indelCost, blosum62Dict, lookupDict)
    (solution2, lookupDict) = getSolution(s, currentSuffixT, indelCost, blosum62Dict, lookupDict)
    (solution3, lookupDict) = getSolution(currentSuffixS, currentSuffixT, indelCost, blosum62Dict, lookupDict)
    answers = [(indelCost+solution1[0],s[0]+solution1[1],'-'+solution1[2]),
               (indelCost+solution2[0],'-'+solution2[1],t[0]+solution2[2]),
               (blosum62Dict[(s[0],t[0])]+solution3[0],s[0]+solution3[1],t[0]+solution3[2])]
    bestAlignmentIndex = 0
    currentBest = -sys.maxint - 1
    for i in range(len(answers)):
        if answers[i][0]>currentBest:
            currentBest = answers[i][0]
            bestAlignmentIndex = i
    return answers[bestAlignmentIndex]

def getSolution(suffix1, suffix2, indelCost, blosum62Dict, lookupDict):
    if (suffix1, suffix2) in lookupDict:
        solution = lookupDict[suffix1,suffix2]
    else:
        solution = globalAlignment(suffix1, suffix2, indelCost, blosum62Dict, lookupDict)
        lookupDict[suffix1, suffix2]=solution
    return (solution, lookupDict)

def globalAlignmentWrapper(fileName):

    (keys, blosum62) = importMatrix("texts\hw6\hblosum62.txt")
    contents = open(fileName).readlines()
    s = contents[0].strip()
    t = contents[1].strip()
    print(s)
    print(t)
    indelCost = -5
    backtrack = []
    '''for i in range(len(s)):
        backtrack += [[]]
        for j in range(len(t)):
            backtrack += [0]'''
    lookupDict = {} #Used for memoization. 
    blosum62Dict = {} #dict mapping amino acid pair to blosum62 score.
                      #seems easier than using lists
    for i in range(len(keys)):
        for j in range(len(keys)):
            blosum62Dict[(keys[i],keys[j])] = blosum62[i][j]
    (score, s1, t1) = globalAlignment(s,t, indelCost, blosum62Dict, lookupDict)
    print(score)
    print("\n")
    print(s1)
    print("\n")
    print(t1)

def getSolution2(suffix1, suffix2, indelCost, pam250Dict, lookupDict, best):
    if (suffix1, suffix2) in lookupDict:
        solution = lookupDict[suffix1,suffix2]
    else:
        (best, solution) = localAlignment(suffix1, suffix2, indelCost, pam250Dict, lookupDict, best)
        lookupDict[suffix1, suffix2]=solution
    return (solution, lookupDict)

def localAlignment(s,t,indelCost, pam250Dict, lookupDict, best):

    if not(len(s))>0:
        lookupDict[s,t]=(0,"",t)
        return (best, (0,"",t))
    if not(len(t))>0:
        lookupDict[s,t]=(0,s,"")
        return (best, (0,s,""))
    currentSuffixS = s[1:]
    currentSuffixT = t[1:]
    (solution1, lookupDict) = getSolution2(currentSuffixS, t, indelCost, pam250Dict, lookupDict, best)
    (solution2, lookupDict) = getSolution2(s, currentSuffixT, indelCost, pam250Dict, lookupDict, best)
    (solution3, lookupDict) = getSolution2(currentSuffixS, currentSuffixT, indelCost, pam250Dict, lookupDict, best)
    answers = [(indelCost+solution1[0],s[0]+solution1[1],'-'+solution1[2]),
               (indelCost+solution2[0],'-'+solution2[1],t[0]+solution2[2]),
               (pam250Dict[(s[0],t[0])]+solution3[0],s[0]+solution3[1],t[0]+solution3[2]),
               (0, s[0]+"-", t[0]+"-")]
    bestAlignmentIndex = 0
    currentBest = -sys.maxint - 1
    for i in range(len(answers)):
        if answers[i][0]>currentBest:
            currentBest = answers[i][0]
            bestAlignmentIndex = i
    if bestAlignmentIndex == 3:
        #lookupDict[currentSuffixS, currentSuffixT] = (0, s[0]+solution3[1], t[0]+solution3[2])
        lookupDict[s, currentSuffixT] = (0, "", "")
        lookupDict[currentSuffixS, t] = (0, "", "")
    best += [answers[bestAlignmentIndex]]
    return (best, answers[bestAlignmentIndex])
        
#localAlignmentWrapper("texts\hw6\localAlignment.txt")    
def localAlignmentWrapper(fileName):

    (keys, pam250) = importMatrix("texts\hw6\hpam250.txt")
    contents = open(fileName).readlines()
    s = contents[0].strip()
    t = contents[1].strip()
    prin(s)
    print(t)
    indelCost = -5
    backtrack = []
    lookupDict = {} #Used for memoization. 
    pam250Dict = {} #dict mapping amino acid pair to blosum62 score.
                      #seems easier than using lists
    for i in range(len(keys)):
        for j in range(len(keys)):
            pam250Dict[(keys[i],keys[j])] = pam250[i][j]
    best = []
    (best, (score, s1, t1)) = localAlignment(s,t, indelCost, pam250Dict, lookupDict, best)
    highestScore = 0
    (bestS, bestT) = ("", "")
    for i in best:
        if i[0] > highestScore:
            highestScore = i[0]
            (bestS, bestT) = (i[1],i[2])
    print(highestScore)
    print("\n")
    print(bestS)
    print("\n")
    print(bestT)

def editDistance(s,t, lookupDict):

    if len(s)==0:
        return len(t)
    elif len(t)==0:
        return len(s)
    else:
        (solution1, lookupDict) = getSolution3(s[:len(s)-1],t[:len(t)-1],lookupDict)
        (solution2, lookupDict) = getSolution3(s[:len(s)-1],t,lookupDict)   
        (solution3, lookupDict) = getSolution3(s,t[:len(t)-1], lookupDict)
        if s[len(s)-1]==t[len(t)-1]:
            ans1 = solution1
        else:
            ans1 = solution1+1
        ans2 = solution2+1
        ans3 = solution3+1
    return min(ans1, ans2, ans3)

def getSolution3(suffix1, suffix2, lookupDict):

    if (suffix1, suffix2) in lookupDict:
        solution = lookupDict[suffix1,suffix2]
    else:
        solution = editDistance(suffix1, suffix2, lookupDict)
        lookupDict[suffix1, suffix2]=solution
    return (solution, lookupDict)

def editDistanceWrapper(fileName):
    contents = open(fileName).readlines()
    s = contents[0].strip()
    t = contents[1].strip()
    print(s)
    print("\n")
    print(t)
    lookupDict = {}
    distance = editDistance(s,t,lookupDict)
    print(distance)

def fittingAlignment(s,t,matchScore, mismatchCost, indelCost, lookupDict,best):
    
    if not(len(s))>0:
        return (best,(indelCost*len(t), '-'*len(t),t))
    if not(len(t))>0:
        return (best,(0 * len(s),s,'-'*len(s)))
    currentSuffixS = s[1:]
    currentSuffixT = t[1:]
    (solution1, lookupDict) = getSolution4(currentSuffixS, t, matchScore, mismatchCost, indelCost, lookupDict,best)
    (solution2, lookupDict) = getSolution4(s,currentSuffixT, matchScore, mismatchCost, indelCost, lookupDict,best)
    (solution3, lookupDict) = getSolution4(currentSuffixS, currentSuffixT, matchScore, mismatchCost, indelCost, lookupDict,best)
    if s[0]==t[0]:
        a3 = matchScore
    else:
        a3 = mismatchCost
    answers = [(indelCost+solution1[1][0],s[0]+solution1[1][1],'-'+solution1[1][2]),
               (indelCost+solution2[1][0],'-'+solution2[1][1],t[0]+solution2[1][2]),
               (a3+solution3[1][0],s[0]+solution3[1][1],t[0]+solution3[1][2])]
    bestAlignmentIndex = 0
    currentBest = -sys.maxint - 1
    for i in range(len(answers)):
        if answers[i][0]>currentBest:
            currentBest = answers[i][0]
            bestAlignmentIndex = i
    best += [answers[bestAlignmentIndex]]
    return (best,answers[bestAlignmentIndex])

def getSolution4(suffix1, suffix2, matchScore, mismatchCost, indelCost, lookupDict,best):
    if (suffix1, suffix2) in lookupDict:
        solution = lookupDict[suffix1,suffix2]
    else:
        solution = fittingAlignment(suffix1, suffix2, matchScore, mismatchCost, indelCost, lookupDict,best)
        lookupDict[suffix1, suffix2]=solution
    return (solution, lookupDict)        

def fittingAlignmentWrapper(fileName):

    contents = open(fileName).readlines()
    s = contents[0].strip()
    t = contents[1].strip()
    indelCost = -1
    matchScore = 1
    mismatchCost = -1
    print(s)
    print("\n")
    print(t)
    lookupDict = {}
    best=[]
    (best,(score, s1, t1)) = fittingAlignment(s,t,matchScore, mismatchCost, indelCost, lookupDict,best)
    highestScore = 0
    (bestS, bestT) = ("", "")
    for i in best:
        if i[0] > highestScore and (checkSeq(i[2]))>=len(t):
            highestScore = i[0]
            (bestS, bestT) = (i[1],i[2])
    (cleanedS, cleanedT) = cleanSeqs(bestS, bestT)
    print(highestScore)#, score
    print("\n")
    print(cleanedS)#, s1
    print("\n")
    print(cleanedT)#, t1

def checkSeq(s):

    eLen = 0
    for i in s:
        if i != '-':
            eLen+=1
    return eLen

def cleanSeqs(s,t):

    cleanedT = t.strip('-')
    sIndex = len(cleanedT)
    cleanedS = s[:sIndex]
    return cleanedS, cleanedT

def getSolution5(suffix1, suffix2, matchScore, mismatchCost, indelCost, lookupDict,best):
    if (suffix1, suffix2) in lookupDict:
        solution = lookupDict[suffix1,suffix2]
    else:
        solution = overlapAlignment(suffix1, suffix2, matchScore, mismatchCost, indelCost, lookupDict,best)
        lookupDict[suffix1, suffix2]=solution
    return (solution, lookupDict)

def overlapAlignment(s,t,matchScore, mismatchCost, indelCost, lookupDict,best):
    
    if not(len(s))>0:
        return (best,(0*len(t), '-'*len(t),t))
    if not(len(t))>0:
        return (best,(indelCost * len(s),s,'-'*len(s)))
    currentSuffixS = s[1:]
    currentSuffixT = t[1:]
    (solution1, lookupDict) = getSolution5(currentSuffixS, t, matchScore, mismatchCost, indelCost, lookupDict,best)
    (solution2, lookupDict) = getSolution5(s,currentSuffixT, matchScore, mismatchCost, indelCost, lookupDict,best)
    (solution3, lookupDict) = getSolution5(currentSuffixS, currentSuffixT, matchScore, mismatchCost, indelCost, lookupDict,best)
    if s[0]==t[0]:
        a3 = matchScore
    else:
        a3 = mismatchCost
    answers = [(indelCost+solution1[1][0],s[0]+solution1[1][1],'-'+solution1[1][2]),
               (indelCost+solution2[1][0],'-'+solution2[1][1],t[0]+solution2[1][2]),
               (a3+solution3[1][0],s[0]+solution3[1][1],t[0]+solution3[1][2])]
    bestAlignmentIndex = 0
    currentBest = -sys.maxint - 1
    for i in range(len(answers)):
        if answers[i][0]>currentBest:
            currentBest = answers[i][0]
            bestAlignmentIndex = i
    best += [answers[bestAlignmentIndex]]
    return (best,answers[bestAlignmentIndex])

def overlapAlignmentWrapper(fileName):

    contents = open(fileName).readlines()
    s = contents[0].strip()
    t = contents[1].strip()
    indelCost = -2
    matchScore = 1
    mismatchCost = -2
    print(s)
    print("\n")
    print(t)
    lookupDict = {}
    best=[]
    (best,(score, s1, t1)) = overlapAlignment(s,t,matchScore, mismatchCost, indelCost, lookupDict,best)
    highestScore = 0
    (bestS, bestT) = ("", "")
    #print best
    for i in best:
        if i[0] > highestScore and (checkSeq(i[2]))>=len(t):
            highestScore = i[0]
            (bestS, bestT) = (i[1],i[2])
    (cleanedS, cleanedT) = cleanSeqs2(bestS, bestT)
    print(highestScore)#, score
    print("\n")
    print(cleanedS)#, s1
    print("\n")
    print(cleanedT)#, t1

def cleanSeqs2(s,t):

    cleanedS = s.strip('-')
    tIndex = len(cleanedS)
    cleanedT = t[:tIndex]
    return cleanedS, cleanedT

def getSolutionL(s,t,blosum62Dict,lowerDict,upperDict,middleDict,gapOpen,gapExtend):
    if (s, t) in lowerDict:
        solution = lowerDict[s,t]
    else:
        solution = affineGapAlignment(s,t,blosum62Dict,lowerDict,upperDict,middleDict,gapOpen,gapExtend)
        lowerDict[s, t]=solution
    return (solution, lowerDict, upperDict, middleDict)

def getSolutionU(s,t,blosum62Dict,lowerDict,upperDict,middleDict,gapOpen,gapExtend):
    if (s, t) in upperDict:
        solution = upperDict[s,t]
    else:
        solution = affineGapAlignment(s,t,blosum62Dict,lowerDict,upperDict,middleDict,gapOpen,gapExtend)
        upperDict[s, t]=solution
    return (solution, lowerDict, upperDict, middleDict)

def getSolutionM(s,t,blosum62Dict,lowerDict,upperDict,middleDict,gapOpen,gapExtend):
    if (s, t) in middleDict:
        solution = middleDict[s,t]
    else:
        solution = affineGapAlignment(s,t,blosum62Dict,lowerDict,upperDict,middleDict,gapOpen,gapExtend)
        middleDict[s, t]=solution
    return (solution, lowerDict, upperDict, middleDict)

#affineGapAlignmentWrapper("texts\hw6\haffineGap.txt")
def affineGapAlignment(s,t,blosum62Dict,lowerDict,upperDict,middleDict,gapOpen,gapExtend):

    if not(len(s))>0:
        return (gapOpen + gapExtend*len(t), '-'*len(t),t)
    if not(len(t))>0:
        return (gapOpen + gapExtend*len(s), s, '-'*len(s))
    
    currentSuffixS = s[1:]
    currentSuffixT = t[1:]
    
    (lower1, lowerDict, upperDict, middleDict)=getSolutionL(currentSuffixS,t,blosum62Dict,lowerDict,upperDict,middleDict,gapOpen,gapExtend)
    #print "lower1: " + str(lower1)
    (lower2, lowerDict, upperDict, middleDict)=getSolutionM(currentSuffixS,t,blosum62Dict,lowerDict,upperDict,middleDict,gapOpen,gapExtend)
    #print "lower2: " + str(lower2)
    (upper1, lowerDict, upperDict, middleDict)=getSolutionU(s,currentSuffixT,blosum62Dict,lowerDict,upperDict,middleDict,gapOpen,gapExtend)
    #print "upper1: " + str(upper1)
    (upper2, lowerDict, upperDict, middleDict)=getSolutionM(s,currentSuffixT,blosum62Dict,lowerDict,upperDict,middleDict,gapOpen,gapExtend)
    #print "upper2: " + str(upper2)
    (middle1, lowerDict, upperDict, middleDict)=getSolutionL(currentSuffixS,currentSuffixT,blosum62Dict,lowerDict,upperDict,middleDict,gapOpen,gapExtend)
    #print "middle1: " + str(middle1)
    (middle2, lowerDict, upperDict, middleDict)=getSolutionM(currentSuffixS,currentSuffixT,blosum62Dict,lowerDict,upperDict,middleDict,gapOpen,gapExtend)
    #print "middle2: " + str(middle2)
    (middle3, lowerDict, upperDict, middleDict)=getSolutionU(currentSuffixS,currentSuffixT,blosum62Dict,lowerDict,upperDict,middleDict,gapOpen,gapExtend)
    #print "middle3: " + str(middle3)    
    leftAnswer1 = (gapExtend+lower1[0],s[0]+lower1[1],'-'+lower1[2])
    leftAnswer2 = (gapOpen + lower2[0],s[0]+lower2[1],'-'+lower2[2])
    
    upperAnswer1 = (gapExtend+upper1[0],'-'+upper1[1],t[0]+upper1[2])
    upperAnswer2 = (gapOpen + upper2[0],'-'+upper2[1],t[0]+upper2[2])
    
    middleAnswer1 = (middle1[0], middle1[1], middle1[2])
    middleAnswer2 = (blosum62Dict[(s[0],t[0])]+middle2[0],s[0]+middle2[1],t[0]+middle2[2])
    middleAnswer3 = (middle3[0], middle3[1], middle3[2])
    
    bestLower = max(leftAnswer1[0],leftAnswer2[0])
    bestUpper = max(upperAnswer1[0],upperAnswer2[0])
    bestMiddle = max(middleAnswer1[0], middleAnswer2[0], middleAnswer3[0])
    
    if bestLower==leftAnswer1[0]:
        bestLowerAnswer=leftAnswer1
    else:
        bestLowerAnswer=leftAnswer2
        
    if bestUpper==upperAnswer1[0]:
        bestUpperAnswer=upperAnswer1
    else:
        bestUpperAnswer=upperAnswer2

    if bestMiddle==middleAnswer1[0]:
        bestMiddleAnswer=middleAnswer1
    elif bestMiddle==middleAnswer2[0]:
        bestMiddleAnswer=middleAnswer2
    else:
        bestMiddleAnswer=middleAnswer3
        
    bestOverallAnswer=max(bestLower,bestUpper,bestMiddle)

    if bestOverallAnswer==bestLower:
        return bestLowerAnswer
    elif bestOverallAnswer==bestUpper:
        return bestUpperAnswer
    else:
        return bestMiddleAnswer

def affineGapAlignmentWrapper(fileName):

    (keys, blosum62) = importMatrix("texts\hw6\hblosum62.txt")
    contents = open(fileName).readlines()
    s = contents[0].strip()
    t = contents[1].strip()
    print(s)
    print(t)
    gapOpen = -11
    gapExtend = -1
    lowerDict = {}
    upperDict = {}
    middleDict = {}
    blosum62Dict = {} #dict mapping amino acid pair to blosum62 score.
                      #seems easier than using lists
    for i in range(len(keys)):
        for j in range(len(keys)):
            blosum62Dict[(keys[i],keys[j])] = blosum62[i][j]
    (score, s1, t1) = affineGapAlignment(s,t,blosum62Dict,lowerDict,upperDict,middleDict,gapOpen,gapExtend)
    print(score)
    print(s1)
    print(t1)
    
sys.setrecursionlimit(5000)
