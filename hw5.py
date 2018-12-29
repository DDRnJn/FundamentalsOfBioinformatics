#Bioinformatics hw5
#Dhruv Ranjan

import sys

def dpChange(money, coins):

    minNumCoins = []
    minNumCoins += [0]
    for m in xrange(money):
        minNumCoins += [0]
    for m in xrange(1,money+1):
        minNumCoins[m] = sys.maxint
        for i in xrange(len(coins)):
            if m >= coins[i]:
                if minNumCoins[m-coins[i]] + 1 < minNumCoins[m]:
                    minNumCoins[m] = minNumCoins[m-coins[i]]+1
    return minNumCoins[money]

def dpChangeWrapper(fileName):

    contents = open(fileName).readlines()
    money = contents[0].rstrip()
    numString = contents[1]
    coins = []
    currentNum = ""
    for s in numString:
        if s == ",":
            coins += [currentNum]
            currentNum = ""
        else:
            currentNum += s
    coins += [currentNum]
    coins[len(coins)-1] = coins[len(coins)-1].rstrip()
    coins = map(int, coins)
    minCoins = dpChange(int(money), coins)
    print minCoins

def manhattanTourWrapper(fileName):

    contents = open(fileName).readlines()
    (n, m) = contents[0].split()
    n = int(n)
    m = int(m)
    contents.pop(0)
    downMatrix = []
    rightMatrix = []
    currentLine = []
    currentNum = ""
    seenDash = False
    for line in contents:
        for num in line:
            if num == "-":
                seenDash = True
                break
            if num == " ":
                currentLine += [currentNum]
                currentNum = "" 
            else:
                currentNum += num
        if seenDash == True and len(currentLine) != 0:
            currentLine += [currentNum.rstrip()]
            rightMatrix += [map(int,currentLine)]
        elif seenDash == False:
            currentLine += [currentNum.rstrip()]
            downMatrix += [map(int,currentLine)]
        currentLine = []
        currentNum = ""
    manhattanTour = manhattanTourist(n,m,downMatrix,rightMatrix)
    print str(manhattanTour)

def manhattanTourist(n,m,down,right):

    s = []
    for i in xrange(n+1):
        s += [[]]
        for j in xrange(m+1):
            s[i] += [0]
    print s
    for i in xrange(1,n+1):
        s[i][0] = s[i-1][0] + down[i-1][0]
    for j in xrange(1,m+1):
        s[0][j] = s[0][j-1] + right[0][j-1]
    for i in xrange(1,n+1):
        for j in xrange(1,m+1):
            s[i][j] = max(s[i-1][j] + down[i-1][j],
                          s[i][j-1] + right[i][j-1])
    print s
    return s[n][m]

def outputLCS(backtrack,v,i,j, LCS):

    if i == 0 or j == 0:
        return LCS
    if backtrack[i][j] == 1:
        outputLCS(backtrack, v, i-1, j, LCS)
    elif backtrack[i][j] == 2:
        outputLCS(backtrack, v, i, j-1, LCS)
    else:
        print v[i-1]
        return outputLCS(backtrack, v, i-1, j-1, v[i-1] + LCS)     

def LCSBacktrack(v,w):
    
    s = []
    backtrack = []
    for i in xrange(len(v)+1):
        s += [[]]
        backtrack += [[]]
        for j in xrange(len(w)+1):
            s[i] += [0]
            backtrack[i] += [0]
    for i in xrange(len(v)):
        for j in xrange(len(w)):
            if v[i] != w[j]:
                s[i][j] = max(s[i-1][j],s[i][j-1],s[i-1][j-1])
            else:
                s[i][j] = max(s[i-1][j],s[i][j-1],s[i-1][j-1]+1)
            if s[i][j]==s[i-1][j]:
                backtrack[i][j]=1
            elif s[i][j]==s[i][j-1]:
                backtrack[i][j]=2
            else:
                backtrack[i][j]=3
    return backtrack

def outputLCSWrapper(fileName):

    contents = open(fileName).readlines()
    s = contents[0].rstrip()
    t = contents[1].rstrip()
    backtrack = LCSBacktrack(s,t)
    LCS = ""
    LCSa = outputLCS(backtrack, s, len(s), len(t),LCS)
    print LCSa
                
sys.setrecursionlimit(10000)        
