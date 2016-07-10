import math,copy,itertools
from sets import Set as set
from random import shuffle


def getSets(N, D,distanceO,K):
    S=[]
    combinations=[]
    for i in range(1,K+1):
        combinations+=list(itertools.combinations(N,i))

    for tour in combinations:
        tour=list(tour)
        d,tour=distanceO.minimumTourDistance(tour,D)
        S.append([tour, d])
        
    return S    

def doExtraGreedy(N,D,distanceO,K):
    numbersCollected=set()
    lengthToHave = len(N) 
    C=[]
    while not len(numbersCollected)>=lengthToHave:
        C,numbersCollected = selectNextBest(numbersCollected,C,distanceO,N,D,K)
        allNodesInC=[node[0] for tour in C for node in tour[0]]
        assert len(allNodesInC)==len(set(allNodesInC))," all nodes: %s \n set: %s" %(allNodesInC, set(allNodesInC))
    return C
        
def selectNextBest(numbersCollected,C,distanceO,N,D,K):
    newNodes = [node for node in N if node[0] not in numbersCollected]
    assert len(newNodes)>0, "Did not find new node in C: \n %s \n which isnt already in \n %s " %(C, numbersCollected)
    if len(newNodes)<=K:
        distanceT, tour=distanceO.minimumTourDistance(newNodes, D)
        newNumbers=[node[0] for node in tour]
        numbersCollected.update(newNumbers)
        C.append([tour,distanceT])
        return C,numbersCollected
           
    newTour = [getClosestNodeToADepot(D,newNodes, distanceO)]
    for i in range(1,K):
        latestNode = newTour[-1]
        options = [node for node in newNodes if node not in newTour]
        newTour.append(getClosestNode(latestNode, options, distanceO))
    options = [node for node in newNodes if node not in newTour] 
    latestNode = newTour[-1]   
    distance = float('inf')
    nodeNew=None
    for node in options:
        d=distanceO.EuclideanNode(latestNode, node)
        d +=distanceO.distanceToDepot(node, D)
        if d<distance:
            distance=d
            nodeNew=node
            
    assert nodeNew, 'No new node found, %s' %len(newNodes)
    distanceT, newTour = distanceO.minimumTourDistance(newTour, D)
    numbersCollected.update([node[0] for node in newTour])
    C.append([newTour,distanceT ])
    return C,numbersCollected

def getClosestNode(node1,N,distanceO):
    distance=float('inf')
    closestNode=None
    for node2 in N:
        if node1 is node2:
            continue
        d=distanceO.EuclideanNode(node1,node2)
        if d<distance:
            distance = d
            closestNode=node2
    assert closestNode, 'Closest node does not exist'
    return closestNode
    
def getClosestNodeToADepot(D,N,distanceO):
    distance=float('inf')
    closestNode=None
    for node in N:
        d=distanceO.distanceToDepot(node, D)
        if d<distance:
            distance = d
            closestNode=node
    assert closestNode, 'Closest node does not exist'
    return closestNode
    
       
def doGreedy(N,D, distanceO,K):
    S = getSets(N, D,distanceO,K)
    Scopy = S
    C=[]
    numbersCollected=set()
    lengthToHave = len(N) 
    while not len(numbersCollected)>=lengthToHave:
        assert len(Scopy)>0, 'O dear'
        s = getBestSubset(Scopy,numbersCollected,lengthToHave)
        C.append(s)
        Scopy.remove(s)
        numbersCollected.update([j for j in s[0]])
    #print C
    #print "With total length: %s" %(sum([i[1] for i in C]))
    return C

def getBestSubset(Scopy,numbersCollected,lengthToHave):
    cost = float("inf")
    bestSet=[]
    for s in Scopy:
        amountnew = amountNew(s,numbersCollected)
        if amountnew==0:
            continue
        setCost  = (s[1]/ float(amountnew))#cost per new element, check for int and floats etc
        if setCost < cost:
            cost=setCost
            bestSet=s
    assert len(bestSet)>0, "getBestSubset got weird result, %s \n amountnew :%s \n numbers: %s \n and length to have %s" %(Scopy, amountnew, sorted(numbersCollected), lengthToHave)
    return bestSet

def amountNew(s, numbersCollected):
    sIDs = [i for i in s[0]]
    return len(set(sIDs)- set(numbersCollected))
