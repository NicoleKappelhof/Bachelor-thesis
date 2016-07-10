import copy, numpy
from sets import Set as set
import matplotlib.pyplot as plt
import Greedy.greedy as greedy
import Exact.ExactDepots as depots
from random import shuffle
from Main import Distances 

problemConfiguration = None
gridSize = 15
stopAdd = 1
generations=10
similarTours=4.0/5
amountOfParents=7
minimumAmountOfNewParents = 5
BestParentsX=[]
BestParentsY=[]
AvParentsX=[]
AvParentsY=[]
BestSolution=[float('inf')]  

def checkNoDoubleNodeTours(tours):
    for tour in tours:
        nodeTour=[node[0] for node in tour]
        assert len(nodeTour)==len(set(nodeTour)), "node tour %s \n not same as \n %s" %(nodeTour,set(nodeTour))    
    allNodes = [node[0] for tour in tours for node in tour]
    assert len(allNodes)==len(set(allNodes)),"node tour %s \n not same as \n %s" %(allNodes,set(allNodes)) 

def toSimpleTours(tours):
    simple=[]
    for tour in tours:
        simple.append([node[0] for node in tour])
    return simple  

def nodesListShareElement(list1, list2):
    nodeNumbers1 = set([node[0] for node in list1])
    nodeNumbers2 = set([node[0] for node in list2])   
    return len(nodeNumbers1.intersection(nodeNumbers2))>0

def nofNodesListHaveInCommon(list1, list2):
    nodeNumbers1 = set([node[0] for node in list1])
    nodeNumbers2 = set([node[0] for node in list2])   
    return len(nodeNumbers1.intersection(nodeNumbers2))

def toursAndNodeListAreEqual(tours,list1):
    nodesSet1 = set(node[0] for tour in tours for node in tour)
    nodeSet2 = set(node[0] for node in list1)
    return (nodesSet1==nodeSet2)

def toursContainExactSameNodes(tours1,tours2):
    nodesSet1 = set(node[0] for tour in tours1 for node in tour)
    nodesSet2 = set(node[0] for tour in tours2 for node in tour)
    return (nodesSet1==nodesSet2)

def toursAndNodeListShareElement(tours,list1):
    nodesSet1 = set(node[0] for tour in tours for node in tour)
    nodeSet2 = set(node[0] for node in list1)
    return len(nodesSet1.intersection(nodeSet2))>0

def toursTooFarRemoved(tour1,tour2,distanceO,D):
    if len(tour1)==0 or len(tour2)==0:
        return False
    d=[]
    d.append(distanceO.distanceToDepot(tour1[0], D))
    d.append(distanceO.distanceToDepot(tour1[-1], D))
    for i in range(1,len(tour1)):
        d.append(distanceO.EuclideanNode(tour1[i], tour1[i-1]))
    maximimDistance1=max(d)
    d=[]
    d.append(distanceO.distanceToDepot(tour2[0], D))
    d.append(distanceO.distanceToDepot(tour2[-1], D))
    for i in range(1,len(tour2)):
        d.append(distanceO.EuclideanNode(tour2[i], tour2[i-1]))
    maximimDistance2=max(d)        
    d=[]
    for node1 in tour1:
        for node2 in tour2:
            d.append(distanceO.EuclideanNode(node1, node2))
    minimumDistanceBetweenTours=min(d)  
    if minimumDistanceBetweenTours > max([maximimDistance1,maximimDistance2]):
        return True        
    return False      

def extraGreedySolutionToTours(greedySolution,N):
    tours = [s[0] for s in greedySolution]
    length = sum([i[1] for i in greedySolution])   
    return tours,length

def extraGreedySolutionToParent(greedySolution,depot1, depot2, N):
    tours,length = extraGreedySolutionToTours(greedySolution,N)          
    solution = [[depot1,depot2 ]] + tours + [length]
    return solution

def tabuChange(D,tabu):
    tabuS = sorted(tabu) 
    for t in tabuS:
        if  t-0.2<D[1] < t+0.2 and D[1]<9.8:
            D[1]+=0.2
    tabuS = sorted(tabu, reverse=True) 
    for t in tabuS:
        if  t-0.2<D[0] < t+0.2 and D[0]>0.2:
            D[0]-=0.2
    if len(tabu)>6:
        tabu = tabu[2:]
    tabu.append(D[0])
    tabu.append(D[1])    
    return D 

def getChildDistanceAndDepots(child,distanceO,tabu=None):
    assert len(child)>2, "Incorrect child: %s" %child
    childTours = child[1:-1]
    nodesConnectedToDepot=[[tour[0], tour[-1]] for tour in childTours]
    nodesConnectedToDepot = [i for subtour in nodesConnectedToDepot for i in subtour]
    solution = depots.exactDepots(nodesConnectedToDepot,2);

    D = [solution[1],solution[2]]
    if tabu:
        D = tabuChange(D,tabu)      
    distance = distanceO.getDistanceTours(childTours,D)

    return distance,D

def replaceDepotsFromChild(child,distanceM,tabu=None):
    distance1, D= getChildDistanceAndDepots(child,distanceM,tabu)
    child[0]=D
    child[-1]=distance1
    return child

def deltaInSolution(parent1,parent2,b):
    distance1=parent1[-1]
    distance2=parent2[-1] 
    Dp1 = parent1[0]    
    d1p1,d2p1=min(Dp1),max(Dp1)
    Dp2=parent2[0]
    d1p2,d2p2=min(Dp2),max(Dp2)  
    if abs(d1p1-d1p2)<b/200.0 and abs(d2p1-d2p2)<b/200.0 and abs(distance1-distance2)<0.5:
        return False
    return True

def parentsAreTooSimilar(parent1,parent2):
    global similarTours
    tours1=parent1[1:-1]
    tours2=parent2[1:-1]
    similar = [tour for tour in tours1 if tour in tours2 ]
    if len(similar)>max(similarTours*len(tours1), similarTours*len(tours2)):
        return True
    return False

def shuffleGeneration(parents):
    for i,parent in enumerate(parents):
        middle = parent[1:-1]
        shuffle(middle)
        parents[i]=[parent[0]]+middle+[parent[-1]]
    return parents

def initialize(N,b,distanceO):
    global amountOfParents
    parents = []
    stepsize = float(b)/amountOfParents
    depotLocations = numpy.arange(0,b,stepsize-0.01)
    for i in range(1,int(amountOfParents/2.0)):
        depot1 = depotLocations[i-1]
        depot2 = depotLocations[-i]
        parents.append(extraGreedySolutionToParent(greedy.doExtraGreedy(N,[depot1, depot2],distanceO,K),depot1,depot2,N))
    newstart= int(amountOfParents/2.0)
    for i in range(newstart,amountOfParents + 1):
        depot1 = depotLocations[i-newstart+1]
        depot2 = depotLocations[-(i-newstart+1)]
        tours = [[node] for node in N]
        d=distanceO.getDistanceTours(tours, [depot1,depot2])
        parents.append([[depot1,depot2]]+ tours + [d])
    
    p1=parents[0][1:-1]
    p2=parents[1][1:-1]
    same=[tour for tour in p1 if tour in p2]
    assert len(same)<(3.0/4)*len(p2), '\n p1: %s \n p2: %s \n same: %s\n %s vs %s, \n %s \n depotLocations: %s \n and %s' %(p1,p2,same,len(same),(3.0/4)*(len(p2)),len(same)<(3.0/4)*len(p2),parents[0][0],parents[1][0])
    assert len(parents)==amountOfParents, "len parents %s, amoutnOfParents %s" %(len(parents),amountOfParents)
    return parents

def interSwapTour(tour,D, distanceO):  
    tour=list(tour) 
    distanceOriginal = distanceO.getDistanceTour(tour,D) 
    bestTour=tour
    for i in range(0,len(tour)-1):
        for j in range(i+1, len(tour)): 
            newTour = copy.deepcopy(tour)
            newTour[j] = copy.deepcopy(tour[i])
            newTour[i]= copy.deepcopy(tour[j])
            distanceSwap = distanceO.getDistanceTour(newTour,D) 
            if distanceSwap < distanceOriginal:
                distanceOriginal=distanceSwap
                bestTour=newTour
    return bestTour
            
def interSwap(child,distanceM):
    D=child[0]
    childTours = child[1:-1]
    for i,tour in enumerate(childTours):
        if len(tour)<=2:
            continue
        childTours[i] = interSwapTour(tour,D,distanceM)
    child = [D] + childTours + [0]
    return child     

def removeLongEdge(tour, distanceO,D):
    for i in range(1,len(tour)):
        dToDepot1 = distanceO.distanceToDepot(tour[i-1], D)
        dToDepot2 = distanceO.distanceToDepot(tour[i], D)
        dOriginal= distanceO.EuclideanNode(tour[i-1],tour[i])
        if dToDepot1 + dToDepot2 < dOriginal:
            return [tour[:i], tour[i:]]
    return [tour]

def removeLongEdges(child,distanceO):
    D=child[0]
    oldTours = child[1:-1]
    newTours = []
    for tour in oldTours:
        newTours += removeLongEdge(tour, distanceO,D)
    child = [D]+ newTours + [0]
    return child

def twoSwapNodes(tour1,tour2,D,distanceO):
    originalDistance = distanceO.getDistanceTour(tour1,D)+distanceO.getDistanceTour(tour2,D)
    bestTour1,bestTour2=tour1,tour2
    bestDelta = 0
    tour1=list(tour1)
    tour2=list(tour2)
    for i,node1 in enumerate(tour1):
        for j,node2 in enumerate(tour2):
            newTour1 = copy.deepcopy(tour1)
            newTour2 = copy.deepcopy(tour2)
            newTour1[i],newTour2[j]=node2,node1
            distance1 = distanceO.getDistanceTour(newTour1,D)
            distanceR = distanceO.getDistanceTour(newTour2,D)
            delta = distance1+distanceR-originalDistance
            if delta < bestDelta:
                assert newTour1 is not tour1 or newTour2 is not tour2, "1 tour1n: %s \n tour1: %s \n tour2N: %s \n, tour2 %s \n distanceO: %s, distance new: %s"%(newTour1,tour1,newTour2,tour2,bestDelta, distance1+distanceR)
                bestDelta=delta
                bestTour1,bestTour2=newTour1,newTour2       
    return bestTour1,bestTour2,bestDelta 

def twoSwapTour(tour1, D, tours, index,distanceM):
    bestDelta = 0
    index2=-1
    bestTour1=tour1
    bestTour2=None
    for i in range(index+1,len(tours)):
        tour2=tours[i]
        if tour2 is tour1 or toursTooFarRemoved(tour1,tour2,distanceM, D):
            continue
        tour1N, tour2N,delta = twoSwapNodes(tour1,tour2,D,distanceM)
        if delta < bestDelta:
            assert tour2N is not tour2
            bestTour1=tour1N
            bestTour2=tour2N
            index2=i
            bestDelta=delta
    if bestTour2:
        tours[index]=bestTour1
        tours[index2]=bestTour2
    return tours    
    

def twoSwap(child,distanceO):
    childTours = child[1:-1]
    D=child[0]
    distance = child[-1]
    for i in range(0,len(childTours)):
        tour = childTours[i]
        childTours = twoSwapTour(tour, D, childTours, i,distanceO)    
    newChild = [D] + childTours + [distance]
    return newChild

def getBestInsertedDistance(tour, node,D,distanceO):
    tour=list(tour)
    distance=float('inf')
    bestTour=None
    for i in range(0,len(tour)+1):
        tour=copy.deepcopy(tour)
        tour.insert(i,node)
        d = distanceO.getDistanceTour(tour, D)
        if d<distance:
            bestTour=tour
    return bestTour,distance
        
def addNewNodeFromTourToTour(tour1, tour2,D,distanceO):
    """
    A node from tour1 is inserted in the best place of tour2
    If this results in an improvement, the new tours are returned
    """
    bestTour1,bestTour2=tour1,tour2
    distance = distanceO.getDistanceTour(tour1, D) + distanceO.getDistanceTour(tour2, D)
    bestDelta=0
    for node in tour1:
        newTour1=copy.deepcopy(tour1)
        newTour1.remove(node)
        distanceNew1 = distanceO.getDistanceTour(newTour1,D)
        newTour2,distanceNew2 = getBestInsertedDistance(tour2, node,D, distanceO)
        delta = distanceNew1+distanceNew2 - distance
        if delta < bestDelta:
            bestDelta = delta
            bestTour1,bestTour2=newTour1,newTour2
            
    return bestTour1,bestTour2,bestDelta

def addNodeFromTourToTour(tour1,tours,D,index,distanceO):
    bestDelta=0
    index2=-1
    bestTour1,bestTour2=tour1,None
    for i in range(0,len(tours)):
        tour2 = tours[i]
        if (tour1 is tour2 or toursTooFarRemoved(tour1, tour2, distanceO, D) 
            or len(tour2)==K):
            continue
        tour1N, tour2N, delta = addNewNodeFromTourToTour(tour1, tour2,D,distanceO)
        if  delta < bestDelta:
            assert tour1N is not tours,'\n %s\n %s \n %s, \n %s' %(tour1N, tour1,bestDelta,delta ) 
            assert len([node for node in tour1N if node in tour2N])==0,'%s, %s' %(tour1N,tour2N)
            index2=i
            bestTour1,bestTour2=tour1N,tour2N
            bestDelta=delta
    
    if bestTour2:
        tours[index]=bestTour1
        tours[index2]=bestTour2
    return tours

def addOne(child,distanceO):
    D = child[0]
    tours = child[1:-1]
    for i in range(0,len(tours)):
        if len(tours[i])==0:
            continue
        tours = addNodeFromTourToTour(tours[i],tours,D,i,distanceO)
    distance = child[-1]
    filledTours = [tour for tour in tours if len(tour)>0]
    child = [D]+ filledTours + [distance]
    return child

def improve(children,distanceO,tabu=None):
    print "Children are educated"
    for i,child in enumerate(children):
        distance=child[-1]
        
        children[i] = removeLongEdges(child,distanceO)
        children[i] = interSwap(children[i],distanceO)
        children[i] = twoSwap(children[i],distanceO) 
        children[i]= addOne(children[i],distanceO) 
        
        children[i] = replaceDepotsFromChild(children[i],distanceO,tabu)
         
        if distance < children[i][-1]: print('No distance improvement found, %s vs %s' %(children[i][-1],distance))
        else: print('Distance improvement found,  %s vs %s' %(children[i][-1],distance)) 
        for tour in children[i][1:-1]:
            assert len(tour)<=K,'len tour:%s, K: %s' %(len(tour),K)
            assert len(tour)>0
    print "All children are educated"
    return children

def missingGenes(tours1,tours2, genesFrom1, genesFrom2,N, lenthOfToursToAdd):       
    extraTours = tours1
    if len(genesFrom2)<len(genesFrom1):
        extraTours = tours2

    extraGenes=[]
    nodesFrom1 = [i for tour in genesFrom1 for i in tour]
    nodesFrom2 = [i for tour in genesFrom2 for i in tour]
    missingNodes =[node for node in N if not (node in nodesFrom1 or node in nodesFrom2) ]
    for tour in extraTours:
        if nofNodesListHaveInCommon(tour,missingNodes) == lenthOfToursToAdd:
            extra = [node for node in tour if node in missingNodes]
            if len(extra)>0:
                extraGenes.append(extra)   
    return extraGenes

def createChild(genesFrom1, genesFrom2,N,b, distanceO):
    nodesFrom1 = [i for tour in genesFrom1 for i in tour]
    nodesFrom2 = [i for tour in genesFrom2 for i in tour]
    missingNodes =[node for node in N if not (node in nodesFrom1 or node in nodesFrom2) ]
    
    nodesConnectedToDepot = [[tour[0], tour[-1]] for tour in genesFrom1]
    nodesConnectedToDepot+= [[tour[0], tour[-1]] for tour in genesFrom2]
    nodesConnectedToDepot = [node for sublist in nodesConnectedToDepot for node in sublist]
    solution = depots.quickDepots(nodesConnectedToDepot,2);
    D = [solution[1],solution[2]]
    newTours = []
    if len(missingNodes) >0:
        greedySolution = greedy.doExtraGreedy(missingNodes,D,distanceO,K)
        newTours, __ = extraGreedySolutionToTours(greedySolution,N)
        
    child = [[0,0]]+genesFrom1 + genesFrom2 + newTours +[0]
    distance, D = getChildDistanceAndDepots(child,distanceO)
    return  [D]+ child[1:-1] + [distance]

def addNewGene(tour1, tours2,forbiddenToursFor2,genesFrom1):
    setToAdd = [node for node in tour1]
    newForbidden = [tour for tour in tours2 if nodesListShareElement(tour, setToAdd)]
    forbiddenToursFor2 += newForbidden
    assert len(tour1)>0
    genesFrom1.append(tour1)
    return genesFrom1, forbiddenToursFor2

def addNewGeneNotFromForbidden(tour, genesFrom1,tours2,forbiddenToursFrom2):         
    assert len(tour)>0
    genesFrom1.append(tour)    
    newForbidden = [tour2 for tour2 in tours2 if nodesListShareElement(tour, tour2) ] 
    forbiddenToursFrom2 += newForbidden
    return genesFrom1,forbiddenToursFrom2
    
def obtainGenesFrom1(genesFrom1, tours1, tours2, forbiddenToursFrom2): 
    for tour in tours1:
        if len(genesFrom1) > (len(tours1)/2.0) or len(forbiddenToursFrom2) > (3.0*len(tours2)/4.0)  :
            return  genesFrom1,forbiddenToursFrom2 
        if tour in genesFrom1:
            continue
        elif toursAndNodeListShareElement(forbiddenToursFrom2,tour):
            genesFrom1, forbiddenToursFrom2 = addNewGene(tour, tours2,forbiddenToursFrom2,genesFrom1)
        elif toursContainExactSameNodes(genesFrom1,forbiddenToursFrom2):
            genesFrom1,forbiddenToursFrom2 = addNewGeneNotFromForbidden(tour, genesFrom1,tours2,forbiddenToursFrom2)
    return  genesFrom1,forbiddenToursFrom2   

def reproduceOnce(parent1,parent2,N,b,distanceO):###children are too similar, maybe scramble parents??
    genesFrom1,forbiddenToursFrom2 =[],[]
    tours1 = parent1[1:-1]
    tours2 = parent2[1:-1]   
    
    amountPreviousGenesFrom1=0
    genesFrom1,forbiddenToursFrom2 = addNewGeneNotFromForbidden(tours1[0], genesFrom1,tours2,forbiddenToursFrom2)
    tours1Copy = [tour for tour in tours1 if tour not in genesFrom1]    
    while len(genesFrom1)> amountPreviousGenesFrom1 and len(genesFrom1) < len(tours1)/2.0:
        assert len([tour for tour in tours1Copy if tour in genesFrom1])==0
        amountPreviousGenesFrom1 = len(genesFrom1)
        genesFrom1, forbiddenToursFrom2 = obtainGenesFrom1(genesFrom1,tours1Copy,tours2,forbiddenToursFrom2)
        tours1Copy = [tour for tour in tours1Copy if tour not in genesFrom1]    
    
    genesFrom2 = [tour for tour in tours2 if tour not in forbiddenToursFrom2]
    missingGenes(tours1,tours2, genesFrom1, genesFrom2,N,len(tours2[0])-1)
    child=createChild(genesFrom1, genesFrom2,N,b,distanceO)
    checkNoDoubleNodeTours(child[1:-1])
    return child

def reproduce(parents,N,b,distanceO):
    parents = shuffleGeneration(parents)
    print "Reproduction phase is started, %s" %(len(parents))
    children=[]
    for j in range(0,len(parents)):
        for i in range(j+1,len(parents)):
            if not deltaInSolution(parents[i],parents[j],b):
                print "combo removed1"
                continue
            if parentsAreTooSimilar(parents[i],parents[j]):
                print "combo removed2"
                continue
            children.append(reproduceOnce(parents[i],parents[j],N,b,distanceO))
    print "Reproduction phase is ended"
    return children

def extractBestSoFar(children, parents):
    global amountOfParents, minimumAmountOfNewParents,BestSolution
    temp= children+parents
    totalScores=sorted([person[-1] for person in temp])[:amountOfParents]
    childrenScores=sorted([child[-1] for child in children])[:minimumAmountOfNewParents]
    
    if len(set(childrenScores).intersection(set(totalScores))) >= minimumAmountOfNewParents:
        newChilds1 =  [ child for child in children if child[-1] in totalScores]
        newChilds2= [ parent for parent in parents if parent[-1] in totalScores]
        newParents=newChilds1+newChilds2
        newParents=sorted(newParents, key=lambda child: child[-1])    
    else:
        oldParents = amountOfParents-minimumAmountOfNewParents
        parentScores = sorted([person[-1] for person in parents])[:oldParents]
        newChilds1 =  [child for child in children if child[-1] in childrenScores]
        newChilds2= [ parent for parent in parents if parent[-1] in parentScores]
        newParents=newChilds1+newChilds2
        newParents=sorted(newParents, key=lambda child: child[-1])      
    
    if BestSolution[-1]>newParents[0][-1]:
        BestSolution=newParents[0]
    setBestAndAvParents(newParents)
    return newParents

def setBestAndAvParents(population):
    global BestSolution
    population=sorted(population, key=lambda child: child[-1]) 
    BestParentsY.append(population[0][-1])
    z = len(BestParentsX)
    BestParentsX.append(z)
    scores=[new[-1] for new in population]
    AvParentsY.append(reduce(lambda x, y: x + y, scores) / len(scores))
    AvParentsX.append(z) 
    if population[0][-1] <BestSolution[-1]:
        BestSolution = population[0]
    

def doGenetic(N,h,b,distanceO,maxNodeTours):
    """
    Parents is a list of solutions. A solution looks like:
    [[float , float],[a,b,c],[e,r,t] ...[g,h,j], float]
    where [float,float] denotes the place of the two depots
    [a,b,c] is a list containing the nodes in the first tour.
    The node 'a' is denoted the same as in N: [number, xcoordinate, ycoordinate]. 
    This represents the tour (d_i, a, b, c, d_j) with d_i and d_j being the nodes closest to a and c respectively
    The float at the end represents to total length of the tour. 
    """
    global similarTours, K, BestParentsX, BestParentsY, AvParentsX, AvParentsY,BestSolution
    BestSolution=[float('inf')]
    BestParentsX=[]
    BestParentsY=[]
    AvParentsX=[]
    AvParentsY=[]   
    
    K=maxNodeTours

    parents = initialize(N,b,distanceO)
    similarTours=1
    parents = improve(parents,distanceO)
    children = reproduce(parents,N,b,distanceO)
    
    similarTours=3.0/4.0
    children = improve(children,distanceO)
    parents = extractBestSoFar(children,parents)

    for i in xrange(generations):
        #parents = selection(children)
        children = reproduce(parents,N,b,distanceO)
        children = improve(children,distanceO)
        parents = improve(parents,distanceO) #improve above or below extractBestSoFar?
        parents = extractBestSoFar(children,parents)
    
    extractBestSoFar(children,parents)
    showSolution()
    global BestSolution
    return BestSolution

def showSolution():
    plt.scatter(BestParentsX, BestParentsY)
    plt.show()
    plt.clf()
    plt.scatter(AvParentsX, AvParentsY)
    plt.show()
    plt.clf()
    
    
def doLocalSearch(N,h,b,distanceO,maxNodeTours):
    """
    Algorithm which only executes the local search procedure
    """
    global K, BestParentsX, BestParentsY, AvParentsX, AvParentsY,BestSolution
    BestSolution=[float('inf')]
    BestParentsX=[]
    BestParentsY=[]
    AvParentsX=[]
    AvParentsY=[]   
    tabu=[]
    
    K=maxNodeTours
    K=maxNodeTours

    parents = initialize(N,b,distanceO)
    parents = extractBestSoFar(parents,parents)

    parents = improve(parents,distanceO,tabu)
    setBestAndAvParents(parents)

    for i in xrange(2):
        parents = improve(parents,distanceO,tabu)
        setBestAndAvParents(parents)
    
    extractBestSoFar(parents,parents)
    showSolution()
    return BestSolution
