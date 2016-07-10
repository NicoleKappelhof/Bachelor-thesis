import itertools

N=[]
S =[]

def bestPossibleTrips(N,K,distanceO):
    allCombinationOfTrips=tripRecursion(N,K)
    distance=float('inf')
    bestCombination=None
    for combination in allCombinationOfTrips:
        d, tours = combinationDistance(combination,distanceO)
        if d<distance:
            bestCombination = tours
    return [[0,0]]+ bestCombination + [distance]
        
def tripRecursion(N,K):#Ik weet niet of dit werkt
    if K>=len(N):
        return [N]
    tripList=[]
    nodeCombo=N[0]
    combinations=list(itertools.combinations(N[1:], K-1))
    for combination in combinations:
        combination=list(combination) + nodeCombo
        NWithout =[node for node in N if node not in combination]
        combosOfcombo = bestPossibleTrips(NWithout,K)
        for comboOfcombo in combosOfcombo:
            trip = comboOfcombo + [combination]
            tripList.append(trip)
    return tripList

def combinationDistance(combination,distanceO):
    distance=0
    tours=[]
    for trip in combination:
        d, bestTrip =distanceO.minimumTourDistanceNoDepots(trip)
        distance+=d
        tours.append(bestTrip)
    return distance, tours

def combinationsDistanceWithDepots(combination,D, distanceO):
    distance=0
    for trip in combination:
        distance+= distanceO.minimumTourDistance(trip, D)[0]
    return distance
    
def exactAlgorithm(Nodes,D,K, distanceO):
    global S
    global N
    N = Nodes
    bestCombination=None
    allCombinationOfTrips=tripRecursion(N,K)
    distance=float('inf')
    for combination in allCombinationOfTrips:
        d = combinationsDistanceWithDepots(combination,distanceO,D)
        if d<distance:
            bestCombination = combination
    
    assert bestCombination
    return bestCombination
    print "The best solution is: %s " %(bestCombination)

    
    
    


