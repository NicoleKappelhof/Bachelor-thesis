import copy, math
import itertools
gamma = 0.002   
eta=0.05

def bestDepotLocation(N,K,numberDepots):
    minAmountOfTours = int(float(len(N))/K + 1)
    combinations=list(itertools.combinations(N, minAmountOfTours))
    distance=float('inf')
    solution=None
    for combo in combinations:
        combo=list(combo)
        s=exactDepots(combo,numberDepots)
        d=s[0]
        if d<distance:
            distance=d
            solution=s
    return exactDepotSolutionToRegSolution(solution)
    
#[([1, 3, 4], [2, 3, 4], [3, 2, 3]), ([1, 3, 4], [2, 3, 4], [4, 6, 7]), ([1, 3, 4], [3, 2, 3], [4, 6, 7]), ([2, 3, 4], [3, 2, 3], [4, 6, 7])]

def exactDepotSolutionToRegSolution(s):
    D = [s[1],s[2]]
    d=s[-1]
    S1 = s[3]
    S2=s[4]
    tours=[]
    for node in S1:
        tours.append(['depot1', node[0]])
    for node in S2:
        tours.append(['depot2', node[0]])
    solution =[D] + tours + [d]
    return solution 

def exactDepots(N, numberDepots, x1=None,x2=None):
    assert numberDepots==2, "number of depots not implemented" 
    NOrder = sorted(N, key=lambda node:node[1])
    S1 =[]
    S2=copy.deepcopy(NOrder)
    d,x1,x2 = findMin(S1,S2,x1,x2)
    minimum = d
    solution = [d, x1,x2,S1,S2]
    for i in range(0,len(NOrder)):
        S1 += [S2[0]]
        S2 = S2[1:]
        d,x1,x2 = findMin(S1,S2,x1,x2)
        if d < minimum:
            minimum =d
            copyS1 = copy.deepcopy(S1)
            copyS2 = copy.deepcopy(S2)
            solution = [minimum, x1,x2,copyS1,copyS2]
    return solution


def quickDepots(N, numberDepots):
    global eta
    global gamma
    etaOri = eta
    gammaOri=gamma
    eta=0.1
    gamma=0.012
    assert numberDepots==2, "number of depots not implemented"
    start=int((1.0/3.0)*len(N))
    end=int((2.0/3.0)*len(N))
    NOrder = sorted(N, key=lambda node:node[1])
    S1=NOrder[:start]
    S2=NOrder[start:]
    d,x1,x2 = findMin(S1,S2)
    min = d
    solution = [d, x1,x2,S1,S2]   
    for n in NOrder[start:end]: 
        S1 += [n]
        S2.remove(n)
        d,x1,x2 = findMin(S1,S2)
        if d < min:
            min =d
            copyS1 = S1
            copyS2 = S2
            solution = [d, x1,x2,copyS1,copyS2]
    print "The depots should be placed at %s and %s" %(solution[1],solution[2])
    print "They service %s and %s with a total distance of %s" %(solution[3],solution[4], solution[0])
    eta=etaOri
    gamma=gammaOri
    return solution


def findMin(S1,S2,x1=None,x2=None):
    if not x1 or x2:
        x1=1.4
        x2=9.8
    if len(S1)>0:
        gamma = 1/float(len(S1))
    while len(S1)>0 and (derivative(x1,S1)> eta or derivative(x1,S1)< -eta):
        x1=x1 -  gamma * derivative(x1,S1) 
        
    if len(S2)>0:
        gamma = 1/float(len(S2))
    while len(S2)>0 and (derivative(x2,S2)> eta or derivative(x2,S2)< -eta):
        x2= x2 - gamma *  derivative(x2,S2); 
        
    y1 = functionValue(x1,S1)
    y2= functionValue(x2,S2)
    return y1+y2,x1,x2

def functionValue(x, S):
    y = 0
    for node in S:
        xn = node[1]
        yn=node[2]
        y += math.sqrt(math.pow((xn-x),2) + math.pow(yn,2))
    return y
    
    
def derivative(x,S):
    derivative = 0.0
    for node in S:
        xn = float(node[1])
        yn=float(node[2])
        derivative += -float(xn-x)/math.sqrt(math.pow((xn-x),2) + math.pow(yn,2))
    return derivative

def secondDerivative(x,S):
    secondDerivative = 0.0
    for node in S:
        xn = float(node[1])
        yn=float(node[2])
        secondDerivative += (math.pow(yn,2))/((math.pow((xn-x),2) + math.pow(yn,2))*math.sqrt(math.pow((xn-x),2) + math.pow(yn,2)))
    return secondDerivative