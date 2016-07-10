"""
parameters:
ImproveStartTime or ImproveStartGap can be set to a specific value. If this value is reached,
the MIP will no longer continue to prove optimality but rather improve the bound

MIPFocus, set to 1 for focus on feasible solution, set to 2 for focus on optimal solution, 
set to 3 for improving the bound

MIPGap, the solver will terminate when the relative gap between the lower and upper
objective bound is less than MIPGap times the upper bound. Default is 1e-4, maximum is 1.
    
Method, set to 0 for primal simplex, 1 for dual simplex, 2 for the parallel barrier and
3 for the concurrent (only for root node)

NodeMethod, set to 0 for primal simplex, 1 for dual simplex (default), 2 for the parallel barrier and
3 for the concurrent (only for non-root nodes)

NodefileStart, if the amount of memory in GB used to store nodes exceeds the specified 
parameter value, nodes are written to disk (a setting of 0.5 is recommended)

Sifting, if you have many more variables than restrictions (100 to 1 or more). The sifting parameter 
is also set by Gurobi according to this ratio, so you prob wont have to touch it 

Thread, each thread in parallel MIP requires a copy of the model. Reducing the
Threads parameter can sometimes significantly reduce memory usage. By default, concurrent and barrier
will use all available cores in your machine.

Heuristics, controls the fraction of runtime spent on feasibility heuristics. Increasing the parameter can 
lead to more and better feasible solutions, but it will also reduce the rate of progress in the best bound.
Default is 0.05, maximum is 1.

Cuts, (finer grain: FlowCoverCuts, MIRCuts) can be set to Aggressive (2), Conservative (1), Automatic (-1), 
or None (0). The more specific parameters override the more general. 
Very easy models can sometimes benefit from turning cuts off, while extremely difficult models can benefit
from turning them to their Aggressive setting.

Presolve, the aggressiveness level of presolve. Options are Aggressive (2), Conservative (1), Automatic (-1),
or None (0). More aggressive application of presolve takes more time, but can sometimes lead to a significantly
tighter model.

PreDepRow, the presolve dependent row reduction, which eliminates linearly dependent constraints from the constraint 
matrix. The default setting (-1) applies the reduction to continuous models but not to MIP models. Setting 0 turns 
the reduction off for all models. Setting 1 turns it on for all models.

Aggregate, controls whether presolve performs constraint aggregation, 1 is default (and maximum), 0 is minimum

SubMIPNodes, maximum number of nodes explored by the RINS heuristic. Exploring more nodes can produce
better solutions, but it generally takes longer. Range 0-MAXINT, 500 is default. 

TuneTimeLimit, maximal total tuning runtime of tuning (in seconds). The default setting (-1) chooses a time limit automatically.

"""

from gurobipy import *
import itertools
from ThreeOrMore import Distances


#N=[[1,1,1],[2,7,1],[3,1,1],[4,5,1],[5,9,1],[6,9,1],[7,9,1],[8,9,1]]
#distanceObject = Distances.distances(N)


def lpK_sets(N, distanceO,D,K):
    """
    Set LP
    """ 
    # Create a new model
    m = Model("lpK_sets")
    
    # Create variables of combinations of nodes
    combinations=[]
    removedCombinations=[]
    for i in range(1,K+1):
        combinations += list(itertools.combinations(N, i))
    
    variables=[]
    for i,combo in enumerate(combinations):
        d, tour=distanceO.minimumTourDistanceWithCheck(list(combo),D)
        if d is False and tour is None:
            removedCombinations.append(combo)
            continue
        comboints=[node[0] for node in tour]
        variablename = 'combo' + '-'.join(map(str,comboints))
        variables.append(m.addVar(vtype=GRB.BINARY,obj=d, name=variablename))
    
    # Integrate new variables
    m.modelSense = GRB.MINIMIZE
    m.update()

    # Every node needs to be serviced once and only once
    combinations = [combo for combo in combinations if combo not in removedCombinations]
    for node in N:
        containingNode = [i for i,combo in enumerate(combinations) if node in combo]
        m.addConstr(quicksum(variables[i] for i in containingNode) == 1,"service constraint %s" %node[0])    
                       
    m.params.tuneResults = 1
    m.params.TuneTimeLimit=max(2*(len(N)-20),0)                 
    m.tune()
       
    if m.tuneResultCount > 0:
        # Load the best tuned parameters into the model
        m.getTuneResult(0)                       
                       
    m.params.NodefileStart = 0.8
    #m.params.LogToConsole=1
    m.optimize()
    m.printAttr('X')
    print('Obj: %g' % m.objVal)
    
    distance = 0
    tours=[]
    for v in m.getVars():
        if v.x==0:
            continue
        assert v.varName.startswith('combo')
        nodes = v.varName.split('combo')[1]
        tour=[]
        for node in nodes.split('-'):
            tour.append(N[int(node)-1])           
        d=distanceO.getDistanceTour(tour,D)
        distance+=d
        tours.append(tour)
    assert distance >= m.objVal-0.01 and distance <= m.objVal + 0.01, ' %s, %s' %(distance, m.objVal)
    solution = [D] +tours +[distance]

    return solution



def lpK_edges(N, distanceO,D,K):
    """
    The obvious LP
    """
    # Create a new model
    m = Model("lpK_edges")

    # Create variables of edges between nodes
    variables=[]
    for node1 in N:
        variables.append([])
        i=node1[0]
        for node2 in [node for node in N if node[0]>i]:
            j=node2[0]
            variable = 'x' + str(node1[0])+'-'+str(node2[0])
            variables[i-1].append(m.addVar(vtype=GRB.BINARY,obj=distanceO.EuclideanNode(node1,node2), name=variable))
    
    # Create variables of edges between depots and nodes
    toDepotVariables=[]
    fromDepotVariables=[]
    for node in N:
        i=node[0]-1
        toDepotVariables.append([])
        fromDepotVariables.append([])
        depotString, d = distanceO.closestDepotToString(D[0],D[1],node)
        variable = depotString + 'n' + str(node[0])
        fromDepotVariables[i].append(m.addVar(vtype=GRB.BINARY,obj=d, name=variable))
        variable = 'n' +str(node[0])+depotString
        toDepotVariables[i].append(m.addVar(vtype=GRB.BINARY,obj=d, name=variable))
    
    # Integrate new variables
    m.modelSense = GRB.MINIMIZE
    m.update()

    # Add in-out constraints nodes      
    for i in range(0,len(N)):
        iVariables=variables[i] + [variables[k][i-k-1] for k in range(0,i)]
        m.addConstr(sum(iVariables) + quicksum(fromDepotVariables[i])+ quicksum(toDepotVariables[i])== 2,"in-out %s" % str(i))

    # Add tour constraints
    combinations = list(itertools.combinations(range(0,len(N)),K))
    for combo in combinations:
        tourConstraints(m,combo,variables,K,toDepotVariables,fromDepotVariables)
          
    m.params.tuneResults = 1
    m.params.TuneTimeLimit=3*len(N)
    m.params.LogToConsole=0  
    m.tune()
#      
    if m.tuneResultCount > 0:
        # Load the best tuned parameters into the model
        m.getTuneResult(0)    

    m.params.MIPGap=0.05
    
    

    
    m.optimize()    
    m.printAttr('X')
    
    edges=[]
    for v in m.getVars():
        if v.x==0:
            continue
        if v.varName.startswith('x'):
            nodes = v.varName.split('x')[1]
            nodes=nodes.split('-')
            assert len(nodes)==2
            edges.append([nodes[0],nodes[1]])
        elif v.varName.startswith('depot'):
            e = v.varName.split('depot')[1]
            node = e.split('n')[1]
            depot = e.split('n')[0]
            s = 'depot' + str(depot)
            edges.append([node, s]) 
        else:
            assert v.varName.startswith('n') 
            e = v.varName.split('n')[1]
            depot = e.split('depot')[1]
            node = e.split('depot')[0]
            s = 'depot' + str(depot)
            edges.append([node, s])     
    print('Obj: %g' % m.objVal)
    solution = [D] + edges + [m.objVal]   
    return solution


def lpK_flow(N, distanceO,D,K):
    """
    The flow LP
    """
    # Create a new model
    m = Model("lpK_flow")

    # Create variables for edges between nodes and depots        
    fromDepot = []
    toDepot=[]
    for node in N:
        i = node[0]-1
        depotString, d = distanceO.closestDepotToString(D[0],D[1], node)
        variable = depotString+'-'+str(node[0])
        fromDepot.append([])
        fromDepot[i].append(m.addVar(vtype=GRB.BINARY,obj=d, name=variable))     
        variable = str(node[0])+ '-' + depotString   
        toDepot.append([]) 
        toDepot[i].append(m.addVar(vtype=GRB.BINARY,obj=d, name=variable))
        
    # Create variables for edges between nodes and nodes
    levels=[0]*K
    for i in range(1,K):
        levels[i-1] = levelK(N,distanceO,m,D, i)
     
     
    levelLast=[]
    for node in N:
        i=node[0]-1
        levelLast.append([])
        variable = 'lvl'+str(K)+ '-' + 'ending' '-' + str(node[0])
        levelLast[i].append(m.addVar(vtype=GRB.BINARY,obj=0, name=variable))  
    levels[K-1] =  levelLast   
    
    # Create variables for edges between nodes and service nodes
    serviced=[0]*K
    for i in range(1,K+1):
        serviced[i-1]= servicedK(N,distanceO,m,i)
            
    # Integrate new variables
    m.modelSense = GRB.MINIMIZE
    m.update()
    
    #in-out constraints
    for node in N:
        i=node[0]-1
        for j in range(0,K):
            level = levels[j]
            m.addConstr(quicksum(level[i])<= 1, 'out level %s- %s'%(j, str(i+1)))
        
        m.addConstr(2*fromDepot[i][0]-(quicksum(levels[0][i]))-quicksum(serviced[0][i])==0, 'in-out level 1- %s'%str(i+1))
        for k in range(1,K):
            m.addConstr(2*quicksum(levels[k-1][j][i] for j in range(0,len(N)))-quicksum(levels[k][i])-serviced[k][i][0] - levels[k-1][i][i]==0, 'in-out level %s-%s'%(k,str(i+1)))
        m.addConstr(quicksum(levels[-1][i])-toDepot[i][0]==0, 'to depot-%s'%str(i+1))
    
    
    # Service constraints
    for node in N:
        i=node[0]-1
        m.addConstr(quicksum(serviced[k][i][0] for k in range(0,K))==1, 'service Constraint- %s'%str(i))
    

    m.params.tuneResults = 1
    m.params.TuneTimeLimit=max(2*(len(N)-20),0)
     
    m.tune()
      
    if m.tuneResultCount > 0:
        # Load the best tuned parameters into the model
        m.getTuneResult(0)    
     
    
    m.params.LogToConsole=1
    
      
    
    m.optimize()    
    m.display()
    m.printAttr('X')
    
    edges=[]
    for v in m.getVars():
        if v.x==0 or v.varName.startswith('serviced'):
            continue
        if v.varName.startswith('lvl'):
            nodes = v.varName.split('-')[1:]
            if nodes[0]=='ending':
                continue
            assert len(nodes)==2
            edges.append([nodes[0],nodes[1]])
        elif v.varName.startswith('depot'):
            e = v.varName.split('depot')[1]
            node = e.split('-')[1]
            depot = e.split('-')[0]
            s = 'depot' + str(depot)
            edges.append([node, s]) 
        else:
            e = v.varName.split('-')
            assert e[1].startswith('depot'), '%s' %e
            node = e[0] 
            depot = e[1].split('depot')[1]
            s = 'depot' + str(depot)
            edges.append([node, s])     
    print('Obj: %g' % m.objVal)
    solution = [D] + edges + [m.objVal]   
    return solution

def levelK(N,distanceO,m,D, level):
    levelK=[]
    for node1 in N:
        i=node1[0]-1
        levelK.append([])
        for node2 in N:
            d=distanceO.EuclideanNodeWithCheck(node1,node2,D)
            if d <0:
                levelK[i].append(0)
                continue
            variable = 'lvl' +str(level) + '-' + str(node1[0])+ '-' + str(node2[0])
            levelK[i].append(m.addVar(vtype=GRB.BINARY,obj=d, name=variable))
    return levelK

def servicedK(N,distanceO,m,level):
    servicedK=[]
    for node in N:
        i=node[0]-1
        servicedK.append([])
        variable = 'serviced'+str(level)+ '-' + str(node[0])
        servicedK[i].append(m.addVar(vtype=GRB.BINARY,obj=0, name=variable))
    return servicedK


def singleTourConstraint(m,combination,variables,K,toDepotVariables,fromDepotVariables):
    constraints=[]
    variableList = list(combination) 
    depotTours=[]
    strartNode=variableList[0]
    endNode=variableList[-1]
    depotTours.append(quicksum(fromDepotVariables[strartNode]))
    depotTours.append(quicksum(toDepotVariables[endNode]))
    for i in range(1, len(variableList)):
        j = variableList[i-1]
        k=variableList[i]
        if k < j:
            constraints.append(variables[k][j-k-1])
        else:
            constraints.append(variables[j][k-j-1])
    tourString = '-'.join(map(str,[variable+1 for variable in variableList]))
    m.addConstr(2*quicksum(constraints) - (quicksum(depotTours))<= 2*(K-2),"New tour constraint %s" %(tourString))    
    return 

def singleTourConstraint2(m,combination,variables,K,depotVariablesIn,depotVariablesOut):
    constraints=[]
    variableList = list(combination) 
    strartNode=variableList[0]
    endNode=variableList[-1]
    depotTours=[]
    depotTours += depotVariablesOut[strartNode-1]
    depotTours+= depotVariablesIn[endNode-1]
    for i in range(1, len(variableList)):
        j = variableList[i-1]
        k=variableList[i]
        constraints.append(variables[j-1][k-1])
    tourString = '-'.join(map(str,[variable+1 for variable in variableList]))
    m.addConstr(sum(constraints) - (sum(depotTours))<= K -3,"New tour constraint %s" %(tourString))    
    return 


def tourConstraints(m,combination,variables,K,toDepotVariables,fromDepotVariables):
    variableList = list(combination) 
    permutations = list(itertools.permutations(variableList, len(variableList)))
    tours=[list(permutations[0])]
    reversedTours = [list(reversed(tours[0]))]
    for p in permutations[1:]:
        p=list(p)
        if p in reversedTours:
            continue
        tours.append(p)
        reversedTours.append(list(reversed(p)))
    for tour in tours:
        singleTourConstraint(m,tour,variables,K,toDepotVariables,fromDepotVariables)
    return 

def tourConstraints2(m,combination,variables,K,depotVariablesIn,depotVariablesOut):
    variableList = list(combination) 
    permutations = list(itertools.permutations(variableList, len(variableList)))
    tours=[list(permutations[0])]
    for p in permutations[1:]:
        p=list(p)
        tours.append(p)
    for tour in tours:
        singleTourConstraint2(m,tour,variables,K,depotVariablesIn,depotVariablesOut)
    return 
    
def variblesOfEdgesToDepots(set1, toDepotVariables,fromDepotVariables):
    variables=[]
    for i in set1:
        variables.append(toDepotVariables[i][0])
        variables.append(toDepotVariables[i][1])
        variables.append(fromDepotVariables[0][i])
        variables.append(fromDepotVariables[1][i])
        
    return variables

def lpK_DoubleEdges(N, distanceO,D,K):
    """
    Variant of the Obvious LP where the edges are directed
    """
    # Create a new model
    m = Model("lpK_DoubleEdges")

    # Create variables of edges between nodes
    variables=[]
    for node1 in N:
        variables.append([])
        i=node1[0]
        for node2 in N:
            j=node2[0]
            variable = 'x' + str(node1[0])+'-'+str(node2[0])
            variables[i-1].append(m.addVar(vtype=GRB.BINARY,obj=distanceO.EuclideanNode(node1,node2), name=variable)) 

    # Create variables of edges between depots and nodes
    depotVariablesOut=[]
    for node1 in N:
        depotVariablesOut.append([])
        i=node1[0]
        for j,depot in enumerate(D):
            variable = 'd' + str(j+1)+'n'+str(node1[0])
            depotVariablesOut[i-1].append(m.addVar(vtype=GRB.BINARY,obj=distanceO.Euclidean(node1[1],node1[2],depot,0), name=variable))
   
    # Create variables of edges between depots and nodes
    depotVariablesIn=[]
    for node1 in N:
        depotVariablesIn.append([])
        i=node1[0]
        for j,depot in enumerate(D):
            variable = 'n' +str(node1[0])+'d'+ str(j+1)
            depotVariablesIn[i-1].append(m.addVar(vtype=GRB.BINARY,obj=distanceO.Euclidean(node1[1],node1[2],depot,0), name=variable))


    # Integrate new variables
    m.modelSense = GRB.MINIMIZE
    m.update()

    combinations = list(itertools.combinations(range(0,len(N)),K))
    for combo in combinations:
        tourConstraints2(m,combo,variables,K,depotVariablesIn,depotVariablesOut)
   
    # Add in-out constraints nodes      
    for i in range(0,len(N)):
        m.addConstr(variables[i][i] == 0,"out %s" % str(i))
    
    for i in range(0,len(N)):
        m.addConstr(sum(variables[i]) + sum(depotVariablesIn[i]) == 1,"out %s" % str(i))
        
    for i in range(0,len(N)):
        inEdges = [variables[j][i] for j in range(0,len(N))]
        m.addConstr(sum(inEdges)+ sum(depotVariablesOut[i]) == 1,"in %s" % str(i))       
    
    
    m.optimize()    
    m.display()
    m.printAttr('X')
    
    edges=[]
    for v in m.getVars():
        if v.x==0:
            continue
        if v.varName.startswith('x'):
            nodes = v.varName.split('x')[1]
            nodes=nodes.split('-')
            assert len(nodes)==2
            edges.append([nodes[0],nodes[1]])
        elif v.varName.startswith('d'):
            e = v.varName.split('d')[1]
            node = e.split('n')[1]
            depot = e.split('n')[0]
            s = 'depot' + str(depot)
            edges.append([node, s]) 
        else:
            assert v.varName.startswith('n') 
            e = v.varName.split('n')[1]
            depot = e.split('d')[1]
            node = e.split('d')[0]
            s = 'depot' + str(depot)
            edges.append([node, s])         
        
    solution = [D] + edges + [m.objVal]   
    return solution

       

def doLP(N, distanceObject, modelname,D,K):
    functions=globals().copy()
    method = functions.get(modelname)
    if not method:
        raise Exception('Method %s not found' %modelname)
    solution =method(N,distanceObject,D,K)
    return solution

def bestPossibleTrips(N,K,distanceO):
    """
    LP to find the cheapest way to construct paths with at most K nodes from N,
    in such a way that every node from N is in one path    
    """
    # Create a new model
    m = Model("bestPossibleTrips")
    
    # Create variables of combinations of nodes
    combinations=[]
    removedCombinations=[]
    for i in range(1,K+1):
        combinations += list(itertools.combinations(N, i))
    
    variables=[]
    for i,combo in enumerate(combinations):
        d, tour=distanceO.minimumTourDistanceNoDepots(combo)
        comboints=[node[0] for node in tour]
        variablename = 'combo' + '-'.join(map(str,comboints))
        variables.append(m.addVar(vtype=GRB.BINARY,obj=d, name=variablename))
    
    # Integrate new variables
    m.modelSense = GRB.MINIMIZE
    m.update()

    # Every node needs to be serviced once and only once
    combinations = [combo for combo in combinations if combo not in removedCombinations]
    for node in N:
        containingNode = [i for i,combo in enumerate(combinations) if node in combo]
        m.addConstr(quicksum(variables[i] for i in containingNode) == 1,"service constraint %s" %node[0])   
        
    # Maximum length of tours
    combinations = [combo for combo in combinations if combo not in removedCombinations]
    oldNum=0
    for i in range(0,K):
        length = float(K-i)
        numTours = int((len(N)-oldNum)/length)
        containingNode = [i for i,combo in enumerate(combinations) if len(combo)==length]
        m.addConstr(quicksum(variables[i] for i in containingNode)>= numTours,"maximumTour length %s" %length) 
        oldNum += numTours*length   
                       
    m.optimize()
    m.display()
    m.printAttr('X')
    print('Obj: %g' % m.objVal)
    
    edges=[]
    for v in m.getVars():
        if v.x==0:
            continue
        assert v.varName.startswith('combo')
        nodes = v.varName.split('combo')[1]
        nodes=nodes.split('-')
        for i in range(1,len(nodes)):
            edges.append([nodes[i-1],nodes[i]])           
    solution = [[-1,-1]] +edges +[m.objVal]

    return solution
