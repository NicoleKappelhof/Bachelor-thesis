import os, shutil
import numpy
import networkx as nx
import matplotlib.pyplot as plt
import re, time
import Exact.ExactDepots as depots
import ThreeOrMore.GeneticMore as geneticMore
import ThreeOrMore.Distances as Distances
import Exact.linearProgram as LinProg

n = 300
l = 10
h = 10
N = []
K=4
kolmo=1
repetitions = 1
normalNumber = False

def randomNormalNumber(minV, maxV):
    mu = (maxV+minV)/ 2.0
    sigma  = (mu-minV)/3.0
    while True:
        number = numpy.random.normal(mu, sigma, 1)[0]
        if number> minV and number < maxV:
            return number

def uniformRandom(minV, maxV):
    while True:
        number = numpy.random.uniform(minV, maxV, 1)[0]
        if number> minV and number < maxV:
            return number

def main():
    string= ""
    for rep in range(1,repetitions+1):
        distanceObject = createNodes(rep)
        string+= doSolution(N,K,distanceObject, rep)
        string+= '\n'
    print string

    
def createNodes(rep):
    global N
    N=[]
    if normalNumber:
        for i in range(1,n+1):
            N.append([i, uniformRandom(0,l),uniformRandom(0,h)]);
            distanceObject = Distances.distances(N)     
            showNodes(N, title='Nodes' + str(len(N)) + '-'+str(rep))
            return distanceObject
    for i in range(1,int((1/6.0)*n)):
        N.append([i, randomNormalNumber(0, l), randomNormalNumber(0, h)]);
    for i in range(int((1/6.0)*n),int((3/6.0)*n)):
        N.append([i, randomNormalNumber(l/4.0, 3*l/4.0), randomNormalNumber(h/4.0, 3*h/4.0)]);
    for i in range(int((3/6.0)*n),int((4/6.0)*n)):
        N.append([i, randomNormalNumber(0, l*(1/4.0)), randomNormalNumber(0, h/2.0)]);
    for i in range(int((4/6.0)*n),int((5/6.0)*n)):
        N.append([i, randomNormalNumber(l*(3/4.0), l), randomNormalNumber(3*h/4.0, h)]);
    for i in range(int((5/6.0)*n),n+1):
        N.append([i, randomNormalNumber(l*(3/4.0), l), randomNormalNumber(0, (1/5.0)*h)]);
    distanceObject = Distances.distances(N)     
    showNodes(N, title='Nodes' + str(len(N)) + '-'+str(rep))
    return distanceObject

def doSolution(N,K,distanceObject,rep):
    sting = ''
    if K==2:
        sting += 'Kolmo'
        lengthFloat = writeToKolmogorov()
        raw_input("Press Enter to continue...")
        print "Continue"
        sting += obtainKolmogorov('solution.txt',distanceObject,lengthFloat,title='KolmogorovRND' + str(len(N))+ '-' + str(rep) )
        sting+= '\n'
    
# Functionality to find the least distance `matching' of length K (or least distance T using the genetic notation) 
#     time1=time.time()
#     solution =  LinProg.bestPossibleTrips( N,K, distanceObject)
#     time2=time.time()
#     dif=time2-time1
#     sting += 'bestPossibleTrips got %s, took %0.3f s' % (solution[-1],(dif))
#     showSolutionEdgesLp(solution,distanceObject,N,title='bestPossibleTrips' + str(len(N)) +' time'+str(dif))
    
    time1=time.time()
    lpiterations = 9
    solution =executeLP(K,N,distanceObject, lpiterations, 'lpK_sets')
    depots = '&'.join(['%.5f' %(depot) for depot in solution[0]])
    time2=time.time()
    dif=time2-time1
    sting += 'lpK_sets & %s & %s & %s & %.6f & %.4f depots %s \n' %(str(len(N)),str(lpiterations),str(rep),solution[-1], dif, depots)
    title = 'lpK_sets%s-%s time %.4f depots %s' %(str(len(N)),str(rep), dif, depots)
    showSolution(solution,N,distanceObject,title)
    
    time1=time.time()
    lpiterations = 1
    solution = executeLP(K,N,distanceObject, lpiterations, 'lpK_edges')
    depots = '&'.join(['%.5f' %(depot) for depot in solution[0]])
    time2=time.time()
    dif=time2-time1
    sting += 'lpK_edges & %s & %s & %s & %.6f & %.4f depots %s \n' %(str(len(N)),str(lpiterations),str(rep),solution[-1], dif, depots)
    title = 'lpK_edges%s-%s time %.4f depots %s' %(str(len(N)),str(rep), dif, depots)
    showSolutionEdgesLp(solution,distanceObject,N, title)
      
    time1=time.time()
    lpiterations = 9
    solution = executeLP(K,N,distanceObject, lpiterations, 'lpK_flow')
    time2=time.time()
    dif=time2-time1
    depots = '&'.join(['%.5f' %(depot) for depot in solution[0]])
    sting += 'lpK_flow & %s & %s & %s & %.6f & %.4f depots %s \n' %(str(len(N)),str(lpiterations),str(rep),solution[-1], dif, depots)  
    title = 'lpK_flow%s-%s time %.4f depots %s' %(str(len(N)),str(rep), dif, depots)
    showSolutionEdgesLp(solution,distanceObject,N, title)
     
    time1=time.time()
    solution = executeLP(K,N,distanceObject, 5, 'lpK_DoubleEdges')
    depots = '-'.join([str(depot) for depot in solution[0]])
    time2=time.time()
    dif=time2-time1
    sting += 'lpK_DoubleEdges got %s, took %0.3f s' % (solution[-1],(dif))
    title = 'lpK_DoubleEdges%s time %.4f depots %s' %(str(len(N)), dif, depots)
    showSolutionEdgesLp(solution,distanceObject,N, title)
     
    time1=time.time()
    solution = geneticMore.doGenetic(N,h,l,distanceObject,K)
    depots = '&'.join(['%.5f' %(depot) for depot in solution[0]])
    time2=time.time()
    dif=time2-time1
    sting += 'genetic & %s & %s  & %.6f & %.4f depots %s \n' %(str(len(N)),str(rep),solution[-1], dif, depots)  
    title = 'genetic%s-%s time %.4f depots %s' %(str(len(N)), str(rep),dif, depots)
    showSolution(solution, N, distanceObject, title)
    
    return sting
    
    
def getNodes():
    assert len(N) != 0, "Something went wrong with initialization"
    return N

def getField():
    return {'Length': l, 'Height': h} 

def executeLP(K,N,distanceObject, iterations, lp):
    depotLocation = depots.exactDepots(N, 2)
    Doud=[depotLocation[1], depotLocation[2]]
    solution = LinProg.doLP(N, distanceObject, lp,Doud,K)
    if iterations==1:
        return solution
    bestSolution=solution
    youd = solution[-1]
    Dnieuw=[(2/6.0)*l, (4/6.0)*l]
    solution = LinProg.doLP(N, distanceObject, lp,Dnieuw,K)
    ynieuw = solution[-1]
    change=(youd-ynieuw)
    i=2
    if ynieuw < youd:
        bestSolution = solution
    while abs(change) > len(N)/100 and i<iterations:
        i+=1
        d1 = Dnieuw[0] + (1.0/len(N))*(float(youd-ynieuw)/float(Doud[0]-Dnieuw[0]))
        d2 = Dnieuw[1] + (1.0/len(N))*(float(youd-ynieuw)/float(Doud[1]-Dnieuw[1]))
        youd=ynieuw
        Doud=Dnieuw
        Dnieuw=[d1,d2]
        solution = LinProg.doLP(N, distanceObject, lp,Dnieuw,K)
        ynieuw = solution[-1]
        change = (youd-ynieuw)
        if change>0:
            bestSolution = solution
    return bestSolution


def writeToKolmogorov():
    global kolmo, N
    file = open('KolmogorovProblem.txt','w')
    file.write("NAME : Kolmogorov" + str(kolmo)+"\n")
    file.write("TYPE : TSP \n")
    file.write("DIMENSION: " +str(n) +"\n")
    file.write("EDGE_WEIGHT_TYPE : EUC_2D \n")
    file.write("NODE_COORD_TYPE : TWOD_COORDS \n")
    file.write("NODE_COORD_SECTION \n")
    kolmo+=1
    lengthFloat = 10000000
    for node in N[:-1]:
        i=[node[0],0,0]
        i[1]=int(node[1]*lengthFloat)
        i[2]=int(node[2]*lengthFloat)
        file.write("  " + ' '.join(map(str, i)) +"\n")
    i=[N[-1][0],0,0]
    i[1]=int(N[-1][1]*lengthFloat)
    i[2]=int(N[-1][2]*lengthFloat)
    file.write("  " + ' '.join(map(str, i)) +"EOF")
    file.close()
    a = shutil.move(file.name,'C:\Users\Admin\workspacecpp\cppScriptie\KolmogorovProblem.txt' )
    return lengthFloat
    

def getInput():
    assert os.path.isfile('C:\Users\Admin\workspacecpp\cppScriptie\KolmogorovProblem.txt')
    file=open('C:\Users\Admin\workspacecpp\cppScriptie\KolmogorovProblem.txt','r')
    line=file.readline()
    match = re.match(r'^\s*?[0-9]', line)
    while match is None:
        line =file.readline()
        match=re.match(r'^\s*?[0-9]', line)
    Nodes=[]
    match = re.match(r'EOF',line)
    while match is None and line is not None:
        coordinates = line.split()
        Nodes.append(coordinates)
        line = file.readline()
        match = re.match(r'.*EOF',line)
    coordinates = line.split()
    coordinates[2]=re.findall('\d+',coordinates[2])[0]
    Nodes.append(coordinates)
    return Nodes
    
def getCoordinates(n, Nodes, rounding):
    for node in Nodes:
        if int(node[0])==(n+1): # Kolmogorov is off by one
            return float(node[1])/rounding,float(node[2])/rounding
    raise Exception('Node %s not found'%(n+1))  

def obtainKolmogorov(filename, distanceO,lengthFloat, title):
    global N
    Nodes = getInput()
    timeD=time.time()
    sol =  depots.exactDepots(N, 2, 2.5, 7.5)
    timeD2=time.time()
    print 'Depot Allocation took %0.6f s' % ( timeD2-timeD)
    depot1 = sol[1]
    depot2=sol[2]
    file = os.path.join('C:\Users\Admin\workspacecpp\cppScriptie',filename)
    assert os.path.isfile(file)
    file = open(str(file), 'r')
    firstline = file.readline()
    parameters = firstline.split()
    assert len(parameters)==2, "File not correct"
    matchings = int(parameters[1]) #Not sure if correct
    matching=[]
    G = nx.Graph()
    G.add_node('depot1', pos=(depot1,0))
    G.add_node('depot2', pos=(depot2,0))
    objValue=sol[0]
    for i in range(0,matchings):
        oneMatching = file.readline().split()
        matching.append(oneMatching)
        a = getCoordinates(int(oneMatching[0]),Nodes, lengthFloat)
        G.add_node(oneMatching[0], pos=(a[0],a[1]))
        b = getCoordinates(int(oneMatching[1]),Nodes,lengthFloat)
        G.add_node(oneMatching[1], pos=(b[0],b[1]))
        G.add_edges_from([(oneMatching[0], oneMatching[1])])
        
        todepot = distanceO.closestDepotToString(depot1,depot2,['', a[0],a[1]])[0]
        G.add_edges_from([(todepot, oneMatching[0])])
        todepot = distanceO.closestDepotToString(depot1,depot2,['', b[0],b[1]])[0]
        G.add_edges_from([(todepot, oneMatching[1])])
        
        objValue += distanceO.Euclidean(b[0],b[1],a[0],a[1])
    objValue = objValue  
    string = 'Depot:%0.6f & %.07f' % (timeD2-timeD, objValue) 
    file.close()
    pos=nx.get_node_attributes(G,'pos')
    plt.title(title + 'obj:' + str(objValue))
    nx.draw(G,pos, node_size=5, with_labels=False)
    plt.savefig(title +'.png')
    plt.clf()
    return string
    
def showNodes(N,title=None):
    G = nx.Graph()
    for node in N:
        G.add_node(str(node[0]), pos=(float(node[1]),float(node[2])))
    pos=nx.get_node_attributes(G,'pos')
    if title:
        plt.title(title)
    nx.draw(G,pos, node_size=5, with_labels=False)
    if title:
        plt.savefig(title + '.png')
    plt.clf()
       
def showSolutionEdgesLp(allEdges,distanceO,N,title=None):
    G = nx.Graph()
    D=allEdges[0]
    objValue=allEdges[-1]
    for node in N:
        G.add_node(str(node[0]), pos=(float(node[1]),float(node[2])))
    for i,depot in enumerate(D):
        if depot>=0:
            s='depot' + str(i+1)
            G.add_node(s, pos=(float(depot),0))
    allEdges=allEdges[1:-1] 
    for edge in allEdges:
        G.add_edges_from([(str(edge[0]), str(edge[1]))]) 
    pos=nx.get_node_attributes(G,'pos')
    if title:
        plt.title(title + 'obj:' + str(objValue))
    nx.draw(G,pos, node_size=5, with_labels=False)
    if title:
        plt.savefig(title + '.png')
    plt.clf()
     
def showSolution(solution,N,distanceO,title=None):
    print title
    G = nx.Graph()
    assert solution is not None
    depot1 = solution[0][0]
    G.add_node('depot1', pos=(float(depot1),0))
    depot2= solution[0][1]
    G.add_node('depot2', pos=(float(depot2),0))
    for tour in solution[1:-1]:
        startingDepot,__ = distanceO.closestDepotToString(depot1,depot2,tour[0])        
        endingDepot,__= distanceO.closestDepotToString(depot1,depot2,tour[-1])
        tourSting = startingDepot
        previousNode=None
        for node in tour:
            tourSting=tourSting+str(node[0])
            G.add_node(str(node[0]), pos=(float(node[1]),float(node[2])))
            if previousNode:
                G.add_edges_from([(previousNode, str(node[0]))]) 
            previousNode = str(node[0])
        G.add_edges_from([(startingDepot, str(tour[0][0]))])
        G.add_edges_from([(endingDepot, str(tour[-1][0]))])   

    pos=nx.get_node_attributes(G,'pos')
    if title:
        plt.title(title + 'obj:' + str(solution[-1]) )
    nx.draw(G,pos, node_size=5, with_labels=False)
    if title:
        plt.savefig(title + '.png')
    plt.clf()
    return 

    
main()
    