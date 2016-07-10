import math, itertools


class distances:

    def __init__(self, N):
        self.DistancesM=[[]]
        self.DistancesM = [0] * (len(N)+1)
        for i,column in enumerate(self.DistancesM):
            self.DistancesM[i]=[0]*(len(N)+1)
        for node1 in N:
            for node2 in N:
                d = self.Euclidean(node1[1], node1[2], node2[1], node2[2])
                self.DistancesM[node1[0]][node2[0]] =d
                self.DistancesM[node2[0]][node1[0]] =d
                
            
    def Euclidean(self,x1,y1,x2,y2):    
        return math.sqrt(math.pow((x1-x2),2)+math.pow((y1-y2),2))  
    
    def EuclideanNode(self,node1, node2):
        return self.DistancesM[node1[0]][node2[0]]      
    
    def EuclideanNodeWithCheck(self, node1,node2,D):
        d = self.DistancesM[node1[0]][node2[0]]
        if  self.distanceToDepot(node1, D) + self.distanceToDepot(node2, D) < d :
            return -1
        return d
             
    def minimumTourDistanceNoDepots(self,trip):
        distance = float('inf')
        bestTour=None
        tourPermutations = itertools.permutations(trip)
        for tour in tourPermutations:
            d = self.getDistanceTourNoDepots(tour)
            if d<distance:
                distance = d
                bestTour=list(tour)
        return distance,bestTour
    
    def checkTour(self,tour,distance,D):
        for i in range(1, len(tour)):
            node1=tour[i-1]
            node2=tour[i]
            if self.distanceToDepot(node1, D) + self.distanceToDepot(node2, D) < self.EuclideanNode(node1, node2):
                return False, None
        return distance,tour
    
    def minimumTourDistance(self,tour,D):
        combinations  =list(itertools.permutations(tour,len(tour)))
        distance=float('inf')
        bestTour = None
        for combi in combinations:
            d=self.getDistanceTour(combi,D)
            if d< distance:
                distance=d
                bestTour = list(combi)
        return distance,bestTour
    
    def minimumTourDistanceWithCheck(self,tour,D):
        combinations  =list(itertools.permutations(tour,len(tour)))
        distance=float('inf')
        bestTour = None
        for combi in combinations:
            d=self.getDistanceTour(combi,D)
            if d< distance:
                distance=d
                bestTour = list(combi)
        return self.checkTour(bestTour,distance,D)
    
    def getDistanceTourNoDepots(self,tour) :
        if len(tour)==0:
            return 0   
        previousNode=tour[0]
        distance=0
        for node in tour[1:]:
            distance+=self.EuclideanNode(previousNode, node)
            previousNode=node    
        return distance   
    
    def distanceToDepot(self,node, D):
        x=node[1]
        y=node[2]
        return min(self.Euclidean(D[0],0,x,y),self.Euclidean(D[1],0,x,y ))
    
    def minimumDistanceToDepot(self,tour, D):
        if len(tour)==0:
            return 0
        distance=float('inf')
        toDepot=None
        for node in tour:
            d=self.distanceToDepot(node, D)
            if d<distance:
                distance=d
                toDepot = node
        return distance,toDepot
    
    def getDistanceTour(self,tour,D) :
        if len(tour)==0:
            return 0   
        distance = self.distanceToDepot(tour[0], D)
        previousNode=tour[0]
        for node in tour[1:]:
            distance+=self.EuclideanNode(previousNode, node)
            previousNode=node    
        distance+= self.distanceToDepot(tour[-1], D)
        return distance  
    
    def getDistanceTours(self,tours,D):
        d=0
        for tour in tours:
            d+=self.getDistanceTour(tour,D)
        return d
    
    def closestDepotToString(self, depot1,depot2,node):
        x=node[1]
        y=node[2]
        d1=self.Euclidean(depot1,0,x,y)
        d2=self.Euclidean(depot2,0,x,y )
        if d1< d2:
            return 'depot1', d1
        return 'depot2', d2
    
    
