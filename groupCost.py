

from group import Group
from orderCost1 import *
import mrca
import ast
from operator import itemgetter
from copy import deepcopy
from itertool import izip

def main():

    y = Group(1)

    familyData = readFamilies('dupDelAll.txt')

    groupDict = initializeGroups(familyData)

    tree = mrca.readTree('testATree')

    f=open('FamSpAdjD.txt','r')
    s=f.readline()
    FamSpAdjD=ast.literal_eval(s)
    f.close()

    costs = {}
    for x in range(0, len(groupDict)):
        for y in range(x+1, len(groupDict)):
            costs[(x,y)] = groupCost(groupDict[x], groupDict[y], mrca.readTree('testATree'), FamSpAdjD)

    with open('orderCosts.txt', 'w+') as f:
        f.write(str(costs))


def initializeMagicalMatrix(groups, tree, famSpAdjD):

    distancesD = {}

    for x in range(0, len(groups)):
        for y in range(x+1, len(groups)):
            distancesD[(x,y)] = calcDist(groups[x], groups[y], tree, famSpAdjD)

    return distancesD

def mergeGroups(distancesD, groups, dCutoff, oCutoff, tree, famSpAdjD):

    # oldGroups = deepcopy(groups)
    startOver = False

    recalculate = []
    newGroups = []
    for x in range(0, len(groups)):
        for y in range(x+1, len(groups)):
            if x not in recalculate and y not in recalculate and distancesD[(x,y)] != None:
                if distancesD[(x,y)][0] <= dCutoff and distancesD[(x,y)][1][0] <= oCutoff:
                    if groups[x] != None and groups[y] != None:
                        groups[x] = deepcopy(groups[x].mergeGroup(groups[y], distancesD[(x,y)][1][1]))
                        groups[y] = None
                        recalculate.extend([x,y])


    xss = recalculate[0::2]
    yss = recalculate[1::2]

    for coord in distancesD.iterkeys():
        if coord[1] in yss:
            del distancesD[coord]
        elif coord[0] in xss:
            distancesD[coord] = calcDist(groups[coord[0]], groups[coord[1]], tree, famSpAdjD)

    # it = pairwise(recalculate)
    # while True:
    #     try:
    #         key = it.next()
    #     except StopIteration:
    #         break

    #     for loc in distancesD.iterkeys():
    #         if loc[0] == key[0]:
    #             distancesD[loc] = calcDist(groups[x], groups[y], tree, famSpAdjD)
    #         elif loc[1] == key[1]:
    #             del distancesD[loc]




    
def calcDist(groupA, groupB, tree, famSpAdjD):
    '''Compares the dupDel model and order cost of two groups'''


    ddDiff = calcDiff(groupA.getDuplications(), groupB.getDuplications())
                +calcDiff(groupA.getDeletions(), groupB.getDeletions())

    ordCost = groupCost(groupA, groupB, tree, famSpAdjD)

    return (ddDiff, ordCost)

def calcDiff(la, lb, n):
    '''Take two duplication or deletions models and compare them.
    Score starts at zero, and increases by 1 for each event difference.'''
    diff = 0
    for i in range(0,2*n-1):
        for x in range(0, len(i)):
            if la[i][x][0] != lb[i][x][0]:
                diff += 1
    return diff


def initializeGroups(familyData):
    '''Make every family into a group'''

    groups = [None]

    for index, value in enumerate(familyData):
        groups.append(Group(value, index, 6))

    return groups

def readFamilies(filename):
    '''Reads in the results of dupDel and stores families as a dictionary'''

    families = [None]

    with open(filename, 'r') as f:
        # famNum = 1
        while True:
            line = f.readline()

            if not line:
                break
            temp = list(ast.literal_eval(line))
            if temp[1][1] != [] or temp[1][2] != []:
                families.append(temp)
                # temp.append(famNum)
                # families.append(temp)
            # famNum+=1

    return families


def groupCost(groupA, groupB, tree, adjInfo):
    '''Find the cost of merging two groups'''

    pairs = [(groupA.getFront()[0], groupB.getFront()[0]), (groupA.getBack()[0], groupB.getFront()[0]), 
                (groupA.getFront()[0], groupB.getBack()[0]), (groupA.getBack()[0], groupB.getBack()[0])] 

    costs = [] #float('inf')
    for index, adj in enumerate(pairs):
        cost = pairOrderCost(tree, 1, adjInfo, groupA.getFamilies(), groupB.getFamilies(), adj, {})
        costs.append((cost, index))

    return sorted(costs, key=itemgetter(0,1))


def mergeLists(order, l1, l2):
    '''Merge l1 and l2 in the correct order based on the pairs in groupCost'''
    lb = l2
    if len(l2) > 1:
        lb = l2.reverse()


    return {0: lb+l1, 1: l1+l2, 2: l2+l1, 3: l1+lb}[order]

def pairwise(iterable):
    "s -> (s0,s1), (s2,s3), (s4, s5), ..."
    a = iter(iterable)
    return izip(a, a)


if __name__ == "__main__":
   main()
