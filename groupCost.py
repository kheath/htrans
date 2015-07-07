

from group import Group
from orderCost1 import *
import mrca
import ast

def main():

    y = Group(1)

    familyData = readFamilies('dupDelAll.txt')

    groupDict = initializeGroups(familyData)

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


def initializeGroups(familyData):
    '''Make every family into a group'''

    groups = {}

    for index, family in enumerate(familyData):
        temp = Group(index)
        groups[index] = temp

    return groups

def readFamilies(filename):

    families = []

    with open(filename, 'r') as f:
        famNum = 1
        while True:
            line = f.readline()

            if not line:
                break
            temp = list(ast.literal_eval(line))
            if temp[1][1] != [] or temp[1][2] != []:
                temp.append(famNum)
                families.append(temp)
            famNum+=1

    return families


def groupCost(groupA, groupB, tree, adjInfo):
    '''Find the cost of merging two groups'''

    pairs = [(groupA.getFamilies()[0], groupB.getFamilies()[0]), (groupA.getFamilies()[-1], groupB.getFamilies()[0]), 
                (groupA.getFamilies()[0], groupB.getFamilies()[-1]), (groupA.getFamilies()[-1], groupB.getFamilies()[-1])] 


    minCost = float('inf')
    bestOrder = 0
    for index, adj in enumerate(pairs):
        cost = pairOrderCost(tree, 1, adjInfo, groupA.getFamilies(), groupB.getFamilies(), adj, {})
        if cost < minCost:
            minCost = cost
            bestOrder = index

    return (minCost, mergeLists(bestOrder, groupA.getFamilies(), groupB.getFamilies()))


def mergeLists(order, l1, l2):
    '''Merge l1 and l2 in the correct order based on the pairs in groupCost'''
    lb = l2
    if len(l2) > 1:
        lb = l2.reverse()


    return {0: lb+l1, 1: l1+l2, 2: l2+l1, 3: l1+lb}[order]

if __name__ == "__main__":
   main()
