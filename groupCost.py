

import group
import orderCost1

def groupCost(groupA, groupB, tree, adjInfo):
    '''Find the cost of merging two groups'''

    pairs = [(groupA.getFamilies()[0], groupB.getFamilies()[0]), (groupA.getFamilies()[-1], groupB.getFamilies()[0]), 
                (groupA.getFamilies()[0], groupB.getFamilies()[-1]), (groupA.getFamilies()[-1], groupB.getFamilies()[-1])] 

    f=open('FamSpAdjD.txt','r')
    s=f.readline()
    FamSpAdjD=ast.literal_eval(s)
    f.close()

    minCost = float(inf)
    bestOrder = 0
    for index, adj in enumerate(pairs):
        cost = pairOrderCost(mrca.readTree('testATree'), 1, FamSpAdjD, groupA.getFamilies(),
                 groupB.getFamilies, adj, {})
        if cost < minCost:
            minCost = cost
            bestOrder = index

    return minCost, mergeLists(bestOrder, groupA.getFamilies(), groupB.getFamilies())


def mergeLists(order, l1, l2):
    '''Merge l1 and l2 in the correct order based on the pairs in groupCost'''

    return {0: l2.reverse()+l1, 1: l1+l2, 2: l2+l1, 3: l1+l2.reverse()}[order]
