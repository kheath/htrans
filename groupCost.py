import group

def groupCost(groupA, groupB, tree, adjInfo):
    '''Find the cost of merging two groups'''

    pairs = [(groupA.getFamilies()[0], groupB.getFamilies()[0]), (groupA.getFamilies()[-1], groupB.getFamilies()[0]), 
                (groupA.getFamilies()[0], groupB.getFamilies()[-1]), (groupA.getFamilies()[-1], groupB.getFamilies()[-1])] 

    minCost = float(inf)
    bestOrder = 0
    for index, adj in enumerate(pairs):
        cost = hOrderCost(tree, adj, 0, adjInfo, groupA, groupB)
        if cost < minCost:
            minCost = cost
            bestOrder = index

    return minCost, mergeLists(bestOrder, groupA.getFamilies(), groupB.getFamilies())


def mergeLists(order, l1, l2):
    return {0: l2.reverse()+l1, 1: l1+l2, 2: l2+l1, 3: l1+l2.reverse()}[order]
