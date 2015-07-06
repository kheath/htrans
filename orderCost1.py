
import mrca,sys,ast


def main(argv):
    tree=mrca.readTree(argv[0])
    cost=int(argv[1])
    FamSpAdjDFile=argv[2]
    f=open(FamSpAdjDFile,'r')
    s=f.readline()
    FamSpAdjD=ast.literal_eval(s)
    f.close()
    groupA=(ast.literal_eval(argv[3]))
    print groupA
    
    groupB=(ast.literal_eval(argv[4]))
    print groupB
    pair=(ast.literal_eval(argv[5]))
    print type(pair[1])
    result=pairOrderCost(tree,cost,FamSpAdjD,groupA,groupB,pair,{})
    print result

def pairOrderCost(tree,cost,FamSpAdjD,groupA,groupB,pair,memo):
    
    if tree[1]==(): #base case: At a leaf:
        if ((pair[1],tree[0]) in FamSpAdjD) == False: 
            return cost
        else:
            if ((pair[0],tree[0]) in FamSpAdjD) == False: #If first family got deleted:
                print "pathD"
                if remainFam(groupA,tree[0],FamSpAdjD) == []:
                    return 0
                else:
                    remain=remainFam(groupA,tree[0],FamSpAdjD)
                    print remain
                    for fam in remain:
                        if fam in FamSpAdjD[(pair[1],tree[0])]: 
                            return 0
                        else: return cost
                    '''
                    ##Following code needs adjacency information within a group.##
                    closest=pair[0]
                    unselectable=[]


                    while True:
                        if closest in remain: break
                        adjacentFam=FamSpAdjD[(closest,tree[0])]
                        for fam in adjacentFam:
                            if (fam in groupA) and (fam not in unselectable): 
                                unselectable.append(closest)
                                closest=fam
                                break
                    if closest in FamSpAdjD[(pair[1],tree[0])]: return 0

                    else: return cost
                    '''
            else: #If both families are present:

                if pair[1] in FamSpAdjD[(pair[0],tree[0])]: 
                    #print "pathE"
                    return 0 #If they are adjacent, no charge.
                else:
                    #print "pathF"
                    return cost #Otherwise, charge.

    elif (tree,pair) in memo: return memo[(tree,pair)]
    else: #General case: At a subtree:
        leftLeaves=leafList(tree[1])
        options = []
        for leaf in leftLeaves:
            if (pair[0],leaf) in FamSpAdjD:
                adjacentFam=FamSpAdjD[(pair[0],leaf)]
                for fam in adjacentFam:
                    if (fam not in options):
                        options.append(fam)
        if (pair[1] not in options): options.append(pair[1])
        
        
        minLeftCost = float('inf')
        for option in options:
            leftSubCost=pairOrderCost(tree[1],cost,FamSpAdjD,groupA,groupB,(pair[0],option),memo)
            if option == pair[1]:
                leftCost=leftSubCost
            else:
                leftCost=cost+leftSubCost
            #print leftCost
            if leftCost < minLeftCost:
                minLeftCost=leftCost
        #print 'minLeftCost:',minLeftCost
        
        rightLeaves=leafList(tree[2])
        options = []
        for leaf in rightLeaves:
            if (pair[0],leaf) in FamSpAdjD:
                adjacentFam=FamSpAdjD[(pair[0],leaf)]
                for fam in adjacentFam:
                    if (fam not in options):
                        options.append(fam)
        
        if (pair[1] not in options): options.append(pair[1])

        minRightCost = float('inf')
        for option in options:
            rightSubCost=pairOrderCost(tree[2],cost,FamSpAdjD,groupA,groupB,(pair[0],option),memo)
            if option == pair[1]:
                rightCost=rightSubCost
            else:
                rightCost=cost+rightSubCost
            if rightCost < minRightCost:
                minRightCost=rightCost
        #print 'minRightCost:', minRightCost
    memo[(tree,pair)]= minLeftCost+minRightCost
    return minLeftCost+minRightCost


def remainFam(group,leaf,FamSpAdjD):
    result=[]
    for fam in group:
        if (fam,leaf) in FamSpAdjD:
            result.append(fam)
    return result
                

                
def descendantNodes(node, Tree):
    '''Returns the descendant nodes of a given node.'''
    return leafList(subtree(node, Tree))

def find(node, Tree):
    ''' Returns True if node is in the Tree and False otherwise. '''
    if Tree == (): return False
    elif Tree[0] == node: return True # found it at the Root!
    else:
        return find(node, Tree[1]) or find(node, Tree[2])

def nodeList(Tree):
    """Returns the list of all nodes in Tree."""
    if Tree[1]==():
        return [Tree[0]]
    else:
        return [Tree[0]]+nodeList(Tree[1])+nodeList(Tree[2])
def leafList(Tree):
    if Tree[1]==():
        return [Tree[0]]
    else:
        return leafList(Tree[1])+leafList(Tree[2])
def subtree(node, Tree):
    '''The function subtree() returns the subtree at the inputted node.'''
    if node == Tree[0]:
        return Tree
    elif find(node, Tree[1]) == True:
        return subtree(node, Tree[1])
    else:
        return subtree(node, Tree[2])


if __name__ == "__main__":
    main(sys.argv[1:])
