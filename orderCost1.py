from htrans import *
import mrca,sys,ast,group


def main(argv):
    tree=mrca.readTree(argv[0])
    originTree=tree
    cost=int(argv[1])
    FamSpAdjDFile=argv[2]
    f=open(FamSpAdjDFile,'r')
    s=f.readline()
    FamSpAdjD=ast.literal_eval(s)
    f.close()

    
        
    pair=(ast.literal_eval(argv[3]))
    groupL=initializeGroups(readFamilies('dupDelAll.txt'))
    groupA=groupL[1221]
    groupB=groupL[1220]
    famGroupD=setFamGroupDict(groupL)
    result=pairOrderCost(tree,originTree,cost,FamSpAdjD,famGroupD,groupL,groupA,groupB,pair,{})
    print result

def pairOrderCost(tree,originTree,cost,FamSpAdjD,famGroupD,groupL,groupA,groupB,pair,memo):
    '''groupA and groupB are actual groups'''
    originHtrans=groupA.getMrcag()
    groupAFam=groupA.getFamilies()
    groupBFam=groupB.getFamilies()
    if tree[1]==(): #base case: At a leaf:
        if (remainFam(groupAFam,tree[0],FamSpAdjD) == []) or (remainFam(groupBFam,tree[0],FamSpAdjD) == []):
            return 0
        else:
            repA=representative(groupAFam,tree[0],pair[0],FamSpAdjD)
            repB=representative(groupBFam,tree[0],pair[1],FamSpAdjD)
            adjFList=FamSpAdjD[(repA,tree[0])]
            rAdjFList={}
            for adjF in adjFList:
                currentGroupID=famGroupD[adjF]
                currentGroup=groupL[currentGroupID]
                currentHtrans=currentGroup.getMrcag()
                if lowerThanOrigin(subtree(originHtrans,originTree),originHtrans,currentHtrans) == True:
                    otherEndFam=otherEnd(currentGroup.getFamilies(),adjF,tree[0],FamSpAdjD)
                    otherEndAdjFL=FamSpAdjD[(otherEndFam,tree[0])]
                    for otherEndAdjF in otherEndAdjFL:
                        if (otherEndAdjF != repA) and (otherEndAdj not in currentGroup.getFamilies()):
                            rAdjFList.append(otherEndAdjF)
                else:
                    rAdjFList.append(adjF)
            if repB in rAdjFList:
                return 0
            else:
                return cost                                                                
    elif (tree,pair) in memo: return memo[(tree,pair)]
    else: #General case: At a subtree:
        leftLeaves=leafList(tree[1])
        options = []
        for leaf in leftLeaves:
            if (pair[0],leaf) in FamSpAdjD:
                adjacentFam=FamSpAdjD[(pair[0],leaf)]
                for fam in adjacentFam:
                    if (fam not in options) and (fam in groupBFam):
                        options.append(fam)
        if (pair[1] not in options): options.append(pair[1])
        
        
        minLeftCost = float('inf')
        for option in options:
            leftSubCost=pairOrderCost(tree[1],originTree,cost,FamSpAdjD,famGroupD,groupL,groupA,groupB,(pair[0],option),memo)
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
                    if (fam not in options) and (fam in groupBFam):
                        options.append(fam)
        
        if (pair[1] not in options): options.append(pair[1])

        minRightCost = float('inf')
        for option in options:
            rightSubCost=pairOrderCost(tree[2],originTree,cost,FamSpAdjD,famGroupD,groupL,groupA,groupB,(pair[0],option),memo)
            if option == pair[1]:
                rightCost=rightSubCost
            else:
                rightCost=cost+rightSubCost
            if rightCost < minRightCost:
                minRightCost=rightCost
        #print 'minRightCost:', minRightCost
    memo[(tree,pair)]= minLeftCost+minRightCost
    return minLeftCost+minRightCost

def otherEnd(group,fam,leaf,FamSpAdjD):
    if len(group)==1:
        return group[0]
    else:
        if fam == group[0]:
            if (group[-1],leaf) in FamSpAdjD:
                return group[-1]
            else:
                return otherEnd(group[:-2],fam,leaf,FamSpAdjD)
        else:
            if (group[0],leaf) in FamSpAdjD:
                return group[0]
            else:
                return otherEnd(group[1:],fam,leaf,FamSpAdjD)
    

def lowerThanOrigin(Tree,originNode,currentNode):
    if Tree[1]==():
        return False
    elif Tree[0]==currentNode and Tree[0]!=originNode:
        return True
    else:
        return lowerThanOrigin(Tree[1],originNode,currentNode) or lowerThanOrigin(Tree[2],originNode,currentNode)
        

def remainFam(group,leaf,FamSpAdjD):
    result=[]
    for fam in group:
        if (fam,leaf) in FamSpAdjD:
            result.append(fam)
    return result
                
def representative(group,leaf,fam,FamSpAdjD):
    if (fam,leaf) in FamSpAdjD:
        return fam
    else:
        if fam==group[0]:
            return representative(group[1:],leaf,group[1],FamSpAdjD)
        elif fam==group[-1]:
            return representative(group[:-1],leaf,group[-2],FamSpAdjD)
        else:
            for family in group:
                if (family,leaf) in FamSpAdjD:
                    return family
            return None
                
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
