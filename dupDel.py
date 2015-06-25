"""
  Author: Kevin Heath, Zunyan Wang
  Date: 9 June 2015
  Generated for summer research in the Bush lab 

  Duplication/Deletion (dupDel)
"""

import mrca,sys

def main(argv):
    '''Do stuff'''
    famT=(7,3,0,104,0,2)
    tree=mrca.readTree(argv[0]) 
    delCost=int(argv[1]) #user input of deletion cost.
    dupCost=int(argv[2]) #user input of duplication cost.
    currentcopynum=int(argv[3]) #user input of initial copy numbers. Normally, 1 should be inputted.
    mrcaA = mrca.mrca(tree, famT) #find the most recent common ancestor (mrca).
    subtreeA = subtree(mrcaA, tree) #find the subtree which has the mrca as its root.
    result=dupDel(subtreeA,famT,delCost,dupCost,currentcopynum,{}) #store the dupDel result.
    print result

def dupDel(tree, famT, delCost, dupCost, currentcopynum,memo):
    '''Calculate the minimal cost of duplication and deletion that could happen and the sequences of 
    duplication and deletion events that correspond to that minimal cost.'''

    if tree[1] == ():
        return (0,[],[]) #base case, we get to a leaf.
    elif (currentcopynum, tree) in memo: return memo[(currentcopynum,tree)] #We memoize this function to make it run faster.
    else:  # This is a subtree
        leftLeaves=descendantFam(tree[1][0],tree,famT) #find all the possible family copy numbers of leaves under the left subtree.
        minLeftCost=float('inf') #this variable stores the minimal cost, initially set to infinity.
        
        for i in range(min(leftLeaves),max(leftLeaves)+1): #loop over all possible copy numbers.
            leftSub=dupDel(tree[1],famT,delCost,dupCost,i,memo) #recursion step.
            leftSubCost=leftSub[0]
            leftDelList=leftSub[1]
            leftDupList=leftSub[2]

            if i >= currentcopynum: #calculate the cost of this particular copy number i.
                leftCost=(i-currentcopynum)*dupCost+leftSubCost                
            else:
                leftCost=(currentcopynum-i)*delCost+leftSubCost
                
            if leftCost<minLeftCost: #if this cost is lower than current minimum, this result should be stored and replace current minimum.
                minLeftCost=leftCost
                realLeftDelList=leftDelList #stores the duplication and deletion lists of lower level operations.
                realLeftDupList=leftDupList
                if i > currentcopynum: #find the required duplication or deletion events on this level.
                    leftDupAction=[tree[1][0]]*(i-currentcopynum)
                    leftDelAction=[]
                elif i < currentcopynum:
                    leftDupAction=[]
                    leftDelAction=[tree[1][0]]*(currentcopynum-i)
                else:
                    leftDelAction=[]
                    leftDupAction=[]
        realLeftDelList+=leftDelAction #add the list of previous events and events on this level together
        realLeftDupList+=leftDupAction

        rightLeaves=descendantFam(tree[2][0],tree,famT) #repeat the first half codes on the right tree
        minRightCost=float('inf')
        for i in range(min(rightLeaves),max(rightLeaves)+1):
            
            rightSub=dupDel(tree[2],famT,delCost,dupCost,i,memo)
            rightSubCost=rightSub[0]
            rightDelList=rightSub[1]
            rightDupList=rightSub[2]
            if i >= currentcopynum:
                rightCost=(i-currentcopynum)*dupCost+rightSubCost
            else:
                rightCost=(currentcopynum-i)*delCost+rightSubCost
            if rightCost<minRightCost:
                minRightCost=rightCost
                realRightDelList=rightDelList
                realRightDupList=rightDupList
                if i > currentcopynum:
                    rightDupAction=[tree[2][0]]*(i-currentcopynum)
                    rightDelAction=[]
                elif i < currentcopynum:
                    rightDupAction=[]
                    rightDelAction=[tree[2][0]]*(currentcopynum-i)
                else:
                    rightDelAction=[]
                    rightDupAction=[]
        realRightDelList+=rightDelAction
        realRightDupList+=rightDupAction

        memo[(currentcopynum,tree)]=(minLeftCost+minRightCost,realLeftDelList+realRightDelList,realLeftDupList+realRightDupList) #memoize calculated result
        return (minLeftCost+minRightCost,realLeftDelList+realRightDelList,realLeftDupList+realRightDupList)
 
def calcSubCost(tree, famT, delCost, dupCost, currentcopynum,memo):

    leaves=descendantFam(tree[0],tree,famT) #find all the possible family copy numbers of leaves under the left subtree.
    minCost=float('inf') #this variable stores the minimal cost, initially set to infinity.
    fullDelList = []
    fullDupList = []

    for i in range(min(leaves),max(leaves)+1): #loop over all possible copy numbers.
            subTree=dupDel(tree,famT,delCost,dupCost,i,memo) #recursion step.
            subCost=subTree[0]
            delList=subTree[1]
            dupList=subTree[2]

            if i >= currentcopynum: #calculate the cost of this particular copy number i.
                subCost=(i-currentcopynum)*dupCost+subCost                
            else:
                subCost=(currentcopynum-i)*delCost+subCost
                
            if subCost<minCost: #if this cost is lower than current minimum, this result should be stored and replace current minimum.
                minCost = subCost
                fullDelList=delList #stores the duplication and deletion lists of lower level operations.
                fullDupList=dupList
                if i > currentcopynum: #find the required duplication or deletion events on this level.
                    fullDelList.extend([tree[0]]*(i-currentcopynum))
                    # delAction=[]
                elif i < currentcopynum:
                    # dupAction=[]
                    fullDelList.extend([tree[0]]*(currentcopynum-i))
                else:
                    # delAction=[]
                    # dupAction=[]
        return minCost, fullDelList, fullDupList


      
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

def descendantNodes(node, Tree):
    '''Returns the descendant nodes of a given node.'''
    return leafList(subtree(node, Tree))

def descendantFam(node,Tree,famT):
    result=[]
    leaves=descendantNodes(node,Tree)
    for leaf in leaves:
        result.append(famT[leaf])
    return result
    


if __name__ == "__main__":
   main(sys.argv[1:])
