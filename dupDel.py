"""
  Author: Kevin Heath
  Date: 29 May 2015
  Generated for summer research in the Bush lab 

  Most Recent Common Ancestor

  Takes tree file as input. 
  Example: python mrca.py treeFile
"""

import mrca,sys


def main(argv):
    '''Do stuff'''
    famT=(3,0,107,0,0,0)
    tree=mrca.readTree(argv[0])
    delCost=int(argv[1])
    dupCost=int(argv[2])
    currentcopynum=int(argv[3])
    
    
    mrcaA = mrca.mrca(tree, famT)

    subtreeA = subtree(mrcaA, tree)
    result=dupDel(subtreeA,famT,delCost,dupCost,currentcopynum,{})
    print result

def dupDel(tree, famT, delCost, dupCost, currentcopynum,memo):
    '''Calculate all possible duplications and deletions for a family to get to curent species'''

    if tree[1] == ():
        return 0
    else:  # This is a subtree
        leftLeaves=descendantFam(tree[1][0],tree,famT)
        minLeftCost=float('inf')
        for i in range(min(leftLeaves),max(leftLeaves)+1):
            if (i,tree[1]) in memo:
                leftSubCost=memo[(i,tree[1])]
            else:
                leftSubCost=dupDel(tree[1],famT,delCost,dupCost,i,memo)
                memo[(i,tree[1])]=leftSubCost
            if i >= currentcopynum:
                leftCost=(i-currentcopynum)*dupCost+leftSubCost
                
            else:
                leftCost=(currentcopynum-i)*delCost+leftSubCost
                
            if leftCost<minLeftCost:
                minLeftCost=leftCost
                bestLeftI=i
                
        print 'bestLeftI:',bestLeftI
        rightLeaves=descendantFam(tree[2][0],tree,famT)
        minRightCost=float('inf')
        for i in range(min(rightLeaves),max(rightLeaves)+1):
            if (i,tree[2]) in memo:
                rightSubCost=memo[(i,tree[2])]
            else:
                rightSubCost=dupDel(tree[2],famT,delCost,dupCost,i,memo)
                memo[(i,tree[2])]=rightSubCost
            if i >= currentcopynum:
                rightCost=(i-currentcopynum)*dupCost+rightSubCost
            else:
                rightCost=(currentcopynum-i)*delCost+rightSubCost
            if rightCost<minRightCost:
                minRightCost=rightCost
                bestRightI=i
        print 'bestRightI:',bestRightI
        return minLeftCost+minRightCost
 
      
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
