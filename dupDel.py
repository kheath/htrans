"""
  Author: Kevin Heath
  Date: 29 May 2015
  Generated for summer research in the Bush lab 

  Most Recent Common Ancestor

  Takes tree file as input. 
  Example: python mrca.py treeFile
"""

import mrca

def main(argv):
    '''Do stuff'''
    mrca = mrca(tree, famT)

    subtree = subtree(mrca, tree)


def dupDel(tree, famT, delCost, dupCost, currentcopynum):
    '''Calculate all possible duplications and deletions for a family to get to curent species'''

    if tree[1] == ():
      if currentcopynum == famT(tree[0]):
        return 0
      else:
        return float('inf')

    else:  # This is a subtree
      r1 = dupDel(tree[1], famT, delCost, dupCost, currentcopynum)





      r2 = dupDel(tree[2], famT, delCost, dupCost, currentcopynum)


def getCounts(tree, famT):

    leaves = descendantNodes(tree[0], tree)
    counts = 
      
def idk(tree, famT, delCost, dupCost, distance):
    leaves = descendantNodes(tree[0], tree)

    for i in range(min(leaves), max(leaves)):




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
    return nodeList(subtree(node, Tree))[1:]




if __name__ == "__main__":
   main(sys.argv[1:])
