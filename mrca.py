"""
  Author: Kevin Heath
  Date: 29 May 2015
  Generated for summer research in the Bush lab 

  Most Recent Common Ancestor
"""

import copy, ast

def main():


    tree = ''
    with open('testATree', 'r') as f:
      temp = f.read()
      tree = ast.literal_eval(temp)#"".join(temp.split()))

    # print tree

    
    family = [1, -1, 1, 1, -1, -1, 0, 0, 0, 0, 0]

    print mrca(tree, family)


def mrca(tree, family):
    """ Identifies the most recent common ancestor of a gene family
    
    General:
        1. If both nodes/tips = -1
            return -1
        2. Neither nodees/tips = -1
            return current node
        3. One value = n, other = -1
            return n

    """
    
    if tree[1] == () and tree[2] == ():
      # print 'Found tip #'+str(tree[0])+'\t'+str([family[tree[0]], tree[0]])
      if family[tree[0]] == -1:
        return -1
      else:
        return tree[0]
    else:
      r1 = mrca(tree[1], family)
      r2 = mrca(tree[2], family)

      # print 'Comparing r1 = '+str(r1)+' and r2 = '+str(r2)

      if r1 == -1 and r2 == -1:
        # print 'Returning '+str(-1)
        return -1
      elif r1 >= 0 and r2 >= 0:
        # print 'Returning '+str(tree[0])
        return tree[0]
      elif r1 >= 0:
        # print 'Returning '+str(r1)
        return r1
      else:
        # print 'Returning '+str(r2)
        return r2




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

def subtree(node, tree):
    '''Returns the subtree rooted at node'''
    if tree == ():
      return ()
    elif tree[0] == node:
      return tree
    elif find(node, tree[1]):
      return subtree(node, tree[1])
    else:
      return subtree(node, tree[2])

def parent(node, Tree):
    '''Returns the immediate parent of a given node.'''
    if node == Tree[0]:
        return Tree[0]
    else:
        if find(node, Tree[1]) == True:
            if node == Tree[1][0]:
                return Tree[0]
            else:
                return parent(node, Tree[1])
        else:
            if node == Tree[2][0]:
                return Tree[0]
            else:
                return parent(node, Tree[2])

def scale(Tree, scaleFactor):
    '''Returns the given tree with number values as nodes denoting years since a common ancestor existed scaled by scaleFactor.'''
    if Tree[1] == ():
        return Tree
    else:
        return (Tree[0] * scaleFactor, scale(Tree[1], scaleFactor), scale(Tree[2], scaleFactor))

def descendantNodes(node, Tree):
    '''Returns the descendant nodes of a given node.'''
    return nodeList(subtree(node, Tree))[1:]

def leafCount(Tree):
    '''Returns the number of leaves in a tree.'''
    if Tree[1] == ():
        return 1
    else:
        return leafCount(Tree[1]) + leafCount(Tree[2])

if __name__ == "__main__":
   main()




