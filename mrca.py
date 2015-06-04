"""
  Author: Kevin Heath
  Date: 29 May 2015
  Generated for summer research in the Bush lab 

  Most Recent Common Ancestor

  Takes tree file as input. 
  Example: python mrca.py treeFile
"""

import copy, ast, sys

def main(argv):


    tree = readTree(argv[0])
    
    family = (1, -1, 1, 1, -1, -1, 0, 0, 0, 0, 0)

    print mrca(tree, family)

def readTree(infile):
    """ Reads a text file with a tree represented as a tuple """
    with open(infile, 'r') as f:
      temp = f.read()
      tree = ast.literal_eval(temp)

    return tree


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




if __name__ == "__main__":
   main(sys.argv[1:])




