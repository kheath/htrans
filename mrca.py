"""
  Author: Kevin Heath
  Date: 29 May 2015
  Generated for summer research in the Bush lab 

  Most Recent Common Ancestor

  Takes tree file as input. 
  Example: python mrca.py treeFile
"""

import copy, ast, sys, argparse

def main(argv):

    # parser = argparse.ArgumentParser()
    # requiredNamed = parser.add_argument_group('required named arguments')
    # requiredNamed.add_argument('-t', '--tree', help='Tree file name ()', required=True)
    # requiredNamed.add_argument('-f', '--family', help='family file name', required=True)
    # requiredNamed.add_argument('-o', '--output', help='Output file name', required=True)
    # requiredNamed.add_argument('-n', help='Number of species', type=int, required=True)
    # requiredNamed.add_argument('-p', help='Partition file name', required=True)  
    # parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    # args = parser.parse_args()

    tree = readTree('testAtree')
    
    family = (3, 0, 1, 0, 0, 0)

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
    if max(family) <= 0:
        return -1

    if tree[1] == () and tree[2] == ():
        # print 'Found tip #'+str(tree[0])+'\t'+str([tree[0], family[tree[0]]])
        if family[tree[0]] == 0:
            return (tree[0], 0)
        else:
            return (tree[0], 1)
    else:
        r1 = mrca(tree[1], family)
        r2 = mrca(tree[2], family)

        # print 'Comparing r1 = '+str(r1)+' and r2 = '+str(r2)

        if r1[1] == 0 and r2[1] == 0:
            # print 'Returning '+str(-1)
            return (tree[0],0)
        elif r1[1] > 0 and r2[1] > 0:
            # print 'Returning '+str(tree[0])
            return (tree[0], 1)
        elif r1[1] > 0:
            # print 'Returning '+str(r1)
            return r1
        else:
            # print 'Returning '+str(r2)
            return r2




if __name__ == "__main__":
   main(sys.argv[1:])




