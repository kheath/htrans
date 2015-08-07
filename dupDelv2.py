"""
  Author: Kevin Heath, Zunyan Wang
  Date: 9 June 2015
  Generated for summer research in the Bush lab 

  Duplication/Deletion (dupDel)
"""

import mrca, sys, ast, argparse
from functools import wraps

def main(argv):

    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-t', '--tree', help='Tree file name (testATree)', required=True)
    requiredNamed.add_argument('-f', '--fams', help='File with family tuple (famInfoResult.txt)', required=True)
    requiredNamed.add_argument('-d', help='Deletion cost', type=int, required=True)
    requiredNamed.add_argument('-c', help='Duplication cost', type=int, required=True) 
    # parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('-n', help='Initial copy number of genes at the mrca', type=int, required=False, default=1)
    args = parser.parse_args()
    
    # verbose = args.verbose


    fullTree = mrca.readTree(args.tree)

    famInfo = readFamInfo(args.fams) #Tuple of family tuples
    results=[]
    it = iter(famInfo)
    while True:
        try:
            fam = it.next()
        except:
            break
        mrcaF = mrca.mrca(fullTree, fam)
        sub = subtree(mrcaF, fullTree)
        results.append((mrcaF, dupDel(sub, fam, args.d, args.c, args.n)))

    # Write out the cost, duplication, deletion, and mrca information
    with open('dupDelAll.txt', 'w+') as f:
        for family in results:
            f.write(str(family)+'\n')


def readFamInfo(infile):
    '''Read in the famInfoResult.txt'''
    with open(infile, 'r') as f:
        line = ast.literal_eval(f.readline())
        if type(line) == tuple:
            return line[1:]
        else:
            print 'Error reading in family info'
            return False


def memoize(func):
    '''Memo function for dupDel to avoid repeated calculations'''
    cache = {}
    @ wraps(func)
    def wrap(*args):
        if args not in cache:
            cache[args] = func(*args)
            print cache
        return cache[args]

    return wrap

@memoize
def dupDel(tree, famT, delCost, dupCost, currentcopynum):
    '''Calculate the minimal cost of duplication and deletion that could happen and the sequences of 
    duplication and deletion events that correspond to that minimal cost.'''

    if tree[1] == ():
        return (0,[],[]) #base case, we get to a leaf.
    else:  # This is a subtree
        
        # Get the cost, duplications, and deletions for both subtrees
        lcost, lDels, lDups = calcSubCost(tree[1], famT, delCost, dupCost, currentcopynum)
        rcost, rDels, rDups = calcSubCost(tree[2], famT, delCost, dupCost, currentcopynum)

        return (lcost+rcost, lDels+rDels, lDups+rDups)

 
def calcSubCost(tree, famT, delCost, dupCost, currentcopynum):
    '''Recursive helper function of dupDel.
    Takes a subtree and the remaining info from dupDel and uses recursion to find
    the min cost and associated duplications and deletions.'''

    leaves=descendantFam(tree[0],tree,famT) #find all the possible family copy numbers of leaves under the left subtree.
    minCost=float('inf') #this variable stores the minimal cost, initially set to infinity.
    fullDelList = []
    fullDupList = []

    for i in range(min(leaves),max(leaves)+1): #loop over all possible copy numbers.
            subTree=dupDel(tree,famT,delCost,dupCost,i) #recursion step.
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
                elif i < currentcopynum:
                    fullDelList.extend([tree[0]]*(currentcopynum-i))

    return minCost, fullDelList, fullDupList


######### Stuff I Shameless Copied from a CS 5 Green Assignment #########

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
   main(sys.argv)
