'''
Author: Kevin Heath & Zunyan Wang
Date: August 1, 2015

Description:  Wrapper file for the htrans project.
    This file every function required for the pipeline
    after preprocessing.  This wrapper is not project
    specific, but requires a large amount of preprocessed
    data to work.  This will also run for a while, depending
    on how many gene families you are analyzing.

Example:  If you have the same files as me...(you probably don't)
    and want to receive a notification:
    python htrans.py -f fam.out -m geneSpeciesMap.txt -d dbList.txt -g geneOrder.txt 
            -t testATree -b 3 -c 5 -s 6 -o orthologs.out &> htrans.out; echo "Process 
            finished" | mail -s "Message" EmailAdress
'''

import sys, io, ast, math, os, argparse, time
from group import Group
from operator import itemgetter
from copy import deepcopy
from collections import Counter
from collections import defaultdict
from itertools import *
from group import Group
from copy import deepcopy
from os import path
from functools import wraps
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import resource
import gc
import cPickle as pickle




def main(argv):

    plt.ioff()  # Turn off interactive mode for plots
    
    ############## ---- Parse Command line Arguments ---- ################
    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-f', '--families', help='Family file name (fam.out)', required=True)
    requiredNamed.add_argument('-m', '--map', help='Gene to species map (geneSpeciesMap.txt)', required=True)
    requiredNamed.add_argument('-d', '--database', help='File with species of interest (dbList.txt)', required=True)
    requiredNamed.add_argument('-g', '--geneOrder', help='Gene order file (geneOrder.txt)', required=True)
    requiredNamed.add_argument('-t', '--tree', help='Phylogenetic tree file (testATree)', required=True) 
    requiredNamed.add_argument('-b', help='Deletion cost', type=int, required=True)
    requiredNamed.add_argument('-c', help='Duplication cost', type=int, required=True) 
    requiredNamed.add_argument('-s', '--species', help='Number of species', type=int, required=True)
    requiredNamed.add_argument('-o', help='Orthologs files (orthologs.out)', required=True) 
    parser.add_argument('-n', help='Initial copy number of genes at the mrca', type=int, required=False, default=1)
    args = parser.parse_args()

    ############## ---- processFamGenes ---- ################
    print 'Processing families'
    gsMap, geneNums = readGeneSpeciesMap(args.map)  # Gene-Species map, and gene Numbers
    species = readSpecies(args.database)            # List of species we're looking at
    famD = sortFamData(args.families)               # Silix family data

    famData = famInfo(famD, gsMap, species)         # Family data in tuple form
    
    adjInfo = adjacencyInfo(args.geneOrder, geneNums)   # Gene adjacency information

    ############## ---- dupDel ---- ################
    # dupDel models currently not fully utilized.  This will probably change in further iterations.
    print 'Calculating dupDel models'

    tree = readTree(args.tree)                  # Phylogenetic tree
    dupDelResults= dupDelAll(tree, famData, args.b, args.c, args.n)

    ############## ---- Group Cost ---- ################
    print 'Making groups'

    familyData = dupDelResults # Make families based on dupDel results.  Really just needs a list of families with mrca info.
    groupsList = initializeGroups(familyData, len(species))  # List of groups.
    famGroupL = setFamGroupDict(groupsList)                  # List of Family # <-> Group #
    
    famSpAdjD = GFSdict(famD, gsMap, speciesDict(args.o, species), geneNums, adjInfo)   # Adjacency information
    leafCache = memoLeafList(tree, {})  # Precaclate the leaves for each node in the tree. Saves a lot of recursive calls.

    print 'Memory usage: %s (mb)' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000)
    ######### Distances Matrix ##########
    print 'Calculating distances'

    try:
        # If htrans has generated the initial distances before, it will try to reuse them.  Delete or rename
        # initialDistances.pickle to redo them.

        print 'Reading in distances from file'
        start = time.clock()
        distancesDict = readPickle('initialDistances.pickle')
        end = time.clock()
        print 'Reading the pickle took '+str(end-start)+' seconds.'
    except:
        # If initial distances are not found, calculate them and save them to a file.
        print 'No file found, calculating distances for the first time.'
        distancesDict = initializeMagicalMatrix(groupsList, tree, famSpAdjD, famGroupL, leafCache)
        print 'Finished initializing distances'
        print 'Memory usage: %s (mb)' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000)
        try:
            writePickle(distancesDict, 'initialDistances.pickle')   # Save as a compressed binary file
        except:
            print 'Failed to write distances to pickle.'
    
    print 'Merging...'
    print 'Starting with '+str(len(groupsList))+' groups.'
    # Stringently merge all groups until no more groups can be merged, or we do 10000 merges.  Will add the limit as a
    # command line argument in the future.

    distancesDict, groupsList, famGroupL, temp = mergeWrapper(distancesDict, groupsList, 0, 0, tree, famSpAdjD, famGroupL, leafCache, 200, 10000)
    print 'Merged '+str(temp)+' times.'
    print 'There are '+str(numGroups(groupsList))+' remaining \n'
    print 'Memory usage: %s (mb)' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000)

    writeGroups('groupsL.txt', groupsList)

    # Code for making heatmaps of the groups.  Do not remove.
    # print 'Making heatmaps'
    # for index, matrixD in enumerate(matricies):
    #     print 'Heatmap #'+str(index)
    #     print 'Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    #     matrix = makeArray(matrixD, len(famGroupL))
    #     # np.savetxt('heatmap'+str(index)+'.txt', matrix)
    #     fig, ax = plt.subplots()
    #     heatmap = ax.pcolor(matrix, cmap=plt.cm.Blues)#, edgecolors='k')

    #     # put the major ticks at the middle of each cell
    #     ax.set_xticks(np.arange(matrix.shape[0])+0.5, minor=False)
    #     ax.set_yticks(np.arange(matrix.shape[1])+0.5, minor=False)

    #     cbar = plt.colorbar(heatmap)

    #     # want a more natural, table-like display
    #     ax.invert_yaxis()
    #     ax.xaxis.tick_top()

    #     ax.set_xticklabels([''], minor=False, fontsize=6)
    #     ax.set_yticklabels([''], minor=False, fontsize=6)
        
    #     plt.savefig('heatmap'+str(index)+'.png', bbox_inches='tight')
        # plt.clf()
 

############## ---- General helper functions ---- ################

def readGroups(infile):
    '''Read groups from a file'''
    groups = []
    with open(infile, 'r') as handle:
        while True:
            line = handle.readline().rstrip()
            if not line:
                break
            if line == 'None':
                groups.append(None)
            else:
                groups.append(Group.withList(ast.literal_eval(line)))
    return groups

def writeGroups(outfile, groups):
    '''Writes groups to a file using the str representation'''

    with open(outfile, 'w+') as f:
        for group in groups:
            if group == None:
                f.write('None\n')
            else:
                f.write(str(group)+'\n')


def numDists(distances, cutoff):
    '''Returns the number of elements in a distance dictionary that are less than or equal to the cut off'''
    count = 0
    for value in distances.itervalues():
        if value[0][0] <= cutoff:
            count+=1
    return count


def numGroups(groups):
    '''Returns the number of groups that are not None'''
    count = 0
    for group in groups:
        if group != None:
            count+=1
    return int(count)

def makeArray(distancesD, size):
    '''Make a numpy heatmap'''
    heat = np.zeros((size, size))

    for (x,y), value in distancesD.iteritems():
        if value != None:
            heat[x,y] = 1.0/(1.0+value[0][0])
            heat[y,x] = 1.0/(1.0+value[0][0])

    return heat

def writeDistances(dists, outfile):
    '''Write out distance matrix to text file'''
    with open(outfile, 'w+') as f:
        for key, value in dists.iteritems():
            f.write(str(key)+'\t'+str(value)+'\n')
    return True

def writePickle(dists, outfile):
    '''Write out dictionary as a compressed binary file'''
    with open(outfile, 'wb') as handle:
        pickle.dump(dists, handle, pickle.HIGHEST_PROTOCOL)
    return True

def readPickle(infile):
    '''Read in compressed binary file as a dictionary'''
    with open(infile, 'rb') as handle:
        return pickle.load(handle)


def getFiles(directory, pattern):
    '''Given a directory, gets the files ending with the pattern'''

    files = filter(lambda x: x.endswith(pattern), os.listdir(directory))
    return [directory + s for s in files]


def createGeneSpeciesMap(broadDir, ncbiDir):
    '''Read in protein file names and makes gene-species  and gene name <-> number maps'''
    
    protFiles = getFiles(broadDir, '.simp.faa') + getFiles(ncbiDir, '.simp.faa')

    geneSpeciesMap = {} # dictionary of {gene : species}
    numbers = {}        # dictionary of gene number <-> gene name
    count = 1

    for fileName in protFiles:
        species = fileName.split('/')[-1].split('.')[0]
            
        strainInfo=loadFasta(fileName)
        for gene in strainInfo:
            gName=gene[0][1:].split()[0]
            geneSpeciesMap[gName] = species
            numbers[count] = gName
            numbers[gName] = count
            count += 1

    return geneSpeciesMap, numbers





############## ---- processFamGenes ---- ################
def readGeneSpeciesMap(infile):
    '''Read in the geneSpeciesMap.
    Returns a dictionary of {gene : species}
    and a dictionary of gene number <-> gene name'''
    data = {}
    numbers = {}
    count = 1
    with open(infile) as f:
        while True:
            line = f.readline()
            if not line:
                break
            line = line.rstrip().split()
            for gene in line[1:]:
                data[gene] = line[0]
                numbers[count] = gene
                numbers[gene] = count
                count+=1
    return data, numbers

############## ---- processFamGenes ---- ################
def readSpecies(infile):
    '''Read in the database of species.
    Returns list of species'''
    species = []
    with open(infile) as f:
        while True:
            line = f.readline()
            if not line:
                break
            species.append(line.strip())
    return species

############## ---- processFamGenes ---- ################
def sortFamData(famData):
    '''Reads fam.out and outputs family# + genes'''

    data = {}
    with open(famData, 'r') as f:
        while True:
            line=f.readline()
            if not line:
                break
            line=line.rstrip().split()

            if int(line[0]) not in data.keys():
                data[int(line[0])] = []
            data[int(line[0])].extend(line[1:])

    return data
     
############## ---- processFamGenes ---- ################
def famInfo(famGeneMap, gsMap, species):
    '''For each family, count the number of times the genes from each species occurs in the family.
    Families index from 1 curtesy of SiLiX.'''
    data = {}
    for family, genes in famGeneMap.iteritems():
        data[family] = [0]*len(species)
        for gene in genes:
            if gsMap[gene] in species:
                data[family][species.index(gsMap[gene])]+=1

    return tuple([0]+[tuple(data[key]) for key in sorted(data.iterkeys())])

############## ---- processFamGenes ---- ################
def adjacencyInfo(geneOrder,geneNumbers):
    '''Generate adjacency info.
    Takes: 
        gene order information
        dictionary of gene name <-> gene number
    '''

    result = [()]*(len(geneNumbers)/2+1)
        
    with open(geneOrder,'r') as f:
        while True:
            s=f.readline()
            if not s: 
                break
            L=s.rstrip().split()[1:]
            for i in range(len(L)):
                if L[i]!='|||':
                    if i==0:
                        result[geneNumbers[L[i]]]=(geneNumbers[L[i+1]],)
                    elif i == len(L)-1:
                        if L[i-1] != '|||':
                            result[geneNumbers[L[i]]]=(geneNumbers[L[i-1]],)
                        else:
                            result[geneNumbers[L[i]]]=()
                    else:
                        if L[i-1] != '|||' and L[i+1] !='|||':
                            result[geneNumbers[L[i]]]=(geneNumbers[L[i-1]],geneNumbers[L[i+1]])
                        elif L[i-1] != '|||':
                            result[geneNumbers[L[i]]]=(geneNumbers[L[i-1]],)
                        elif L[i+1] != '|||':
                            result[geneNumbers[L[i]]]=(geneNumbers[L[i+1]],)
                        else:
                            result[geneNumbers[L[i]]]=()
        
    return result

############## ---- MRCA ---- ################

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
      #print 'Found tip #'+str(tree[0])+'\t'+str([family[tree[0]], tree[0]])
      if family[tree[0]] == 0:
        return -1
      else:
        return tree[0]
    else:
      r1 = mrca(tree[1], family)
      r2 = mrca(tree[2], family)

      #print 'Comparing r1 = '+str(r1)+' and r2 = '+str(r2)

      if r1 == -1 and r2 == -1:
        #print 'Returning '+str(-1)
        return -1
      elif r1 >= 0 and r2 >= 0:
        #print 'Returning '+str(tree[0])
        return tree[0]
      elif r1 >= 0:
        #print 'Returning '+str(r1)
        return r1
      else:
        #print 'Returning '+str(r2)
        return r2




############## ---- dupDel ---- ################

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
    '''Memo function to avoid repeated calculations.  Needs more
        testing to determine how much time this saves.
    '''
    cache = {}
    @ wraps(func)
    def wrap(*args):
        if args not in cache:
            cache[args] = func(*args)
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
    '''Find and return all the leaves in the given tree'''
    if Tree[1]==():
        return [Tree[0]]
    else:
        return leafList(Tree[1])+leafList(Tree[2])

def memoLeafList(tree, cache):
    '''Builds a cache of leafList results for every subtree of the given tree'''
    if tree[1] != () and tree[1] not in cache:
        cache = memoLeafList(tree[1], cache)
    if tree[2] != () and tree[2] not in cache:
        cache = memoLeafList(tree[2], cache)
    if tree not in cache:
        cache[tree] = leafList(tree)
    return cache

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

def dupDelAll(tree, familyTuples, delCost, dupCost, currentcopynum):
    '''Runs dupDel for every gne family'''
    results=[None]
    it = iter(familyTuples[1:])
    while True:
        try:
            fam = it.next()
        except:
            break
        mrcaF = mrca(tree, fam)
        sub = subtree(mrcaF, tree)
        results.append((mrcaF, dupDel(sub, fam, delCost, dupCost, currentcopynum)))

    return results


############## ---- Group Cost ---- ################

def initializeMagicalMatrix(groups, tree, famSpAdjD, famGroupL, leafCache):
    '''Generate a distances dictionary for pairwise group costs'''

    print 'Initializing distances matrix'
    numG = numGroups(groups)
    total = int((numG*(numG-1))/2)
    fiveP = int(total/20)
    oneP = int(total/100)
    distancesD = {}
    progress = 0
    part = 0
    print 'Updates every '+str(oneP)+' calculations.'
    start = time.clock()

    # Iterate over every group that isn't None
    for x in xrange(1, len(groups)):
        for y in xrange(x+1, len(groups)):
            if groups[x] != None and groups[y] != None:
                distancesD[(x,y)] = calcDist(groups[x], groups[y], tree, famSpAdjD, famGroupL, groups, leafCache)
                progress += 1
            if progress == oneP:
                part+=1
                print str(part)+'%'
                # sys.stdout.write(str(part)+'%...')
                progress = 0
    end = time.clock()
    print 'Initializing matrix took '+str(end-start)+' seconds'

    return distancesD

def mergeWrapper(distancesD, groups, dCutoff, oCutoff, tree, famSpAdjD, famGroupL, leafCache, itsPerRecalc, limit):
    '''Wrapper function for mergeOneByOne.
    Takes care of repeated iterations and recalculating the full matrix every so often.
    '''
    begin = time.clock()
    mergeMore = True
    iterations = 0
    count = 0
    count2 = 0
    onePer = limit/100
    progress = 0
    while mergeMore and count < limit:
        start = time.clock()
        distancesD, groups, famGroupL, mergeMore = mergeOneByOne(distancesD, groups, dCutoff, oCutoff, tree, famSpAdjD, famGroupL, leafCache)
        end = time.clock()
        print 'Merge '+str(count)+' took '+str(end-start)+'seconds.'
        print 'Running for '+str((time.clock()-begin)/60.0) + ' minutes.\n'

        iterations += 1
        count += 1
        count2 += 1
        if iterations == itsPerRecalc:
            print 'Writing out groups'
            with open('groupsL.txt', 'w+') as f:
                for group in groups:
                    if group == None:
                        f.write('None\n')
                    else:
                        f.write(group.printG()+'\n')
            iterations = 0
            print 'Recalculating the full matrix'
            distancesD = initializeMagicalMatrix(groups, tree, famSpAdjD, famGroupL, leafCache)
            
        if count2 == onePer:
            progress += 1
            print str(progress)+'%'
            count2 = 0
            print 'Number of groups left = '+str(numGroups(groups))


    return distancesD, groups, famGroupL, count





def mergeOneByOne(distancesD, groups, dCutoff, oCutoff, tree, famSpAdjD, famGroupL, leafCache):
    '''Merge two groups greedy style, then recalculate.  Since distancesD is a dictionary,
       the order of groups it analyzes is psudo-random.  We consider this a feature.
    '''
    reX = 0
    reY = 0
    recalc = False

    start = time.clock()
    for key, value in distancesD.iteritems():
        if value[0][0] <= oCutoff:
            x = key[0]
            y = key[1]
            print 'Merging groups '+str(x)+' and '+str(y)
            groups[x].mergeGroup(groups[y], distancesD[(key[0],key[1])][0][1])
            for familyNum in groups[y].getFamilies():
                famGroupL[familyNum] = x
            groups[y] = None
            reX = x
            reY = y
            recalc = True
            break
    end = time.clock()
    print 'Finding and merging took '+str(end-start)+' seconds'

    if recalc:
        print 'Recalculating...'
        newDists = list(set([(reX, y) for y in xrange(reX+1, len(groups)) if y != reY and groups[y] != None]+[(x, reX) for x in xrange(1, reX) if groups[x] != None]))
        delDists = list(set([(x, reY) for x in xrange(1, reY) if groups[x] != None]+[(reY, y) for y in xrange(reY+1,  len(groups)) if groups[y] != None]))
        
        for key in delDists:
            del distancesD[key]
        
        for x,y in newDists:
            distancesD[(x,y)] = calcDist(groups[x], groups[y], tree, famSpAdjD, famGroupL, groups, leafCache)
        
        return distancesD, groups, famGroupL, recalc
    else:
        return distancesD, groups, famGroupL, recalc


    
def calcDist(groupA, groupB, tree, famSpAdjD, famGroupL, groupD, leafCache):
    '''Compares the pairOrderCost value of two groups.
       Currently set up for stringent merging and only calculates
       cost for groups with the same mrca.
    '''

    defaultCost = float('inf')

    # mrcaCost = moveMRCAcost(tree, groupA.getMrcag(), groupB.getMrcag())
    if groupA.getMrcag() == groupB.getMrcag() and groupA.getMrcag() != tree[0]:
        return groupCost(groupA, groupB, tree, famSpAdjD, famGroupL, groupD, leafCache)
    else:
        return [(defaultCost,0)]


def initializeGroups(familyData, numSpecies):
    '''Make every family into a group.  Assumes the list of families starts with 'None'
    as the first element.'''

    groups = [None]

    # index is both the group idNum and family number
    for index, value in enumerate(familyData[1:], 1):
        groups.append(Group(index, value[0], [index], numSpecies, index, index))    #groups.append(Group(value, index, 6))

    return groups

##### - Family # -> Group # -- #####
def setFamGroupDict(groups):
    '''Initializes the initial map of family number to group number. Now a List!
    This can also be used to update the map based on current groups
    '''

    # famGroupDict = {family: g.getIdNum() for g in groups if g != None for family in g.getFamilies()}
    famGroupDict = [famNum for famNum in range(0, len(groups))]
    famGroupDict[0] = None
    return famGroupDict



def readFamilies(filename):
    '''Reads in the results of dupDel and stores families as a dictionary.
    Note, the first element is None to index families from 1.'''

    families = [None]

    with open(filename, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            temp = list(ast.literal_eval(line))
            families.append(temp)

    return families


def groupCost(groupA, groupB, tree, famSpAdjD, famGroupL, groupD, leafCache):
    '''Find the cost of merging two groups'''
    pairs = [(groupA.getFront(), groupB.getFront()), (groupA.getBack(), groupB.getFront()), 
                (groupA.getFront(), groupB.getBack()), (groupA.getBack(), groupB.getBack())] 

    costs = [] #float('inf')
    memo = {}
    for index, adj in enumerate(pairs):
        if adj in memo:
            costs.append((memo[adj], index))
        else:
            cost = pairOrderCost(tree, tree, 1, famSpAdjD, famGroupL, groupD, groupA, groupB, adj, {}, leafCache)
            # if cost > 0:
                # print (cost, index, groupA.getIdNum(), groupB.getIdNum())
            costs.append((cost, index))
            memo[adj] = cost
    results = sorted(costs, key=itemgetter(0,1))
    # if results[0][0] > 0:
    #     print results
    return results


def mergeLists(order, l1, l2):
    '''Merge l1 and l2 in the correct order based on the pairs in groupCost'''
    lb = l2
    if len(l2) > 1:
        lb = l2.reverse()


    return {0: lb+l1, 1: l1+l2, 2: l2+l1, 3: l1+lb}[order]

def pairwise(iterable):
    "s -> (s0,s1), (s2,s3), (s4, s5), ..."
    a = iter(iterable)
    return izip(a, a)
############## ---- END ---- ################

############## ---- Order Cost Preprocessing ---- ################

def geneSpeciesDict(data):
    '''Takes genesSpeciesMap.txt -> gene-species dictionary'''
    result={}
    with open(data,'r') as f:
        while True:
            s=f.readline()
            if s=='':
                break
            s=s.rstrip('\n').split()
            result[s[0]]=s[1:]
            for i in range(1,len(s)):
                result[s[i]]=s[0]
    return result

def speciesDict(orthologs, speciesList):
    '''Reads in orthologs.out and number of species and creates a dictionary of species name <-> number'''

    specD={}
    with open(orthologs, 'r') as f:
        while True:
            try:
                line = f.readline()
            except:
                break
            dat = line.rstrip('\n').split()
            if len(dat) > 2:
                break
            if dat[1] in speciesList:
                specD[int(dat[0])] = dat[1]
                specD[dat[1]] = int(dat[0])

    return specD

def GFSdict(famNumGeneMap, geneSpeciesDict, speciesDict, geneDict, adjInfo):
    '''
    Takes 5 data structures:
        1) famNumGeneMap -- sortFamData()
        2) geneSpeciesDict -- geneSpeciesDict()
        3) speciesDict -- speciesDict()
        4) geneDict -- readGeneSpeciesMap() 2nd outputs
        5) adjInfo -- adjacencyInfo()
    
    Outputs FamSpAdjD
    '''

    results=defaultdict(list)
    temp={}
    famSpAdjD = {}

    for key, value in famNumGeneMap.iteritems():
        for gene in value:
            speciesName = geneSpeciesDict[gene]
            speciesNum = speciesDict[speciesName]
            geneNum = geneDict[gene]
            famNum = key
            temp[geneNum] = famNum
            results[(famNum, speciesNum)].append(geneNum)
        
    for key, value in results.iteritems():
        adjacentFam = []
        for gene in value:
            adjacent = adjInfo[gene]
            for adj in adjacent:
                if adj in temp:
                    adjacentFam.append(temp[adj])
                else:
                    adjacentFam.append('*')
        famSpAdjD[key] = adjacentFam
                    
    return famSpAdjD

def pairOrderCost(tree,originTree,cost,FamSpAdjD,famGroupD,groupL,groupA,groupB,pair,memo, leafCache):
    '''Calculates the cost of putting two groups together the mrca's genome'''

    originHtrans=groupA.getMrcag()
    groupAFam=groupA.getFamilies()
    groupBFam=groupB.getFamilies()
    if tree[1]==(): #base case: At a leaf:
        if (remainFam(groupAFam,tree[0],FamSpAdjD) == []) or (remainFam(groupBFam,tree[0],FamSpAdjD) == []):

            return 0
        else:
            # print "!"
            repA=representative(groupAFam,tree[0],pair[0],FamSpAdjD)
            repB=representative(groupBFam,tree[0],pair[1],FamSpAdjD)
            adjFList=FamSpAdjD[(repA,tree[0])]
            
            rAdjFList=[] # might have called this originallyAdjacent
                         # or something like that. Its intended to be
                         # a list for things we figure out were really
                         # next to repA, but have been interrupted by
                         # later htrans events etc.
            
            for adjF in adjFList:
                currentGroupID=famGroupD[adjF]
                
                currentGroup=groupL[currentGroupID] # a potentially intervening group
                currentHtrans=currentGroup.getMrcag()

                # Its possible that a later horizontal transfer event
                # might have interrupted a pre-existing group (the
                # group we want to merge here). We could recognize
                # this if the intervening, unrelated group has an
                # mrcag lower in the tree than the groups we're
                # considering. (we just check groupA here).

                if lowerThanOrigin(subtree(originHtrans,originTree),originHtrans,currentHtrans):
                    # print "!",adjF,tree[0]
                    otherEndFam=otherEnd(currentGroup.getFamilies(),adjF,tree[0],FamSpAdjD)
                    
                    otherEndAdjFL=FamSpAdjD[(otherEndFam,tree[0])]
                    for otherEndAdjF in otherEndAdjFL:
                        if (otherEndAdjF != repA) and (otherEndAdjF not in currentGroup.getFamilies()):
                            rAdjFList.append(otherEndAdjF)
                else:
                    rAdjFList.append(adjF)
            if repB in rAdjFList:
                return 0
            else:
                # print "node:",tree[0]
                return cost
                
    elif (tree,pair) in memo: return memo[(tree,pair)]
    else: #General case: At a subtree:

        leftLeaves = leafCache[tree[1]] #leafList(tree[1])

        # We're imagining different order change events on this
        # branch. options represents the possibilities we'll
        # consider. Of course, it only makes sense to try orderchange
        # events putting pair[0] next to things it actually occurs
        # next to in the data for the tips. (also we'll only consider
        # things in group b)

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
            leftSubCost=pairOrderCost(tree[1],originTree,cost,FamSpAdjD,famGroupD,groupL,groupA,groupB,(pair[0],option),memo, leafCache)
            if option == pair[1]:
                leftCost=leftSubCost
            else:
                leftCost=cost+leftSubCost
            #print leftCost
            if leftCost < minLeftCost:
                minLeftCost=leftCost
        #print 'minLeftCost:',minLeftCost
        
        rightLeaves = leafCache[tree[2]] #leafList(tree[2])
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
            rightSubCost=pairOrderCost(tree[2],originTree,cost,FamSpAdjD,famGroupD,groupL,groupA,groupB,(pair[0],option),memo, leafCache)
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
    '''Go to the other end of an intervening group.'''
    if len(group)==1:
        return group[0]
    else:
        if fam == group[0]:
            if (group[-1],leaf) in FamSpAdjD:
                return group[-1]
            else:
                return otherEnd(group[:-1],fam,leaf,FamSpAdjD)
        else:
            if (group[0],leaf) in FamSpAdjD:
                return group[0]
            else:
                return otherEnd(group[1:],fam,leaf,FamSpAdjD)
    

def lowerThanOrigin(Tree,originNode,currentNode):
    
    if Tree[0]==currentNode and Tree[0]!=originNode:
        return True
    elif Tree[1]==():
        return False
    else:
        return lowerThanOrigin(Tree[1],originNode,currentNode) or lowerThanOrigin(Tree[2],originNode,currentNode)
        

def remainFam(group,leaf,FamSpAdjD):
    '''Returns a list of the families in group which have not been deleted.'''
    result=[]
    for fam in group:
        if (fam,leaf) in FamSpAdjD:
            result.append(fam)
    return result
                
def representative(group,leaf,fam,FamSpAdjD):
    '''In trying to merge this group, we've specified one family which
will be merged with another group. That is given in the fam argument
here. If that family is not deleted, we're fine and we just return
that. If the family is deleted, we need to determine which family
should act as its surrogate, based on adjacency.
    '''
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

def findMRCA(tree,mrcaA,mrcaB):
    if (mrcaA == tree[0]) or (mrcaB==tree[0]):
        return tree[0]
    elif ((mrcaA in descendantNodes(tree[1][0],tree)) and (mrcaB in descendantNodes(tree[2][0],tree))) or ((mrcaB in descendantNodes(tree[1][0],tree)) and (mrcaA in descendantNodes(tree[2][0],tree))):
        return tree[0]
    elif (mrcaA in descendantNodes(tree[1][0],tree)) and (mrcaB in descendantNodes(tree[1][0],tree)):
        return findMRCA(tree[1],mrcaA,mrcaB)
    else:
        return findMRCA(tree[2],mrcaA,mrcaB)
        

def moveMRCAcost(tree,mrcaA,mrcaB):
    commonMRCA=findMRCA(tree,mrcaA,mrcaB)
    mrcaAcost=subMoveMRCAcost(tree,commonMRCA,mrcaA)
    mrcaBcost=subMoveMRCAcost(tree,commonMRCA,mrcaB)
    return (mrcaAcost+mrcaBcost,commonMRCA)

def subMoveMRCAcost(tree,commonMRCA,mrcaX):
    if commonMRCA==mrcaX:
        return 0
    elif mrcaX in descendantNodes(subtree(commonMRCA,tree)[1][0],tree):
        return 1+subMoveMRCAcost(tree,subtree(commonMRCA,tree)[1][0],mrcaX)
    else:
        return 1+subMoveMRCAcost(tree,subtree(commonMRCA,tree)[2][0],mrcaX)

############## ---- You made it! ---- ################

if __name__ == "__main__":
   main(sys.argv)
