'''
Wrapper for htrans
'''

import mrca, sys, io, ast, math, os, argparse
from group import Group
# from orderCost1 import *
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
# import matplotlib




def main(argv):
    # matplotlib.use()
    plt.ioff()
    # matplotlib.use('Agg')
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
    print 'Calculating dupDel models'

    tree = readTree(args.tree)                  # Phylogenetic tree
    dupDelResults= dupDelAll(tree, famData, args.b, args.c, args.n)

    ############## ---- Group Cost ---- ################
    print 'Making groups'

    familyData = readFamilies('dupDelSmaller.txt') #dupDelResults
    groupsD = initializeGroups(familyData)
    famGroupL = setFamGroupDict(groupsD)
    # tree = readTree('testATree')
    famSpAdjD = GFSdict(famD, gsMap, speciesDict(args.o, species), geneNums, adjInfo)
    leafCache = memoLeafList(tree, {})
    ######### Distances Matrix ##########
    print 'Calculating distances'

    matricies = []
    distancesDict = initializeMagicalMatrix(groupsD, tree, famSpAdjD, famGroupL, leafCache)
    print 'Distances #1 done'
    # matricies.append(makeArray(distancesDict, len(famGroupL)))
    matricies.append(distancesDict)

    print 'Merging...'
    count = 1
    mergeMore = True
    while mergeMore and count < 6:
        print 'Merge #'+str(count)
        distancesDict, groupsD, temp = mergeGroups(distancesDict, groupsD, 0, 0, tree, famSpAdjD, famGroupL, leafCache)
        # matricies.append(makeArray(distancesDict, len(famGroupL)))
        matricies.append(distancesDict)
        if temp <= 1:
            mergeMore = False
        count+=1

    # print 'Making heatmaps'
    # for index, matrix in enumerate(matricies):
    #     print 'Heatmap #'+str(index)
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
    #     plt.clf()
 
    
def makeArray(distancesD, size):
    heat = np.zeros((size, size))

    for (x,y), value in distancesD.iteritems():
        if value == None:
            heat[x,y] = 0.0
            heat[y,x] = 0.0
        else:
            heat[x,y] = 1.0/(1.0+value[0][0])
            heat[y,x] = 1.0/(1.0+value[0][0])

    return heat


############## ---- geneSpecies Map -- IGNORE ---- ################
def loadFasta(filename):
    """Load fasta or multifasta, return list of tuples (header,seq)."""
    with open(filename,'r') as f:
        header=None
        outL=[]
        tempSeqL=[]
        while True:
            Str=f.readline()
            if Str=="":
                outSeq="".join(tempSeqL)
                tempSeqL=[]
                outSeq="".join(outSeq.split()) # remove all whitespace
                outL.append((header,outSeq))
                break
            if Str[0]==">":
                # this is new header, put together previous header,seq, and move one
                if header!=None:
                    outSeq="".join(tempSeqL)
                    tempSeqL=[]
                    outSeq="".join(outSeq.split()) # remove all whitespace
                    outL.append((header,outSeq))
                header=Str[:-1]
            else:
                # this is a seq line,
                tempSeqL.append(Str)
    return(outL)    

def getFiles(directory, pattern):
    '''Given a directory, gets the files ending with the pattern'''

    files = filter(lambda x: x.endswith(pattern), os.listdir(directory))

    return [directory + s for s in files]

    # broadFiles = filter(lambda x: x.endswith('.simp.faa'), os.listdir(broadDir))
    # ncbiFiles = filter(lambda x: x.endswith('.simp.faa'), os.listdir(ncbiDir))

    # return [broadDir + s for s in broadFiles] + [ncbiDir + t for t in ncbiFiles]

def createGeneSpeciesMap(broadDir, ncbiDir):
    '''Read in protein file names and makes gene-species  and gene name <-> number maps'''
    
    protFiles = getFiles(broadDir, '.simp.faa') + getFiles(ncbiDir, '.simp.faa')

    geneSpeciesMap = {} # dictionary of {gene : species}
    numbers = {}        # dictionary of gene number <-> gene name
    count = 1

    for fileName in protFiles:
        species = fileName.split('/')[-1].split('.')[0]
     # with open(protFile,'r') as f:
     #    while True:
     #        s=f.readline()
     #        if not s:
     #            break
     #        s=s.rstrip('\n')
     #        # L=s.rstrip().split('/')
     #        species=s.rstrip().split('/')[-1].split('.')[0]
            
        strainInfo=loadFasta(fileName)
        for gene in strainInfo:
            gName=gene[0][1:].split()[0]
            geneSpeciesMap[gName] = species
            numbers[count] = gName
            numbers[gName] = count
            count += 1

    return geneSpeciesMap, numbers

############## ---- Gene Order -- IGNORE ---- ################
'''
ls seq/broad/*genome_summary_per_gene.txt > geneFilesBroadTemp.txt
ls seq/broad/*simp.faa > protFilesBroadTemp.txt # needed for names
ls seq/ncbi/*simp.faa > protFilesNCBITemp.txt # needed for names
python code/extractGeneOrder.py geneFilesBroadTemp.txt protFilesBroadTemp.txt protFilesNCBITemp.txt > full/geneOrder.txt
'''
def extractGeneOrder(geneBroadDir,protBroadDir,NCBIDir):
    
    geneBroadFile = getFiles(geneBroadDir, 'genome_summary_per_gene.txt')
    protBroadFile = getFiles(protBroadDir, 'simp.faa')
    ncbiFile = getFiles(NCBIDir, 'simp.faa')    

    f=open(geneBroadFile,'r')
    q=open(protBroadFile,'r')

    while True:
        tempDict={} 
        tempOrder=[]
        allNames=[]
        s=f.readline()
        t=q.readline()
        if s=='':
            break
        s=s.rstrip('\n')
        L=s.rstrip().split('/')
        k=L[2].rstrip('genome_summary_per_gene_.txt')
        print k+'_',
        t=t.rstrip('\n')
        
      
        
        
        r=open(t,'r')
        while True:
            u=r.readline()
            if u=='': break
            if u[0]=='>':
                u=u.split()[0][1:]
                allNames.append(u)
        
        r.close()

        g=open(s,'r')
        h=g.readline()
        while True:
            h=g.readline()
            if h=='':break
            geneInfo=h.rstrip().split('\t')
            if geneInfo[0] in allNames:
                chromosomeNum=int(geneInfo[8])
                if chromosomeNum in tempDict:
                    tempDict[chromosomeNum].append((geneInfo[0],int(geneInfo[4])))
                else:
                    tempDict[chromosomeNum]=[(geneInfo[0],int(geneInfo[4]))]
                    tempOrder.append(chromosomeNum)
        for chromosome in tempOrder:
            tempDictGene={}
            tempOrderGene=[]
            tempList=tempDict[chromosome]
            for gene in tempList:
                tempDictGene[gene[1]]=gene[0]
                tempOrderGene.append(gene[1])
            tempOrderGene=sorted(tempOrderGene)
            for startLocation in tempOrderGene:
                print tempDictGene[startLocation],
            print '|||',
        print
    f.close()
    
    f=open(NCBIFile,'r')
     
    while True:
        s=f.readline()
        if s=='':
            break
        s=s.rstrip('\n')
        L=s.rstrip().split('/')
        k=L[2].split('.')[0]
        print k,
            
        g=open(s,'r')
        tempDict={} 
        tempOrder=[]
        while True:
            h=g.readline()
            if h=='':break
            if h[0]=='>':
                
                j=h.split()[1]
                j=j.split('-')[0]
                if j[0]=='c':
                    j=j[1:]
                j=int(j)
                name=h.split()[0][1:]
                tempDict[j]=name
                tempOrder.append(j)
        tempOrder=sorted(tempOrder)
        for startLocation in tempOrder:
            print tempDict[startLocation],
        print
    return




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
    Families index from 1, not sure why.'''
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
    '''Memo function for dupDel to avoid repeated calculations'''
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
############## ---- END ---- ################


############## ---- Group Cost ---- ################

def initializeMagicalMatrix(groups, tree, famSpAdjD, famGroupL, leafCache):
    print 'Initializing distances matrix'

    distancesD = {}

    for x in range(1, len(groups)):
        for y in range(x+1, len(groups)):
            if groups[x] != None and groups[y] != None:
                distancesD[(x,y)] = calcDist(groups[x], groups[y], tree, famSpAdjD, famGroupL, groups, leafCache)

    return distancesD

def mergeGroups(distancesD, groups, dCutoff, oCutoff, tree, famSpAdjD, famGroupL, leafCache):

    # print 'Merging groups'
    # print len(groups), ' total groups'

    recalculate = []
    for x in range(1, len(groups)):
        for y in range(x+1, len(groups)):
            if x not in recalculate and y not in recalculate and (x,y) in distancesD:
                if distancesD[(x,y)][0][0] <= dCutoff and distancesD[(x,y)][0][0] <= oCutoff:
                    if groups[x] != None and groups[y] != None:
                        groups[x].mergeGroup(groups[y], distancesD[(x,y)][0][1])
                        groups[y] = None
                        famGroupL[y] = x
                        recalculate.extend([x,y])

    if len(recalculate) > int(len(distancesD)/3):
        newDists = initializeMagicalMatrix(groups, tree, famSpAdjD, famGroupL, leafCache)

    else:
        xss = recalculate[0::2]
        yss = recalculate[1::2]

        newDists = deepcopy(distancesD)
        for (x,y) in distancesD.iterkeys():
            if y in yss:
                del newDists[(x,y)]
            elif x in xss and groups[y] != None:
                newDists[(x,y)] = calcDist(groups[x], groups[y], tree, famSpAdjD, famGroupL, groups, leafCache)


    return newDists, groups, len(recalculate)/2



    
def calcDist(groupA, groupB, tree, famSpAdjD, famGroupL, groupD, leafCache):
    '''Compares the dupDel model and order cost of two groups'''

    # mrcaDif = 
    # ddDiff = calcDiff(groupA.getDuplications(), groupB.getDuplications())+calcDiff(groupA.getDeletions(), groupB.getDeletions())

    ordCost = groupCost(groupA, groupB, tree, famSpAdjD, famGroupL, groupD, leafCache)
    return ordCost
    # return (ddDiff, ordCost)

# def calcDiff(la, lb, n):
#     '''Take two duplication or deletions models and compare them.
#     Score starts at zero, and increases by 1 for each event difference.'''
#     diff = 0
#     for i in range(0,2*n-1):
#         for x in range(0, len(i)):
#             if la[i][x][0] != lb[i][x][0] or set(la[i][x][1]) != set(lb[i][x][1]):
#                 diff += 1
#     return diff


def initializeGroups(familyData):
    '''Make every family into a group.  Assumes the list of families starts with 'None'
    as the first element.'''

    groups = [None]

    for index, value in enumerate(familyData[1:], 1):
        groups.append(Group(value, index, 6))

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


def groupCost(groupA, groupB, tree, adjInfo, famGroupL, groupD, leafCache):
    '''Find the cost of merging two groups'''
    pairs = [(groupA.getFront(), groupB.getFront()), (groupA.getBack(), groupB.getFront()), 
                (groupA.getFront(), groupB.getBack()), (groupA.getBack(), groupB.getBack())] 

    costs = [] #float('inf')
    for index, adj in enumerate(pairs):
        cost = pairOrderCost(tree, tree, 1, adjInfo, famGroupL, groupD, groupA, groupB, adj, {}, leafCache)
        costs.append((cost, index))

    return sorted(costs, key=itemgetter(0,1))


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
                specD[dat[0]] = dat[1]
                specD[dat[1]] = dat[0]

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
            
            rAdjFList=[]
            
            for adjF in adjFList:
                currentGroupID=famGroupD[adjF]
                
                currentGroup=groupL[currentGroupID]
                currentHtrans=currentGroup.getMrcag()
                
                if lowerThanOrigin(subtree(originHtrans,originTree),originHtrans,currentHtrans) == True:
                    print "!",adjF,tree[0]
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
                print "node:",tree[0]
                return cost
                
    elif (tree,pair) in memo: return memo[(tree,pair)]
    else: #General case: At a subtree:
        leftLeaves=leafCache(tree[1])
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
        
        rightLeaves=leafCache(tree[2])
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
    
    if Tree[0]==currentNode and Tree[0]!=originNode:
        return True
    elif Tree[1]==():
        return False
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


if __name__ == "__main__":
   main(sys.argv)
