'''
Wrapper for htrans
'''

from group import Group
from orderCost1 import *
import mrca
from operator import itemgetter
from copy import deepcopy
import sys
import io, ast, math
from collections import Counter
from itertools import *
from group import Group
from copy import deepcopy
import os
from os import path



def main(argv):







    ############## ---- processFamGenes ---- ################
    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-f', '--families', help='Family file name (fam.out)', required=True)
    requiredNamed.add_argument('-m', '--map', help='Gene to species map (geneSpeciesMap.txt)', required=True)
    requiredNamed.add_argument('-d', '--database', help='File with species of interest (dbList.txt)', required=True)
    requiredNamed.add_argument('-g', '--geneOrder', help='Gene order file (geneOrder.txt)', required=True) 
    
    args = parser.parse_args()

    gsMap, geneNums = readGeneSpeciesMap(args.map)
    species = readSpecies(args.database)
    famData = sortFamData(args.families)

    data = famInfo(famData, gsMap, species)
    
    adjInfo = adjacencyInfo(args.geneOrder, geneNums)
    ############## ---- END ---- ################


############## ---- geneSpecies Map ---- ################
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
    
############## ---- END ---- ################

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






if __name__ == "__main__":
   main(sys.argv[1:])
