"""
    Author: Kevin Heath
    Date: 18 June 2015
    Generated for summer research in the Bush lab 

    Find inversion costs

    'I like trains'
"""

import copy, ast, sys

def main(argv):

    familyGeneMap = readFamData('result3.txt')

    geneNameToNumber = readGeneData('fjdkal;fjadkl;.txt')


def hInvCost(tree, cost, adjPair, adjTup):
    '''Calculate the cost of inversions

    Cases:

    Edge Case - Leaf node
        Cost = 0 if matches, inf if not

    General Case
        Recurse on both branches
            Add costs together

        Enumerate all possibilies '''

    # Base case
    if tree[1] == ():


def speciesToGeneNum(familyNum, geneNameToNum):
    ''' Get adjacency info for leaf '''



def readFamData(infile):
    """ Reads a text file with Family number to gene names.
    Returns dictionary with key=Family Number, value = ["""

    dataStruct = {}
    with open(infile, 'r') as f:
      temp = tuple(f.readline().split())
      dataStruct[int(temp[0])] = temp[1:]

    return dataStruct

def readGeneData(infile):
    """ Reads a text file with Family number to gene names.
    Returns dictionary with key=Family Number, value = genes"""

    dataStruct = {}
    with open(infile, 'r') as f:
      temp = f.read()
      dataStruct = ast.literal_eval(temp)

    return dataStruct

def readSpeciesToGenes(infile, speciesNameToNumber):
    ''' Takes text file with: speciesName \t geneName1 \t geneName2
    speciesNameToNumber is a dictionary with key = name, value = number?

    Returns tuple of ((geneName1, geneName2,...), (...,...),...)
    Where each index corresponds to a species number '''

    with open(infile, 'r') as f:
        temp = f.readline().split()

    data = []

    for i in range(0, len(speciesNameToNumber)):
        data.append()


if __name__ == "__main__":
   main(sys.argv[1:])
