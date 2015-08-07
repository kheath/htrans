"""
    Author: Kevin Heath
    Date: 22 June 2015
    Generated for summer research in the Bush lab 
    
    Given Family# and Species# return gene name

    'Just in time to explore dank memes'
"""

def main():

    fams = readFamilies('result3.txt')

    # for k, v in fams.iteritems():
    #     print k, v

    spec = readSpecies('geneSpeciesMap.txt')

    # print spec[0]

    print findGeneName(187, 0, fams, spec)

def readFamilies(infile):
    '''Read in result3.txt and turn into dictionary'''

    data = {}
    with open(infile, 'r') as f:
        while True:
            ln = f.readline()
            if not ln:
                break
            temp = tuple(ln.strip().split())
            data[int(temp[0])] = tuple(temp[1:])

    return data

def readSpecies(infile):
    '''Read in species names & gene names -> species # and gene names dictionary'''

    data={}
    with open(infile, 'r') as f:
        for i, l in enumerate(f):
            data[i] = tuple(l.split()[1:])

    return data

def findGeneName(fam, species, famData, specData):
    return [gene for gene in famData[fam] if gene in specData[species]]


if __name__ == "__main__":
   main()

