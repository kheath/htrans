import sys
geneSpeciesMap=sys.argv[1]

def numberGenes(geneSpeciesMap):
    result={}
    counter=1
    f=open(geneSpeciesMap,'r')
    while True:
        
        s=f.readline()
        if s=='':break
        L=s.rstrip().split()
        for i in range(1,len(L)):
            result[counter]=L[i]
            result[L[i]]=counter
            counter+=1
    print result
    return  

numberGenes(geneSpeciesMap)
