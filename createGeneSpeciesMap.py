import sys,fasta
protFile=sys.argv[1]

def createGeneSpeciesMap(protFile):
    f=open(protFile,'r')
     
    while True:
        s=f.readline()
        if s=='':
            break
        s=s.rstrip('\n')
        L=s.rstrip().split('/')
        k=L[2].split('.')[0]
        print k,
        
        strainInfo=fasta.load(s)
        for gene in strainInfo:
            q=gene[0]
            q=q[1:]
            q=q.split()
            q=q[0]
            print q,
        print

    return

createGeneSpeciesMap(protFile)
             
