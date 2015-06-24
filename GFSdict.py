import sys,ast
result3=sys.argv[1]
geneSpeciesDictF=sys.argv[2]
speciesDictF=sys.argv[3]
geneDictF=sys.argv[4]
adjInfoF=sys.argv[5]


def GFSdict(result3,geneSpeciesDictF,speciesDictF,geneDictF,adjInfoF):
    f=open(geneSpeciesDictF,'r')
    s=f.readline()
    geneSpeciesDict=ast.literal_eval(s)
    f.close()

    f=open(speciesDictF,'r')
    s=f.readline()
    speciesDict=ast.literal_eval(s)
    f.close()

    f=open(geneDictF,'r')
    s=f.readline()
    geneDict=ast.literal_eval(s)
    f.close()
    
    f=open(adjInfoF,'r')
    s=f.readline()
    adjInfo=ast.literal_eval(s)
    f.close()
    
    f=open(result3,'r')
    result={}
    temp={}
    while True:
        s=f.readline()
        if s=="": break
        s=s.rstrip('\n').split()
        for i in range(1,len(s)):
            speciesName=geneSpeciesDict[s[i]]
            speciesNum=speciesDict[speciesName]
            geneNum=geneDict[s[i]]
            famNum=int(s[0])
            temp[geneNum]=famNum
            if (famNum,speciesNum) in result:
                result[(famNum,speciesNum)].append(geneNum)
            else:
                result[(famNum,speciesNum)]=[geneNum]
        
    for key in result:
        adjacentFam=[]
        for gene in result[key]:
            adjacent=adjInfo[gene]
            for adj in adjacent:
                if adj in temp:
                    adjacentFam.append(temp[adj])
                else:
                    adjacentFam.append('*')
        result[key]=adjacentFam
                    
    print result
    return

GFSdict(result3,geneSpeciesDictF,speciesDictF,geneDictF,adjInfoF)
