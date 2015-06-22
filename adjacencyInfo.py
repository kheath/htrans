import sys,ast

geneOrder=sys.argv[1]
geneNumbers=sys.argv[2]


def adjacencyInfo(geneOrder,geneNumbers):
    f=open(geneNumbers,'r')
    s=f.readline()
    geneNum=ast.literal_eval(s)
    f.close()
    result=[]
    for i in range(len(geneNum)):
        result.append(())
        
    f=open(geneOrder,'r')
    
    while True:
        s=f.readline()
        if s=='': break
        L=s.rstrip().split()[1:]
        for i in range(len(L)):
            if L[i]!='|||':
                if i==0:
                    result[geneNum[L[i]]]=(geneNum[L[i+1]],)
                elif i == len(L)-1:
                    if L[i-1] != '|||':
                        result[geneNum[L[i]]]=(geneNum[L[i-1]],)
                    else:
                        result[geneNum[L[i]]]=()
                else:
                    if L[i-1] != '|||' and L[i+1] !='|||':
                        result[geneNum[L[i]]]=(geneNum[L[i-1]],geneNum[L[i+1]])
                    elif L[i-1] != '|||':
                        result[geneNum[L[i]]]=(geneNum[L[i-1]],)
                    elif L[i+1] != '|||':
                        result[geneNum[L[i]]]=(geneNum[L[i+1]],)
                    else:
                        result[geneNum[L[i]]]=()
        
    print result
    return

adjacencyInfo(geneOrder,geneNumbers)
            
