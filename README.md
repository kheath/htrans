# Horizontal Gene Transfer project
## Collaborative project for Eliot Bush

## Goals:  
+ Reconstruct gene family lineages in Escherichia Coli
+ Identify horizontal gene transfer events
+ ???
+ Profit

### superalignv2.py & phylibalign.py  

These two python scripts are very similar.  They take a file with blocks of genes (.afa) and convert it to
fasta format.

### fasta2phylip.sh  

This is a shell script that converts fasta files to phylip format.


### mrca.py

mrca(tree, family) returns the most recent common ancestor for a given gene family.

### testATree

The correct phylogenetic tree for our 5 sample species.  I typed this out by hand, there was no script directly involved in making this.

### processFamGenes.py

This is a collection of scripts put together to generate the needed information
to run dupDel.  It's a good idea to run this once and keep the files around.
It is important to note that the input files are made once with other scripts I
haven't looked into.

Input:
+ fam.out (silix results)
+ geneSpeciesMap.txt
+ dbList.txt (file with sample species)
+ geneOrder.txt 

Output:
+ famGenes.txt
+ famInfoResult.txt
+ adjacencyInfo.txt

Sample command: python processFamGenes.py -f fam.out -m geneSpeciesMap.txt -d dbList.txt -g geneOrder.txt

### dupDel.py

This script calculates the minimum cost and associated duplications/deletions for every gene family.

Input:
+ testATree
+ famInfoResult.txt

Output:
+ dupDelAll.txt

Sample command: python dupDel.py -t testATree -f famInfoResult.txt -d 3 -c 5 -n 1


### Authors & Contributors

John Snow knows nothing.

### Support or Contact

For additional help, pray to Cthulhu
