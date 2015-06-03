"""
  Author: Kevin Heath
  Date: 29 May 2015
  Generated for summer research in the Bush lab """


import io, sys, getopt
from itertools import islice


def main(argv):
    try:
      opts, args = getopt.getopt(argv,"hs:n:i:o:p:v")
    except getopt.GetoptError:
      print 'test.py -s <orthologsfile> -n <# of species> -i <inputfile> -o <outputfile> -p <partition output file name>'
      sys.exit(2)
    
    global verbose
    verbose = False

    # Set the variables
    for opt, arg in opts:
      if opt == '-h':
          print 'test.py -s <orthologsfile> -n <# of species> -i <inputfile> -o <outputfile> -p <partition output file name>'
      elif opt == '-o':
          output_filename = arg
      elif opt == '-i':
          input_filename = arg
      elif opt == '-s':
          ortholog_filename = arg
      elif opt == '-n':
          numSpecies = int(arg)
      elif opt == '-p':
          partoutput_filename = arg
      elif opt == '-v':
          verbose = True

    # Generate a list of species from the ortholog file
    species = getSpecies(ortholog_filename, numSpecies)

    # if verbose:
    #   print 'Species found:'
    #   for sp in species:
    #     print sp

    # Create a dictionary with each species and the corresponding genes
    souperAlignment = superAlign(species, input_filename)

    # Output the super alignment with all the genes concatenated together as a .fasta file
    writeSuperAlignment(souperAlignment, species, output_filename)

    # Output each gene as a partition.
    makePartition(souperAlignment.itervalues().next(), partoutput_filename)

    if verbose:
      print "\nWhen cryptography is outlawed, bayl bhgynjf jvyy unir cevinpl."

def getSpecies(filename, numSpecies):
    """Extracts the species names from the orthologs.txt file"""
    if verbose:
      print 'Extracting species...'


    with open(filename, 'r') as f:
      species = []

      while numSpecies > 0:
          line = f.readline()
          temp = line.split()
          spName = temp[1].replace('_', '')
          spName = spName.replace('.', '')


          if len(spName) > 10:
            shortname = spName[-10:]
          elif len(spName) == 10:
            shortname = spName
          else:
            shortname = '{:0<10}'.format(spName)
          print spName+' > '+shortname
          # name = temp[1][-9:].replace('_', '')
          species.append(shortname)
          numSpecies-=1

    return species


def superAlign(species, alignments):
    """Concatenate gene blocks together

    Assumes genes are stored in blocks,
    Each species is listed in the same order"""

    if verbose:
      print 'Generating super alignment...'
    

    numSpecies = len(species)
    
    collection = {sp:[] for sp in species}
    numBlocks = 0
    with open(alignments, 'r') as f:
      while True:
        next_n_lines = list(islice(f, 1, numSpecies*2+1, 2))
        
        if not next_n_lines:
          break
        numBlocks+=1
        for index,key in enumerate(species):
          collection[key].append(next_n_lines[index].strip())

    if verbose:
      print 'Processed '+str(numBlocks)+' blocks'

    return collection



def writeSuperAlignment(superalignment, species, outputfile):
    """Write out super alignment to file"""

    if verbose:
      print 'Writing super alignment, prepare for awesomeness...'

    with open(outputfile, 'w+') as f:
      for key in species:
        f.writelines(['>'+key+'\n', ''.join(superalignment.get(key))+'\n'])
      f.write('\n')

def makePartition(geneList, output):
    """Write out the partition.

    Takes a single entry in the superalignment dictionary
    and uses the length of each gene to make the partition/

    Assumes the sequence starts at 1, based on examples"""

    if verbose:
      print 'Writing partition file...'

    start,end = 1,1
    count = 1
    with open(output, 'w+') as f:
      it = iter(geneList)
      gene=''
      while True:
        try:
          gene = it.next()
            
        except StopIteration:
          if verbose:
            print "Exausted list of partitions"
          break
        end=start+len(gene)-1
        temp = "DNA, gene{0}={1}-{2}\n".format(count,start,end)
        f.write(temp)
        start=end+1
        count+=1
        



if __name__ == "__main__":
   main(sys.argv[1:])



    

