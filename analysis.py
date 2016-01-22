'''
File to automate some simple analysis of htrans results.
    Groups per Node
    Groups per Group Size
    Genes in GAD AFI
    Size of distance matrix
        Average/mode/median value
'''

from htrans import *
mpl.use('Agg')


def main(argv):

    plt.ioff()  # Turn off interactive mode for plots

    ############## ---- Parse Command line Arguments ---- ################
    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-d', '--distances', help='Pickle file with the distance matrix (initialDistances.pickle)', required=True)
    requiredNamed.add_argument('-g', '--groups', help='File with remaining groups (groupsL.txt)', required=True)
    args = parser.parse_args()


    groupsL = readGroups(args.groups)
    distances = readPickle(args.distances)

    groupsBySize, groupsByNode = analyzeGroups(groupsL)

    famMap, groupMap = famGroupMap(groupsL)

    print 'Groups by Size:'
    sumValues = 0
    for key, value in groupsBySize.iteritems():
        print str(key)+': '+str(len(value))
        sumValues += len(value)

    print str(sumValues)+' total groups.'
    print str(sumValues - len(groupsBySize[1]))+' groups of 2 or more families.\n'

    print 'Groups by Node:'
    for key, value in groupsByNode.iteritems():
        print str(key)+': '+str(len(value))


    print 'Statistics of distances:'
    dValues = distances.values()
    dMean = np.mean(dValues)
    dStdDev = np.std(dValues)
    dVar = np.var(dValues)
    dMin = np.amin(dValues)
    dMax = np.amax(dValues)

    print 'Min: '+str(dMin)
    print 'Max: '+str(dMax)
    print 'Mean: '+str(dMean)
    print 'Std. Deviation: '+str(dStdDev)
    print 'Variance: '+str(dVar)

    fig, ax = plt.subplots()
    counts, bins, patches = ax.hist(dValues, facecolor='blue')
    #ax.set_xticks(bins)
    plt.subplots_adjust(bottom=0.15)
    plt.savefig('histogram.png', bbox_inches='tight')
    plt.clf()

    counts = [0,0,0,0,0,0]
    for value in distances.itervalues():
        if value == -10:
            counts[0]+=1
        elif value == -2:
            counts[1]+=1
        elif value == -1:
            counts[2]+= 1
        elif value == 0:
            counts[3]+=1
        elif value == 1:
            counts[4]+=1
        elif value == 2:
            counts[5]+=1
        else:
            print 'Unexpected value:', value

    print 'Done!'


def analyzeGroups(groups):
    '''Count the number of groups of each size, and number of groups at each node'''

    sortedBySize = defaultdict(list)
    sortedByNode = defaultdict(list)

    for group in groups:
        if group != None:
            sortedBySize[len(group.getFamilies())].append(group)
            sortedByNode[group.getMrcag()].append(group)

    return sortedBySize, sortedByNode

def famGroupMap(groups):
    '''Family -> Group Dictionary and Group -> Family Dictionary for ease of use'''

    famMap = {}
    groupMap = {}

    for group in groups:
        if group != None:
            groupMap[group.getIdNum()] = group.getFamilies()
            for family in group.getFamilies():
                famMap[family] = group.getIdNum()

    return famMap, groupMap

def findOtherFams(family, famMap, groupMap):
    '''Finds which other families are in the same group as the specified family'''

    return groupMap[famMap[family]]


def makeHistogram(data):

    fig, ax = plt.subplots()
    counts, bins, patches = ax.hist(data, facecolor='blue')
    ax.set_xticks(bins)
    plt.subplots_adjust(bottom=0.15)
    plt.savefig('histogram.png', bbox_inches='tight')
    plt.clf()



if __name__ == "__main__":
    main(sys.argv)