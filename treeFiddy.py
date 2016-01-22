import sys, io, ast, math, os, argparse, time, string
from operator import itemgetter
from copy import deepcopy
from collections import Counter
from collections import defaultdict
from itertools import *
from group import Group
from copy import deepcopy
from os import path
from functools import wraps
import numpy as np

def main(argv):

    readEvents(argv[0])


def readEvents(infile):
    '''Read events from ALF log file'''
    lgtList = []
    with open(infile, 'r') as handle:
        for line in islice(handle, 13, None):
            words = line.translate(string.maketrans("",""), string.punctuation)
            info = words.rstrip().split()
            if not info:
                break
            # print info
            if info[0] == 'time':
                if info[2] == 'lgt':
                    lgtList.append(info)

    # print lgtList

    # Dict to sort events
    # key: gene ID, data: [donor organism #, donor organism gene ID, recipient organism #, recipient gene ID]
    geneSort = {}

    for event in lgtList:
        if int(event[-1]) in geneSort.iterkeys():
            print 'We already have this gene?'
        else:
            geneSort[int(event[-1])] = [int(event[5]), int(event[8]), int(event[11]), int(event[-1])]

    # Organisms in clade of interest
    print 'HGT events'
    for key, val in geneSort.iteritems():
        print key, val

    print 'Novel HGT events'
    for key, val in geneSort.iteritems():
        # Change the ranges to reflect which species to include
        # or exclude in looking for HGT events
        if val[2] in range(0,12) and val[0] not in range(0,12):
            print key, val



if __name__ == "__main__":
    main(sys.argv[1:])