"""
    Author: Kevin Heath
    Date: 9 June 2015
    Generated for summer research in the Bush lab 

    Visualize family data

    'Where is Waldo?'
"""

import sys, matplotlib.pyplot as plt
import io, ast, math
import numpy as np
from collections import Counter
from itertools import *
from group import Group
from copy import deepcopy


def main(argv):

    familyData = readFamilies(argv[0])
    # print familyData

    uData = uniq(familyData)

    print 'Original length = ', len(familyData)
    print 'Uniqued length = ', len(uData)

    data, row_labels, column_labels = evaluate(uData, 1.0, 2.0)
    
    # gs = initializeGroups(familyData)
    # for key, value in gs.iteritems():
        # print key, ' : ', value.getFamilies(), '--', value.getDuplications(), '--', value.getDeletions()


    idstuff = idGroups(familyData, data, 1.0)

    for key, value in idstuff.iteritems():
        if len(value.getFamilies()) > 1:
            print key, ' : ', value.getFamilies(), '--', value.getDuplications(), '--', value.getDeletions()
            print ''

    

    # x,y = makeXY(familyData)

    # plt.plot(x, y, 'ro')
    # plt.axis([4, max(x)+1, 0, max(y)])
    # plt.show()

    # print(max(y))

    # plt.hist(y)
    # plt.show()

    # plt.scatter(rand_jitter(x),rand_jitter(y))
    # plt.show()

    fig, ax = plt.subplots()
    heatmap = ax.pcolor(data, cmap=plt.cm.Blues)#, edgecolors='k')

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)

    cbar = plt.colorbar(heatmap)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticklabels([''], minor=False, fontsize=6)
    ax.set_yticklabels([''], minor=False, fontsize=6)
    
    plt.savefig('heatmap.png', bbox_inches='tight')
    plt.close(fig)

def rand_jitter(arr):
    stdev = .01*(max(arr)-min(arr))
    return arr + np.random.randn(len(arr)) * stdev


def readFamilies(filename):

    families = []

    with open(filename, 'r') as f:
        famNum = 1
        while True:
            line = f.readline()

            if not line:
                break
            temp = list(ast.literal_eval(line))
            if temp[1][1] != [] or temp[1][2] != []:
                temp.append(famNum)
                families.append(temp)
            famNum+=1

    return families

# def makeXY(familyData):

#     x = []
#     y = []

#     for family in familyData:
#         heuristic = len(family[1][1])+len(family[1][2])
#         if heuristic> 0 and heuristic < 80:
#             x.append(family[0])
#             y.append(len(family[1][1])+len(family[1][2]))

#     return (x,y)

# def something(heatmap, familyData):
#     for x in range(0, heatmap.shape):
#         for y in range(x+1, heatmap.shape):

def initializeGroups(familyData):
    '''Make every family into a group'''

    groups = {}

    for index, family in enumerate(familyData):
        groups[index] = Group(family, 6)

    return groups





def idGroups(familyData, heatmap, cutoff):
    '''Merge families into groups based on a cutoff'''

    groups = []
    bigG = {}
    newMap = deepcopy(heatmap)

    for x in range(0, len(familyData)):
        for y in range(x+1, len(familyData)):
            if newMap[x,y] >= cutoff:
                groups.append((x,y))

    print 'List of groups: '
    print groups

    count = 0
    while len(groups) > 0:

        commonSp = [item for item in groups if count in item]
        if len(commonSp) > 0:
            flat = uniq(list(sum(commonSp,())))
            print 'Merging families: ', flat
            bigG[count] = Group(familyData[flat[0]], 6)
            
            if len(flat[1:]) > 0:
                for i in range(1, len(flat)):
                    bigG[count].addFamily(familyData[flat[i]])
            temp = len(groups)
            for num in flat:
                for tup in groups:
                    if num in tup:
                        groups.remove(tup)
                # groups.remove(tup)
            print 'Removed ', temp-len(groups), ' tuples' 
        count += 1

    return bigG



def evaluate(familyData, cut, maxVal):
    '''Creates a heatmap of familyData vs familyData'''

    column_labels = [str(family[0]) for family in familyData]
    row_labels = column_labels


    heat = np.zeros((len(column_labels), len(row_labels)))

    for (x,y), value in np.ndenumerate(heat):
        # print familyData[x][0], familyData[y][0]
        # if familyData[x][0] == familyData[y][0]:
        #     heat[x,y] = 1.0
        # elif familyData[x]
        # diff = len(set(familyData[x][1][1]).symmetric_difference(familyData[y][1][1]))+len(set(familyData[x][1][2]).symmetric_difference(familyData[y][1][2]))
        diff = calcDiff(familyData[x][1][1], familyData[y][1][1], 6) + calcDiff(familyData[x][1][2], familyData[y][1][2], 6)
        # if familyData[x][1][1] != [] and familyData[y][1][1] != [] and familyData[x][1][2] != [] and familyData[y][1][2] != []:
        heat[x,y] = 1.0/(diff+1.0)
        if familyData[x][0] == familyData[y][0]:
            heat[x,y] += 1.0

        if heat[x,y] >= cut:
            heat[x,y] = maxVal
        
        
    return heat, row_labels, column_labels

def calcDiff(l1, l2, n):
    # Calculate the difference in the number of times each node appears in two list
    diff = 0
    for i in range(0,2*n-1):
        diff += abs(l1.count(i)-l2.count(i))
    return diff

    # Sum of unique events / total number of events

def uniq(input):
  output = []
  for x in input:
    if x not in output:
      output.append(x)
  return output

def unique(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if not (x in seen or seen_add(x))]


if __name__ == "__main__":
   main(sys.argv[1:])
