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



def main(argv):

    familyData = readFamilies(argv[0])

    data, row_labels, column_labels = evaluate(familyData)
    
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
    heatmap = ax.pcolor(data, cmap=plt.cm.Blues)

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)

    cbar = plt.colorbar(heatmap)
    
    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticklabels(row_labels, minor=False)
    ax.set_yticklabels(column_labels, minor=False)
    
    plt.savefig('heatmap.png', bbox_inches='tight')
    plt.close(fig)

def rand_jitter(arr):
    stdev = .01*(max(arr)-min(arr))
    return arr + np.random.randn(len(arr)) * stdev


def readFamilies(filename):

    families = []

    with open(filename, 'r') as f:
        while True:
            line = f.readline()

            if not line:
                break
            temp = ast.literal_eval(line)
            if temp[1][1] != [] or temp[1][2] != []:
                families.append(temp)

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

def evaluate(familyData):
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
        
        
    return heat, row_labels, column_labels

def calcDiff(l1, l2, n):
    # Calculate the difference in the number of times each node appears in two list
    diff = 0
    for i in range(0,2*n-1):
        diff += abs(l1.count(i)-l2.count(i))
    return diff



if __name__ == "__main__":
   main(sys.argv[1:])
