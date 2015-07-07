"""
    Author: Kevin Heath
    Date: 8 June 2015
    Generated for summer research in the Bush lab 

    Group class
    'Holds gene families and their point of origin'
    
    'I like trains'
"""
from copy import deepcopy
from collections import Counter

class Group:

    def __init__(self, family):
        # In this simplified version, a family is simply a number
        # self.mrcag = family[0]
        self.families = [family]            # List of family number that retain order
        # self.numSpecies = numSpecies
        self.id = family
        # self.deletions = []
        # self.duplications = []

        # initData()

        # self.mergeFamily(family)
        


    def getSize(self):
        return len(self.families)

    # def getMrcag(self):
        # return self.mrcag

    def getFamilies(self):
        return self.families

    # def getDeletions(self):
        # return self.deletions

    # def getDuplications(self):
        # return self.duplications

    def addFamily(self, family):
        self.families.append(family)
        # self.mergeFamily(family)

    # def setMrcag(self, node):
        # self.mrcag = node

    # def initData(self):

        # Initialze self.duplications and self.deletions if they're empty.  They should have an empty
        # List for every node in the tree

        # if self.duplications == [] and self.deletions == []:
            # for i in range(0,2*self.numSpecies-1):
                # self.duplications.append([])
                # self.deletions.append([])

    # def mergeFamily(self, family):
    #     ''' Gene history dictionary:
    #     Duplications = [[[value,[families]], [value, [families]]], [families], families] where each index
    #     represents a node.  This way we get node -- number of dups/del -- family number '''


    #     dups = family[1][1]
    #     dels = family[1][2]
    #     familyNum = family[2]

    #     # Add events to duplications structure
    #     dupsD = Counter(dups)
    #     for key, value in dupsD.iteritems():
    #         isCoEvent = False
    #         for event in self.duplications[key]:
    #             if len(event) > 0:
    #                 if event[0] == value:
    #                     isCoEvent = True
    #                     if familyNum not in event[1]:
    #                         event[1].append(familyNum)
    #                         break
    #         if not isCoEvent:
    #             self.duplications[key].append([value, [familyNum]])

    #     # Add events to deletions structure
    #     delsD = Counter(dels)
    #     for key, value in delsD.iteritems():
    #         isCoEvent = False
    #         for event in self.deletions[key]:
    #             if len(event) > 0:
    #                 if event[0] == value:
    #                     isCoEvent = True
    #                     if familyNum not in event[1]:
    #                         event[1].append(familyNum)
    #         if not isCoEvent:
    #             self.deletions[key].append([value, [familyNum]])


    def mergeGroups(self, famOrder):
        '''Merges two groups with 1 or more families each together.
        This should currently fail if they have the same origin.'''

        self.families = famOrder

        # Check if both groups have the same origin
        # if self.mrcag != group.getMrcag():
        #     return False

        # Takes in given family order - need to figure this out outside the class
        # self.families = famOrder

        # # Merge duplications
        # newDups = group.getDuplications()
        # for nodeNum, node in enumerate(newDups):           
        #     for eventNum, event in enumerate(node):
        #         mergeFound = False
        #         for a, b in enumerate(self.duplications[nodeNum]):
        #             if event[0] == b[a][0]:
        #                 self.duplications[nodeNum][a][1] = list(set(self.duplications[nodeNum][eventNum][1]) 
        #                                                         |set(deepcopy(event[1])))
        #                 mergeFound = True
        #                 break
        #         if not mergeFound:
        #             self.duplications[nodeNum].append(event)

        # # Merge deletions
        # newDels = group.getDeletions()
        # for nodeNum, node in enumerate(newDels):
        #     for eventNum, event in enumerate(node):
        #         mergeFound = False
        #         for a,b in enumerate(self.deletions[nodeNum]):
        #             if event[0] == b[a][0]:
        #                 self.deletions[nodeNum][a][1] = list(set(self.deletions[nodeNum][eventNum][1]) 
        #                                                         |set(deepcopy(event[1])))
        #                 mergeFound = True
        #                 break
        #         if not mergeFound:
        #             self.deletions[nodeNum].append(event)

        # return True






        # Add duplications
        # it = iter(family[1][1])
        # while True:
        #     try:
        #         next = it.next()
        #     except StopIteration:
        #         break
        #     self.history[it][0] += 1
        #     self.history[it][1].add(familyNum)

        # Add deletions
        # it = iter(family[1][2])
        # while True:
        #     try:
        #         next = it.next()
        #     except StopIteration:
        #         break
        #     self.history[it][2] += 1
        #     self.history[it][3].add(familyNum)