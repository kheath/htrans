"""
    Author: Kevin Heath
    Date: 8 June 2015
    Generated for summer research in the Bush lab 

    Group class
    'Holds gene families and their point of origin'

    'I like trains'
"""
import copy import deepcopy

class Group:

    def __init__(self, node, family, familyName, numSpecies):
        self.origin = node
        self.families = [family]
        self.numSpecies = numSpecies
        self.deletions = []
        self.duplications = []

        mergeFamily(family, familyName)
        


    def size(self):
        return len(self.families)

    def origin(self):
        return self.origin

    def families(self):
        return self.families

    def deletitions(self):
        return self.deletions

    def duplications(self):
        return self.duplications

    def addFamily(self, family):
        self.families.append(family)

    def setOrigin(self, node):
        self.origin = node

    def mergeFamily(self, family, familyNum):
        ''' Gene history dictionary:
        Duplications = [[[value,[families]], [value, [families]]], [families], families] where each index
        represents a node '''

        if self.duplications == [] and self.deletions == []:
            for i in range(0,2*self.numSpecies-1):
                self.duplications[].append([])
                self.deletions[].append([])

        dups = family[1][1]
        dels = family[1][2]

        if familyNum in self.families:
            return False

        # Add events to duplications structure
        dupsD = Counter(dups)
        for key, value in dupsD.iteritems():
            isCoEvent = False
            for event in self.duplications[key]:
                if event[0] == value:
                    isCoEvent = True
                    if familyNum not in event[1]:
                        self.duplications[key][1].append(familyNum)
            if not isCoEvent:
                self.duplications[key].append([value, familyNum])

        # Add events to deletions structure
        delsD = Counter(dels)
        for key, value in delsD.iteritems():
            isCoEvent = False
            for event in self.deletions[key]:
                if event[0] == value:
                    isCoEvent = True
                    if familyNum not in event[1]:
                        self.deletions[key][1].append(familyNum)
            if not isCoEvent:
                self.deletions[key].append([value, familyNum])


    def mergeGroups(self, group):
        if self.origin != group.origin():
            return False
        self.families = list(set(self.families)|set(deepcopy(group.families())))
        



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





