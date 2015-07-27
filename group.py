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

    def __init__(self, family, idnum, numSpecies):
        self.id = idnum
        self.mrcag = family[0]
        self.families = [idnum]            # Need to retain order
        self.numSpecies = numSpecies
        

        self.front = idnum
        self.back = idnum

    def __str__(self):
        return [self.id, self.mrcag, self.families, self.numSpecies, self.front, self.back]

    def readIn(self, params):
        self.id = params[0]
        self.mrcag = params[1]
        self.families = params[2]
        self.numSpecies = params[3]
        self.front = params[4]
        self.back = params[5]

    def printG(self):
        return str(self.id)+str('\t')+str(self.mrcag)+str('\t')+str(self.families)+str('\t')+str(self.numSpecies)+str('\t')+str(self.front)+str('\t')+str(self.back)

    def getSize(self):
        return len(self.families)

    def getMrcag(self):
        return self.mrcag

    def getFamilies(self):
        return self.families

    def getNumSpecies(self):
        return self.numSpecies

    def getIdNum(self):
        return self.id

    def getDeletions(self):
        return self.deletions

    def getDuplications(self):
        return self.duplications

    def getFront(self):
        return self.front

    def getBack(self):
        return self.back

    def setMrcag(self, node):
        self.mrcag = node
        return True

    def updateFront(self, famNum):
        self.families.insert(0, famNum)
        self.front = famNum
        return True

    def updateBack(self, famNum):
        self.families.append(famNum)
        self.back = famNum
        return True


    def mergeGroup(self, group, direction):
        '''Merges two groups with 1 or more families each together.
        This should currently fail if they have the same origin.'''

        # Check if both groups have the same origin
        # if self.mrcag != group.getMrcag():
            # return False

        # Combine family lists in correct order
        if direction == 0:
            if len(group.getFamilies()) > 1:
                self.families = list(reversed(group.getFamilies()))+self.families
            else:
                self.families = group.getFamilies()+self.families
            self.front = group.getBack()
        elif direction == 1:
            self.families = self.families+group.getFamilies()
            self.back = group.getBack()
        elif direction == 2:
            self.families = group.getFamilies()+self.families
            self.front = group.getFront()
        elif directions == 3:
            if len(group.getFamilies()) > 1:
                self.families = self.families+list(reversed(group.getFamilies()))
            else:
                self.families = self.families+group.getFamilies()
            self.back = group.getFront()



        return True





