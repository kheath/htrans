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
        self.mrcag = family[0]
        self.families = [idnum]            # Need to retain order
        self.numSpecies = numSpecies
        self.deletions = []
        self.duplications = []
        self.id = idnum
        self.duplications, self.deletions = self.makeModel(family, idnum)

        self.front = (idnum, self.duplications, self.deletions)
        self.back = (idnum, self.duplications, self.deletions)


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

    def updateFront(self, family, idnum):
        dups, dels = self.makeModel(family)
        self.front = (idnum, dups, dels)
        return True

    def updateBack(self, family, idnum):
        dups, dels = self.makeModel(family)
        self.back = (idnum, dups, dels)
        return True

    def initData(self):

        # Initialze self.duplications and self.deletions if they're empty.  They should have an empty
        # List for every node in the tree
        dups =  []
        dels = []
        for i in range(0,2*self.numSpecies-1):
            dups.append([])
            dels.append([])

        return dups, dels

    def makeModel(self, family, familyNum):
        '''Take a family and generate the dup/del structures for it'''

        dups = family[1][1]
        dels = family[1][2]
        

        duplications, deletions = self.initData()

        # Add events to duplications structure
        dupsD = Counter(dups)
        for key, value in dupsD.iteritems():
            isCoEvent = False
            for event in duplications[key]:
                if len(event) > 0:
                    if event[0] == value:
                        isCoEvent = True
                        if familyNum not in event[1]:
                            event[1].append(familyNum)
                            break
            if not isCoEvent:
                duplications[key].append([value, [familyNum]])

        # Add events to deletions structure
        delsD = Counter(dels)
        for key, value in delsD.iteritems():
            isCoEvent = False
            for event in deletions[key]:
                if len(event) > 0:
                    if event[0] == value:
                        isCoEvent = True
                        if familyNum not in event[1]:
                            event[1].append(familyNum)
            if not isCoEvent:
                deletions[key].append([value, [familyNum]])

        return duplications, deletions


    # def mergeModel(self, m1, m2):

    #     # Merge duplications
    #     newDups = group.getDuplications()
    #     for nodeNum, node in enumerate(newDups):           
    #         for eventNum, event in enumerate(node):
    #             mergeFound = False
    #             for a, b in enumerate(self.duplications[nodeNum]):
    #                 if event[0] == b[a][0]:
    #                     self.duplications[nodeNum][a][1] = list(set(self.duplications[nodeNum][eventNum][1]) 
    #                                                             |set(deepcopy(event[1])))
    #                     mergeFound = True
    #                     break
    #             if not mergeFound:
    #                 self.duplications[nodeNum].append(event)

    #     # Merge deletions
    #     newDels = group.getDeletions()
    #     for nodeNum, node in enumerate(newDels):
    #         for eventNum, event in enumerate(node):
    #             mergeFound = False
    #             for a,b in enumerate(self.deletions[nodeNum]):
    #                 if event[0] == b[a][0]:
    #                     self.deletions[nodeNum][a][1] = list(set(self.deletions[nodeNum][eventNum][1]) 
    #                                                             |set(deepcopy(event[1])))
    #                     mergeFound = True
    #                     break
    #             if not mergeFound:
    #                 self.deletions[nodeNum].append(event)

    #     return True


    def mergeFamily(self, family, familyNum):
        ''' Gene history dictionary:
        Duplications = [[[value,[families]], [value, [families]]], [families], families] where each index
        represents a node.  This way we get node -- number of dups/del -- family number '''


        dups = family[1][1]
        dels = family[1][2]
        # familyNum = family[2]

        # Add events to duplications structure
        dupsD = Counter(dups)
        for key, value in dupsD.iteritems():
            isCoEvent = False
            for event in self.duplications[key]:
                if len(event) > 0:
                    if event[0] == value:
                        isCoEvent = True
                        if familyNum not in event[1]:
                            event[1].append(familyNum)
                            break
            if not isCoEvent:
                self.duplications[key].append([value, [familyNum]])

        # Add events to deletions structure
        delsD = Counter(dels)
        for key, value in delsD.iteritems():
            isCoEvent = False
            for event in self.deletions[key]:
                if len(event) > 0:
                    if event[0] == value:
                        isCoEvent = True
                        if familyNum not in event[1]:
                            event[1].append(familyNum)
            if not isCoEvent:
                self.deletions[key].append([value, [familyNum]])

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
                self.families = group.getFamilies().reverse()+self.families
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
                self.families = self.families+group.getFamilies().reverse()
            else:
                self.families = self.families+group.getFamilies()
            self.back = group.getFront()




        # Merge duplications
        newDups = group.getDuplications()
        for nodeNum, node in enumerate(newDups):           
            for eventNum, event in enumerate(node):
                mergeFound = False
                for a, b in enumerate(self.duplications[nodeNum]):
                    if event[0] == b[a][0]:
                        self.duplications[nodeNum][a][1] = list(set(self.duplications[nodeNum][eventNum][1]) 
                                                                |set(deepcopy(event[1])))
                        mergeFound = True
                        break
                if not mergeFound:
                    self.duplications[nodeNum].append(event)

        # Merge deletions
        newDels = group.getDeletions()
        for nodeNum, node in enumerate(newDels):
            for eventNum, event in enumerate(node):
                mergeFound = False
                for a,b in enumerate(self.deletions[nodeNum]):
                    if event[0] == b[a][0]:
                        self.deletions[nodeNum][a][1] = list(set(self.deletions[nodeNum][eventNum][1]) 
                                                                |set(deepcopy(event[1])))
                        mergeFound = True
                        break
                if not mergeFound:
                    self.deletions[nodeNum].append(event)

        return True




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





