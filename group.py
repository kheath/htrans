"""
    Author: Kevin Heath
    Date: 8 June 2015
    Generated for summer research in the Bush lab 

    Group class
    'Holds gene families and their point of origin'

    'I like trains'
"""

class Group:

    def __init__(self, node, families, deletions, duplications):
        self.origin = node
        self.families = families
        self.deletions = deletions
        self.duplications = duplications

    def size(self):
        return len(self.families)

    def origin(self):
        return self.origin

    def families(self):
        return self.families

    def addFamily(self, family):
        self.families.append(family)

    def setOrigin(self, node):
        self.origin = node



