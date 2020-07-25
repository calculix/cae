#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Joël Cugnoni, September 2006 """


# A simple FEM object structure
class FEM:
    def __init__(self):
        self.nodes = []
        self.elements = []
        self.nsets = []
        self.esets = []


# A single node object
class Node:
    def __init__(self, num, coords):
        self.num = num
        self.coords = coords


# A single finite element object
class Element:
    def __init__(self, num, etype, nodes):
        self.num = num
        self.type = etype
        self.nodes = nodes


# FE group object
class Group:
    def __init__(self, name, gtype, items):
        self.name = name
        self.type = gtype
        self.items = items
        self.nitems = len(items)

