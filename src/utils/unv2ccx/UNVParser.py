#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
© Joël Cugnoni, September 2006 - original code, www.caelinux.com
© Ihor Mirzov, August 2019 - refactoring
Distributed under GNU General Public License v3.0

This is a set of objects & functions to read a Universal File
into a simple FEM object structure in order to simplify the
conversion of Mesh definitions from UNV to any other format.

This code is based on two main objects:
1) FEM object structure to store nodes, elements and groups.
2) UNVParser which provides a simple & modular solution to read
some datasets from UNV file and store them in a FEM object structure.

UNV format documentation:
http://sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse/universal-file-datasets-summary
"""

import logging

from utils.unv2ccx import FEM
FLAG = '    -1'


# Universal file parser class
class UNVParser:


    def __init__(self, filename):
        self.file = None
        self.filename = filename
        self.fem = FEM.FEM()
        self.sections = []

        # List of supported datasets and corresponding dataset handler functions
        self.datasetsIds = [2411, 2412, 2467, 2477]
        self.datasetsHandlers = [UNV2411Reader, UNV2412Reader, UNV2467Reader, UNV2467Reader]


    # Read file & fill the section list
    def scanfile(self):
        while True:
            line = self.file.readline()
            if len(line) > 0:
                if line.startswith(FLAG):

                    # Identify section & save offset
                    gid = int(self.file.readline())
                    offset = self.file.tell()
                    self.sections.append([gid, offset])

                    # Ignore data until end of section
                    while not self.file.readline().startswith(FLAG):
                        continue
            else:
                break

        # Rewind file
        self.file.seek(0)


    # Parse UNV file to fill the FEM data structure
    def parse(self):
        self.file = open(self.filename, 'r')
        self.scanfile()
        for sectionId, offset in self.sections:
            if (sectionId in self.datasetsIds):
                self.file.seek(offset)
                func = self.datasetsHandlers[self.datasetsIds.index(sectionId)]
                self.fem = func(self.file, self.fem)
        self.file.close()
        return self.fem


# Reads an UNV2411 dataset (nodes) from file and store data in FEM object
# http://sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse/universal-dataset-number-2411
def UNV2411Reader(f, fem):
    while True:
        line1 = f.readline()
        line2 = f.readline().strip()
        if len(line2) and not line1.startswith(FLAG):
            dataline = Line2Int(line1)
            line2 = line2.replace('D', 'E') # replacement is inserted by Prool
            coords = Line2Float(line2)
            n = FEM.Node(dataline[0], coords)
            fem.nodes.append(n)
        else:
            break
    return fem


# Reads an UNV2412 dataset (elements) from file and store data in FEM object
# http://sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse/universal-dataset-number-2412
def UNV2412Reader(f, fem):
    SpecialElemTypes = [11] # types of elements which are defined on 3 lines
    while True:
        line1 = f.readline()
        line2 = f.readline().strip()
        if len(line2) and not line1.startswith(FLAG):
            dataline = Line2Int(line1)
            etype = dataline[1]
            nnodes = dataline[-1]
            if etype < 33:
                # 1D elements have an additionnal line in definition
                nodes = Line2Int(f.readline())
            else:
                # Standard elements have connectivities on secnd line
                nodes = Line2Int(line2)
                while nnodes > 8:
                    nodes.extend(Line2Int(f.readline()))
                    nnodes -= 8
            e = FEM.Element(dataline[0], etype, nodes)
            fem.elements.append(e)
        else:
            break
    return fem


# Reads an UNV2467 dataset (groups) from file and store data in FEM object
# http://sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse/universal-dataset-number-2467
def UNV2467Reader(f, fem):
    while True:
        # Read Record 1:
        line1 = f.readline()
        # Read Record 2
        line2 = f.readline().strip()
        if len(line2) and not line1.startswith(FLAG):
            # Read group
            dataline = Line2Int(line1)
            group_name = line2.replace(' ', '_')
            nitems = dataline[-1]
            nlines = int((nitems + 1) / 2)

            # check if there are to be any items in the group. 
            if nlines == 0:
                fpos = f.tell()
                # If not, skip row - or if not empty: let row be
                # and start next read of Record 1
                if f.readline().strip():
                    f.seek(fpos)

            else:
                # Read group items
                group_items = []
                for i in range(nlines):
                    dat = Line2Int(f.readline())
                    group_items.append(dat[0:3])
                    if len(dat) > 4:
                        group_items.append(dat[4:7])

                # Split group in node and element sets
                nset = []; eset = []
                for item in group_items:
                    if item[0] == 7:
                        nset.append(item[1])
                    if item[0] == 8:
                        eset.append(item[1])

                # Store non empty groups
                if len(nset):
                    nset = FEM.Group(group_name, 7, nset)
                    fem.nsets.append(nset)
                if len(eset):
                    eset = FEM.Group(group_name, 8, eset)
                    fem.esets.append(eset)
        else:
            break
    return fem


# Convert a string into a list of Float
def Line2Float(line):
    return [float(x) for x in line.split()]


# Convert a string into a list of Int
def Line2Int(line):
    return [int(x) for x in line.split()]
