#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
© Lukas Bante, Sep 2017 - original code https://gitlab.lrz.de/snippets/238
© Ihor Mirzov, Jan 2020 - bugfix, refactoring and improvement
Distributed under GNU General Public License v3.0

This module contains classes for parsing CalculiX .frd files """

import os
import re
import logging
import math

import numpy as np


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


# Nodal Point Coordinate Block
# cgx_2.15 Manual, § 11.3
class NodalPointCoordinateBlock:


    # Read nodal coordinates
    def __init__(self, in_file):
        line = readByteLine(in_file)
        self.nodes = {} # dictionary with nodes {num:Node}
        while True:
            line = readByteLine(in_file)

            # End of block
            if line == '-3':
                break

            regex = '^-1(.{10})' + '(.{12})'*3
            match = parseLine(regex, line)
            node_number = int(match.group(1))
            node_coords = [ float(match.group(2)),
                            float(match.group(3)),
                            float(match.group(4)), ]
            self.nodes[node_number] = Node(node_number, node_coords)
            # logging.debug('Node {}: {}'.format(node_number, node_coords))

        self.numnod = len(self.nodes) # number of nodes in this block
        logging.info('{} nodes'.format(self.numnod)) # total number of nodes


# Element Definition Block
# cgx_2.15 Manual, § 11.4
class ElementDefinitionBlock:

    # Parse elements
    def __init__(self, in_file):
        self.in_file = in_file
        line = readByteLine(in_file)
        self.elements = [] # list of Element objects

        while True:
            line = readByteLine(in_file)

            # End of block
            if line == '-3':
                break

            self.parseElement(line)

        self.numelem = len(self.elements) # number of elements in this block
        logging.info('{} cells'.format(self.numelem)) # total number of elements


    # Read element composition
    def parseElement(self, line):
        """
            -1         1    1    0AIR
            -2         1         2         3         4         5         6         7         8
            -1         1   10    0    1
            -2         1         2         3         4         5         6         7         8
            -1         2   11    0    2
            -2         9        10
            -1         3   12    0    2
            -2        10        12        11
         """
        element_num = int(line.split()[1])
        element_type = int(line.split()[2])
        element_nodes = []
        num_nodes = self.amount_of_nodes_in_frd_element(element_type)
        for j in range(self.num_lines(element_type)):
            line = readByteLine(self.in_file)
            nodes = [int(n) for n in line.split()[1:]]
            element_nodes.extend(nodes)

        elem = Element(element_num, element_type, element_nodes)
        self.elements.append(elem)
        # logging.debug('Element {}: {}'.format(element_num, element_nodes))


    # Amount of nodes in frd element
    def amount_of_nodes_in_frd_element(self, etype):
        # First value is meaningless, since elements are 1-based
        return (0, 8, 6, 4, 20, 15, 10, 3, 6, 4, 8, 2, 3)[etype]


    # Amount of lines in element connectivity definition
    def num_lines(self, etype):
        # First value is meaningless, since elements are 1-based
        return (0, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1)[etype]


# Nodal Results Block
# cgx_2.15 Manual, § 11.6
class NodalResultsBlock:


    # Read calculated values
    def __init__(self, in_file, node_block):
        self.in_file = in_file
        self.node_block = node_block
        self.components = [] # component names
        self.results = {} # dictionary with nodal result {node:data}

        self.readStepInfo()
        self.readVarsInfo()
        self.readComponentsInfo()
        results_counter = self.readNodalResults()
        self.appendStresses() # append Mises and principal stresses
        self.appendStrains() # append principal strains

        if self.value < 1:
            time_str = 'time {:.2e}, '.format(self.value)
        else:
            time_str = 'time {:.1f}, '.format(self.value)
        logging.info('Step {}, '.format(self.numstep) +\
                    time_str +\
                    '{}, '.format(self.name) +\
                    '{} components, '.format(len(self.components)) +\
                    '{} values'.format(results_counter))


    # Read step information
    def readStepInfo(self):
        """
            CL  101 0.36028E+01         320                     3    1           1
            CL  101 1.000000000         803                     0    1           1
            CL  101 1.000000000          32                     0    1           1
            CL  102 117547.9305          90                     2    2MODAL      1
        """
        line = readByteLine(self.in_file)[7:]
        regex = '^(.{12})\s+\d+\s+\d+\s+(\d+)'
        match = parseLine(regex, line)
        self.value = float(match.group(1)) # could be frequency, time or any numerical value
        self.numstep = int(match.group(2)) # step number


    # Read variables information
    def readVarsInfo(self):
        """
            -4  V3DF        4    1
            -4  DISP        4    1
            -4  STRESS      6    1
            -4  DOR1  Rx    4    1
        """
        line = readByteLine(self.in_file)[4:]
        regex = '^(\w+)' + '\D+(\d+)'*2
        match = parseLine(regex, line)
        self.ncomps = int(match.group(2)) # amount of components

        # Rename result block to the name from .inp-file
        inpname = {
            'DISP':'U',
            'NDTEMP':'NT',
            'STRESS':'S',
            'TOSTRAIN':'E',
            'FORC':'RF',
            'PE':'PEEQ',
            }
        self.name = match.group(1) # dataset name
        if self.name in inpname:
            self.name = inpname[self.name]


    # Iterate over components
    def readComponentsInfo(self):
        """
            -5  D1          1    2    1    0
            -5  D2          1    2    2    0
            -5  D3          1    2    3    0
            -5  ALL         1    2    0    0    1ALL

            -5  DFDN        1    1    1    0
            -5  DFDNFIL     1    1    2    0

            -5  V1          1    2    1    0
            -5  V2          1    2    2    0
            -5  V3          1    2    3    0
            -5  ALL         1    2    0    0    1ALL

            -5  SXX         1    4    1    1
            -5  SYY         1    4    2    2
            -5  SZZ         1    4    3    3
            -5  SXY         1    4    1    2
            -5  SYZ         1    4    2    3
            -5  SZX         1    4    3    1
        """
        for i in range(self.ncomps):
            line = readByteLine(self.in_file)[4:]
            regex = '^\w+'
            match = parseLine(regex, line)

            # Exclude variable name from the component name: SXX->xx, EYZ->yz
            component_name = match.group(0)
            if component_name.startswith(self.name):
                component_name = component_name[len(self.name):].lower()

            if 'ALL' in component_name:
                self.ncomps -= 1
            else:
                self.components.append(component_name)


    # Iterate over nodal results
    def readNodalResults(self):
        """
            -1         1-7.97316E+10-3.75220E-01
            -1         2-8.19094E+10-3.85469E-01

            -1         1-6.93889E-18-9.95185E-01-4.66908E-34
            -1         2-1.94151E-01-9.76063E-01 6.46011E-30

            -1         1 1.47281E+04 1.39140E+04 2.80480E+04 5.35318E+04 6.36642E+03 1.82617E+03
            -2           5.31719E+01 6.69780E+01 2.76244E+01 2.47686E+01 1.99930E+02 2.14517E+02
        """
        # Fill data with zeroes - sometimes FRD result block has only non zero values
        for node_num in self.node_block.nodes.keys():
            self.results[node_num] = [0]*self.ncomps

        # Some warnings repeat too much time - mark them
        before = ''
        after = None
        emitted_warning_types = {'NaNInf':0, 'WrongFormat':0}

        results_counter = 0 # independent results counter
        while True:
            line = readByteLine(self.in_file)

            # End of block
            if line == '-3':
                break

            row_comps = min(6, self.ncomps) # amount of values written in row
            regex = '^-1\s+(\d+)' + '(.{12})' * row_comps
            match = parseLine(regex, line)
            node = int(match.group(1))
            data = []
            for c in range(row_comps):
                m = match.group(c + 2)
                try:
                    # NaN/Inf values will be parsed 
                    num = float(m)
                    if ('NaN' in m or 'Inf' in m):
                        emitted_warning_types['NaNInf'] += 1
                except Exception as e:
                    # Too big number is written without 'E'
                    num = float(re.sub(r'(.+).([+-])(\d{3})', r'\1e\2\3', m))
                    emitted_warning_types['WrongFormat'] += 1
                    before = m
                    after = num
                data.append(num)

            results_counter += 1
            self.results[node] = data

            # Result could be multiline
            for j in range((self.ncomps-1)//6):
                row_comps = min(6, self.ncomps-6*(j+1)) # amount of values written in row
                line = readByteLine(self.in_file)
                regex = '^-2\s+' + '(.{12})' * row_comps
                match = parseLine(regex, line)
                data = [float(match.group(c+1)) for c in range(row_comps)]
                self.results[node].extend(data)

            # logging.debug('Node {}: {}'.format(node, self.results[node]))

        if emitted_warning_types['NaNInf']:
            logging.warning('NaN and Inf are not supported in Paraview ({} warnings).'\
                .format(emitted_warning_types['NaNInf']))
        if emitted_warning_types['WrongFormat']:
            logging.warning('Wrong format, {} -> {} ({} warnings).'\
                .format(before.strip(), after, emitted_warning_types['WrongFormat']))
        return results_counter


    # Append Mises and principal stresses
    def appendStresses(self):
        if self.name == 'S':
            try:
                # component_names = (
                #     'Mises',
                #     'Max Principal',
                #     'Mid Principal',
                #     'Min Principal',
                #     'Tresca',
                #     'Pressure',
                #     'Third Invariant'
                # )
                component_names = (
                    'Mises',
                    'Min Principal',
                    'Mid Principal',
                    'Max Principal',
                    )
                for i in range(len(component_names)):
                    # c = Component()
                    # c.ictype = 1; c.name = component_names[i]
                    self.components.append(component_names[i])
                    self.ncomps += 1

                # Iterate over nodes
                for node_num in self.node_block.nodes.keys():
                    data = self.results[node_num] # list with results for current node
                    Sxx = data[0]; Syy = data[1]; Szz = data[2]
                    Sxy = data[3]; Syz = data[4]; Szx = data[5]
                    tensor = np.array([[Sxx, Sxy, Szx], [Sxy, Syy, Syz], [Szx, Syz, Szz]])

                    # Calculate Mises stress for current node
                    mises = 1 / math.sqrt(2) *\
                        math.sqrt(  (Sxx - Syy)**2 +\
                                    (Syy - Szz)**2 +\
                                    (Szz - Sxx)**2 +\
                                    6 * Syz**2 +\
                                    6 * Szx**2 +\
                                    6 * Sxy**2)
                    self.results[node_num].append(mises)

                    # Calculate principal stresses for current node
                    for ps in np.linalg.eigvalsh(tensor).tolist():
                        self.results[node_num].append(ps)

            except:
                logging.error('Additional stresses will not be appended.')


    # Append principal strains
    def appendStrains(self):
        if self.name == 'E':
            try:
                component_names = (
                    'Mises',
                    'Min Principal',
                    'Mid Principal',
                    'Max Principal',
                    )
                for i in range(len(component_names)):
                    self.components.append(component_names[i])
                    self.ncomps += 1

                # Iterate over nodes
                for node_num in self.node_block.nodes.keys():
                    data = self.results[node_num] # list with results for current node
                    Exx = data[0]; Eyy = data[1]; Ezz = data[2]
                    Exy = data[3]; Eyz = data[4]; Ezx = data[5]
                    tensor = np.array([[Exx, Exy, Ezx], [Exy, Eyy, Eyz], [Ezx, Eyz, Ezz]])

                    # Calculate Mises strain for current node
                    mises = math.sqrt(2)/3 *\
                        math.sqrt(   (Exx - Eyy)**2 +\
                                (Eyy - Ezz)**2 +\
                                (Ezz - Exx)**2 +\
                                6 * Eyz**2 +\
                                6 * Ezx**2 +\
                                6 * Exy**2)
                    self.results[node_num].append(mises)

                    # Calculate principal strains for current node
                    for ps in np.linalg.eigvalsh(tensor).tolist():
                        self.results[node_num].append(ps)

            except:
                logging.error('Additional strains will not be appended.')


# Main class
class Parse:


    # Read contents of the .frd file
    def __init__(self, file_name=None):
        self.file_name = None   # path to the .frd-file to be read
        self.node_block = None  # node block
        self.elem_block = None  # elements block
        self.result_blocks = [] # all result blocks in order of appearance
        if file_name:
            self.file_name = file_name
            with open(file_name, 'rb') as in_file:
                key = in_file.read(5).decode().strip()
                while key:
                    # Header
                    if key == '1' or key == '1P':
                        readByteLine(in_file)

                    # Nodes
                    elif key == '2':
                        block = NodalPointCoordinateBlock(in_file)
                        self.node_block = block

                    # Elements
                    elif key == '3':
                        block = ElementDefinitionBlock(in_file)
                        self.elem_block = block

                    # Results
                    elif key == '100':
                        block = NodalResultsBlock(in_file, self.node_block)
                        self.result_blocks.append(block)

                    # End
                    elif key == '9999':
                        break

                    key = in_file.read(5).decode().strip()

            # # Exclude zero nodes added by ccx due to *TRANSFORM
            # nn = sorted(set([len(b.results) for b in self.result_blocks if len(b.results)>0]))
            # if len(nn) == 3:
            #     self.node_block.numnod = nn[1]
            # elif len(nn) == 2:
            #     self.node_block.numnod = nn[0]


# Read byte line and decode
def readByteLine(f):
    try:
        byte = f.read(1).decode()
    except UnicodeDecodeError:
        byte = ' '
    line = byte
    while byte != '\n':
        try:
            byte = f.read(1).decode()
        except UnicodeDecodeError:
            byte = ' '
        line += byte
    return line.strip()
# TODO Use this readByteLine:
# Read byte line and decode: return None after EOF
def read_byte_line(f):

    # Check EOF
    byte = f.read(1)
    if not byte:
        return None

    # Convert first byte
    try:
        line = byte.decode()
    except UnicodeDecodeError:
        line = ' ' # replace endecoded symbols with space

    # Continue reading until EOF or new line
    while byte != b'\n':
        byte = f.read(1)
        if not byte:
            return line.strip() # EOF
        try:
            line += byte.decode()
        except UnicodeDecodeError:
            line += ' ' # replace endecoded symbols with space

    return line.strip()


# Parse regex in line and report problems
def parseLine(regex, line):
    match = re.search(regex, line)
    if match:
        return match
    else:
        logging.error('Can\'t parse line:\n{}\nwith regex:\n{}'\
                .format(line, regex))
        raise Exception
