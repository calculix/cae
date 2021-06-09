#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Parses finite element mesh from the CalculiX .inp-file.
Reads nodes coordinates, elements composition,
node and element sets and surfaces """

# TODO Possible comments and empty lines
# TODO Parse mesh from INCLUDEs
# TODO Newly created elset is not highlighted

# Standard modules
import os
import re
import sys
import time
import logging
import textwrap
import traceback

# My modules
sys_path = os.path.abspath(__file__)
sys_path = os.path.dirname(sys_path)
sys_path = os.path.join(sys_path, '..', '..')
sys_path = os.path.normpath(sys_path)
sys_path = os.path.realpath(sys_path)
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
import importer
import tests
import log

"""
ifile - path to input file
icode - piece of input code
blocks - keyword blocks
old - mesh to be reparsed
"""

class Mesh:

    # Initialization
    def __init__(self, ifile=None,
        icode=None, blocks=None, old=None):

        self.nodes = {} # all mesh nodes with coordinates
        self.nsets = {} # node sets
        self.elements = {} # all mesh elements composition
        self.elsets = {} # element sets
        self.surfaces = {} # with corresponding nodes and element faces

        self.duplicated_warnings = {
            'NSET':[],
            'ELSET':[]} 
        self.warning_counter = 0

        # Mesh being reparsed
        self.old = old
        if not old:
            self.old = self

        # Get lines from INP source
        lines = []
        if ifile is not None: # whole the .inp-file
            lines = importer.read_lines(ifile)
        elif icode is not None:
            lines = icode # some piece of INP code
        elif blocks is not None:
            for b in blocks:
                lines.extend(b.get_inp_code()) # keyword blocks
        if not len(lines):
            self.warn('Nothing to parse!')
            return

        # Call parse methods for everything
        msg_text = 'Mesh parser:'
        for mname in ['nodes', 'nsets', 'elements', 'elsets', 'surfaces']:
            method = self.__dict__[mname]
            try:
                getattr(self, 'parse_' + mname)(lines)
                msg_text += '\n{} {}'.format(len(method), mname)
                # msg_text += str([v.name for v in method.values()])
                # for k,v in method.items():
                #     msg_text += '<br/>\n{0}: {1}'.format(k, v)
            except:
                msg = 'Can\'t parse {}\n'.format(mname) \
                    + traceback.format_exc()
                logging.error(msg)
        logging.info(msg_text)

    # Log some warning
    def warn(self, msg):
        if self.warning_counter < 10:
            logging.warning(msg)
        if self.warning_counter == 10:
            logging.warning('Message limit reached!')
        self.warning_counter += 1

    # # Implement case insensetivity for sets dictionaries
    # def get_set_by_name(self, sets, name):
    #     for s in sets.keys():
    #         if s.upper() == name.upper():
    #             return sets[s]
    #     logging.error('There is no set {}.'.format(name))

    # Parse nodes with coordinates - *NODE keyword
    def parse_nodes(self, lines):
        regex1 = r'^\*[\w\s-]+'
        regex2 = r'NSET\s*=\s*([\w\!\#\%\$\&\"\'\(\)\*\=\+\-\.\/\:\;\<\>\?\@\[\]\^\_\`\{\\\|\}\~]*)'
        for i in range(len(lines)):
            match = re.search(regex1, lines[i])
            if match is None:
                continue
            keyword_name = match.group(0)
            if keyword_name.upper() != '*NODE':
                continue

            nodes = []
            duplicated_nodes = []
            no_coordinates_nodes = []
            lead_line = lines[i]
            match = re.search(regex2, lead_line.upper())

            # Read the whole block
            while i+1<len(lines) and not lines[i+1].startswith('*'):
                if lines[i+1].startswith('**'):
                    i += 1
                    continue
                if not len(lines[i+1]):
                    i += 1
                    continue

                a = lines[i+1].replace(',', ' ').split() # to avoid redundant commas in the end of line
                num = int(a[0]) # node number

                coords = [float(coord) for coord in a[1:]] # node coordinates
                if len(coords) > 3:
                    coords = coords[:3]
                    # TODO Support Abaqus syntax
                    self.warn('Direction cosines are not supported.')

                # if len(coords) == 2: # in 2D case add Z coord equal to zero
                #     coords.append(0)

                # Create NODE
                node = NODE(num, coords)

                # Check duplicates and nodes without coordinates
                if num in self.nodes:
                    duplicated_nodes.append(num)
                    del node
                else:
                    nodes.append(node)
                    self.nodes[num] = node
                    # if len(coords):
                    #     nodes.append(node)
                    #     self.nodes[num] = node
                    # else:
                    #     no_coordinates_nodes.append(num)
                    #     del node

                i += 1

            # Warn about duplicated nodes
            if len(duplicated_nodes):
                msg_text = 'Duplicated nodes: {}.'.format(duplicated_nodes)
                for msg_line in textwrap.wrap(msg_text, width=40):
                    self.warn(msg_line)

            # Warn about nodes without coordinates
            if len(no_coordinates_nodes):
                msg_text = 'Nodes without coordinates: {}.'.format(no_coordinates_nodes)
                for msg_line in textwrap.wrap(msg_text, width=40):
                    self.warn(msg_line)

            # If all nodes are named as a set
            if match is not None:
                name = lead_line[match.start(1):match.end(1)]
                self.create_or_extend_set(self.nsets, name, nodes, NSET)

            # do not return to parse few *NODE sections

    # Parse node sets - *NSET keyword
    def parse_nsets(self, lines):
        rex = r'(\*NSET)\s*,.*NSET\s*=\s*' \
            + r'([\w\!\#\%\$\&\"\'\(\)\*\=\+\-\.\/\:\;\<\>\?\@\[\]\^\_\`\{\\\|\}\~]*)'
        for i in range(len(lines)):
            match = re.search(rex, lines[i].upper())
            if match is None:
                continue

            # name = lines[i][match.start(2):match.end(2)] # node set name
            name = match.group(2)
            nodes = []
            sets_with_non_existent_nodes = {}
            duplicated_nodes = []
            if not 'GENERATE' in lines[i].upper():
                while i+1<len(lines) and not lines[i+1].startswith('*'):
                    if lines[i+1].startswith('**'):
                        i += 1
                        continue
                    if not len(lines[i+1]):
                        i += 1
                        continue

                    a = lines[i+1].replace(',', ' ').split()
                    for n in a:
                        try:
                            # Single node number
                            node = self.old.nodes[int(n)]

                            # Check duplicates
                            if node in nodes:
                                duplicated_nodes.append(int(n))
                            else:
                                nodes.append(node)
                        except ValueError:
                            # Node set name
                            if n.upper() in self.old.nsets:
                                nodes.extend(self.old.nsets[n.upper()].items)
                            else:
                                msg = 'There is no NSET {}.'.format(n)
                                self.warn(msg)
                                # logging.error(msg)
                        except KeyError:
                            # Collect non-existent nodes by sets
                            if not name in sets_with_non_existent_nodes: 
                                sets_with_non_existent_nodes[name] = []
                            sets_with_non_existent_nodes[name].append(int(n))

                    i += 1
            else:
                try:
                    start, stop, step = re.split(',\s*', lines[i+1])
                except:
                    start, stop = re.split(',\s*', lines[i+1])
                    step = 1
                for n in list(range(int(start), int(stop)+1, int(step))):
                    try:
                        node = self.old.nodes[n]
                        nodes.append(node)
                    except KeyError:
                        # Collect non-existent nodes by sets
                        if not name in sets_with_non_existent_nodes: 
                            sets_with_non_existent_nodes[name] = []
                        sets_with_non_existent_nodes[name].append(n)

            # Warn about duplicated nodes in the mesh
            if len(duplicated_nodes):
                msg_text = 'Duplicated nodes {}.'.format(duplicated_nodes)
                for msg_line in textwrap.wrap(msg_text, width=40):
                    self.warn(msg_line)

            # Warn about non-existent nodes in the mesh
            for s, nodes in sets_with_non_existent_nodes.items():
                msg_text = 'NSET {} - mesh hasn\'t nodes {}.'.format(s, nodes)
                for msg_line in textwrap.wrap(msg_text, width=40):
                    self.warn(msg_line)

            self.create_or_extend_set(self.nsets, name, nodes, NSET)
            # do not return to parse few *NSET sections

    # Parse elements composition - *ELEMENT keyword
    def parse_elements(self, lines):
        regex1 = r'^\*[\w\s-]+'
        regex2 = r'TYPE\s*=\s*(\w+)'
        regex3 = r'ELSET\s*=\s*([\w\!\#\%\$\&\"\'\(\)\*\=\+\-\.\/\:\;\<\>\?\@\[\]\^\_\`\{\\\|\}\~]*)'
        for i in range(len(lines)):
            lead_line = lines[i].upper()

            match = re.search(regex1, lead_line)
            if not match:
                continue
            if match.group(0) != '*ELEMENT':
                continue

            match = re.search(regex2, lead_line) # element type
            etype = match.group(1)
            amount = self.amount_of_nodes(etype)
            elements = []
            duplicated_elements = []

            # Read the whole block
            while i+1<len(lines) and not lines[i+1].startswith('*'):
                if lines[i+1].startswith('**'):
                    i += 1
                    continue
                if not len(lines[i+1]):
                    i += 1
                    continue

                # Element nodes could be splitted into few lines
                a = lines[i+1].replace(',', ' ').split()
                while len(a) < amount + 1: # +1 for element number
                    a.extend(lines[i+2].replace(',', ' ').split())
                    i += 1

                num = int(a[0]) # element number

                nodes = []
                create_element = True
                for n in a[1:]: # iterate over element node numbers
                    if int(n) == 0: # it is possible in network element, type=D
                        n = int(a[2]) # take middle node to display in VTK
                    try:
                        node = self.old.nodes[int(n)]
                        nodes.append(node)
                    except KeyError:
                        msg_text = 'Element {} has no node {} and will be removed.'.format(num, n)
                        self.warn(msg_text)
                        create_element = False

                # Create ELEMENT
                if create_element:
                    element = ELEMENT(num, etype, nodes)

                    # Check duplicates
                    if num in self.elements:
                        duplicated_elements.append(num)
                        del element
                    else:
                        elements.append(element)
                        self.elements[num] = element

                i += 1

            # Warn about duplicated nodes
            if len(duplicated_elements):
                msg_text = 'Duplicated elements {}.'.format(duplicated_elements)
                for msg_line in textwrap.wrap(msg_text, width=40):
                    self.warn(msg_line)

            # If all elements are named as a set
            match = re.search(regex3, lead_line) # if all elements are united in a set
            if match is not None:
                name = lead_line[match.start(1):match.end(1)]
                self.create_or_extend_set(self.elsets, name, elements, ELSET)

            # do not return to parse few *ELEMENT sections

    # Parse element sets - *ELSET keyword
    def parse_elsets(self, lines):
        regex = r'(\*ELSET)\s*,.*ELSET\s*=\s*([\w\!\#\%\$\&\"\'\(\)\*\=\+\-\.\/\:\;\<\>\?\@\[\]\^\_\`\{\\\|\}\~]*)'
        for i in range(len(lines)):
            match = re.search(regex, lines[i].upper())
            if match is None:
                continue

            name = lines[i][match.start(2):match.end(2)] # element set name
            elements = []
            sets_with_non_existent_elements = {}

            if not 'GENERATE' in lines[i].upper():
                while i+1<len(lines) and not lines[i+1].startswith('*'):
                    if lines[i+1].startswith('**'):
                        i += 1
                        continue
                    if not len(lines[i+1]):
                        i += 1
                        continue

                    a = lines[i+1].replace(',', ' ').split()
                    for e in a:
                        try:
                            # Single element number
                            element = self.old.elements[int(e)]

                            # Check duplicates
                            if element in elements:
                                msg_text = 'Duplicated element {}.'.format(e)
                                self.warn(msg_text)
                            else:
                                elements.append(element)
                        except ValueError:
                            # Element set name
                            if e.upper() in self.old.elsets:
                                elements.extend(self.old.elsets[e.upper()].items)
                            else:
                                msg = 'There is no ELSET {}.'.format(e)
                                self.warn(msg)
                                # logging.error(msg)
                        except KeyError:
                            # Collect non-existent elements by sets
                            if not name in sets_with_non_existent_elements: 
                                sets_with_non_existent_elements[name] = []
                            sets_with_non_existent_elements[name].append(int(e))
                    i += 1
            else:
                try:
                    start, stop, step = re.split(',\s*', lines[i+1])
                except:
                    start, stop = re.split(',\s*', lines[i+1])
                    step = 1
                for e in list(range(int(start), int(stop)+1, int(step))):
                    try:
                        element = self.old.elements[e]
                        elements.append(element)
                    except KeyError:
                        # Collect non-existent elements by sets
                        if not name in sets_with_non_existent_elements: 
                            sets_with_non_existent_elements[name] = []
                        sets_with_non_existent_elements[name].append(e)

            # Warn about non-existent elements in the mesh
            for s, elements in sets_with_non_existent_elements.items():
                msg_text = 'ELSET {} - mesh hasn\'t elements {}.'.format(s, elements)
                for msg_line in textwrap.wrap(msg_text, width=40):
                    self.warn(msg_line)

            self.create_or_extend_set(self.elsets, name, elements, ELSET)
            # do not return to parse few *ELSET sections

    # Parse surfaces - *SURFACE keyword
    def parse_surfaces(self, lines):
        for i in range(len(lines)):
            skip = True

            # Surface name - required attribute
            name = ''
            match = re.search('\*SURFACE\s*,.*NAME\s*=\s*([\w\!\#\%\$\&\"\'\(\)\*\=\+\-\.\/\:\;\<\>\?\@\[\]\^\_\`\{\\\|\}\~]*)', lines[i].upper())
            if match is not None:
                name = lines[i][match.start(1):match.end(1)]
                skip = False

            # Surface type - optional attribute
            stype = 'ELEMENT' # 'ELEMENT' or 'NODE'
            match = re.search('\*SURFACE\s*,.*TYPE\s*=\s*(\w*)', lines[i].upper())
            if match is not None:
                stype = lines[i][match.start(1):match.end(1)]

            if not skip:
                if name + stype in self.surfaces:
                    msg_text = 'Duplicated surface name {}.'.format(name)
                    self.warn(msg_text)
                items = []

                while i+1<len(lines) and not lines[i+1].startswith('*'):
                    if lines[i+1].startswith('**'):
                        i += 1
                        continue
                    if not len(lines[i+1]):
                        i += 1
                        continue

                    _list = re.split(',\s*', lines[i+1])

                    if stype == 'ELEMENT':
                        """
                            TYPE=ELEMENT:
                            'surf3: [(1, S1), (2, S1), ...]'
                            'surf4: [(elset1, S2), (elset2, S2), ...]'
                        """

                        # Surface with element and face numbers
                        #   1, S1
                        #   2, S1
                        if re.match('^\d+,\s*S\d', lines[i+1]):
                            elem_num = int(_list[0])
                            surf_name = _list[1]
                            items.append((elem_num, surf_name))

                        # Surface with elset and face number
                        #   elset1, S1
                        #   elset2, S2
                        elif re.match('^[\w\-]+,\s*S\d', lines[i+1]):
                            elset_name = _list[0]
                            surf_name = _list[1]
                            if elset_name in self.old.elsets: 
                                for elem_num in self.old.elsets[elset_name].items:
                                    items.append((elem_num, surf_name))
                            else:
                                self.warn('In *SURFACE {} set {} not defined.'\
                                    .format(surf_name, elset_name))

                    elif stype == 'NODE':
                        """
                            TYPE=NODE:
                            'surf1: [1, 2, 3, ...]'
                            'surf2: [nset1, nset2, ...]
                        """
                        for n in _list:
                            if len(n):
                                try:
                                    # Single node number
                                    node = self.old.nodes[int(float(n))]
                                    items.append(node)
                                except ValueError:
                                    # Node set name
                                    if n.upper() in self.old.nsets:
                                        items.extend(self.old.nsets[n.upper()].items)
                                    else:
                                        msg = 'There is no NSET {}.'.format(n)
                                        self.warn(msg)
                                        # logging.error(msg)

                    i += 1

                # Create new SURFACE and append to list
                self.surfaces[name + stype] = SURFACE(name, items, stype)

    # Get amount of nodes by CalculiX element type
    def amount_of_nodes(self, etype):
        # regex = r'[A-Z]+(\dD)*(\d+)*[AEHIMOPRST]*\d?'
        if etype.startswith('Z'): # substructures
            return 0
        if etype.startswith(('U', 'VU')): # user element
            return 0
        if etype not in ['MASS', 'DASHPOTA', 'ITSUNI',
            'ITSCYL', 'GAPSPHER', 'DGAP', 'JOINTC']:
            while etype.endswith(tuple('ABCEHILMNOPRSTVW')):
                etype = etype[:-1]
        if not len(etype):
            raise SystemExit
        if etype.startswith('CAXA'):
            etype = etype[:5]
        try:
            return {
                  'AC1D2': 2,
                  'AC1D3': 3,
                  'AC2D3': 3,
                  'AC2D4': 4,
                  'AC2D6': 6,
                  'AC2D8': 8,
                  'ACAX3': 3,
                  'ACAX4': 4,
                  'ACAX6': 6,
                  'ACAX8': 8,
                  'AC3D4': 4,
                  'AC3D5': 5,
                  'AC3D6': 6,
                  'AC3D8': 8,
                 'AC3D10': 10,
                 'AC3D15': 15,
                 'AC3D20': 20,
                'ACIN2D2': 2,
                'ACIN2D3': 3,
                'ACIN3D3': 3,
                'ACIN3D4': 4,
                'ACIN3D6': 6,
                'ACIN3D8': 8,
                'ACINAX2': 2,
                'ACINAX3': 3,
                   'ASI1': 1,
                   'ASI2': 2,
                 'ASI2D2': 2,
                 'ASI2D3': 3,
                   'ASI3': 3,
                 'ASI3D3': 3,
                 'ASI3D4': 4,
                 'ASI3D6': 6,
                 'ASI3D8': 8,
                   'ASI4': 4,
                   'ASI8': 8,

                    'B21': 2,
                    'B22': 3,
                    'B23': 2,
                    'B31': 2,
                    'B32': 3,
                    'B33': 2,

                   'C3D4': 4,
                   'C3D5': 5,
                   'C3D6': 6,
                   'C3D8': 8,
                  'C3D10': 10,
                  'C3D15': 15,
                  'C3D20': 20,
                  'C3D27': 27,
                   'CCL9': 9,
                  'CCL12': 12,
                  'CCL18': 18,
                  'CCL24': 24,
                'CIN3D12': 12,
                'CIN3D18': 18,
                 'CIN3D8': 8,
                   'CAX3': 3,
                   'CAX4': 4,
                   'CAX6': 6,
                   'CAX8': 8,
                  'CAXA4': 4,
                  'CAXA8': 8,
                 'CINPE4': 4,
                 'CINPE5': 5,
                 'CINPS4': 4,
                 'CINPS5': 5,
                 'CINAX4': 4,
                 'CINAX5': 5,
                   'CPE3': 3,
                   'CPE4': 4,
                   'CPE6': 6,
                   'CPE8': 8,
                  'CPEG3': 3,
                  'CPEG4': 4,
                  'CPEG6': 6,
                  'CPEG8': 8,
                   'CPS3': 3,
                   'CPS4': 4,
                   'CPS6': 6,
                   'CPS8': 8,
                'CONN3D2': 2,
                'CONN2D2': 2,
                 'COH2D4': 4,
                 'COH3D6': 6,
                 'COH3D8': 8,
                 'COHAX4': 4,
                  'CGAX3': 3,
                  'CGAX4': 4,
                  'CGAX6': 6,
                  'CGAX8': 8,

                      'D': 3,
                    'DS3': 3,
                    'DS4': 4,
                    'DS6': 6,
                    'DS8': 8,
               'DASHPOTA': 2,
               'DASHPOT1': 1,
               'DASHPOT2': 2,
                  'DC1D2': 2,
                  'DC1D3': 3,
                  'DC2D3': 3,
                  'DC2D4': 4,
                  'DC2D6': 6,
                  'DC2D8': 8,
                  'DC3D4': 4,
                  'DC3D6': 6,
                  'DC3D8': 8,
                 'DC3D10': 10,
                 'DC3D15': 15,
                 'DC3D20': 20,
                  'DSAX1': 2,
                  'DSAX2': 3,
                  'DCAX3': 3,
                  'DCAX4': 4,
                  'DCAX6': 6,
                  'DCAX8': 8,
                'DCOUP2D': 1,
                'DCOUP3D': 1,
                 'DCC1D2': 2,
                'DCC1D2D': 2,
                 'DCC2D4': 4,
                'DCC2D4D': 4,
                 'DCC3D8': 8,
                'DCC3D8D': 8,
                 'DCCAX4': 4,
                'DCCAX4D': 4,
                 'DCCAX2': 2,
                'DCCAX2D': 2,
                   'DGAP': 2,
                 'DRAG2D': 1,
                 'DRAG3D': 2,

                  'EC3D8': 8,
                'ELBOW31': 2,
                'ELBOW32': 3,
                 'EMC2D3': 3,
                 'EMC2D4': 4,
                 'EMC3D4': 4,
                 'EMC3D8': 8,

                   'F3D8': 8,
                   'F3D6': 6,
                   'F3D4': 4,
                  'FP2D2': 2,
                 'FPC2D2': 2,
                 'FPC3D2': 2,
                  'FP3D2': 2,
                'FRAME2D': 2,
                'FRAME3D': 3,

                  'GAPCY': 2,
               'GAPSPHER': 2,
                   'GAPU': 2,
                  'GK2D2': 2,
                  'GK3D2': 2,
                  'GK3D4': 4,
                  'GK3D6': 6,
                  'GK3D8': 8,
                 'GK3D12': 12,
                 'GK3D18': 18,
                  'GKAX2': 2,
                  'GKAX4': 4,
                  'GKAX6': 6,
                  'GKPS4': 4,
                  'GKPS6': 6,
                  'GKPE4': 4,
                  'GKPE6': 6,

                  'HEATC': 1,

                  'IRS21': 3,
                  'ISL21': 2,
                  'ISL22': 3,
                  'ITT21': 1,
                  'ITT31': 1,
                 'ITSUNI': 2,
                 'ITSCYL': 2,

                 'JOINTC': 2,
                'JOINT2D': 2,
                'JOINT3D': 2,

                    'LS3': 3,
                    'LS6': 6,

                   'M3D3': 3,
                   'M3D4': 4,
                   'M3D6': 6,
                   'M3D8': 8,
                   'M3D9': 9,
                   'MASS': 1,
                   'MAX1': 2,
                   'MAX2': 3,
                   'MCL6': 6,
                   'MCL9': 9,
                  'MGAX1': 2,
                  'MGAX2': 3,

                   'PC3D': 1,
                 'PIPE21': 2,
                 'PIPE22': 3,
                 'PIPE31': 2,
                 'PIPE32': 3,
                  'PSI24': 4,
                  'PSI34': 4,
                  'PSI26': 6,
                  'PSI36': 6,

                   'Q3D4': 4,
                   'Q3D6': 6,
                   'Q3D8': 8,
                  'Q3D10': 10,
                  'Q3D20': 20,

                   'PD3D': 1,

                   'RAX2': 2,
                 'ROTARY': 1,
                   'R2D2': 2,
                   'R3D3': 3,
                   'R3D4': 4,
                  'RB2D2': 2,
                  'RB3D2': 2,

                     'S3': 3,
                     'S4': 4,
                   'S4R5': 4,
                     'S6': 6,
                     'S8': 8,
                   'S8R5': 8,
                   'S9R5': 9,
                   'SAX1': 2,
                   'SAX2': 3,
                 'SAXA11': 2,
                 'SAXA12': 2,
                 'SAXA13': 2,
                 'SAXA14': 2,
                 'SAXA21': 3,
                 'SAXA22': 3,
                 'SAXA23': 3,
                 'SAXA24': 3,
                    'SC6': 6,
                    'SC8': 8,
                'SFMGAX1': 2,
                'SFMGAX2': 2,
                 'SFM3D3': 3,
                 'SFM3D4': 4,
                 'SFM3D6': 6,
                 'SFM3D8': 8,
                 'SFMAX1': 2,
                 'SFMAX2': 3,
                  'STRI3': 3,
                 'STRI65': 6,
                 'SPRING': 2,
                'SPRING1': 1,
                'SPRING2': 2,
                 'SFMCL6': 6,
                 'SFMCL9': 9,

                   'T2D2': 2,
                   'T2D3': 3,
                   'T3D2': 2,
                   'T3D3': 3,

                'WARP2D3': 3,
                'WARP2D4': 4,
            }[etype]
        except KeyError:
            logging.error('Unknown element type - {}.'.format(etype))
            return 1 # minimum possible

    # Replace current mesh attributes with reparsed mesh ones
    """ Delete/add nodes
    Update NSETs, elements, ELSETs, surfaces
    Rebuild ugrid """
    def updateWith(self, reparsedMesh):
        for attrName, attrValue in reparsedMesh.__dict__.items():
            if type(attrValue) == dict and len(attrValue):
                # print('Another mesh:', attrName, attrValue)
                for _setName, _setValue in attrValue.items():
                    # print('Nodes:', _setValue.items)
                    getattr(self, attrName)[_setName] = _setValue

    # Parse input code and update current Mesh
    def reparse(self, icode):
        pass

    # Before modification checks if sets have set with the same name
    def create_or_extend_set(self, sets, name, items, klass):
        # logging.debug('Class: ' + klass.__name__)
        if name.upper() in sets: # check duplicates
            sets[name.upper()].items.extend(items) # append to existing set
            if name.upper() not in self.duplicated_warnings[klass.__name__]:
                self.warn('Duplicated set name {}!'.format(name))
                self.duplicated_warnings[klass.__name__].append(name.upper())
        else:
            sets[name.upper()] = klass(name, items) # create new set


class NODE:
    """
        1: [ 0.0, -1742.5, 0.0],
        2: [74.8, -1663.7, 0.0],
        ...
    """
    def __init__(self, num, coords):
        self.num = num
        self.name = str(num)
        
        # Default coordinates (to support Abaqus)
        self.coords = [0]*3
        for i in range(len(coords)):
            self.coords[i] = coords[i]


class NSET:
    """
        'nset1': [1, 2, 3, 4],
        'nset2': [5, 6, 7, 8],
        ...
    """
    def __init__(self, name, nodes):
        self.name = name
        self.items = nodes


class ELEMENT:
    """
        1: [1, 2],
        2: [3, 4],
        ...
        11: [21, 22, 23],
        12: [24, 25, 26],
        ...
    """
    def __init__(self, num, etype, nodes):
        self.num = num
        self.name = str(num)
        self.type = etype
        self.nodes = nodes

        # Calculate centroid
        x = sum([node.coords[0] for node in self.nodes]) / len(nodes)
        y = sum([node.coords[1] for node in self.nodes]) / len(nodes)
        z = sum([node.coords[2] for node in self.nodes]) / len(nodes)
        self.centroid = [x, y, z] # coordinates of element center


class ELSET:
    """
        'elset1': [1, 2, 3, 4],
        'elset2': [5, 6, 7, 8],
        ...
    """
    def __init__(self, name, elements):
        self.name = name
        self.items = elements


class SURFACE:
    """
        TYPE=NODE:
        'surf1: [1, 2, 3, ...]'
        'surf2: [nset1, nset2, ...]

        TYPE=ELEMENT:
        'surf3: [(1, S1), (2, S1), ...]'
        'surf4: [(elset1, S2), (elset2, S2), ...]'
    """
    def __init__(self, name, items, stype=None):
        self.name = name
        self.items = items
        if stype:
            self.type = stype
        else:
            self.type = 'ELEMENT'


# Run test
# Test mesh parser on all CalculiX examples
@tests.test_wrapper()
def test():

    # Prepare logging
    log_file = __file__[:-3] + '.log'
    log.stop_logging()
    log.add_my_handler(logging.DEBUG)
    log.print(log_file, 'MESH PARSER TEST')

    limit = 50000 # how many files to process
    # examples_dir = '../../../../examples/ccx/test'
    examples_dir = '../../../../examples'
    counter = 0

    examples = tests.scan_all_files_in(examples_dir, '.inp', limit)
    for file_name in examples:
        counter += 1
        relpath = os.path.relpath(file_name, start=os.getcwd())
        log.print(log_file, '\n{}\n{}: {}'.format('='*50, counter, relpath))

        # Parse mesh
        m = Mesh(ifile=file_name)

if __name__ == '__main__':
    test()
