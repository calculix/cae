#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, June 2020
Distributed under GNU General Public License v3.0

Parses finite element mesh from the CalculiX .inp-file.
Reads nodes coordinates, elements composition,
node and element sets and surfaces """

# Standard modules
import os
import sys
import re
import time
import logging
import textwrap
import traceback

# My modules
sys_path = os.path.dirname(__file__)
sys_path = os.path.join(sys_path, '..', '..')
sys_path = os.path.normpath(sys_path)
sys.path.append(sys_path)
import file_tools
import tests
import clean


class Mesh:

    # Initialization
    def __init__(self, INP_file=None, inp_code=None, old=None):
        self.nodes = {} # all mesh nodes with coordinates
        self.nsets = {} # node sets
        self.elements = {} # all mesh elements composition
        self.elsets = {} # element sets
        self.surfaces = {} # with corresponding nodes and element faces

        # Mesh being reparsed
        self.old = old
        if not old:
            self.old = self

        # Get lines from INP source
        if INP_file and not inp_code:
            # Open and read whole the .inp-file
            lines = file_tools.read_lines(INP_file)
        elif inp_code and not INP_file:
            # Parse some piece of INP code
            lines = inp_code
        else:
            logging.warning('Nothing to parse!')
            return

        # Call parse methods for everything
        msg_text = 'Mesh parser:'
        for attrName, attrValue in self.__dict__.items():
            if type(attrValue) == dict:
                try:
                    getattr(self, 'parse_' + attrName)(lines)
                    msg_text += '\n{} {}'.format(len(attrValue), attrName)
                    # msg_text += str([v.name for v in attrValue.values()])
                    # for k,v in attrValue.items():
                    #     msg_text += '<br/>\n{0}: {1}'.format(k, v)
                except:
                    msg = 'Can\'t parse {}\n'.format(attrName) \
                        + traceback.format_exc()
                    logging.error(msg)
        logging.info(msg_text)

    # Parse nodes with coordinates - *NODE keyword
    def parse_nodes(self, lines):
        for i in range(len(lines)):

            # Distinguish 'NODE' and 'NODE PRINT'
            if ',' in lines[i]:
                keyword_name = lines[i].split(',')[0]
            else:
                keyword_name = lines[i]

            if keyword_name.upper() == '*NODE':
                nodes = []
                duplicated_nodes = []
                no_coordinates_nodes = []
                lead_line = lines[i]
                match = re.search('NSET\s*=\s*([\w\!\#\%\$\&\"\'\(\)\*\=\+\-\.\/\:\;\<\>\?\@\[\]\^\_\`\{\\\|\}\~]*)', lead_line.upper())

                # Read the whole block - there will be no comments
                while i+1<len(lines) and not lines[i+1].startswith('*'):
                    a = lines[i+1].replace(',', ' ').split() # to avoid redundant commas in the end of line
                    num = int(a[0]) # node number
                    coords = [float(coord) for coord in a[1:]] # node coordinates
                    if len(coords) == 2: # in 2D case add Z coord equal to zero
                        coords.append(0)

                    # Create NODE
                    node = NODE(num, coords)

                    # Check duplicates and nodes without coordinates
                    if num in self.nodes:
                        duplicated_nodes.append(num)
                        del node
                    else:
                        if len(coords):
                            nodes.append(node)
                            self.nodes[num] = node
                        else:
                            no_coordinates_nodes.append(num)
                            del node

                    i += 1

                # Warn about duplicated nodes
                if len(duplicated_nodes):
                    msg_text = 'Duplicated nodes: {}.'.format(duplicated_nodes)
                    for msg_line in textwrap.wrap(msg_text, width=40):
                        logging.warning(msg_line)

                # Warn about nodes without coordinates
                if len(no_coordinates_nodes):
                    msg_text = 'Nodes without coordinates: {}.'.format(no_coordinates_nodes)
                    for msg_line in textwrap.wrap(msg_text, width=40):
                        logging.warning(msg_line)

                # If all nodes are named as a set
                if match:
                    name = lead_line[match.start(1):match.end(1)]
                    create_or_extend_set(self.nsets, name, nodes, NSET)

                # do not return to parse few *NODE sections

    # Parse node sets - *NSET keyword
    def parse_nsets(self, lines):
        for i in range(len(lines)):
            match = re.search('(\*NSET)\s*,.*NSET\s*=\s*([\w\!\#\%\$\&\"\'\(\)\*\=\+\-\.\/\:\;\<\>\?\@\[\]\^\_\`\{\\\|\}\~]*)', lines[i].upper())
            if match:
                name = lines[i][match.start(2):match.end(2)] # node set name
                nodes = []
                sets_with_non_existent_nodes = {}
                duplicated_nodes = []
                if not 'GENERATE' in lines[i].upper():
                    while i+1<len(lines) and not lines[i+1].startswith('*'):
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
                                nodes.extend(self.old.nsets[n].items)
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
                        logging.warning(msg_line)

                # Warn about non-existent nodes in the mesh
                for s, nodes in sets_with_non_existent_nodes.items():
                    msg_text = 'NSET {} - mesh hasn\'t nodes {}.'.format(s, nodes)
                    for msg_line in textwrap.wrap(msg_text, width=40):
                        logging.warning(msg_line)

                create_or_extend_set(self.nsets, name, nodes, NSET)
                # do not return to parse few *NSET sections

    # Parse elements composition - *ELEMENT keyword
    def parse_elements(self, lines):
        for i in range(len(lines)):
            if lines[i].upper().startswith('*ELEMENT'):
                lead_line = lines[i]
                match = re.search('TYPE\s*=\s*(\w*)', lead_line.upper())
                etype = lead_line[match.start(1):match.end(1)] # element type
                amount = self.amount_of_nodes(etype)
                elements = []
                duplicated_elements = []
                match = re.search('ELSET\s*=\s*([\w\!\#\%\$\&\"\'\(\)\*\=\+\-\.\/\:\;\<\>\?\@\[\]\^\_\`\{\\\|\}\~]*)', lead_line.upper()) # if all elements are united in a set

                # Read the whole block - there will be no comments
                while i+1<len(lines) and not lines[i+1].startswith('*'):

                    # Element nodes could be splitted into few lines
                    a = lines[i+1].replace(',', ' ').split()
                    while len(a) < amount + 1: # +1 for element number
                        a.extend( lines[i+2].replace(',', ' ').split() )
                        i += 1

                    num = int(a[0]) # element number

                    nodes = []
                    create_element = True
                    for n in a[1:]: # iterate over element's node numbers
                        if int(n) == 0: # it is possible in network element, type=D
                            n = int(a[2]) # take middle node to display in VTK
                        try:
                            node = self.old.nodes[int(n)]
                            nodes.append(node)
                        except KeyError:
                            msg_text = 'Element {} has no node {} and will be removed.'.format(num, n)
                            logging.warning(msg_text)
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
                        logging.warning(msg_line)

                # If all elements are named as a set
                if match:
                    name = lead_line[match.start(1):match.end(1)]
                    create_or_extend_set(self.elsets, name, elements, ELSET)

                # do not return to parse few *ELEMENT sections

    # Parse element sets - *ELSET keyword
    def parse_elsets(self, lines):
        for i in range(len(lines)):
            match = re.search('(\*ELSET)\s*,.*ELSET\s*=\s*([\w\!\#\%\$\&\"\'\(\)\*\=\+\-\.\/\:\;\<\>\?\@\[\]\^\_\`\{\\\|\}\~]*)', lines[i].upper())
            if match:
                name = lines[i][match.start(2):match.end(2)] # element set name
                elements = []
                sets_with_non_existent_elements = {}

                if not 'GENERATE' in lines[i].upper():
                    while i+1<len(lines) and not lines[i+1].startswith('*'):
                        a = lines[i+1].replace(',', ' ').split()
                        for e in a:
                            try:
                                # Single element number
                                element = self.old.elements[int(e)]

                                # Check duplicates
                                if element in elements:
                                    msg_text = 'Duplicated element {}.'.format(e)
                                    logging.warning(msg_text)
                                else:
                                    elements.append(element)
                            except ValueError:
                                # Element set name
                                elements.extend(self.old.elsets[e].items)
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
                        logging.warning(msg_line)

                create_or_extend_set(self.elsets, name, elements, ELSET)
                # do not return to parse few *ELSET sections

    # Parse surfaces - *SURFACE keyword
    def parse_surfaces(self, lines):
        for i in range(len(lines)):
            skip = True

            # Surface name - required attribute
            name = ''
            match = re.search('\*SURFACE\s*,.*NAME\s*=\s*([\w\!\#\%\$\&\"\'\(\)\*\=\+\-\.\/\:\;\<\>\?\@\[\]\^\_\`\{\\\|\}\~]*)', lines[i].upper())
            if match:
                name = lines[i][match.start(1):match.end(1)]
                skip = False

            # Surface type - optional attribute
            stype = 'ELEMENT' # 'ELEMENT' or 'NODE'
            match = re.search('\*SURFACE\s*,.*TYPE\s*=\s*(\w*)', lines[i].upper())
            if match:
                stype = lines[i][match.start(1):match.end(1)]

            if not skip:
                if name + stype in self.surfaces:
                    msg_text = 'Duplicated surface name {}.'.format(name)
                    logging.warning(msg_text)
                items = []

                while i+1<len(lines) and not lines[i+1].startswith('*'):
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
                                logging.warning('In *SURFACE {} set {} not defined.'\
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
                                    node = self.old.nodes[int(n)]
                                    items.append(node)
                                except ValueError:
                                    # Node set name
                                    items.extend(self.old.nsets[n].items)

                    i += 1

                # Create new SURFACE and append to list
                self.surfaces[name + stype] = SURFACE(name, items, stype)

    # Get amount of nodes by CalculiX element type
    def amount_of_nodes(self, etype):
        try:
            return {
                   'C3D8': 8,
                   'F3D8': 8,
                  'C3D8R': 8,
                  'C3D8I': 8,
                   'C3D6': 6,
                   'F3D6': 6,
                   'C3D4': 4,
                   'F3D4': 4,
                  'C3D20': 20,
                 'C3D20R': 20,
                  'C3D15': 15,
                  'C3D10': 10,
                 'C3D10T': 10,
                     'S3': 3,
                   'M3D3': 3,
                   'CPS3': 3,
                   'CPE3': 3,
                   'CAX3': 3,
                     'S6': 6,
                   'M3D6': 6,
                   'CPS6': 6,
                   'CPE6': 6,
                   'CAX6': 6,
                     'S4': 4,
                    'S4R': 4,
                   'M3D4': 4,
                  'M3D4R': 4,
                   'CPS4': 4,
                  'CPS4R': 4,
                   'CPE4': 4,
                  'CPE4R': 4,
                   'CAX4': 4,
                  'CAX4R': 4,
                     'S8': 8,
                    'S8R': 8,
                   'M3D8': 8,
                  'M3D8R': 8,
                   'CPS8': 8,
                  'CPS8R': 8,
                   'CPE8': 8,
                  'CPE8R': 8,
                   'CAX8': 8,
                  'CAX8R': 8,
                    'B21': 2,
                    'B31': 2,
                   'B31R': 2,
                   'T2D2': 2,
                   'T3D2': 2,
                 'GAPUNI': 2,
               'DASHPOTA': 2,
                'SPRING2': 2,
                'SPRINGA': 2,
                    'B32': 3,
                   'B32R': 3,
                   'T3D3': 3,
                      'D': 3,
                'SPRING1': 1,
                'DCOUP3D': 1,
                   'MASS': 1,
            } [etype]
        except:
            return 2 # minimum possible

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

    # Parse inp_code and update current Mesh
    def reparse(self, inp_code):
        pass


# Before modification checks if sets have set with the same name
def create_or_extend_set(sets, name, items, klass):
    if name in sets: # check duplicates
        sets[name].items.extend(items) # append to existing set
        logging.warning('Duplicated set name {}!'.format(name))
    else:
        sets[name] = klass(name, items) # create new set


class NODE:
    """
        1: [ 0.0, -1742.5, 0.0],
        2: [74.8, -1663.7, 0.0],
        ...
    """
    def __init__(self, num, coords):
        self.num = num
        self.name = str(num)
        self.coords = coords


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
        self.centroid = [x, y, z] # coordinates of element's center


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


# Test mesh parser on all CalculiX examples
if __name__ == '__main__':
    clean.screen()
    os.chdir(os.path.dirname(__file__))
    start_time = time.perf_counter()
    print = tests.print

    # Prepare logging
    log_file = __file__[:-3] + '.log'
    h = tests.myHandler(log_file)
    logging.getLogger().addHandler(h)
    logging.getLogger().setLevel(logging.INFO)

    limit = 3000 # how many files to process
    examples_dir = '../../../../examples/ccx_2.16.test'
    counter = 0

    print(log_file, 'MESH PARSER TEST\n\n')
    examples = tests.scan_all_files_in(examples_dir, '.inp', limit)
    for file_name in examples:
        counter += 1
        relpath = os.path.relpath(file_name, start=os.getcwd())
        print(log_file, '\n{}\n{}: {}'.format('='*50, counter, relpath))

        # Parse mesh
        m = Mesh(INP_file=file_name)

    print(log_file, '\nTotal {:.1f} seconds.'
        .format(time.perf_counter() - start_time))
    clean.cache()
