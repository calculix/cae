# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    Parses finite element mesh from the CalculiX .inp-file.
    Reads nodes coordinates, elements composition, node and element sets, surfaces.
    It's case insensitive and translates all text uppercase.
"""


import os, re, logging


class Parse:


    # Initialization
    def __init__(self, inp_file):
        self.nodes = {} # all mesh nodes with coordinates
        self.elements = {} # all mesh elements composition
        self.nsets = {} # node sets
        self.elsets = {} # element sets
        self.surfaces = {} # surfaces with corresponding nodes and element faces

        # Mesh bounds to avoid camera flying to infinity
        self.bounds = [1e+6,-1e+6]*3 # Xmin,Xmax, Ymin,Ymax, Zmin,Zmax

        # Open and read whole the .inp-file
        lines = self.parse_lines(inp_file)

        # Parse nodes
        try:
            self.parse_nodes(lines)

            msg_text = '{} nodes'.format(len(self.nodes))
            # msg_text += ': ' + str(list(self.nodes.keys()))
            logging.info(msg_text)
        except:
            logging.error('Can\'t parse nodes')

        # Parse node sets
        try:
            self.parse_nsets(lines)

            msg_text = '{} nsets'.format(len(self.nsets))
            # msg_text += ': ' + str(list(self.nsets.keys()))
            # for k,v in self.nsets.items():
            #     msg_text += '<br/>\n{0}: {1}'.format(k, v)
            logging.info(msg_text)
        except:
            logging.error('Can\'t parse nsets')

        # Parse elements
        try:
            self.parse_elements(lines)

            msg_text = '{} elements'.format(len(self.elements))
            # msg_text += ': ' + str(list(self.elements.keys()))
            logging.info(msg_text)
        except:
            logging.error('Can\'t parse elements')

        # Parse element sets
        try:
            self.parse_elsets(lines)

            msg_text = '{} elsets'.format(len(self.elsets))
            # msg_text += ': ' + str(list(self.elsets.keys()))
            # for k,v in self.elsets.items():
            #     msg_text += '<br/>\n{0}: {1}'.format(k, v)
            logging.info(msg_text)
        except:
            logging.error('Can\'t parse elsets')

        # Parse surfaces
        try:
            self.parse_surfaces(lines)

            msg_text = '{} surfaces'.format(len(self.surfaces))
            # msg_text += ': ' + str(self.surfaces)
            logging.info(msg_text)
        except:
            logging.error('Can\'t parse surfaces')


    # Parse nodes with coordinates - *NODE keyword
    def parse_nodes(self, lines):
        for i in range(len(lines)): # lines are uppercase

            # Distinguish 'NODE' and 'NODE PRINT'
            if ',' in lines[i]:
                keyword_name = lines[i].split(',')[0]
            else:
                keyword_name = lines[i]

            if keyword_name == '*NODE':
                nodes = []
                match = re.search('NSET\s*=\s*(\w*)', lines[i]) # if all nodes are united in a set

                while i+1<len(lines) and not lines[i+1].startswith('*'): # read the whole block and return
                    a = lines[i+1].replace(',', ' ').split() # to avoid redundant commas in the end of line
                    num = int(a[0]) # node number
                    coords = [float(coord) for coord in a[1:]] # node coordinates
                    if len(coords) == 2: # in 2D case add Z coord equal to zero
                        coords.append(0)

                    # Create NODE
                    node = NODE(num, coords)

                    # Check duplicates
                    if num in self.nodes:
                        msg_text = 'Duplicated node {}.'.format(num)
                        logging.warning(msg_text)
                        del node
                    else:
                        if len(coords):
                            nodes.append(node)
                            self.nodes[num] = node

                            # Bounding box
                            for j,coord in enumerate(coords):
                                if coord < self.bounds[j*2]:
                                    self.bounds[j*2] = coord # update min coords values
                                if coord > self.bounds[j*2+1]:
                                    self.bounds[j*2+1] = coord # update max coords values
                        else:
                            msg_text = 'Node {} has no coordinates and will be removed.'.format(num)
                            logging.warning(msg_text)
                            del node

                    i += 1
 
                # If there is node set name for all nodes
                if match:
                    name = match.group(1)
                    self.nsets[name] = NSET(name, nodes)

                # do not return to parse few *NODE sections


    # Parse node sets - *NSET keyword
    def parse_nsets(self, lines):
        for i in range(len(lines)): # lines are uppercase
            match = re.search('(\*NSET)\s*,.*NSET\s*=\s*(\w*)', lines[i])
            if match:
                name = match.group(2) # node set name
                nodes = []

                if not 'GENERATE' in lines[i]:
                    while i+1<len(lines) and not lines[i+1].startswith('*'):
                        a = lines[i+1].replace(',', ' ').split()
                        for n in a:
                            try:
                                # Single node number
                                node = self.nodes[int(n)]
                                nodes.append(node)
                            except ValueError:
                                # Node set name
                                nodes.extend(self.nsets[n].nodes)
                            except KeyError:
                                msg_text = 'NSET {} - there is no node {} in the mesh.'.format(name, n)
                                logging.warning(msg_text)
                        i += 1
                else:
                    try:
                        start, stop, step = re.split(',\s*', lines[i+1])
                    except:
                        start, stop = re.split(',\s*', lines[i+1])
                        step = 1
                    for n in list(range(int(start), int(stop)+1, int(step))):
                        try:
                            node = self.nodes[n]
                            nodes.append(node)
                        except KeyError:
                            msg_text = 'NSET {} - there is no node {} in the mesh.'.format(name, n)
                            logging.warning(msg_text)

                # Create node set
                node_set = NSET(name, nodes)

                # Check duplicates
                if name in self.nsets:
                    msg_text = 'Duplicated nset name {}.'.format(name)
                    logging.warning(msg_text)
                else:
                    self.nsets[name] = node_set

                # do not return to parse few *NSET sections


    # Parse elements composition - *ELEMENT keyword
    def parse_elements(self, lines):
        for i in range(len(lines)): # lines are uppercase
            if lines[i].startswith('*ELEMENT'):
                match = re.search('TYPE\s*=\s*(\w*)', lines[i])
                etype = match.group(1) # element type
                amount = self.amount_of_nodes(etype)
                elements = []
                match = re.search('ELSET\s*=\s*(\w*)', lines[i]) # if all elements are united in a set

                while i+1<len(lines) and not lines[i+1].startswith('*'): # there will be no comments

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
                            node = self.nodes[int(n)]
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
                            msg_text = 'Duplicated element {}.'.format(num)
                            logging.warning(msg_text)
                            del element
                        else:
                            elements.append(element)
                            self.elements[num] = element

                    i += 1
 
                # If there is element set name for all elements
                if match:
                    name = match.group(1)
                    self.elsets[name] = ELSET(name, elements)

                # do not return to parse few *ELEMENT sections


    # Parse element sets - *ELSET keyword
    def parse_elsets(self, lines):
        for i in range(len(lines)): # lines are uppercase
            match = re.search('(\*ELSET)\s*,.*ELSET\s*=\s*(\w*)', lines[i])
            if match:
                name = match.group(2) # element set name
                elements = []

                if not 'GENERATE' in lines[i]:
                    while i+1<len(lines) and not lines[i+1].startswith('*'):
                        a = lines[i+1].replace(',', ' ').split()
                        for e in a:
                            try:
                                # Single element number
                                element = self.elements[int(e)]
                                elements.append(element)
                            except ValueError:
                                # Element set name
                                elements.extend(self.elsets[e].elements)
                            except KeyError:
                                msg_text = 'ELSET {} - there is no element {} in the mesh.'.format(name, e)
                                logging.warning(msg_text)
                        i += 1
                else:
                    try:
                        start, stop, step = re.split(',\s*', lines[i+1])
                    except:
                        start, stop = re.split(',\s*', lines[i+1])
                        step = 1
                    for e in list(range(int(start), int(stop)+1, int(step))):
                        try:
                            element = self.elements[e]
                            elements.append(element)
                        except KeyError:
                            msg_text = 'ELSET {} - there is no element {} in the mesh.'.format(name, e)
                            logging.warning(msg_text)

                # Create element set
                element_set = ELSET(name, elements)

                # Check duplicates
                if name in self.elsets:
                    msg_text = 'Duplicated elset name {}.'.format(name)
                    logging.warning(msg_text)
                else:
                    self.elsets[name] = element_set


    # Parse surfaces - *SURFACE keyword
    def parse_surfaces(self, lines):
        for i in range(len(lines)): # lines are uppercase
            skip = True

            # Surface name - required attribute
            name = ''
            match = re.search('\*SURFACE\s*,.*NAME\s*=\s*(\w*)', lines[i])
            if match:
                name = match.group(1)
                skip = False

            # Surface type - optional attribute
            stype = 'ELEMENT' # 'ELEMENT' or 'NODE'
            match = re.search('\*SURFACE\s*,.*TYPE\s*=\s*(\w*)', lines[i])
            if match:
                stype = match.group(1)

            if not skip:
                if name + stype in self.surfaces:
                    msg_text = 'Duplicated surface name {}.'.format(name)
                    logging.warning(msg_text)
                _set = []

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
                        if re.match('^\d+,\s*S\d+', lines[i+1]):
                            _set.append(tuple(_list))

                        # Surface with elset and face number
                        #   elset1, S1
                        #   elset2, S2
                        elif re.match('^\w+,\s*S\d+', lines[i+1]):
                            elset_name = _list[0]
                            surf_name = _list[1]
                            for element in self.elsets[elset_name].elements:
                                _set.append((element, surf_name))

                    elif stype == 'NODE':
                        """
                            TYPE=NODE:
                            'surf1: [1, 2, 3, ...]'
                            'surf2: [nset1, nset2, ...]
                        """
                        for n in _list:
                            if len(n):
                                try:
                                    _set.append(int(n))
                                except ValueError:
                                    _set.extend(self.nsets[n].nodes)

                    i += 1

                # Create new SURFACE and append to list
                self.surfaces[name + stype] = SURFACE(name, _set, stype)


    # Recurcively read all the lines of the file and its includes.
    # For simplicity convert all lines to uppercase.
    def parse_lines(self, inp_file):
        lines = []
        try:
            inp_file = os.path.abspath(inp_file)
            with open(inp_file, 'r') as f:
                for line in f.readlines():
                    line = line.strip()
                    if (not line.startswith('**')) and len(line): # skip comments and empty lines
                        lines.append(line.upper())

                        # Append lines from include file
                        if line.upper().startswith('*INCLUDE'):
                            path = os.path.dirname(inp_file)
                            inp_file2 = line.split('=')[1].strip()
                            lines.extend(self.parse_lines(path + '/' + inp_file2))
        except:
            msg_text = 'There is no file {}.'.format(inp_file)
            logging.error(msg_text)
        return lines


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


class NODE:
    """
        1: [ 0.0, -1742.5, 0.0],
        2: [74.8, -1663.7, 0.0],
        ...
    """
    def __init__(self, num, coords):
        self.num = num
        self.coords = coords


class NSET:
    """
        'nset1': [1, 2, 3, 4],
        'nset2': [5, 6, 7, 8],
        ...
    """
    def __init__(self, name, nodes):
        self.name = name
        self.nodes = nodes


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
        self.elements = elements


class SURFACE:
    """
        TYPE=NODE:
        'surf1: [1, 2, 3, ...]'
        'surf2: [nset1, nset2, ...]

        TYPE=ELEMENT:
        'surf3: [(1, S1), (2, S1), ...]'
        'surf4: [(elset1, S2), (elset2, S2), ...]'
    """
    def __init__(self, name, _set, stype=None):
        self.name = name
        self.set = _set
        if stype:
            self.type = stype
        else:
            self.type = 'ELEMENT'


# Convert Calculix element type to VTK
def frd2vtk(frd_elem_type):
    """
        Keep in mind that CalculiX expands shell elements
        In vtk elements nodes are numbered starting from 0, not from 1

        For frd see http://www.dhondt.de/cgx_2.15.pdf pages 117-123 (chapter 10)
        For vtk see https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf pages 9-10

            _________________________________________________________________
        |                               |                                 |
        | №№      CalculiX type         |  №№         VTK type            |
        |_______________________________|_________________________________|
        |    |          |               |      |                          |
        |  1 | C3D8     |  8 node brick | = 12 | VTK_HEXAHEDRON           |
        |    | F3D8     |               |      |                          |
        |    | C3D8R    |               |      |                          |
        |    | C3D8I    |               |      |                          |
        |____|__________|_______________|______|__________________________|
        |    |          |               |      |                          |
        |  2 | C3D6     |  6 node wedge | = 13 | VTK_WEDGE                |
        |    | F3D6     |               |      |                          |
        |____|__________|_______________|______|__________________________|
        |    |          |               |      |                          |
        |  3 | C3D4     |  4 node tet   | = 10 | VTK_TETRA                |
        |    | F3D4     |               |      |                          |
        |____|__________|_______________|______|__________________________|
        |    |          |               |      |                          |
        |  4 | C3D20    | 20 node brick | = 25 | VTK_QUADRATIC_HEXAHEDRON |
        |    | C3D20R   |               |      |                          |
        |____|__________|_______________|______|__________________________|
        |    |          |               |      |                          |
        |  5 | C3D15    | 15 node wedge | ~ 13 | VTK_WEDGE                |
        |____|__________|_______________|______|__________________________|
        |    |          |               |      |                          |
        |  6 | C3D10    | 10 node tet   | = 24 | VTK_QUADRATIC_TETRA      |
        |    | C3D10T   |               |      |                          |
        |____|__________|_______________|______|__________________________|
        |    |          |               |      |                          |
        |  7 | S3       |  3 node shell | =  5 | VTK_TRIANGLE             |
        |    | M3D3     |               |      |                          |
        |    | CPS3     |               |      |                          |
        |    | CPE3     |               |      |                          |
        |    | CAX3     |               |      |                          |
        |____|__________|_______________|______|__________________________|
        |    |          |               |      |                          |
        |  8 | S6       |  6 node shell | = 22 | VTK_QUADRATIC_TRIANGLE   |
        |    | M3D6     |               |      |                          |
        |    | CPS6     |               |      |                          |
        |    | CPE6     |               |      |                          |
        |    | CAX6     |               |      |                          |
        |____|__________|_______________|______|__________________________|
        |    |          |               |      |                          |
        |  9 | S4       |  4 node shell | =  9 | VTK_QUAD                 |
        |    | S4R      |               |      |                          |
        |    | M3D4     |               |      |                          |
        |    | M3D4R    |               |      |                          |
        |    | CPS4     |               |      |                          |
        |    | CPS4R    |               |      |                          |
        |    | CPE4     |               |      |                          |
        |    | CPE4R    |               |      |                          |
        |    | CAX4     |               |      |                          |
        |    | CAX4R    |               |      |                          |
        |____|__________|_______________|______|__________________________|
        |    |          |               |      |                          |
        | 10 | S8       |  8 node shell | = 23 | VTK_QUADRATIC_QUAD       |
        |    | S8R      |               |      |                          |
        |    | M3D8     |               |      |                          |
        |    | M3D8R    |               |      |                          |
        |    | CPS8     |               |      |                          |
        |    | CPS8R    |               |      |                          |
        |    | CPE8     |               |      |                          |
        |    | CPE8R    |               |      |                          |
        |    | CAX8     |               |      |                          |
        |    | CAX8R    |               |      |                          |
        |____|__________|_______________|______|__________________________|
        |    |          |               |      |                          |
        | 11 | B21      |  2 node beam  | =  3 | VTK_LINE                 |
        |    | B31      |               |      |                          |
        |    | B31R     |               |      |                          |
        |    | T2D2     |               |      |                          |
        |    | T3D2     |               |      |                          |
        |    | GAPUNI   |               |      |                          |
        |    | DASHPOTA |               |      |                          |
        |    | SPRING2  |               |      |                          |
        |    | SPRINGA  |               |      |                          |
        |____|__________|_______________|______|__________________________|
        |    |          |               |      |                          |
        | 12 | B32      |  3 node beam  | = 21 | VTK_QUADRATIC_EDGE       |
        |    | B32R     |               |      |                          |
        |    | T3D3     |               |      |                          |
        |    | D        |               |      |                          |
        |____|__________|_______________|______|__________________________|
        |    |          |               |      |                          |
        | ?? | SPRING1  |  1 node       | =  1 | VTK_VERTEX               |
        |    | DCOUP3D  |               |      |                          |
        |    | MASS     |               |      |                          |
        |____|__________|_______________|______|__________________________|
    """
    # frd_elem_type : vtk_elem_type
    frd2vtk_num = {
            1: 12,
            2: 13,
            3: 10,
            4: 25,
            5: 13,
            6: 24,
            7:  5,
            8: 22,
            9:  9,
        10: 23,
        11:  3,
        12: 21,
        }
    frd2vtk_txt = {
            'C3D8':12,
            'F3D8':12,
            'C3D8R':12,
            'C3D8I':12,
            'C3D6':13,
            'F3D6':13,
            'C3D4':10,
            'F3D4':10,
            'C3D20':25,
            'C3D20R':25,
            'C3D15':13,
            'C3D10':24,
            'C3D10T':24,
                'S3':5,
            'M3D3':5,
            'CPS3':5,
            'CPE3':5,
            'CAX3':5,
                'S6':22,
            'M3D6':22,
            'CPS6':22,
            'CPE6':22,
            'CAX6':22,
                'S4':9,
                'S4R':9,
            'M3D4':9,
            'M3D4R':9,
            'CPS4':9,
            'CPS4R':9,
            'CPE4':9,
            'CPE4R':9,
            'CAX4':9,
            'CAX4R':9,
                'S8':23,
                'S8R':23,
            'M3D8':23,
            'M3D8R':23,
            'CPS8':23,
            'CPS8R':23,
            'CPE8':23,
            'CPE8R':23,
            'CAX8':23,
            'CAX8R':23,
                'B21':3,
                'B31':3,
            'B31R':3,
            'T2D2':3,
            'T3D2':3,
            'GAPUNI':3,
        'DASHPOTA':3,
            'SPRING2':3,
            'SPRINGA':3,
                'B32':21,
            'B32R':21,
            'T3D3':21,
                'D':21,
            'SPRING1':1,
            'DCOUP3D':1,
            'MASS':1,
        }
    if frd_elem_type in frd2vtk_num:
        return frd2vtk_num[frd_elem_type]
    else:
        if frd_elem_type in frd2vtk_txt:
            return frd2vtk_txt[frd_elem_type]
        else:
            return 0
