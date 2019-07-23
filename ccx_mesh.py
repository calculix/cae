# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, 2019.  
    Distributed under GNU General Public License, version 2.

    Parses finite element mesh from the CalculiX .inp-file.
    Reads nodes coordinates, elements composition, node and element sets, surfaces.
    It's case insensitive and translates all text uppercase.
"""


import ccx_log, os, re


class Parse:

    # Parse nodes with coordinates
    # *NODE keyword
    def get_nodes(self, lines):
        for i in range(len(lines)): # lines are uppercase

            # Distinguish 'NODE' and 'NODE PRINT'
            if ',' in lines[i]:
                keyword_name = lines[i].split(',')[0]
            else:
                keyword_name = lines[i]

            if keyword_name == '*NODE':

                # Node set with all nodes
                # if 'NSET' in lines[i]:
                match = re.search('NSET\s*=\s*\w*', lines[i])
                try:
                    name = match.group(0).split('=')[1].strip()
                except:
                    name = 'ALL'
                self.nsets[name] = ()

                while i+1<len(lines) and not lines[i+1].startswith('*'): # read the whole block and return
                    lines[i+1] = lines[i+1].replace(',', ' ') # to avoid redundant commas in the end of line
                    a = lines[i+1].split()
                    if len(a) == 3: # in 2D case add Z coord equal to zero
                        a.append(0)
                    num = int(a[0]) # node number
                    self.nsets[name] += (num, )
                    self.nodes[num] = () # tuple with node coordinates
                    for j,coord in enumerate(a[1:]):
                        coord = float(coord)
                        self.nodes[num] += (coord, ) # add coordinate to tuple

                        # Bounding box
                        if coord < self.bounds[j*2]:
                            self.bounds[j*2] = coord # update min coords values
                        if coord > self.bounds[j*2+1]:
                            self.bounds[j*2+1] = coord # update max coords values
                    i += 1
                    if not len(self.nodes[num]):
                        msg_text = 'Node {} has no coordinates and will be removed.'.format(num)
                        msg = ccx_log.msg(ccx_log.msgType.WARNING, msg_text)
                        self.msg_list.append(msg)
                        del self.nodes[num]
                # do not return to parse few *NODE sections


    # Parse node sets
    # *NSET keyword
    def get_nsets(self, lines):
        for i in range(len(lines)): # lines are uppercase
            if lines[i].startswith('*NSET'):

                match = re.search('NSET\s*=\s*\w*', lines[i])
                name = match.group(0).split('=')[1].strip()
                self.nsets[name] = ()

                if not 'GENERATE' in lines[i]:
                    while i+1<len(lines) and not lines[i+1].startswith('*'):
                        a = lines[i+1].split(',')
                        for n in a:
                            n = n.strip()
                            if len(n):
                                try:
                                    self.nsets[name] += (int(n), )
                                except ValueError as err:
                                    for node in self.nsets[n]:
                                        self.nsets[name] += (node, )
                        i += 1
                else:
                    try:
                        start, stop, step = re.split(',\s*', lines[i+1].strip())
                    except:
                        start, stop = re.split(',\s*', lines[i+1].strip())
                        step = 1
                    self.nsets[name] += tuple(range(int(start), int(stop)+1, int(step)))


    # Parse elements composition and calculate centroid
    # *ELEMENT keyword
    def get_elements(self, lines):
        for i in range(len(lines)): # lines are uppercase
            if lines[i].startswith('*ELEMENT'):

                # Element set with all elements
                match = re.search('ELSET\s*=\s*\w*', lines[i])
                try:
                    name = match.group(0).split('=')[1].strip()
                except:
                    name = 'ALL'
                self.esets[name] = ()

                match = re.search('TYPE\s*=\s*\w*', lines[i])
                etype = match.group(0).split('=')[1].strip()
                amount = self.get_amount_of_nodes(etype)

                while i+1<len(lines) and not lines[i+1].startswith('*'): # there will be no comments

                    # Element nodes could be splitted into 2 lines
                    a = lines[i+1].replace(',', ' ').split()
                    if len(a) < amount + 1: # +1 for element number
                        a.extend( lines[i+2].replace(',', ' ').split() )
                        i += 1

                    num = int(a[0].strip()) # element number
                    self.esets[name] += (num, )
                    self.types[num] = etype.strip() # save element type
                    self.elements[num] = () # tuple with element nodes
                    for n in a[1:]:
                        self.elements[num] += (int(n.strip()), ) # add node to tuple
                    x=0; y=0; z=0
                    for n in a[1:]: # iterate over element's node numbers
                        node = int(n.strip()) # node number
                        if node == 0: # it is possible in network element, type=D
                            node = int(a[2].strip()) # take middle node to display in VTK
                        if node in self.nodes:
                            x += self.nodes[node][0] # sum up x-coordinates of all nodes of the element
                            y += self.nodes[node][1] # sum up y-coordinates of all nodes of the element
                            try: # 3D case
                                z += self.nodes[node][2] # sum up z-coordinates of all nodes of the element
                            except:
                                pass
                        else:
                            msg_text = 'Theree is no node {} in element {}. '.format(n.strip(), num) + \
                                        'This element will be removed.'
                            msg = ccx_log.msg(ccx_log.msgType.WARNING, msg_text)
                            self.msg_list.append(msg)
                            del self.elements[num]
                    amount = len(a[1:]) # amount of nodes in element
                    x /= amount; y /= amount; z /= amount
                    self.centroids[num] = (x, y, z) # centroid coordinates 3D
                    i += 1
                # do not return to parse few *ELEMENT sections


    # Parse element sets
    # *ELSET keyword
    def get_esets(self, lines):
        for i in range(len(lines)): # lines are uppercase
            if lines[i].startswith('*ELSET'):

                match = re.search('ELSET\s*=\s*\w*', lines[i])
                name = match.group(0).split('=')[1].strip()
                self.esets[name] = ()

                if not 'GENERATE' in lines[i]:
                    while i+1<len(lines) and not lines[i+1].startswith('*'):
                        a = lines[i+1].split(',')
                        for e in a:
                            e = e.strip()
                            if len(e):
                                try:
                                    self.esets[name] += (int(e), )
                                except ValueError as err:
                                    for element in self.esets[e]:
                                        self.esets[name] += (element, )
                        i += 1
                else:
                    try:
                        start, stop, step = re.split(',\s*', lines[i+1].strip())
                    except:
                        start, stop = re.split(',\s*', lines[i+1].strip())
                        step = 1
                    self.esets[name] += tuple(range(int(start), int(stop)+1, int(step)))


    # Parse surfaces
    # *SURFACE keyword
    def get_surfaces(self, lines):
        for line in lines: # lines are uppercase
            if line.startswith('*SURFACE') and \
                    not '*SURFACE INTERACTION' in line and \
                    not '*SURFACE BEHAVIOR' in line:
                match = re.search('NAME\s*=\s*\w*', line)
                name = match.group(0).split('=')[1].strip()
                self.surfaces += (name, )


    # Recurcively read all the lines of the file and its includes
    def get_lines(self, inp_file):
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
                            lines.extend(self.get_lines(path + '/' + inp_file2))
        except:
            msg_text = 'There is no file {}.'.format(inp_file)
            msg = ccx_log.msg(ccx_log.msgType.ERROR, msg_text)
            self.msg_list.append(msg)
        return lines


    # Initialization
    def __init__(self, inp_file):

        self.msg_list = [] # list of messages for logger

        # All mesh nodes with coordinates
        """
            1: ( 0.0, -1742.5, 0.0),
            2: (74.8, -1663.7, 0.0),
            ...
        """
        self.nodes = {}

        # All mesh elements composition
        """
            1: (1, 2),
            2: (3, 4),
            ...
            11: (21, 22, 23),
            12: (24, 25, 26),
            ...
        """
        self.elements = {}
        
        # Element types
        """
            1: 'type1',
            2: 'type2',
            ...
        """
        self.types = {}

        # Coordinates of all elements centroids
        """
            1: ( 0.0, -1742.5, 0.0),
            2: (74.8, -1663.7, 0.0),
            ...
        """
        self.centroids = {}
        
        # Node sets
        """
            'nset1': [1, 2, 3, 4],
            'nset2': [5, 6, 7, 8],
            ...
        """
        self.nsets = {}

        # Element sets
        """
            'eset1': [1, 2, 3, 4],
            'eset2': [5, 6, 7, 8],
            ...
        """
        self.esets = {}

        # Surface names
        """
            'surf1', 'surf2', 'surf3', 
        """
        self.surfaces = ()

        # Mesh bounds to avoid camera flying to infinity
        self.bounds = [1e+6,-1e+6]*3 # Xmin,Xmax, Ymin,Ymax, Zmin,Zmax

        # Open and read whole the .inp-file
        lines = self.get_lines(inp_file)

        # Parse nodes
        self.get_nodes(lines)
        try:

            msg_text = '{} nodes'.format(len(self.nodes))
            msg = ccx_log.msg(ccx_log.msgType.INFO, msg_text)
            self.msg_list.append(msg)
        except:
            msg_text = 'Can\'t parse nodes'
            msg = ccx_log.msg(ccx_log.msgType.ERROR, msg_text)
            self.msg_list.append(msg)

        # Parse node sets
        try:
            self.get_nsets(lines)

            msg_text = '{} nsets'.format(len(self.nsets))
            # for k,v in self.nsets.items():
            #     msg_text += '<br/>\n{0}: {1}'.format(k, v)
            msg = ccx_log.msg(ccx_log.msgType.INFO, msg_text)
            self.msg_list.append(msg)
        except:
            msg_text = 'Can\'t parse nsets'
            msg = ccx_log.msg(ccx_log.msgType.ERROR, msg_text)
            self.msg_list.append(msg)

        # Parse elements
        try:
            self.get_elements(lines)

            msg_text = '{} elements'.format(len(self.elements))
            msg = ccx_log.msg(ccx_log.msgType.INFO, msg_text)
            self.msg_list.append(msg)
        except:
            msg_text = 'Can\'t parse elements'
            msg = ccx_log.msg(ccx_log.msgType.ERROR, msg_text)
            self.msg_list.append(msg)

        # Parse element sets
        try:
            self.get_esets(lines)

            msg_text = '{} esets'.format(len(self.esets))
            # for k,v in self.esets.items():
            #     msg_text += '<br/>\n{0}: {1}'.format(k, v)
            msg = ccx_log.msg(ccx_log.msgType.INFO, msg_text)
            self.msg_list.append(msg)
        except:
            msg_text = 'Can\'t parse esets'
            msg = ccx_log.msg(ccx_log.msgType.ERROR, msg_text)
            self.msg_list.append(msg)

        # Parse surfaces
        try:
            self.get_surfaces(lines)

            msg_text = '{} surfaces'.format(len(self.surfaces))
            msg = ccx_log.msg(ccx_log.msgType.INFO, msg_text)
            self.msg_list.append(msg)
        except:
            msg_text = 'Can\'t parse surfaces'
            msg = ccx_log.msg(ccx_log.msgType.ERROR, msg_text)
            self.msg_list.append(msg)


    # Convert Calculix element type to VTK
    @staticmethod
    def convert_elem_type(frd_elem_type):
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


    # Get amount of nodex by element type
    def get_amount_of_nodes(self, etype):
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
