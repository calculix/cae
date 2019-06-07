# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, UJV Rez, May 2019.  
    Distributed under GNU General Public License, version 2.

    Parses finite element mesh in the Abaqus, Gmsh or CalculiX .inp-file.
    Tested on 2D quadrilateral and triangular first order elements.
    Reads nodes coordinates, elements composition, node and element sets, surfaces.
    Calculates elements cendroid coordinates.
    Generates triangles or quadrangles list to use with matplotlib.
    'project_field_on_centroids' method interpolates node field to elements centroids.
"""

import numpy as np
import matplotlib.tri as tri

# Mesh object, contains methods for .inp-file parsing
class Mesh:

    # All mesh nodes with coordinates
    """
        1: ( 0.0, -1742.5, 0.0),
        2: (74.8, -1663.7, 0.0),
        ...
    """
    nodes = {}

    # All mesh elements composition
    """
        1: (1, 2),
        2: (3, 4),
        ...
        11: (21, 22, 23),
        12: (24, 25, 26),
        ...
    """
    elements = {}
    
    # Element types
    """
        1: 'type1',
        2: 'type2',
        ...
    """
    types = {}

    # Coordinates of all elements centroids
    """
        1: ( 0.0, -1742.5, 0.0),
        2: (74.8, -1663.7, 0.0),
        ...
    """
    centroids = {}
    
    # Node sets
    """
        'nset1': [1, 2, 3, 4],
        'nset2': [5, 6, 7, 8],
        ...
    """
    nsets = {}

    # Element sets
    """
        'eset1': [1, 2, 3, 4],
        'eset2': [5, 6, 7, 8],
        ...
    """
    esets = {}

    # Surface names
    """
        'surf1', 'surf2', 'surf3', 
    """
    surfaces = ()

    # Additional mesh variables
    cx = []; cy = [] # centroid coordinates as numpy array
    nx = []; ny = [] # nodes coordinates as numpy array
    triangles = () # triangles list to use with matplotlib
    quadrangles = [] # quadrangles to use with matplotlib

    # Some parameters
    initialized = False


    # Parse nodes with coordinates
    # *NODE keyword
    def get_nodes(self, lines):
        for i in range(len(lines)):
            if lines[i].startswith('*NODE'):
                while i+1<len(lines) and not lines[i+1].startswith('*'): # read the whole block and return
                    a = lines[i+1].split(',')
                    num = int(a[0].strip()) # node number
                    self.nodes[num] = () # tuple with node coordinates
                    for coord in a[1:]:
                        self.nodes[num] += (float(coord.strip()), ) # add coordinate to tuple
                    i += 1
                return

    # Parse node sets
    # *NSET keyword
    def get_nsets(self, lines):
        for i in range(len(lines)):
            if lines[i].startswith('*NSET'):
                name = lines[i].split('=')[1]
                self.nsets[name] = ()
                while i+1<len(lines) and not lines[i+1].startswith('*'):
                    a = lines[i+1].split(',')
                    for n in a:
                        if len(n.strip()):
                            self.nsets[name] += (int(n), )
                    i += 1

    # Parse elements composition and calculate centroid
    # *ELEMENT keyword
    def get_elements(self, lines):
        for i in range(len(lines)):
            if lines[i].startswith('*ELEMENT'):
                # etype = lines[i].split('=')[1].split(',')[0] # element type
                etype = lines[i].upper().split('TYPE=')[1].split(',')[0] # element type
                while i+1<len(lines) and not lines[i+1].startswith('*'): # there will be no comments
                    a = lines[i+1].split(', ')
                    num = int(a[0]) # element number
                    self.types[num] = etype # save element type
                    self.elements[num] = () # tuple with element nodes
                    # TODO element nodes could be splitted into few lines
                    for n in a[1:]:
                        self.elements[num] += (int(n), ) # add node to tuple
                    x=0; y=0; z=0
                    for n in a[1:]: # iterate over element's node numbers
                        x += self.nodes[int(n)][0] # sum up x-coordinates of all nodes of the element
                        y += self.nodes[int(n)][1] # sum up y-coordinates of all nodes of the element
                        try: # 3D case
                            z += self.nodes[int(n)][2] # sum up z-coordinates of all nodes of the element
                        except:
                            pass
                    amount = len(a[1:]) # amount of nodes in element
                    x /= amount; y /= amount; z /= amount
                    self.centroids[num] = (x, y, z) # centroid coordinates 3D
                    i += 1

    # Parse element sets
    # *ELSET keyword
    def get_esets(self, lines):
        for i in range(len(lines)):
            if lines[i].startswith('*ELSET'):
                name = lines[i].split('=')[1]
                self.esets[name] = ()
                while i+1<len(lines) and not lines[i+1].startswith('*'):
                    a = lines[i+1].split(',')
                    for e in a:
                        try:
                            self.esets[name] += (int(e.strip()), )
                        except:
                            pass
                    i += 1

    # Parse surfaces
    # *SURFACE keyword
    def get_surfaces(self, lines):
        for line in lines:
            if line.startswith('*SURFACE'):
                name = line.upper().split('NAME=')[1].split(',')[0]
                self.surfaces += (name, )

    # Set additional variables
    def set_additional_vars(self):
        self.cx = np.array( [v[0] for k,v in sorted(self.centroids.items())] ) # centroids x-coords sorted by element number
        self.cy = np.array( [v[1] for k,v in sorted(self.centroids.items())] ) # centroids y-coords sorted by element number
        self.nx = np.array( [v[0] for k,v in sorted(self.nodes.items())] ) # list of x-coords sorted by node number
        self.ny = np.array( [v[1] for k,v in sorted(self.nodes.items())] ) # list of y-coords sorted by node number
        first_node_num = list(self.nodes.keys())[0]
        for elem, nodes in sorted(self.elements.items()): # tuples of node numbers sorted by element number
            if len(nodes)==3:
                # Triangles consist of nodes indexes (not
                # numbers), so we may subtract first_node_num
                self.triangles += ((
                                    nodes[0] - first_node_num,
                                    nodes[1] - first_node_num,
                                    nodes[2] - first_node_num), )
            if len(nodes)==4:
                quad = [] # one quadrangle - array of node coordinates couples
                for n in nodes:
                    x = self.nodes[n][0]; y = self.nodes[n][1]
                    coords = np.array([[x, y]]) # coords of one node of the quadrangle
                    if len(quad):
                        quad = np.append(quad, coords, axis=0)
                    else:
                        quad = coords
                if len(self.quadrangles):
                    self.quadrangles = np.append(self.quadrangles, [quad], axis=0)
                else:
                    self.quadrangles = [quad]

    # Initialization
    def __init__(self, inp_file):
        # Open and read all the .inp-file
        lines = []
        with open(inp_file, 'r') as f:
            for i, line in enumerate(f):
                if not '**' in line: # skip comments
                    lines.append(line.strip().upper())
        self.get_nodes(lines) # parse nodes
        self.get_nsets(lines) # parse node sets
        self.get_elements(lines) # parse elements
        self.get_esets(lines) # parse node sets
        self.get_surfaces(lines) # parse surfaces
        print('Total:')
        print('\t{0} nodes'.format(len(self.nodes)))
        print('\t{0} elements'.format(len(self.elements)))
        print('\t{0} centroids'.format(len(self.centroids)))
        print('\t{0} nsets'.format(len(self.nsets)))
        print('\t{0} esets'.format(len(self.esets)))
        self.set_additional_vars()
        self.initialized = True


    # Project field on mesh centroids
    def project_field_on_centroids(self, fx, fy, field_values):
        """
            fx - column with x-coordinates of field points
            fy - column with y-coordinates of field points
            field_values - column with field values in points (fx, fy)
            fx, fy, field_values should be the same length
        """
        triang = tri.Triangulation(fx, fy) # Delaunay triangulation from field points
        interp = tri.LinearTriInterpolator(triang, field_values) # interpolation object
        res = interp(self.cx, self.cy) # perform linear interpolation on centroids
        res = [x if type(x)==np.float64 else np.float64(0) for x in res] # zero values outside field triangles
        res = res / max(res) * max(field_values) # diminish interpolation error
        return res

    # Convert Calculix element type to VTK
    @staticmethod
    def convert_elem_type(frd_elem_type):
        """
            Keep in mind that CalculiX expands shell elements
            In vtk elements nodes are numbered starting from 0, not 1

            For frd see http://www.dhondt.de/cgx_2.15.pdf pages 117-123 (chapter 10)
            For vtk see https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf pages 9-10

             _________________________________________________________________
            |                               |                                 |
            | №№      CalculiX type         |   №№        VTK type            |
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
