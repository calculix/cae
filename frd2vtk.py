# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, August 2019.
    Distributed under GNU General Public License, version 2.

    Convert Calculix element type to VTK
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


def convert_elem_type(frd_elem_type):

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
