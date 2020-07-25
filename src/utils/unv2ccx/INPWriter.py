#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
© Joël Cugnoni, September 2006 - original code, www.caelinux.com
© Ihor Mirzov, August 2019 - refactoring
Distributed under GNU General Public License v3.0

Writes FEM nodes, elements and groups
(node and element sets) into INP file. """

import os
import logging
import re

PADDING = ' '*4 # four spaces

# Main function
def write(FEM, file_name):
    with open(file_name, 'w') as f:

        # Nodes
        f.write('*NODE, NSET=NALL\n')
        for node in FEM.nodes:
            f.write(PADDING + '{:d}, {:.10e}, {:.10e}, {:.10e}\n'\
                .format(node.num, node.coords[0], node.coords[1], node.coords[2]))

        # Split element list by types
        elements = {}
        for e in FEM.elements:
            ccx_element_type = convert_element(e.type) # can return None
            if ccx_element_type:
                if not ccx_element_type in elements:
                    elements[ccx_element_type] = [e, ] # initialize
                else:
                    elements[ccx_element_type].append(e) # append to initialized

        # Write element groups by type
        for etype in elements.keys():
            f.write('*ELEMENT, TYPE={0}, ELSET={0}\n'.format(etype))

            # Select the appropriate map between the nodes
            themap = element_connectivity(etype)

            # Write elements connectivity
            for e in elements[etype]:
                f.write(PADDING + '{:d}, '.format(e.num))
                for i in range(len(e.nodes)):
                    if i > 0 and i % 10 == 0:
                        f.write('\n' + PADDING)
                    node_number = themap[i] - 1
                    f.write('{:d}, '.format(e.nodes[node_number]))
                f.write('\n')

        # Write node sets
        for group in FEM.nsets:
            f.write('*NSET, NSET=' + group.name)
            writeGroup(f, group)

        # Write element sets
        for group in FEM.esets:
            f.write('*ELSET, ELSET=' + group.name)
            writeGroup(f, group)

# Write node or element set
def writeGroup(f, group):
    for i in range(group.nitems):
        if i % 8 == 0:
            f.write('\n' + PADDING)
        f.write('{:d}, '.format(group.items[i]))
    f.write('\n')

# Convert UNV element type to CalculiX
def convert_element(unv_element_type):

    # unv_element_type : ccx_element_type
    dic = {
        11:'B31',       # Rod
        21:'B31',       # Linear beam
        22:'B32',       # Tapered beam
        23:'B32',       # Curved beam
        24:'B32',       # Parabolic beam
        31:'B31',       # Straight pipe
        32:'B32',       # Curved pipe
        41:'CPS3',      # Plane Stress Linear Triangle
        42:'CPS6',      # Plane Stress Parabolic Triangle
        43: None,       # Plane Stress Cubic Triangle
        44:'CPS4',      # Plane Stress Linear Quadrilateral
        45:'CPS8',      # Plane Stress Parabolic Quadrilateral
        46: None,       # Plane Strain Cubic Quadrilateral
        51:'CPE3',      # Plane Strain Linear Triangle
        52:'CPE6',      # Plane Strain Parabolic Triangle
        53: None,       # Plane Strain Cubic Triangle
        54:'CPE4',      # Plane Strain Linear Quadrilateral
        55:'CPE8',      # Plane Strain Parabolic Quadrilateral
        56: None,       # Plane Strain Cubic Quadrilateral
        61:'M3D3',      # Plate Linear Triangle
        62:'M3D6',      # Plate Parabolic Triangle
        63: None,       # Plate Cubic Triangle
        64:'M3D4',      # Plate Linear Quadrilateral
        65:'M3D8',      # Plate Parabolic Quadrilateral
        66: None,       # Plate Cubic Quadrilateral
        71:'M3D4',      # Membrane Linear Quadrilateral
        72:'M3D6',      # Membrane Parabolic Triangle
        73: None,       # Membrane Cubic Triangle
        74:'M3D3',      # Membrane Linear Triangle
        75:'M3D8',      # Membrane Parabolic Quadrilateral
        76: None,       # Membrane Cubic Quadrilateral
        81:'CAX3',      # Axisymetric Solid Linear Triangle
        82:'CAX6',      # Axisymetric Solid Parabolic Triangle
        84:'CAX4',      # Axisymetric Solid Linear Quadrilateral
        85:'CAX8',      # Axisymetric Solid Parabolic Quadrilateral
        91:'S3',        # Thin Shell Linear Triangle
        92:'S6',        # Thin Shell Parabolic Triangle
        93: None,       # Thin Shell Cubic Triangle
        94:'S4',        # Thin Shell Linear Quadrilateral
        95:'S8',        # Thin Shell Parabolic Quadrilateral
        96: None,       # Thin Shell Cubic Quadrilateral
        101:'C3D6',     # Thick Shell Linear Wedge
        102:'C3D15',    # Thick Shell Parabolic Wedge
        103: None,      # Thick Shell Cubic Wedge
        104:'C3D8',     # Thick Shell Linear Brick
        105:'C3D20',    # Thick Shell Parabolic Brick
        106: None,      # Thick Shell Cubic Brick
        111:'C3D4',     # Solid Linear Tetrahedron
        112:'C3D6',     # Solid Linear Wedge
        113:'C3D15',    # Solid Parabolic Wedge
        114: None,      # Solid Cubic Wedge
        115:'C3D8',     # Solid Linear Brick
        116:'C3D20',    # Solid Parabolic Brick
        117: None,      # Solid Cubic Brick
        118:'C3D10',    # Solid Parabolic Tetrahedron
        121: None,      # Rigid Bar
        122: None,      # Rigid Element
        136:'SPRINGA',  # Node To Node Translational Spring
        137:'SPRINGA',  # Node To Node Rotational Spring
        138:'SPRINGA',  # Node To Ground Translational Spring
        139:'SPRINGA',  # Node To Ground Rotational Spring
        141:'DASHPOTA', # Node To Node Damper
        142:'DASHPOTA', # Node To Gound Damper
        151:'GAPUNI',   # Node To Node Gap
        152:'GAPUNI',   # Node To Ground Gap
        161:'MASS',     # Lumped Mass
        171: None,      # Axisymetric Linear Shell
        172: None,      # Axisymetric Parabolic Shell
        181: None,      # Constraint
        191: None,      # Plastic Cold Runner
        192: None,      # Plastic Hot Runner
        193: None,      # Plastic Water Line
        194: None,      # Plastic Fountain
        195: None,      # Plastic Baffle
        196: None,      # Plastic Rod Heater
        201: None,      # Linear node-to-node interface
        202: None,      # Linear edge-to-edge interface
        203: None,      # Parabolic edge-to-edge interface
        204: None,      # Linear face-to-face interface
        208: None,      # Parabolic face-to-face interface
        212: None,      # Linear axisymmetric interface
        213: None,      # Parabolic axisymmetric interface
        221: None,      # Linear rigid surface
        222: None,      # Parabolic rigin surface
        231: None,      # Axisymetric linear rigid surface
        232: None,      # Axisymentric parabolic rigid surface
        }
    if unv_element_type in dic:
        return dic[unv_element_type]
    else:
        return None

# Map of the nodes between Universal and Calculix elements
def element_connectivity(ccx_element_type):

    # Some elements have the same connectivity
    etype = re.sub('CPS|CPE|M3D|CAX', 'S', ccx_element_type)
    if ccx_element_type != etype:
        logging.debug('Using {} connectivity for {}.'.format(etype, ccx_element_type))

    dic = {
        'S6': [1, 3, 5, 2, 4, 6],
        'S8': [1, 3, 5, 7, 2, 4, 6, 8],
        'C3D10': [1, 3, 5, 10, 2, 4, 6, 7, 8, 9],
        'C3D15': [1, 3, 5, 10, 12, 14, 2, 4, 6, 11, 13, 15, 7, 8, 9],
        'C3D20': [1, 3, 5, 7, 13, 15, 17, 19, 2, 4, 6, 8, 14, 16, 18, 20, 9, 10, 11, 12],
        }

    if etype in dic:
        return dic[etype]

    # If no map is available translate as it is
    else:
        return list(range(1,21))
