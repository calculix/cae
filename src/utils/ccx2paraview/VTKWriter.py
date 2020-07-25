#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, January 2020
Distributed under GNU General Public License v3.0

Inspired by C# converter written by Maciek Hawryłkiewicz in 2015.

About the format read:
https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf

Remember that the frd file is node based, so element results are also
stored at the nodes after extrapolation from the integration points:
http://www.dhondt.de/ccx_2.15.pdf """


import logging
import math

import frd2vtk


# Write element connectivity with renumbered nodes
def write_element_connectivity(renumbered_nodes, e, f):
    # frd: 20 node brick element
    if e.type == 4:
        element_string = '20 '
        # Last eight nodes have to be repositioned
        r1 = tuple(range(12)) # 8,9,10,11
        r2 = tuple(range(12, 16)) # 12,13,14,15
        r3 = tuple(range(16, 20)) # 16,17,18,19
        node_num_list = r1 + r3 + r2
        for i in node_num_list:
            node = renumbered_nodes[e.nodes[i]] # node after renumbering
            element_string += '{:d} '.format(node)

    # frd: 15 node penta element
    elif e.type==5 or e.type==2:
        """
            CalculiX elements type 5 are not supported in VTK and
            has to be processed as CalculiX type 2 (6 node wedge,
            VTK type 13). Additional nodes are omitted.
        """
        element_string = '6 '
        for i in [0,2,1,3,5,4]: # repositioning nodes
            node = renumbered_nodes[e.nodes[i]] # node after renumbering
            element_string += '{:d} '.format(node)

    # All other elements
    else:
        n = len(e.nodes)
        element_string = '{} '.format(n)
        for i in range(n):
            node = renumbered_nodes[e.nodes[i]] # node after renumbering
            element_string += '{:d} '.format(node)

    f.write('\t' + element_string + '\n')


# Write data
def write_data(f, b, numnod):
    f.write('FIELD {} 1\n'.format(b.name))
    f.write('\t{} {} {} double\n'.format(b.name, len(b.components), numnod))
    nodes = sorted(b.results.keys())

    # Some warnings repeat too much time - mark them
    emitted_warning_types = {'Inf':0, 'NaN':0}

    for n in range(numnod): # iterate over nodes
        node = nodes[n]
        data = b.results[node]
        f.write('\t')
        for d in data:
            # # Filter small values for smooth zero fields
            # if abs(d) < 1e-9: d = 0
            if math.isinf(d):
                d = 0.0
                emitted_warning_types['Inf'] += 1
            if math.isnan(d):
                d = 0.0
                emitted_warning_types['NaN'] += 1
            f.write('\t{: .8E}'.format(d))
        f.write('\n')

    for k, v in emitted_warning_types.items():
        if v > 0:
            logging.warning('{} {} values are converted to 0.0'.format(v, k))


# Main function
def writeVTK(p, file_name, time): # p is FRDParser object

    with open(file_name, 'w') as f:
        # Header
        f.write('# vtk DataFile Version 3.0\n\n')
        f.write('ASCII\n')
        f.write('DATASET UNSTRUCTURED_GRID\n\n')

        # POINTS section - coordinates of all nodes
        f.write('POINTS ' + str(p.node_block.numnod) + ' double\n')
        new_node_number = 0 # node numbers should start from 0
        renumbered_nodes = {} # old_number : new_number
        for n in p.node_block.nodes.keys():

            # Write nodes coordinates
            coordinates = ''.join('\t{: .8E}'.format(coord) \
                for coord in p.node_block.nodes[n].coords)
            f.write(coordinates + '\n')

            # For VTK nodes should go consequently starting from 0
            renumbered_nodes[n] = new_node_number
            new_node_number += 1

            if new_node_number == p.node_block.numnod:
                break

        f.write('\n')

        # CELLS section - elements connectyvity
        totn = 0 # total number of nodes
        for e in p.elem_block.elements:
            if e.type == 5: totn += 6
            else: totn += len(e.nodes)
        f.write('CELLS {} {}\n'.format(p.elem_block.numelem, p.elem_block.numelem + totn)) # number of cells and size of the cell list
        for e in p.elem_block.elements:
            write_element_connectivity(renumbered_nodes, e, f)
        f.write('\n')

        # CELL TYPES section - write element types:
        f.write('CELL_TYPES {}\n'.format(p.elem_block.numelem))
        for e in p.elem_block.elements:
            vtk_elem_type = frd2vtk.convert_elem_type(e.type)
            f.write('\t{}\n'.format(vtk_elem_type))
        f.write('\n')

        # POINT DATA - from here start all the results
        f.write('POINT_DATA {}\n'.format(p.node_block.numnod))
        for b in p.result_blocks: # iterate over NodalResultsBlock
            if b.value != time: # write results for one time increment only
                continue
            if len(b.results) and len(b.components):
                if b.value < 1:
                    time_str = 'time {:.2e}, '.format(b.value)
                else:
                    time_str = 'time {:.1f}, '.format(b.value)
                logging.info('Step {}, '.format(b.numstep) +\
                            time_str +\
                            '{}, '.format(b.name) +\
                            '{} components, '.format(len(b.components)) +\
                            '{} values'.format(len(b.results)))
                write_data(f, b, p.node_block.numnod)
            else:
                logging.warning(b.name, '- no data for this increment')


"""
    TODO learn and use it for future code improvement:
    https://vtk.org/doc/nightly/html/classvtkUnstructuredGridWriter.html
    https://vtk.org/doc/nightly/html/c2_vtk_e_5.html#c2_vtk_e_vtkUnstructuredGrid
    https://vtk.org/gitweb?p=VTK.git;a=blob;f=Examples/DataManipulation/Python/pointToCellData.py

    TODO Use it to write mesh
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(file_name)
    writer.SetInputData(unstructured_grid)
    writer.Write()
"""