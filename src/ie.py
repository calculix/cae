#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, December 2019
    Distributed under GNU General Public License v3.0

    INP importer:
        Enrich KOM with implementations from parsed file.
        Generate new tree with keyword implementations.
        Parse mesh and build ugrid
        Display ugrid in the VTK widget

    INP exporter:
        Recursively write implementation's INP_code to output .inp-file
"""


from PyQt5.QtWidgets import QFileDialog
import os, logging
from model.kom import item_type, implementation, KOM
from model.parsers.mesh import Mesh
import file_tools


# Menu File -> Import
"""
s - Settings
w - Window
m - Model
t - Tree
j - Job
"""
def importFile(s, w, m, t, j, file_name=None):
    clear_textEdit = False
    if not file_name:
        clear_textEdit = True
        file_name = QFileDialog.getOpenFileName(None, \
            'Import INP/UNV file', j.dir, \
            'INP (*.inp);;UNV (*.unv)')[0]

    if file_name:

        # Clear log window before new import
        if clear_textEdit:
            for h in logging.getLogger().handlers:
                h.flush()

        # Clear selection before import new model
        if s.show_vtk:
            w.VTK.actionSelectionClear()

        # Generate new KOM without implementations
        m.KOM = KOM()

        # Rename job before tree regeneration
        j.rename(file_name[:-4] + '.inp')

        # Convert UNV to INP
        if file_name.lower().endswith('.unv'):
            j.convertUNV()
            if not os.path.isfile(j.inp):
                logging.error('Error converting ' + j.unv)
                return

        # Show model name in window's title
        w.setWindowTitle('CalculiX CAE - ' + j.name)


        # Parse INP and enrich KOM with parsed objects
        logging.info('Loading ' + j.inp + '.')
        lines = file_tools.readLines(j.inp)
        def importer(INP_doc, KOM):
            keyword_chain = []
            impl_counter = {}
            for i in range(len(INP_doc)):
                line = INP_doc[i]

                # Parse keyword
                if line.startswith('*'):

                    # Distinguish 'NODE' and 'NODE PRINT'
                    if ',' in line:
                        keyword_name = line.split(',')[0]
                    else:
                        keyword_name = line

                    # Find KOM keyword path corresponding to keyword_chain
                    keyword_chain.append(keyword_name)
                    path = KOM.getPath(keyword_chain)
                    if path:
                        logging.debug('path found: ' + str([item.name for item in path]))

                        # Read INP_code for the current keyword
                        INP_code = [line] # line is stripped in mesh.py
                        while i+1 < len(INP_doc) and \
                            not INP_doc[i+1].startswith('*'): # here will be no comments - they are removed in mesh.py
                            INP_code.append(INP_doc[i+1])
                            i += 1

                        # Create keyword implementations
                        impl = None
                        path_as_string = '' # string representation of 'path' accounting for implementations
                        for j in range(len(path)):
                            # Choose where to create implementation
                            if impl:
                                # Implementation will be created inside another implementation
                                item = impl.getItemByName(path[j].name)
                            else:
                                # Implementation will be created inside keyword or group
                                item = path[j]

                            path_as_string += '/' + item.name
                            if j == len(path) - 1: # last item is always keyword
                                # Create implementation (for example, MATERIAL-1)
                                impl = implementation(item, INP_code)
                                logging.debug('1')
                            elif item.item_type == item_type.KEYWORD:
                                # If for this keyword implementation was created previously
                                counter = impl_counter[path_as_string] - 1
                                impl = item.items[counter] # first implementation, for example, STEP-1
                                path_as_string += '/' + impl.name
                                logging.debug('2')
                            else:
                                impl = item
                                logging.debug('3')

                        # Count implementation
                        if path_as_string in impl_counter:
                            # If current keyword already has implementations
                            impl_counter[path_as_string] += 1
                        else:
                            # If first implementation was created for current keyword
                            impl_counter[path_as_string] = 1

                    else:
                        logging.warning('Wrong keyword {}.'.format(keyword_name))
        importer(lines, m.KOM) # pass whole INP-file to the parser


        # Add parsed implementations to the tree
        t.generateTreeView(m)

        # Parse mesh and plot it
        m.Mesh = Mesh(INP_file=j.inp) # parse mesh
        if s.show_vtk:
            w.VTK.plotMesh(m.Mesh)
        
        return True
    else:
        return False


# Menu File -> Write INP file
def writeInput(m, j, file_name=None):
    if not file_name:
        file_name = QFileDialog.getSaveFileName(None, \
            'Write INP file', j.dir, \
            'Input files (*.inp)')[0]
    if file_name:

        # Recursively write implementation's INP_code to output .inp-file
        def writer(parent, f, level):

            # Level is used for padding
            if parent.item_type == item_type.IMPLEMENTATION:
                level += 1

            # For each group/keyword from KOM
            for item in parent.items:

                if item.item_type == item_type.ARGUMENT:
                    continue

                if item.item_type == item_type.IMPLEMENTATION:
                    # INP_code is stripped
                    f.write(' '*4*level + item.INP_code[0] + '\n')
                    for line in item.INP_code[1:]:
                        f.write(' '*4*(level+1) + line + '\n')

                # Continue call iterator until dig to implementation
                writer(item, f, level)

        with open(file_name, 'w') as f:
            writer(m.KOM.root, f, 0)
        logging.info('Input written to ' + file_name)

        # Rename job and update treeView if new INP-file name was selected
        j.rename(file_name)
