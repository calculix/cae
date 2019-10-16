# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, September 2019
    Distributed under GNU General Public License v3.0

    INP importer:
        Enrich KOM with implementations from parsed file.
        Generate new tree with keyword implementations.
        Parse mesh and build ugrid
        Display ugrid in the VTK widget

    INP exporter:
        Recursively write implementation's INP_code to output .inp-file

    Depends on cae.model.
"""


from PyQt5.QtWidgets import QFileDialog
import os, logging
from kom import item_type, implementation, KOM
from mesh import Parse, readLines


# Menu File -> Import
def importFile(model, file_name=None):
    clear_textEdit = False
    if not file_name:
        clear_textEdit = True
        file_name = QFileDialog.getOpenFileName(None, \
            'Import INP/UNV file', model.job.dir, \
            'INP (*.inp);;UNV (*.unv)')[0]

    if file_name:

        # Clear log window before new import
        # if clear_textEdit:
        #     model.textEdit.setText('')

        # Clear selection before import new model
        # if model.settings.show_vtk:
        #     model.VTK.actionSelectionClear()

        # Generate new KOM without implementations
        model.KOM = KOM()

        # Rename job before tree regeneration
        model.job.rename(file_name[:-4] + '.inp')

        # Convert UNV to INP
        if file_name.lower().endswith('.unv'):
            model.job.convertUNV()
            if not os.path.isfile(model.job.inp):
                logging.error('Error converting ' + model.job.unv)
                return

        # Show model name in window's title
        # MainWindow.setWindowTitle('CalculiX CAE - ' + model.job.name)


        # Parse INP and enrich KOM with parsed objects
        logging.info('Loading ' + model.job.inp + '.')
        lines = readLines(model.job.inp)
        importer(lines, model.KOM) # pass whole INP-file to the parser


        # Add parsed implementations to the tree
        # model.tree.generateTreeView()

        # Parse mesh
        # model.mesh = Parse(model.job.inp) # parse mesh

        # Create ugrid from mesh
        # if model.settings.show_vtk:
        #     ugrid = model.VTK.mesh2ugrid(model.mesh)

        #     # Plot ugrid in VTK
        #     if ugrid:
        #         model.VTK.mapper.SetInputData(ugrid)
        #         model.VTK.actionViewIso() # iso view after import


# Enrich KOM with keywords from INP_doc
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
            logging.debug('path found: ' + str([item.name for item in path]))

            if path:
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


# Menu File -> Write INP file
def writeInput(file_name=None):
    if not file_name:
        file_name = QFileDialog.getSaveFileName(None, \
            'Write INP file', model.job.dir, \
            'Input files (*.inp)')[0]
    if file_name:
        with open(file_name, 'w') as f:
            writer(model.KOM.root, f, 0)
        logging.info('Input written to ' + file_name)

        # Rename job and update treeView if new INP-file name was selected
        old_job_name = model.job.name
        model.job.rename(file_name)
        if model.job.name != old_job_name:
            model.tree.appendJobName()


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
