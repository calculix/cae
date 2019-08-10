# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    INP importer:
        Enrich DOM with implementations from parsed file.
        Generate new tree with keyword implementations.
        Parse mesh and build ugrid
        Display ugrid in the VTK widget

    INP exporter:
        Recursively write implementation's INP_code to output .inp-file

    Depends on ccx_cae.CAE.
"""


from PyQt5 import QtWidgets
import vtk, os, logging
import ccx_dom, ccx_mesh, ccx_vtk


class IE:


    def __init__(self, CAE):
        self.CAE = CAE

        # Actions
        self.CAE.actionFileImportINP.triggered.connect(self.importINP)
        self.CAE.actionFileExportINP.triggered.connect(self.exportINP)


    # Menu File -> Import INP file
    def importINP(self, file_name=None):

        if not file_name:
            file_name = QtWidgets.QFileDialog.getOpenFileName(None, \
                'Import INP file', '', 'Input files (*.inp);;All Files (*)')[0]

        if file_name:

            # Show model name in window's title
            self.CAE.setWindowTitle('Calculix CAE - ' + os.path.basename(file_name))

            # Clear log window before new import
            self.CAE.textEdit.setText('')
            logging.info('Loading ' + file_name + '.')

            # Clear selection before import new model
            self.CAE.VTK.actionSelectionClear()

            # Generate new DOM without implementations
            self.CAE.DOM = ccx_dom.DOM()

            # Parse INP and enrich DOM with parsed objects
            lines = ccx_mesh.read_lines(file_name)
            self.importer(lines) # pass whole INP-file to the parser

            # Add parsed implementations to the tree
            self.CAE.tree.generateTreeView()

            # Parse mesh
            self.CAE.mesh = ccx_mesh.Parse(file_name) # parse mesh

            # Create ugrid from mesh
            ugrid = self.CAE.VTK.mesh2ugrid(self.CAE.mesh)

            # Plot ugrid in VTK
            if ugrid:
                self.CAE.VTK.mapper.SetInputData(ugrid)
                self.CAE.VTK.actionViewIso() # iso view after import


    # Enrich DOM with keywords from INP_doc
    def importer(self, INP_doc):
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

                keyword_chain.append(keyword_name)

                # Find DOM keyword path corresponding to keyword_chain
                path = self.CAE.DOM.getPath(keyword_chain)
                logging.debug('Path found: ' + str([item.name for item in path]))
                if path:

                    # Read INP_code for the current keyword 
                    INP_code = [line] # here line is only r-stripped
                    while i+1 < len(INP_doc) and \
                        not INP_doc[i+1].lstrip().startswith('*'): # there will be no comments
                        INP_code.append(INP_doc[i+1].rstrip()) # cut '\n'
                        i += 1

                    # Create keyword implementations
                    impl = None
                    path_as_string = '' # string representation of 'path' accounting for implementations
                    for j in range(len(path)):
                        if impl:
                            item = impl.getItemByName(path[j].name)
                        else:
                            item = path[j] # keyword or group
                        path_as_string += '/' + item.name
                        if j == len(path) - 1: # last item is always keyword
                            impl = ccx_dom.implementation(item, INP_code) # create, for example, MATERIAL-1
                        elif item.item_type == ccx_dom.item_type.KEYWORD:
                            # If we are here, then for this keyword implementation was created previously
                            counter = impl_counter[path_as_string]-1
                            impl = item.items[counter] # first implementation, for example, STEP-1
                            path_as_string += '/' + impl.name
                        else:
                            impl = item

                    # Count implementation
                    if path_as_string in impl_counter:
                        impl_counter[path_as_string] += 1
                    else:
                        impl_counter[path_as_string] = 1

                else:
                    logging.info('Wrong keyword {}.'.format(keyword_name))


    # Menu File -> Write INP file
    def exportINP(self):

        file_name = QtWidgets.QFileDialog.getSaveFileName(None, \
            'Write INP file', '', 'Input files (*.inp);;All Files (*)')[0]

        # Recursively iterate over DOM items, write INP_code for each implementation
        if file_name:
            with open(file_name, 'w') as f:
                self.exporter(self.CAE.DOM.root, f, 0)
            logging.info('Input written!')


    # Recursively write implementation's INP_code to output .inp-file
    def exporter(self, parent, f, level):

        # Level is used for padding
        if parent.item_type == ccx_dom.item_type.IMPLEMENTATION:
            level += 1

        # For each group/keyword from DOM
        for item in parent.items:

            if item.item_type == ccx_dom.item_type.ARGUMENT:
                continue

            if item.item_type == ccx_dom.item_type.IMPLEMENTATION:
                # INP_code is stripped
                f.write(' '*4*level + item.INP_code[0] + '\n')
                for line in item.INP_code[1:]:
                    f.write(' '*4*(level+1) + line + '\n')

            # Continue call iterator until dig to implementation
            self.exporter(item, f, level)
