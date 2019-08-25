# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    INP importer:
        Enrich KOM with implementations from parsed file.
        Generate new tree with keyword implementations.
        Parse mesh and build ugrid
        Display ugrid in the VTK widget

    INP exporter:
        Recursively write implementation's INP_code to output .inp-file

    Depends on ccx_cae.CAE.
"""


from PyQt5.QtWidgets import QFileDialog
import os, logging
import ccx_kom, ccx_mesh, ccx_log


class IE:


    def __init__(self, CAE, settings):
        self.CAE = CAE
        self.settings = settings

        # Actions
        self.CAE.actionFileImport.triggered.connect(self.importFile)
        self.CAE.actionFileWriteInput.triggered.connect(self.writeInput)


    # Menu File -> Import
    def importFile(self, file_name=None):

        if not file_name:
            file_name = QFileDialog.getOpenFileName(None, \
                'Import INP file', self.CAE.job.dir, \
                'INP (*.inp);;UNV (*.unv);;All Files (*)')[0]

        if file_name:

            # Clear log window before new import
            self.CAE.textEdit.setText('')

            # Clear selection before import new model
            if self.settings.show_vtk:
                self.CAE.VTK.actionSelectionClear()

            # Generate new KOM without implementations
            self.CAE.KOM = ccx_kom.KOM()

            # Rename job before tree regeneration
            self.CAE.job.rename(file_name[:-4] + '.inp')

            # Convert UNV to INP
            if file_name.lower().endswith('.unv'):
                self.CAE.job.importUNV()
                if not os.path.isfile(self.CAE.job.inp):
                    logging.error('Error converting ' + self.CAE.job.unv)
                    return

            # Show model name in window's title
            self.CAE.setWindowTitle('CalculiX CAE - ' + os.path.basename(self.CAE.job.inp))

            # Parse INP and enrich KOM with parsed objects
            logging.info('Loading ' + self.CAE.job.inp + '.')
            lines = ccx_mesh.readLines(self.CAE.job.inp)
            self.importer(lines) # pass whole INP-file to the parser

            # Add parsed implementations to the tree
            self.CAE.tree.generateTreeView()

            # Parse mesh
            self.CAE.mesh = ccx_mesh.Parse(self.CAE.job.inp) # parse mesh

            # Create ugrid from mesh
            if self.settings.show_vtk:
                ugrid = self.CAE.VTK.mesh2ugrid(self.CAE.mesh)

                # Plot ugrid in VTK
                if ugrid:
                    self.CAE.VTK.mapper.SetInputData(ugrid)
                    self.CAE.VTK.actionViewIso() # iso view after import


    # Enrich KOM with keywords from INP_doc
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

                # Find KOM keyword path corresponding to keyword_chain
                keyword_chain.append(keyword_name)
                path = self.CAE.KOM.getPath(keyword_chain)
                logging.debug('path found: ' + str([item.name for item in path]))

                if path:
                    # Read INP_code for the current keyword
                    INP_code = [line] # line is stripped in ccx_mesh
                    while i+1 < len(INP_doc) and \
                        not INP_doc[i+1].startswith('*'): # here will be no comments - they are removed in ccx_mesh
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
                            impl = ccx_kom.implementation(item, INP_code)
                            logging.info(impl.name)
                            logging.debug('1')
                        elif item.item_type == ccx_kom.item_type.KEYWORD:
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
    def writeInput(self, file_name=None):
        if not file_name:
            file_name = QFileDialog.getSaveFileName(None, \
                'Write INP file', self.CAE.job.dir, \
                'Input files (*.inp);;All Files (*)')[0]
        if file_name:
            with open(file_name, 'w') as f:
                self.writer(self.CAE.KOM.root, f, 0)
            logging.info('Input written to ' + file_name)

            # Rename job and update treeView if new INP-file name was selected
            old_job_name = self.CAE.job.name
            self.CAE.job.rename(file_name)
            if self.CAE.job.name != old_job_name:
                self.CAE.tree.appendJobName()


    # Recursively write implementation's INP_code to output .inp-file
    def writer(self, parent, f, level):

        # Level is used for padding
        if parent.item_type == ccx_kom.item_type.IMPLEMENTATION:
            level += 1

        # For each group/keyword from KOM
        for item in parent.items:

            if item.item_type == ccx_kom.item_type.ARGUMENT:
                continue

            if item.item_type == ccx_kom.item_type.IMPLEMENTATION:
                # INP_code is stripped
                f.write(' '*4*level + item.INP_code[0] + '\n')
                for line in item.INP_code[1:]:
                    f.write(' '*4*(level+1) + line + '\n')

            # Continue call iterator until dig to implementation
            self.writer(item, f, level)
