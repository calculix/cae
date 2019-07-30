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
import vtk, os
import ccx_dom, ccx_mesh, ccx_cae_log, ccx_vtk


class IE:


    def __init__(self, CAE):
        self.CAE = CAE

        # Actions
        self.CAE.actionFileImportINP.triggered.connect(lambda:
            self.CAE.logger.messages(self.importINP()))
        self.CAE.actionFileExportINP.triggered.connect(lambda:
            self.CAE.logger.messages(self.exportINP()))


    # Menu File -> Import INP file
    def importINP(self, file_name=None):
        msg_list = [] # list of messages for logger

        if not file_name:
            file_name = QtWidgets.QFileDialog.getOpenFileName(None, \
                'Import INP file', '', 'Input files (*.inp);;All Files (*)')[0]

        if file_name:

            # Show model name in window's title
            self.CAE.setWindowTitle('Calculix CAE - ' + os.path.basename(file_name))

            # Clear log window before new import
            self.CAE.textEdit.setText('')

            # Clear selection before import new model
            self.CAE.VTK.actionSelectionClear()

            msg_text = 'Loading ' + file_name + '.'
            msg = ccx_cae_log.msg(ccx_cae_log.msgType.INFO, msg_text)
            msg_list.append(msg)

            # Generate new DOM without implementations
            self.CAE.DOM = ccx_dom.DOM()
            msg_list.extend(self.CAE.DOM.msg_list) # 'CalculiX object model generated.'

            # Parse INP and enrich DOM with parsed objects
            with open(file_name, 'r') as f:
                self.importer(f.readlines(), msg_list) # pass whole INP-file to the parser

            # Add parsed implementations to the tree
            self.CAE.tree.generateTreeView()

            # Parse mesh
            self.CAE.mesh = ccx_mesh.Parse(file_name) # parse mesh
            msg_list.extend(self.CAE.mesh.msg_list) # show info about nodes, elements etc.

            # Create ugrid from mesh
            msgs, ugrid = self.CAE.VTK.mesh2ugrid(self.CAE.mesh)
            msg_list.extend(msgs)

            # Plot ugrid in VTK
            if ugrid:
                self.CAE.VTK.mapper.SetInputData(ugrid)
                self.CAE.VTK.actionViewIso() # iso view after import
                # msg_list.extend(self.CAE.VTK.msg_list)

        return msg_list


    # Enrich DOM with keywords from INP_doc
    def importer(self, INP_doc, msg_list):
        keyword_chain = []
        impl_counter = {}
        for i in range(len(INP_doc)):
            line = INP_doc[i].rstrip() # cut '\n'

            # Skip comments and empty lines
            if line.strip() == '':
                continue
            elif line.strip().startswith('**'):
                continue

            # Parse keyword
            elif line.startswith('*'):

                # Distinguish 'NODE' and 'NODE PRINT'
                if ',' in line:
                    keyword_name = line.split(',')[0]
                else:
                    keyword_name = line

                keyword_name = keyword_name.lower().strip()
                keyword_chain.append(keyword_name)

                # Find DOM keyword path corresponding to keyword_chain
                path = self.CAE.DOM.getPath(keyword_chain)
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
                    msg_text = 'Wrong keyword {}.'.format(keyword_name)
                    msg = ccx_cae_log.msg(ccx_cae_log.msgType.INFO, msg_text)
                    msg_list.append(msg)


    # Menu File -> Write INP file
    def exportINP(self):
        msg_list = [] # list of messages for logger

        file_name = QtWidgets.QFileDialog.getSaveFileName(None, \
            'Write INP file', '', 'Input files (*.inp);;All Files (*)')[0]

        if file_name:
            with open(file_name, 'w') as f:
                # Recursively iterate over DOM items, write INP_code for each implementation
                self.exporter(self.CAE.DOM.root, f)

            # Log message
            # self.CAE.logger.info('Input written!')
            msg_text = 'Input written!'
            msg = ccx_cae_log.msg(ccx_cae_log.msgType.INFO, msg_text)
            msg_list.append(msg)

        return msg_list


    # Recursively write implementation's INP_code to output .inp-file
    def exporter(self, parent, f):
        for item in parent.items: # for each group/keyword from DOM

            if item.item_type == ccx_dom.item_type.ARGUMENT:
                continue

            if item.item_type == ccx_dom.item_type.IMPLEMENTATION:
                for line in item.INP_code:
                    f.write(line + '\n')

            self.exporter(item, f) # continue call iterator until dig to implementation
