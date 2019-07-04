# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    INP importer/parser and exporter/writer
"""


from PyQt5 import QtWidgets
import vtk
import ccx_dom, ccx_mesh


class inp:


    def __init__(self, CAE):
        self.CAE = CAE

        # Actions
        self.CAE.actionImportINP.triggered.connect(self.importINP)
        self.CAE.actionExportINP.triggered.connect(self.exportINP)


    # Menu File -> Import INP file
    # Import mesh and display it in the VTK widget
    def importINP(self, file_name=None):

        if not file_name:
            file_name = QtWidgets.QFileDialog.getOpenFileName(None, \
                'Import INP file', '', 'Input files (*.inp);;All Files (*)')[0]

        if file_name:
            self.CAE.logger.info('Loading ' + file_name + '.')

            # Generate CalculiX DOM based on keywords hierarchy from ccx_dom.txt
            self.CAE.DOM = ccx_dom.DOM(self.CAE)

            # Parse model from INP - it will modify DOM
            with open(file_name, 'r') as f:
                INP_doc = f.readlines()
                # Parse CalculiX'es keywords in the INP_doc lines
                self.parser(INP_doc) # enrich DOM with parsed objects

            # Update DOM in ccx_tree class
            # self.CAE.tree.DOM = self.CAE.DOM

            # Regenerate treeView items to account for modifications in DOM
            self.CAE.tree.generateTreeView()

            # Parse mesh and transfer it to VTK
            mesh = ccx_mesh.Parse(file_name, self.CAE) # parse mesh, textEdit passed for logging
            # try:
            points = vtk.vtkPoints()
            for n in mesh.nodes.keys(): # create VTK points from mesh nodes
                points.InsertPoint(n-1, mesh.nodes[n]) # node numbers should start from 0!

            ugrid = vtk.vtkUnstructuredGrid() # create empty grid in VTK
            ugrid.Allocate(len(mesh.elements)) # allocate memory fo all elements
            ugrid.SetPoints(points) # insert all points to the grid

            for e in mesh.elements.keys():
                ccx_element_type = mesh.types[e]
                vtk_element_type = ccx_mesh.Parse.convert_elem_type(ccx_element_type)
                node_numbers = [n-1 for n in mesh.elements[e]] # list of nodes in the element: node numbers should start from 0!
                ugrid.InsertNextCell(vtk_element_type, len(node_numbers), node_numbers) # create VTK element
                # print(ccx_element_type, 'to', vtk_element_type, ':', e, node_numbers)

            self.CAE.VTK.mapper.SetInputData(ugrid) # ugrid is our mesh data
            self.CAE.VTK.actionViewIso() # iso view after import
            self.CAE.logger.info('Rendering OK.')
            # except:
            #     self.CAE.logger.error('Can\'t render INP mesh.')
    def parser(self, INP_doc):
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
                    if j == len(path) - 1:
                        # Last item is always keyword
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
                # print(path_as_string)


    # Menu File -> Write INP file
    # Write input file for CalculiX
    def exportINP(self):

        file_name = QtWidgets.QFileDialog.getSaveFileName(None, \
            'Write INP file', '', 'Input files (*.inp);;All Files (*)')[0]

        if file_name:
            with open(file_name, 'w') as f:
                # Recursively iterate over DOM items, write INP_code for each implementation
                self.writer(self.CAE.DOM.root, f)

            self.CAE.logger.info('Input written!')
    def writer(self, parent, f):
        for item in parent.items: # for each group/keyword from DOM

            if item.item_type == ccx_dom.item_type.ARGUMENT:
                continue

            if item.item_type == ccx_dom.item_type.IMPLEMENTATION:
                for line in item.INP_code:
                    f.write(line + '\n')

            self.writer(item, f) # continue call iterator until dig to implementation
