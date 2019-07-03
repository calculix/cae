# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, June 2019.
    Distributed under GNU General Public License, version 2.

    CalculiX CAE
    Main window application
"""


import sys, os, argparse, vtk
from PyQt5 import QtWidgets, uic, QtCore, QtGui, Qt
import ccx_mesh, ccx_tree, ccx_log, ccx_vtk, ccx_dom


class CAE(QtWidgets.QMainWindow):

    # Create main window
    def __init__(self):

        # Create main window
        QtWidgets.QMainWindow.__init__(self)

        # Load form
        uic.loadUi('ccx_cae.ui', self)

        # Configure logging
        self.logger = ccx_log.logger(self)

        # Create VTK widget
        self.VTK = ccx_vtk.VTK(self) # create everything for model visualization
        self.vl.addWidget(self.VTK.widget) # add vtk_widget to the form
        self.frame.setLayout(self.vl) # apply layout: it will expand vtk_widget to the frame size

        # Generate CalculiX DOM based on keywords hierarchy from ccx_dom.txt
        self.DOM = ccx_dom.DOM(self) # empty DOM w/o implementations

        # Generate empty CalculiX model and treeView items
        self.tree = ccx_tree.tree(self)

        # Default start model could be chosen with command line parameter
        parser = argparse.ArgumentParser()
        parser.add_argument("--mesh", "-mesh",
                            help="Mesh .inp file",
                            type=str, default='./examples/acou1.inp')
        args = parser.parse_args()

        # Creates default ugrid, adds it to the VTK mapper
        self.importINP(args.mesh)

        # Actions
        self.actionImportINP.triggered.connect(self.importINP)
        self.actionWriteINP.triggered.connect(self.writeINP)
        self.treeView.keyPressEvent = self.keyPressEvent


    # Delete keyword's implementation in the treeView by pressing 'Delete' button
    def keyPressEvent(self, e):
        if e.key() == QtCore.Qt.Key_Delete:
            self.tree.actionDeleteImplementation()


    # Menu File -> Import INP file
    # Import mesh and display it in the VTK widget
    def importINP(self, file_name=None):

        if not file_name:
            file_name = QtWidgets.QFileDialog.getOpenFileName(None, \
                'Import INP file', '', 'Input files (*.inp);;All Files (*)')[0]

        self.logger.info('Loading ' + file_name + '.')

        if file_name:
            # Generate CalculiX DOM based on keywords hierarchy from ccx_dom.txt
            self.DOM = ccx_dom.DOM(self)

            # Parse model from INP - it will modify DOM
            with open(file_name, 'r') as f:
                INP_doc = f.readlines()
                self.parseINP(INP_doc) # enrich DOM with parsed objects

            # Update DOM in ccx_tree class
            self.tree.DOM = self.DOM

            # Regenerate treeView items to account for modifications in DOM
            self.tree.generateTreeView()

            # Parse mesh and transfer it to VTK
            mesh = ccx_mesh.Parse(file_name, self) # parse mesh, textEdit passed for logging
            try:
                points = vtk.vtkPoints()
                for n in mesh.nodes.keys(): # create VTK points from mesh nodes
                    points.InsertPoint(n-1, mesh.nodes[n]) # node numbers should start from 0!

                self.ugrid = vtk.vtkUnstructuredGrid() # create empty grid in VTK
                self.ugrid.Allocate(len(mesh.elements)) # allocate memory fo all elements
                self.ugrid.SetPoints(points) # insert all points to the grid

                for e in mesh.elements.keys():
                    vtk_element_type = ccx_mesh.Parse.convert_elem_type(mesh.types[e])
                    node_numbers = [n-1 for n in mesh.elements[e]] # list of nodes in the element: node numbers should start from 0!
                    self.ugrid.InsertNextCell(vtk_element_type, len(node_numbers), node_numbers) # create VTK element

                self.VTK.mapper.SetInputData(self.ugrid) # ugrid is our mesh data
                self.VTK.actionViewIso() # iso view after import
                self.logger.info('Rendering OK.')
            except:
                self.logger.error('Can\'t render INP mesh.')


    # Parse CalculiX'es keywords in the INP_doc lines
    def parseINP(self, INP_doc):
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
                path = self.DOM.getPath(keyword_chain)

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


    # Menu File -> Write INP file
    # Write input file for CalculiX
    def writeINP(self):

        file_name = QtWidgets.QFileDialog.getSaveFileName(self, \
            'Write INP file', '', 'Input files (*.inp);;All Files (*)')[0]

        if file_name:
            with open(file_name, 'w') as self.f:
                # Start recursive DOM iterator
                self.iterator(self.tree.DOM.root)

            self.logger.info('Input written!')


    # Recursively iterate over DOM items, write INP_code for each implementation
    def iterator(self, parent):
        for item in parent.items: # for each group/keyword from DOM

            if item.item_type == item_type.ARGUMENT:
                continue

            if item.item_type == item_type.IMPLEMENTATION:
                for line in item.INP_code:
                    self.f.write(line + '\n')

            self.iterator(item) # continue call iterator until dig to implementation


if __name__ == '__main__':

    # Clean cached files before start
    os.system('py3clean .')

    app = QtWidgets.QApplication(sys.argv)
    window = CAE()
    window.show() # window.showMaximized() or window.show()
    sys.exit(app.exec_())
