#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, UJV Rez, June 2019.
    Distributed under GNU General Public License, version 2.

    CalculiX CAE
    Main window application
"""


import sys, vtk
from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import QFileDialog, QTreeWidgetItem, QVBoxLayout
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from INPParser import Mesh


# Read the form with GUI
Ui_MainWindow, QtBaseClass = uic.loadUiType('ccx_form.ui')


class ccx_cae(QtWidgets.QMainWindow, Ui_MainWindow):

    # Create main window
    def __init__(self):
        # Create main window
        QtWidgets.QMainWindow.__init__(self)

        # Load form
        Ui_MainWindow.__init__(self)
        self.setupUi(self)

        # Write log
        self.log_widget.append('Started.')
        # TODO put here all messages from command prompt

        # Read CalculiX object model
        self.readObjectModel()

        # Creates VTK widget with default ugrid, adds it to the form
        self.drawVTKWidget(self.defaultUGrid())

        # Actions
        self.actionImportInpMesh.triggered.connect(self.importInpMesh)


    # Create VTK widget and draw some object
    def drawVTKWidget(self, ugrid):
        # Replace or create VTK widget
        try:
            self.vl.removeWidget(self.vtk_widget)
        except:
            pass
        self.vtk_widget = QVTKRenderWindowInteractor()
        self.vl.addWidget(self.vtk_widget)

        # Create the graphics structure
        renderer = vtk.vtkRenderer() # renderer renders into the window
        window = self.vtk_widget.GetRenderWindow()
        window.AddRenderer(renderer)
        interactor = window.GetInteractor() # interactor captures mouse events

        # Mapper pushes the geometry into the graphics library
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(ugrid)

        # The actor is a grouping mechanism
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(0.2, 0.6, 0.8)
        actor.RotateX(30.0)
        actor.RotateY(-45.0)

        # Add the actors to the renderer, set the background
        renderer.AddActor(actor)
        renderer.SetBackground(1, 1, 1)
        renderer.ResetCamera()

        # Apply layout: it will expand vtk_widget to the frame size
        self.frame_widget.setLayout(self.vl)

        # Initalize the interactor before an event loop
        interactor.Initialize()

        # Start the event loop
        interactor.Start()

        self.log_widget.append('Rendering OK.')


    # Read CalculiX keywords hierarchy
    def readObjectModel(self):
        # TODO catch click action on keyword
        parent_dict = {}
        try:
            with open('ccx_dom.txt', 'r') as f:
                lines = f.readlines() # read the whole file
                for line in lines: # iterate over keywords

                    # Skip comments and empty lines
                    if line.strip() == '': continue 
                    if line.strip().startswith('**'): continue
                    if line.strip().startswith('-'): continue

                    level = 0
                    keyword = line.rstrip() # cut '\n'
                    while keyword.startswith('    '):
                        level += 1
                        keyword = keyword[4:]

                    parent_dict[level] = QTreeWidgetItem([keyword])
                    if level == 0:
                        self.tree_widget.addTopLevelItem(parent_dict[0])
                    else:
                        parent_dict[level-1].addChild(parent_dict[level])
            self.log_widget.append('CalculiX object model generated.')
        except:
            self.log_widget.append('Error reading keywords hierarchy!')


    # Import mesh and display it in the VTK widget
    def importInpMesh(self):
        file_name = QFileDialog.getOpenFileName(self,\
            'Import .inp mesh', '', 'Input files (*.inp);;All Files (*)')[0]
        if file_name:
            # Parce mesh and transfer it to VTK
            mesh = Mesh(file_name) # parse mesh
            points = vtk.vtkPoints()
            for n in mesh.nodes.keys(): # create VTK points from mesh nodes
                points.InsertPoint(n, mesh.nodes[n])
            ugrid = vtk.vtkUnstructuredGrid() # create empty grid in VTK
            ugrid.Allocate(len(mesh.elements)) # allocate memory fo all elements
            ugrid.SetPoints(points) # insert all points to the grid
            for e in mesh.elements.keys():
                vtk_element_type = Mesh.convert_elem_type(mesh.types[e])
                node_numbers = mesh.elements[e] # list of nodes in the element
                ugrid.InsertNextCell(vtk_element_type, len(node_numbers), node_numbers) # create VTK element

            self.log_widget.append(file_name + '\nSuccessfully loaded.')
            self.drawVTKWidget(ugrid)


    # Default ugrid show on application load
    def defaultUGrid(self):
        nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ]
        element = [0, 1, 2, 3, 4, 5, 6, 7]

        points = vtk.vtkPoints()
        for i in range(len(nodes)):
            points.InsertPoint(i, nodes[i])

        ugrid = vtk.vtkUnstructuredGrid()
        ugrid.Allocate(1)
        ugrid.SetPoints(points)
        ugrid.InsertNextCell(vtk.VTK_HEXAHEDRON, 8, element)
        return ugrid


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = ccx_cae()
    window.showMaximized()
    sys.exit(app.exec_())
    os.system('py3clean .') # doesn't work
