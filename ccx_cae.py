#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, UJV Rez, June 2019.
    Distributed under GNU General Public License, version 2.

    CalculiX CAE
    Main window application
"""


import sys, os, vtk, argparse, logging
from PyQt5 import QtWidgets, uic, QtGui
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from INPParser import Mesh
from ccx_select_style import *
from ccx_dom import *
import ccx_dialog


# Read the form with GUI
Ui_MainWindow = uic.loadUiType('ccx_form.ui')[0]


class ccx_cae(QtWidgets.QMainWindow, Ui_MainWindow):

    # Create main window
    def __init__(self):

        # Create main window
        QtWidgets.QMainWindow.__init__(self)

        # Load form
        Ui_MainWindow.__init__(self)
        self.setupUi(self)

        # Configure logging
        logging.info = self.logging_info
        logging.error = self.logging_error

        # Create VTK widget
        self.createVTKWidget()

        # Read CalculiX object model
        self.dom = None
        self.readObjectModel()

        # Default start model could be chosen with command line parameter
        parser = argparse.ArgumentParser()
        parser.add_argument("--mesh", "-mesh",
                            help="Mesh .inp file",
                            type=str, default='mesh_baffle.inp')
        args = parser.parse_args()

        # Creates VTK widget with default ugrid, adds it to the form
        self.importMeshFromInp(args.mesh)

        # Actions
        self.tree_view.doubleClicked.connect(self.treeViewDoubleClicked)
        self.actionImportMeshFromInp.triggered.connect(self.importMeshFromInp)
        self.actionCameraFitView.triggered.connect(self.cameraFitView)
        self.actionSelectNodes.triggered.connect(self.selectNodes)
        self.actionSelectElements.triggered.connect(self.selectElements)
        self.actionClearSelection.triggered.connect(self.clearSelection)


    # Create VTK widget and draw some object: called once during startup
    def createVTKWidget(self):
        # Create VTK widget
        self.vtk_widget = QVTKRenderWindowInteractor()
        self.vl.addWidget(self.vtk_widget)
        self.frame_widget.setLayout(self.vl) # apply layout: it will expand vtk_widget to the frame size

        # Create the graphics structure
        self.renderer = vtk.vtkRenderer() # renderer renders into the window
        self.window = self.vtk_widget.GetRenderWindow()
        self.window.AddRenderer(self.renderer)
        self.interactor = self.window.GetInteractor() # captures mouse events

        # Add the custom style
        # style = MouseInteractorHighLightActor()
        # style.SetDefaultRenderer(self.renderer)
        # self.interactor.SetInteractorStyle(style)
        self.interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())

        # Mapper pushes the mesh data into the graphics library
        self.mapper = vtk.vtkDataSetMapper()

        # The actor is a grouping mechanism
        self.actor = vtk.vtkActor()
        self.actor.SetMapper(self.mapper)
        self.actor.GetProperty().SetColor(0.2, 0.6, 0.8)

        # Add the actors to the renderer, set the background
        self.renderer.AddActor(self.actor)
        self.renderer.SetBackground(1, 1, 1)

        # Tune camera
        """
            self.camera = self.renderer.GetActiveCamera()
            self.logging_info('Camera\'s focal point is ' + str(self.camera.GetFocalPoint()))
            self.logging_info('Camera\'s position is ' + str(self.camera.GetPosition()))
            self.logging_info('Camera\'s ViewUp is ' + str(self.camera.GetViewUp()))
            self.logging_info('Camera\'s distance is ' + str(self.camera.GetDistance()))
            self.logging_info('Camera\'s Roll is ' + str(self.camera.GetRoll()))
            self.logging_info('Camera\'s ViewAngle is ' + str(self.camera.GetViewAngle()))
            self.logging_info('Camera\'s ParallelScale is ' + str(self.camera.GetParallelScale()))
            self.logging_info('Camera\'s ClippingRange is ' + str(self.camera.GetClippingRange()))
            self.logging_info('Camera\'s WindowCenter is ' + str(self.camera.GetWindowCenter()))
            self.logging_info('Camera\'s orientation is ' + str(self.camera.GetOrientation())) # TODO camera is rotated
        """

        # Add orientation axes
        self.axes = vtk.vtkAxesActor()
        self.axes.SetShaftTypeToCylinder()
        self.axes.SetXAxisLabelText('X')
        self.axes.SetYAxisLabelText('Y')
        self.axes.SetZAxisLabelText('Z')
        self.axes.SetTotalLength(1.0, 1.0, 1.0)
        self.axes.SetCylinderRadius(0.5 * self.axes.GetCylinderRadius())
        tprop = vtk.vtkTextProperty()
        tprop.SetColor(0, 0, 0) # font color
        self.axes.GetXAxisCaptionActor2D().SetCaptionTextProperty(tprop)
        self.axes.GetYAxisCaptionActor2D().SetCaptionTextProperty(tprop)
        self.axes.GetZAxisCaptionActor2D().SetCaptionTextProperty(tprop)
        self.omw = vtk.vtkOrientationMarkerWidget()
        self.omw.SetOrientationMarker(self.axes)
        self.omw.SetViewport(0, 0, 0.15, 0.15) # position lower left in the viewport
        self.omw.SetInteractor(self.interactor)
        self.omw.EnabledOn()
        self.omw.InteractiveOn()

        # Render mesh
        self.window.Render()

        # Start the event loop
        self.interactor.Start()


    # Logging functions
    def logging_info(self, msg):
        self.log_widget.append('<span style=\'color:Black;\'>' + msg + '</span>')
        self.log_widget.moveCursor(QtGui.QTextCursor.End) # scroll text to the end
    def logging_error(self, msg):
        self.log_widget.append('<span style=\'color:Red;\'>' + msg + '</span>')
        self.log_widget.moveCursor(QtGui.QTextCursor.End) # scroll text to the end


    # Read CalculiX keywords hierarchy
    def readObjectModel(self):
        try:
            self.dom = ccx_dom() # generate DOM based on keywords hierarchy from ccx_dom.txt
            self.model = QtGui.QStandardItemModel()
            parent = self.model.invisibleRootItem() # top element in QTreeView
            self.addToTree(parent, self.dom.root.items) # pass root - group 'Model'
            self.tree_view.setModel(self.model)
            self.tree_view.expandAll()
            self.logging_info('CalculiX object model generated.')
        except:
            self.logging_error('Can\'t generate keywords hierarchy!')
    def addToTree(self, parent, children):
        """
            parent is QtGui.QStandardItem
            children are items of ccx_dom object
        """
        for item in children:
            # print('\t'*level + item.name)
            if ('keyword' in item.item_type) or ('group' in item.item_type):
                tree_element = QtGui.QStandardItem(item.name)
                tree_element.setData(item)
                parent.appendRow(tree_element)
                self.addToTree(tree_element, item.items)


    # Import mesh and display it in the VTK widget
    def importMeshFromInp(self, file_name=None):
        if not file_name:
            file_name = QtWidgets.QFileDialog.getOpenFileName(self,\
                'Import .inp mesh', '', 'Input files (*.inp);;All Files (*)')[0]

        self.logging_info('Loading ' + file_name + '.')

        if file_name:
            # Parse mesh and transfer it to VTK
            mesh = Mesh(file_name) # parse mesh
            points = vtk.vtkPoints()
            for n in mesh.nodes.keys(): # create VTK points from mesh nodes
                points.InsertPoint(n-1, mesh.nodes[n]) # node numbers should start from 0!
            self.ugrid = vtk.vtkUnstructuredGrid() # create empty grid in VTK
            self.ugrid.Allocate(len(mesh.elements)) # allocate memory fo all elements
            self.ugrid.SetPoints(points) # insert all points to the grid
            for e in mesh.elements.keys():
                vtk_element_type = Mesh.convert_elem_type(mesh.types[e])
                node_numbers = [n-1 for n in mesh.elements[e]] # list of nodes in the element: node numbers should start from 0!
                self.ugrid.InsertNextCell(vtk_element_type, len(node_numbers), node_numbers) # create VTK element

            self.mapper.SetInputData(self.ugrid) # ugrid is our mesh data
            self.cameraFitView() # reset camera

            self.logging_info('Rendering OK.')


    # Automatically set up the camera
    def cameraFitView(self):
        # TODO rotate camera to initial position
        self.renderer.ResetCamera()
        self.window.Render()


    # Double click on tree_view item: create keyword in DOM
    def treeViewDoubleClicked(self, index):
        item = self.tree_view.model().itemFromIndex(index)
        if item.text().startswith('*'):
            item.data().printAll()
            # self.logging_info(item.text())

            # TODO dialog window with iten edit

            # Staring Functions for Execution
            dinput = ['LastName','Country','Age']
            # Call the UI and get the inputs
            dialog = ccx_dialog.Dialog(dinput)
            if dialog.exec_() == ccx_dialog.Dialog.Accepted:
                name, item, value = dialog.get_output()
                print(name, item, value)


    # Menu Select->Nodes
    def selectNodes(self):
        self.clearSelection() # clear selection before new call
        style = ccx_select_nodes_style(self.renderer, self.window, self.log_widget)
        self.interactor.SetInteractorStyle(style)


    # Menu Select->Elements
    def selectElements(self):
        self.clearSelection() # clear selection before new call
        style = ccx_select_elements_style(self.renderer, self.window, self.ugrid, self.log_widget)
        self.interactor.SetInteractorStyle(style)


    # Menu Select->Clear selection
    def clearSelection(self):
        self.renderer.RemoveAllViewProps()
        self.renderer.AddActor(self.actor)
        self.window.Render()

        self.logging_info('Clear selection.')


    """
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
    """


if __name__ == '__main__':

    # Clean cached files before start
    os.system('py3clean .')

    app = QtWidgets.QApplication(sys.argv)
    window = ccx_cae()
    window.show() # window.showMaximized() or window.show()
    sys.exit(app.exec_())
