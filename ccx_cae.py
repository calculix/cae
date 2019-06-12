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
import INPParser, ccx_dom, ccx_dialog, ccx_select_style


class CAE(QtWidgets.QMainWindow):

    # Create main window
    def __init__(self):

        # Create main window
        QtWidgets.QMainWindow.__init__(self)

        # Load form
        uic.loadUi('ccx_cae.ui', self)

        # Configure logging
        logging.info = self.logging_info
        logging.error = self.logging_error

        # Create VTK widget
        self.createVTKWidget()

        # Generate CalculiX DOM based on keywords hierarchy from ccx_dom.txt
        try:
            self.dom = ccx_dom.DOM()
            self.logging_info('CalculiX object model generated.')
            self.generateTreeView()
        except:
            self.logging_error('Can\'t generate keywords hierarchy!')

        # Default start model could be chosen with command line parameter
        parser = argparse.ArgumentParser()
        parser.add_argument("--mesh", "-mesh",
                            help="Mesh .inp file",
                            type=str, default='mesh_baffle.inp')
        args = parser.parse_args()

        # Creates VTK widget with default ugrid, adds it to the form
        self.importMeshFromInp(args.mesh)

        # Actions
        self.treeView.doubleClicked.connect(self.treeViewDoubleClicked)
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
        self.frame.setLayout(self.vl) # apply layout: it will expand vtk_widget to the frame size

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
        self.textEdit.append('<p style=\'color:Black; margin:0px;\'>' + msg + '</p>')
        self.textEdit.moveCursor(QtGui.QTextCursor.End) # scroll text to the end
    def logging_error(self, msg):
        self.textEdit.append('<p style=\'color:Red; margin:0px;\'>ERROR! ' + msg + '</p>')
        self.textEdit.moveCursor(QtGui.QTextCursor.End) # scroll text to the end


    # Recursively generate treeView widget items based on DOM
    def generateTreeView(self):
        self.model = QtGui.QStandardItemModel()
        parent = self.model.invisibleRootItem() # top element in QTreeView
        self.addToTree(parent, self.dom.root.items) # pass root - group 'Model'
        self.treeView.setModel(self.model)
        self.treeView.expandAll() # expanded looks better
    def addToTree(self, parent, children):
        """
            parent is QtGui.QStandardItem
            children are items of ccx_dom.DOM object
        """
        for item in children:
            if (item.item_type == 'keyword') or (item.item_type == 'group'):
                tree_element = QtGui.QStandardItem(item.name)
                tree_element.setData(item)
                parent.appendRow(tree_element)

                # Draw keyword's children for its implementation
                for i in range(len(item.implementations)):
                    impl = item.implementations[i] # keyword implementation object
                    e = QtGui.QStandardItem(impl.name)
                    e.setData(impl)
                    tree_element.setText(item.name + ' (' + str(len(item.implementations)) + ')')
                    tree_element.appendRow(e)
                    self.addToTree(e, item.items)

                # Do not draw keyword's children if it doesn't have implementations
                if item.item_type == 'group':
                    self.addToTree(tree_element, item.items)


    # Import mesh and display it in the VTK widget
    def importMeshFromInp(self, file_name=None):
        if not file_name:
            file_name = QtWidgets.QFileDialog.getOpenFileName(self,\
                'Import .inp mesh', '', 'Input files (*.inp);;All Files (*)')[0]

        self.logging_info('Loading ' + file_name + '.')

        if file_name:
            # Parse mesh and transfer it to VTK
            mesh = INPParser.Mesh(file_name) # parse mesh
            points = vtk.vtkPoints()
            for n in mesh.nodes.keys(): # create VTK points from mesh nodes
                points.InsertPoint(n-1, mesh.nodes[n]) # node numbers should start from 0!
            self.ugrid = vtk.vtkUnstructuredGrid() # create empty grid in VTK
            self.ugrid.Allocate(len(mesh.elements)) # allocate memory fo all elements
            self.ugrid.SetPoints(points) # insert all points to the grid
            for e in mesh.elements.keys():
                vtk_element_type = INPParser.Mesh.convert_elem_type(mesh.types[e])
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


    # Double click on treeView item: edit the keyword via dialog
    def treeViewDoubleClicked(self, index):
        item = self.treeView.model().itemFromIndex(index) # treeView item, we obtain it from 'index'
        item = item.data() # now it is ccx_dom.group, ccx_dom.keyword or ccx_dom.implementation 

        # Only double clicking on ccx_dom.keyword creates dialog, not on ccx_dom.group
        if item.item_type == 'keyword':
            dialog = ccx_dialog.Dialog(item) # create dialog window and and pass ccx_dom.keyword object
            if dialog.exec_() == ccx_dialog.Dialog.Accepted: # if user pressed 'OK'

                # The generated piece of .inp code for the CalculiX input file
                INP_code = dialog.onOk()
                for line in INP_code.split('\n'):
                    self.logging_info(line) # show it

                # Create implementation object
                ccx_dom.implementation(item, INP_code)

                # Update treeView widget
                self.generateTreeView()

        # Click on keyword's implementation: show piece of INP_code
        elif item.item_type == 'implementation':
            for line in item.INP_code.split('\n'):
                self.logging_info(line) # show it


    # Menu Select->Nodes
    def selectNodes(self):
        self.clearSelection() # clear selection before new call
        style = ccx_select_style.nodes(self.renderer, self.window, self.textEdit)
        self.interactor.SetInteractorStyle(style)


    # Menu Select->Elements
    def selectElements(self):
        self.clearSelection() # clear selection before new call
        style = ccx_select_style.elements(self.renderer, self.window, self.ugrid, self.textEdit)
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
    window = CAE()
    window.show() # window.showMaximized() or window.show()
    sys.exit(app.exec_())
