#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, UJV Rez, June 2019.
    Distributed under GNU General Public License, version 2.

    CalculiX CAE
    Main window application
"""


import sys, os, vtk, logging, argparse
from PyQt5 import QtWidgets, uic, QtGui
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from INPParser import Mesh
from ccx_style import MouseInteractorHighLightActor


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

        # Configure logging and write log
        logging.info = self.log_widget.append
        logging.error = self.log_widget.append
        logging.info('Started.')

        # Read CalculiX object model
        self.readObjectModel()

        # Default start model could be chosen with command line parameter
        parser = argparse.ArgumentParser()
        parser.add_argument("--mesh", "-mesh",
                            help="Mesh .inp file",
                            type=str, default='ccx_hex.inp')
        args = parser.parse_args()

        # Creates VTK widget with default ugrid, adds it to the form
        self.importMeshFromInp(args.mesh)

        # Actions
        self.actionImportMeshFromInp.triggered.connect(self.importMeshFromInp)
        self.actionResetCamera.triggered.connect(self.resetCamera)
        self.tree_view.doubleClicked.connect(self.treeViewDoubleClicked)


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
        self.renderer = vtk.vtkRenderer() # renderer renders into the window
        self.window = self.vtk_widget.GetRenderWindow()
        self.window.AddRenderer(self.renderer)
        interactor = self.window.GetInteractor() # captures mouse events

        # Add the custom style
        style = MouseInteractorHighLightActor()
        style.SetDefaultRenderer(self.renderer)
        interactor.SetInteractorStyle(style)

        # Add orientation axes
        # TODO crashes app after new mesh import
        # axes = vtk.vtkAxesActor()
        # axes.SetShaftTypeToCylinder()
        # axes.SetXAxisLabelText('X')
        # axes.SetYAxisLabelText('Y')
        # axes.SetZAxisLabelText('Z')
        # axes.SetTotalLength(1.0, 1.0, 1.0)
        # axes.SetCylinderRadius(0.5 * axes.GetCylinderRadius())
        # axes.SetConeRadius(1.025 * axes.GetConeRadius())
        # axes.SetSphereRadius(1.5 * axes.GetSphereRadius())
        # tprop = vtk.vtkTextProperty()
        # tprop.SetColor(0, 0, 0) # font color
        # axes.GetXAxisCaptionActor2D().SetCaptionTextProperty(tprop)
        # axes.GetYAxisCaptionActor2D().SetCaptionTextProperty(tprop)
        # axes.GetZAxisCaptionActor2D().SetCaptionTextProperty(tprop)
        # self.om1 = vtk.vtkOrientationMarkerWidget()
        # self.om1.SetOrientationMarker(axes)
        # self.om1.SetViewport(0, 0, 0.2, 0.2) # position lower left in the viewport
        # self.om1.SetInteractor(interactor)
        # self.om1.EnabledOn()
        # self.om1.InteractiveOn()

        # Mapper pushes the mesh data into the graphics library
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(ugrid) # ugrid is our mesh data

        # The actor is a grouping mechanism
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(0.2, 0.6, 0.8)

        # TODO rotate camera to initial position
        # actor.RotateX(30.0) # doesn't rotate axes
        # actor.RotateY(-45.0) # doesn't rotate axes

        # Add the actors to the renderer, set the background
        self.renderer.AddActor(actor)
        self.renderer.SetBackground(1, 1, 1)

        # Apply layout: it will expand vtk_widget to the frame size
        self.frame_widget.setLayout(self.vl)

        # Tune camera
        # camera = self.renderer.GetActiveCamera()
        # (x,y,z) = ugrid.GetCenter()
        # if abs(x) > 1e+5: x = 0
        # if abs(y) > 1e+5: y = 0
        # if abs(z) > 1e+5: z = 0
        # camera.SetFocalPoint(x,y,z) # model's center
        # logging.info('Setting focal point to ' + str((x,y,z)))
        # logging.info('Data bounds are ' + str(ugrid.GetBounds()))
        # logging.info('Camera\'s position is ' + str(camera.GetPosition()))
        # logging.info('Camera\'s ViewUp is ' + str(camera.GetViewUp()))
        # logging.info('Camera\'s distance is ' + str(camera.GetDistance()))
        # logging.info('Camera\'s Roll is ' + str(camera.GetRoll()))
        # logging.info('Camera\'s ViewAngle is ' + str(camera.GetViewAngle()))
        # logging.info('Camera\'s ParallelScale is ' + str(camera.GetParallelScale()))
        # logging.info('Camera\'s ClippingRange is ' + str(camera.GetClippingRange()))
        # logging.info('Camera\'s WindowCenter is ' + str(camera.GetWindowCenter()))
        # logging.info('Camera\'s orientation is ' + str(camera.GetOrientation())) # TODO camera is rotated

        # Render mesh
        self.renderer.ResetCamera() # automatically set up the camera based on the visible actors
        self.window.Render()

        # Initalize before an event loop
        # interactor.Initialize()

        # Start the event loop
        interactor.Start()

        logging.info('Rendering OK.')
        self.log_widget.ensureCursorVisible() # scroll text to the end


    # Read CalculiX keywords hierarchy
    def readObjectModel(self):
        self.model = QtGui.QStandardItemModel()
        tree_dict = {}
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

                    if level == 0:
                        tree_dict[0] = self.model.invisibleRootItem()
                    else:
                        tree_dict[level] = QtGui.QStandardItem(line.strip())
                        tree_dict[level-1].setChild(tree_dict[level-1].rowCount(), tree_dict[level])

                self.tree_view.setModel(self.model)
                self.tree_view.expandAll()
            logging.info('CalculiX object model generated.')
        except:
            logging.error('<span style=\'color: Red\'>Error reading keywords hierarchy!</span>')
        self.log_widget.ensureCursorVisible() # scroll text to the end


    # Import mesh and display it in the VTK widget
    def importMeshFromInp(self, file_name=None):
        if not file_name:
            file_name = QtWidgets.QFileDialog.getOpenFileName(self,\
                'Import .inp mesh', '', 'Input files (*.inp);;All Files (*)')[0]

        logging.info('Loading ' + file_name + '.')

        if file_name:
            # Parse mesh and transfer it to VTK
            mesh = Mesh(file_name) # parse mesh
            points = vtk.vtkPoints()
            for n in mesh.nodes.keys(): # create VTK points from mesh nodes
                points.InsertPoint(n-1, mesh.nodes[n]) # node numbers should start from 0!
            ugrid = vtk.vtkUnstructuredGrid() # create empty grid in VTK
            ugrid.Allocate(len(mesh.elements)) # allocate memory fo all elements
            ugrid.SetPoints(points) # insert all points to the grid
            for e in mesh.elements.keys():
                vtk_element_type = Mesh.convert_elem_type(mesh.types[e])
                node_numbers = [n-1 for n in mesh.elements[e]] # list of nodes in the element: node numbers should start from 0!
                ugrid.InsertNextCell(vtk_element_type, len(node_numbers), node_numbers) # create VTK element

            logging.info('Successfully loaded ' + file_name + '.')
            self.log_widget.ensureCursorVisible() # scroll text to the end
            self.drawVTKWidget(ugrid)


    # Automatically set up the camera based on the visible actors + render
    def resetCamera(self):
        self.renderer.ResetCamera()
        self.window.Render()


    # Double click on tree_view item
    def treeViewDoubleClicked(self, index):
        item = index.model().itemFromIndex(index)
        if item.text().startswith('*'):
            logging.info(item.text())
            self.log_widget.ensureCursorVisible() # scroll text to the end


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
    window.showMaximized() # or it could be window.show()
    sys.exit(app.exec_())
