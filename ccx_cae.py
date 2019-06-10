#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, UJV Rez, June 2019.
    Distributed under GNU General Public License, version 2.

    CalculiX CAE
    Main window application
"""


import sys, os, vtk, logging, argparse, math
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

        # Configure logging
        logging.info = self.log_widget.append
        logging.error = self.log_widget.append
        logging.info('Window created.')

        # Create VTK widget
        self.createVTKWidget()

        # Read CalculiX object model
        self.readObjectModel()

        # Default start model could be chosen with command line parameter
        parser = argparse.ArgumentParser()
        parser.add_argument("--mesh", "-mesh",
                            help="Mesh .inp file",
                            type=str, default='ccx_mesh_baffle.inp')
        args = parser.parse_args()

        # Creates VTK widget with default ugrid, adds it to the form
        self.importMeshFromInp(args.mesh)

        # Actions
        self.tree_view.doubleClicked.connect(self.treeViewDoubleClicked)
        self.actionImportMeshFromInp.triggered.connect(self.importMeshFromInp)
        self.actionCameraFitView.triggered.connect(self.cameraFitView)
        self.actionSelectNodes.triggered.connect(self.selectNodes)
        self.actionSelectElements.triggered.connect(self.selectElements)


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
        self.camera = self.renderer.GetActiveCamera()
        # (x,y,z) = ugrid.GetCenter()
        # if abs(x) > 1e+5: x = 0
        # if abs(y) > 1e+5: y = 0
        # if abs(z) > 1e+5: z = 0
        # self.camera.SetFocalPoint(x,y,z) # model's center
        # logging.info('Setting focal point to ' + str((x,y,z)))
        # logging.info('Data bounds are ' + str(ugrid.GetBounds()))
        logging.info('Camera\'s position is ' + str(self.camera.GetPosition()))
        logging.info('Camera\'s ViewUp is ' + str(self.camera.GetViewUp()))
        logging.info('Camera\'s distance is ' + str(self.camera.GetDistance()))
        logging.info('Camera\'s Roll is ' + str(self.camera.GetRoll()))
        logging.info('Camera\'s ViewAngle is ' + str(self.camera.GetViewAngle()))
        logging.info('Camera\'s ParallelScale is ' + str(self.camera.GetParallelScale()))
        logging.info('Camera\'s ClippingRange is ' + str(self.camera.GetClippingRange()))
        logging.info('Camera\'s WindowCenter is ' + str(self.camera.GetWindowCenter()))
        logging.info('Camera\'s orientation is ' + str(self.camera.GetOrientation())) # TODO camera is rotated

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
            self.ugrid = vtk.vtkUnstructuredGrid() # create empty grid in VTK
            self.ugrid.Allocate(len(mesh.elements)) # allocate memory fo all elements
            self.ugrid.SetPoints(points) # insert all points to the grid
            for e in mesh.elements.keys():
                vtk_element_type = Mesh.convert_elem_type(mesh.types[e])
                node_numbers = [n-1 for n in mesh.elements[e]] # list of nodes in the element: node numbers should start from 0!
                self.ugrid.InsertNextCell(vtk_element_type, len(node_numbers), node_numbers) # create VTK element

            self.mapper.SetInputData(self.ugrid) # ugrid is our mesh data
            self.cameraFitView() # reset camera

            logging.info('Rendering OK.')
            self.log_widget.ensureCursorVisible() # scroll text to the end


    # Automatically set up the camera
    def cameraFitView(self):
        # TODO rotate camera to initial position
        self.renderer.ResetCamera()
        self.window.Render()


    # Double click on tree_view item: create keyword in DOM
    def treeViewDoubleClicked(self, index):
        item = index.model().itemFromIndex(index)
        if item.text().startswith('*'):
            logging.info(item.text())
            self.log_widget.ensureCursorVisible() # scroll text to the end


    # Menu Select -> Nodes
    def selectNodes(self):
        point_picker = vtk.vtkPointPicker()
        self.interactor.SetPicker(point_picker)
        point_picker.AddObserver("StartPickEvent", self.pickPoint)
        point_picker.AddObserver("PickEvent", self.pickPoint)
        point_picker.AddObserver("EndPickEvent", self.pickPoint)


    def pickPoint(self, object, event):
        node_id = object.GetPointId()
        logging.info('Node ' + str(node_id))

        points = vtk.vtkPoints()
        points.InsertPoint(0, object.GetMapperPosition())

        vertex = vtk.vtkVertex()
        vertex.GetPointIds().SetId(0, 0)

        vertices = vtk.vtkCellArray()
        vertices.InsertNextCell(vertex)

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetVerts(vertices)

        selected_mapper = vtk.vtkPolyDataMapper()
        selected_mapper.SetInputData(polydata)

        actor = vtk.vtkActor()
        actor.SetMapper(selected_mapper)
        actor.GetProperty().SetPointSize(10)
        actor.GetProperty().SetColor(1.0, 0.0, 0.0)

        self.renderer.AddActor(actor)
        self.window.Render()


    # Menu Select -> Elements
    def selectElements(self):
        cell_picker = vtk.vtkCellPicker()
        self.interactor.SetPicker(cell_picker)
        cell_picker.AddObserver("PickEvent", self.pickCell)


    def pickCell(self, object, event):
        cell_id = object.GetCellId() # object is cell_picker
        logging.info('Element ' + str(cell_id))

        ids_to_hightlight = [cell_id, ]
        ids = vtk.vtkIdTypeArray() 
        ids.SetNumberOfComponents(1) 
        ids.Allocate(len(ids_to_hightlight)) 
        for i in ids_to_hightlight: 
            ids.InsertNextValue(i)

        selection_node = vtk.vtkSelectionNode() 
        selection_node.SetFieldType(vtk.vtkSelectionNode.CELL) # CELL POINT VERTEX
        selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
        selection_node.SetSelectionList(ids)

        selection = vtk.vtkSelection() 
        selection.AddNode(selection_node) 

        extract_selection = vtk.vtkExtractSelection() 
        extract_selection.SetInputData(0, self.ugrid) 
        extract_selection.SetInputData(1, selection) 
        extract_selection.Update() 

        grid_selected = vtk.vtkUnstructuredGrid() 
        grid_selected.ShallowCopy(extract_selection.GetOutput()) 

        selected_mapper = vtk.vtkDataSetMapper()
        selected_mapper.SetInputData(grid_selected) 

        actor = vtk.vtkActor()
        actor.SetMapper(selected_mapper)
        actor.GetProperty().SetColor(1.0, 0.0, 0.0)

        self.renderer.AddActor(actor)
        self.window.Render()


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
