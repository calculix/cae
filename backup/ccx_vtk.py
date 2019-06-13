# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, UJV Rez, June 2019.
    Distributed under GNU General Public License, version 2.

    Methods to work with VTK widget
"""


from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import vtk, ccx_select_style, ccx_log


class VTK:


    # Create empty VTK widget: called once during startup
    def __init__(self, textEdit):

        # Configure logging
        self.textEdit = textEdit
        self.logger = ccx_log.logger(self.textEdit)

        # Create the graphics structure
        self.widget = QVTKRenderWindowInteractor()
        self.renderer = vtk.vtkRenderer() # renderer renders into the window
        self.window = self.widget.GetRenderWindow()
        self.window.AddRenderer(self.renderer)
        self.interactor = self.window.GetInteractor() # captures mouse events

        # Add the custom style
        self.interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())

        # Draw a point to avoid error logging on startup
        p = vtk.vtkPointSource()
        p.SetNumberOfPoints(1)

        # Mapper pushes the mesh data into the graphics library
        self.mapper = vtk.vtkDataSetMapper()
        self.mapper.SetInputConnection(p.GetOutputPort())

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
            self.logger.info('Camera\'s focal point is ' + str(self.camera.GetFocalPoint()))
            self.logger.info('Camera\'s position is ' + str(self.camera.GetPosition()))
            self.logger.info('Camera\'s ViewUp is ' + str(self.camera.GetViewUp()))
            self.logger.info('Camera\'s distance is ' + str(self.camera.GetDistance()))
            self.logger.info('Camera\'s Roll is ' + str(self.camera.GetRoll()))
            self.logger.info('Camera\'s ViewAngle is ' + str(self.camera.GetViewAngle()))
            self.logger.info('Camera\'s ParallelScale is ' + str(self.camera.GetParallelScale()))
            self.logger.info('Camera\'s ClippingRange is ' + str(self.camera.GetClippingRange()))
            self.logger.info('Camera\'s WindowCenter is ' + str(self.camera.GetWindowCenter()))
            self.logger.info('Camera\'s orientation is ' + str(self.camera.GetOrientation())) # TODO camera is rotated
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


    # Automatically set up the camera
    def cameraFitView(self):
        # TODO rotate camera to initial position
        self.renderer.ResetCamera()
        self.window.Render()


    # Menu Select->Nodes
    # TODO Should be active after mesh import
    def selectNodes(self):
        self.clearSelection() # clear selection before new call
        style = ccx_select_style.nodes(self.renderer, self.window, self.textEdit)
        self.interactor.SetInteractorStyle(style)


    # Menu Select->Elements
    # TODO Should be active after mesh import
    def selectElements(self):
        self.clearSelection() # clear selection before new call
        style = ccx_select_style.elements(self.renderer, self.window, self.textEdit)
        self.interactor.SetInteractorStyle(style)


    # Menu Select->Clear selection
    # TODO Should be active after mesh import
    def clearSelection(self):
        self.renderer.RemoveAllViewProps()
        self.renderer.AddActor(self.actor)
        self.window.Render()

        self.logger.info('Clear selection.')


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
