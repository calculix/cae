# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    Methods to work with VTK widget.
"""

from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import vtk, ccx_select_style


class VTK:


    # Create empty VTK widget: called once during startup
    def __init__(self, CAE):

        self.CAE = CAE

        # Create the graphics structure
        self.widget = QVTKRenderWindowInteractor()
        self.renderer = vtk.vtkRenderer() # renderer renders into the window
        self.window = self.widget.GetRenderWindow()
        self.window.AddRenderer(self.renderer)
        self.interactor = self.window.GetInteractor() # captures mouse events

        # Add the custom style
        self.interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())

        # Draw a point to avoid error logs on startup
        p = vtk.vtkPointSource()
        p.SetNumberOfPoints(1)

        # Mapper pushes the mesh data into the graphics library
        self.mapper = vtk.vtkDataSetMapper()
        self.mapper.SetInputConnection(p.GetOutputPort())

        # The actor is a grouping mechanism
        self.actor = vtk.vtkActor()
        self.actor.SetMapper(self.mapper)
        self.actor.GetProperty().SetColor(0.25, 0.60, 0.25) # mesh color

        # Add the actors to the renderer, set the background
        self.renderer.AddActor(self.actor)
        self.renderer.SetBackground(1, 1, 1)

        # Add orientation axes
        self.axes = vtk.vtkAxesActor()
        self.axes.SetShaftTypeToCylinder()
        self.axes.SetXAxisLabelText('X')
        self.axes.SetYAxisLabelText('Y')
        self.axes.SetZAxisLabelText('Z')
        self.axes.SetTotalLength(1.0, 1.0, 1.0)
        self.axes.SetCylinderRadius(1.0 * self.axes.GetCylinderRadius())
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

        # Actions
        CAE.actionSelectionNodes.triggered.connect(self.actionSelectionNodes)
        CAE.actionSelectionElements.triggered.connect(self.actionSelectionElements)
        CAE.actionSelectionClear.triggered.connect(self.actionSelectionClear)
        CAE.actionViewParallel.triggered.connect(self.actionViewParallel)
        CAE.actionViewPerspective.triggered.connect(self.actionViewPerspective)
        CAE.actionViewFront.triggered.connect(self.actionViewFront)
        CAE.actionViewBack.triggered.connect(self.actionViewBack)
        CAE.actionViewTop.triggered.connect(self.actionViewTop)
        CAE.actionViewBottom.triggered.connect(self.actionViewBottom)
        CAE.actionViewLeft.triggered.connect(self.actionViewLeft)
        CAE.actionViewRight.triggered.connect(self.actionViewRight)
        CAE.actionViewIso.triggered.connect(self.actionViewIso)
        CAE.actionViewFit.triggered.connect(self.actionViewFit)
        CAE.actionViewWireframe.triggered.connect(self.actionViewWireframe)
        CAE.actionViewSurface.triggered.connect(self.actionViewSurface)
        CAE.actionViewSurfaceWithEdges.triggered.connect(self.actionViewSurfaceWithEdges)

        # Apply some action functions by default
        self.actionViewParallel()
        self.actionViewSurfaceWithEdges()


    # View toolbar

    def actionViewParallel(self):
        camera = self.renderer.GetActiveCamera() # get camera
        camera.SetParallelProjection(True)
        self.window.Render() # render updated view

    def actionViewPerspective(self):
        camera = self.renderer.GetActiveCamera() # get camera
        camera.SetParallelProjection(False)
        self.window.Render() # render updated view

    def actionViewFront(self):
        camera = self.renderer.GetActiveCamera() # get camera
        camera.SetPosition(0, 0, 1) # camera's positions
        camera.SetViewUp(0, 1, 0) # which axis will be set as vertical
        camera.SetFocalPoint(0, 0, 0) # set focal foint to center
        self.renderer.ResetCamera() # fit view
        self.window.Render() # render updated view

    def actionViewBack(self):
        camera = self.renderer.GetActiveCamera() # get camera
        camera.SetPosition(0, 0, -1) # camera's positions
        camera.SetViewUp(0, 1, 0) # which axis will be set as vertical
        camera.SetFocalPoint(0, 0, 0) # set focal foint to center
        self.renderer.ResetCamera() # fit view
        self.window.Render() # render updated view

    def actionViewTop(self):
        camera = self.renderer.GetActiveCamera() # get camera
        camera.SetPosition(0, 1, 0) # camera's positions
        camera.SetViewUp(0, 0, -1) # which axis will be set as vertical
        camera.SetFocalPoint(0, 0, 0) # set focal foint to center
        self.renderer.ResetCamera() # fit view
        self.window.Render() # render updated view

    def actionViewBottom(self):
        camera = self.renderer.GetActiveCamera() # get camera
        camera.SetPosition(0, -1, 0) # camera's positions
        camera.SetViewUp(0, 0, 1) # which axis will be set as vertical
        camera.SetFocalPoint(0, 0, 0) # set focal foint to center
        self.renderer.ResetCamera() # fit view
        self.window.Render() # render updated view

    def actionViewLeft(self):
        camera = self.renderer.GetActiveCamera() # get camera
        camera.SetPosition(-1, 0, 0) # camera's positions
        camera.SetViewUp(0, 1, 0) # which axis will be set as vertical
        camera.SetFocalPoint(0, 0, 0) # set focal foint to center
        self.renderer.ResetCamera() # fit view
        self.window.Render() # render updated view

    def actionViewRight(self):
        camera = self.renderer.GetActiveCamera() # get camera
        camera.SetPosition(1, 0, 0) # camera's positions
        camera.SetViewUp(0, 1, 0) # which axis will be set as vertical
        camera.SetFocalPoint(0, 0, 0) # set focal foint to center
        self.renderer.ResetCamera() # fit view
        self.window.Render() # render updated view

    def actionViewIso(self):
        camera = self.renderer.GetActiveCamera() # get camera
        camera.SetPosition(1.3, 1.0, 1.3) # camera's positions
        camera.SetViewUp(0, 1, 0) # which axis will be set as vertical
        camera.SetFocalPoint(0, 0, 0) # set focal foint to center
        self.renderer.ResetCamera() # fit view
        self.window.Render() # render updated view

    def actionViewFit(self):
        camera = self.renderer.GetActiveCamera() # get camera
        camera.SetFocalPoint(self.mapper.GetCenter()) # set focal foint to center
        self.renderer.ResetCamera() # fit view
        self.window.Render() # render updated view

        # Some logs for debugging
        # self.CAE.logger.info('Camera\'s focal point is ' + str(camera.GetFocalPoint()))
        # self.CAE.logger.info('Camera\'s position is ' + str(camera.GetPosition()))
        # self.CAE.logger.info('Camera\'s ViewUp is ' + str(camera.GetViewUp()))
        # self.CAE.logger.info('Camera\'s distance is ' + str(camera.GetDistance()))
        # self.CAE.logger.info('Camera\'s Roll is ' + str(camera.GetRoll()))
        # self.CAE.logger.info('Camera\'s ViewAngle is ' + str(camera.GetViewAngle()))
        # self.CAE.logger.info('Camera\'s ParallelScale is ' + str(camera.GetParallelScale()))
        # self.CAE.logger.info('Camera\'s ClippingRange is ' + str(camera.GetClippingRange()))
        # self.CAE.logger.info('Camera\'s WindowCenter is ' + str(camera.GetWindowCenter()))
        # self.CAE.logger.info('Camera\'s orientation is ' + str(camera.GetOrientation()))

    def actionViewWireframe(self):
        self.actor.GetProperty().SetRepresentationToWireframe()
        self.window.Render() # render updated view

    def actionViewSurface(self):
        self.actor.GetProperty().EdgeVisibilityOff()
        self.actor.GetProperty().SetRepresentationToSurface()
        self.window.Render() # render updated view

    def actionViewSurfaceWithEdges(self):
        self.actor.GetProperty().EdgeVisibilityOn()
        # meshFreeSurfaceFilter = vtk.vtkDataSetSurfaceFilter() 
        # meshFreeSurfaceFilter.SetInputData(self.mesh) 
        # meshFreeSurfaceFilter.Update()
        self.actor.GetProperty().SetRepresentationToSurface()
        self.window.Render() # render updated view


    # Menu Select

    def actionSelectionNodes(self):
        self.actionSelectionClear() # clear selection before new call
        style = ccx_select_style.nodes(self.renderer, self.window, self.CAE.textEdit)
        self.interactor.SetInteractorStyle(style)

    def actionSelectionElements(self):
        self.actionSelectionClear() # clear selection before new call
        style = ccx_select_style.elements(self.renderer, self.window, self.CAE.textEdit)
        self.interactor.SetInteractorStyle(style)

    def actionSelectionClear(self):
        self.renderer.RemoveAllViewProps()
        self.renderer.AddActor(self.actor)
        self.window.Render()
        # self.CAE.logger.info('Clear selection.')


    # Highlight node sets, element sets and surfaces

    def highlight(self, _set, field_type):
        # Get our mesh
        ugrid = vtk.vtkUnstructuredGrid()
        ugrid.ShallowCopy(self.mapper.GetInput())

        # From _set get IDs of the selected items
        IDs = vtk.vtkIdTypeArray()
        IDs.SetNumberOfComponents(1)
        for i in _set:
            IDs.InsertNextValue(i-1) # in VTK nodes start from 0

        selectionNode = vtk.vtkSelectionNode()
        selectionNode.SetFieldType(field_type) # select nodes/cells by index
        selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES)
        selectionNode.SetSelectionList(IDs)

        selection = vtk.vtkSelection()
        selection.AddNode(selectionNode)

        extractSelection = vtk.vtkExtractSelection()
        extractSelection.SetInputData(0, ugrid)
        extractSelection.SetInputData(1, selection)
        extractSelection.Update()

        # Construct new mesh from selected nodes/cells
        selected = vtk.vtkUnstructuredGrid()
        selected.ShallowCopy(extractSelection.GetOutput())

        selectedMapper = vtk.vtkDataSetMapper()
        selectedMapper.SetInputData(selected)

        selectedActor = vtk.vtkActor()
        selectedActor.SetMapper(selectedMapper)
        selectedActor.GetProperty().SetColor(1.0, 0.6, 0.0) # color of selected nodes/cells
        selectedActor.GetProperty().SetPointSize(7) # size of selected nodes (if any)

        self.renderer.AddActor(selectedActor)
        self.window.Render() # update window

    def highlightSURFACE(self, _set):
        # Get our mesh
        ugrid = vtk.vtkUnstructuredGrid()
        ugrid.ShallowCopy(self.mapper.GetInput())

        # From _set get IDs of the selected items
        types = []
        faces = vtk.vtkCellArray()
        points = vtk.vtkPoints() # all points of the ugrid
        for cell, F in _set:
            cell = ugrid.GetCell(int(cell) - 1)
            F = cell.GetFace(int(F[1:]) - 1)
            pp = cell.GetPoints()
            for i in range(pp.GetNumberOfPoints()):
                points.InsertNextPoint(pp.GetPoint(i))
            types.append(F.GetCellType())
            faces.InsertNextCell(F)

        # Construct new mesh from selected nodes/cells
        selected = vtk.vtkUnstructuredGrid()
        selected.SetPoints(points)
        selected.SetCells(types, faces)

        selectedMapper = vtk.vtkDataSetMapper()
        selectedMapper.SetInputData(selected)
        selectedActor = vtk.vtkActor()
        selectedActor.SetMapper(selectedMapper)
        selectedActor.GetProperty().SetColor(1.0, 0.6, 0.0) # mesh color

        self.renderer.AddActor(selectedActor)
        self.window.Render() # update window

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
