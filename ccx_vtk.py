# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    Methods to work with VTK widget.
"""

from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import vtk, ccx_cae_log, ccx_mesh
# import ccx_vtk_select_style



class VTK:


    # Create empty VTK widget: called once during startup
    def __init__(self):
        self.bounds = [-1,1]*3

        # Create the graphics structure
        self.widget = QVTKRenderWindowInteractor()
        self.renderer = vtk.vtkRenderer() # renderer renders into the window
        self.window = self.widget.GetRenderWindow()
        self.window.AddRenderer(self.renderer)
        self.interactor = self.window.GetInteractor() # captures mouse events
        self.camera = vtk.vtkCamera()
        self.renderer.SetActiveCamera(self.camera)

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
        self.renderer.SetBackground(1, 1, 1) # white color

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

        # Apply some action functions by default
        self.actionViewParallel()
        self.actionViewSurfaceWithEdges()


    # Generate VTK unstructured grid from ccx_mesh object
    def mesh2ugrid(self, mesh):
        msg_list = [] # list of messages for logger
        ugrid = None

        self.bounds = mesh.bounds

        try:
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

            msg_text = 'VTK ugrid is built.'
            msg = ccx_cae_log.msg(ccx_cae_log.msgType.INFO, msg_text)
            msg_list.append(msg)

        except:
            msg_text = 'Can\'t build VTK ugrid.'
            msg = ccx_cae_log.msg(ccx_cae_log.msgType.INFO, msg_text)
            msg_list.append(msg)

        return msg_list, ugrid


    # View toolbar

    def actionViewParallel(self):
        self.camera.SetParallelProjection(True)
        self.window.Render() # render updated view

    def actionViewPerspective(self):
        self.camera.SetParallelProjection(False)
        self.window.Render() # render updated view

    def actionViewFront(self):
        self.camera.SetPosition(0, 0, 1) # camera's positions
        self.camera.SetViewUp(0, 1, 0) # which axis will be set as vertical
        self.camera.SetFocalPoint(0, 0, 0) # set focal foint to center
        self.actionViewFit()

    def actionViewBack(self):
        self.camera.SetPosition(0, 0, -1) # camera's positions
        self.camera.SetViewUp(0, 1, 0) # which axis will be set as vertical
        self.camera.SetFocalPoint(0, 0, 0) # set focal foint to center
        self.actionViewFit()

    def actionViewTop(self):
        self.camera.SetPosition(0, 1, 0) # camera's positions
        self.camera.SetViewUp(0, 0, -1) # which axis will be set as vertical
        self.camera.SetFocalPoint(0, 0, 0) # set focal foint to center
        self.actionViewFit()

    def actionViewBottom(self):
        self.camera.SetPosition(0, -1, 0) # camera's positions
        self.camera.SetViewUp(0, 0, 1) # which axis will be set as vertical
        self.camera.SetFocalPoint(0, 0, 0) # set focal foint to center
        self.actionViewFit()

    def actionViewLeft(self):
        self.camera.SetPosition(-1, 0, 0) # camera's positions
        self.camera.SetViewUp(0, 1, 0) # which axis will be set as vertical
        self.camera.SetFocalPoint(0, 0, 0) # set focal foint to center
        self.actionViewFit()

    def actionViewRight(self):
        self.camera.SetPosition(1, 0, 0) # camera's positions
        self.camera.SetViewUp(0, 1, 0) # which axis will be set as vertical
        self.camera.SetFocalPoint(0, 0, 0) # set focal foint to center
        self.actionViewFit()

    def actionViewIso(self):
        self.camera.SetPosition(1.3, 1.0, 1.3) # camera's positions
        self.camera.SetViewUp(0, 1, 0) # which axis will be set as vertical
        self.camera.SetFocalPoint(0, 0, 0) # set focal foint to center
        self.actionViewFit()

    def actionViewFit(self):
        self.renderer.ResetCamera(self.bounds) # avoid camera flying to infinity
        self.window.Render() # render updated view

        # Some logs for debugging
        self.msg_list = [] # list of messages for logger
        messages = [
            'Camera\'s focal point is ' + str(self.camera.GetFocalPoint()),
            'Camera\'s position is ' + str(self.camera.GetPosition()),
            # 'Camera\'s ViewUp is ' + str(self.camera.GetViewUp()),
            'Camera\'s distance is {:.1f}'.format(self.camera.GetDistance()),
            # 'Camera\'s Roll is {:.1f}'.format(self.camera.GetRoll()),
            # 'Camera\'s ViewAngle is {:.1f}'.format(self.camera.GetViewAngle()),
            'Camera\'s ParallelScale is {:.1f}'.format(self.camera.GetParallelScale()),
            'Camera\'s ClippingRange is ' + str(self.camera.GetClippingRange()),
            # 'Camera\'s WindowCenter is ' + str(self.camera.GetWindowCenter()),
            # 'Camera\'s orientation is ' + str(self.camera.GetOrientation())
        ]
        for msg_text in messages:
            msg = ccx_cae_log.msg(ccx_cae_log.msgType.INFO, msg_text)
            self.msg_list.append(msg)

    def actionViewWireframe(self):
        self.actor.GetProperty().SetRepresentationToWireframe()
        self.window.Render() # render updated view

    def actionViewSurface(self):
        self.actor.GetProperty().EdgeVisibilityOff()
        self.actor.GetProperty().SetRepresentationToSurface()
        self.window.Render() # render updated view

    def actionViewSurfaceWithEdges(self):
        self.actor.GetProperty().EdgeVisibilityOn()
        self.actor.GetProperty().SetRepresentationToSurface()
        self.window.Render() # render updated view


    # Menu Select

    # def actionSelectionNodes(self):
    #     self.actionSelectionClear() # clear selection before new call
    #     style = ccx_vtk_select_style.nodes(self.renderer, self.window, self.CAE.textEdit)
    #     self.interactor.SetInteractorStyle(style)

    # def actionSelectionElements(self):
    #     self.actionSelectionClear() # clear selection before new call
    #     style = ccx_vtk_select_style.elements(self.renderer, self.window, self.CAE.textEdit)
    #     self.interactor.SetInteractorStyle(style)

    def actionSelectionClear(self):
        self.renderer.RemoveAllViewProps()
        self.renderer.AddActor(self.actor)
        # self.actor.SetVisibility(True)
        self.window.Render()


    # Highlight node sets, element sets and surfaces
    # TODO error: create new set and highlight it
    def highlight(self, _set, field_type):
        # Clear selection
        self.actionSelectionClear()

        # Get our mesh
        ugrid = self.mapper.GetInput()

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
        # Clear selection
        self.actionSelectionClear()

        # Get our mesh
        ugrid = self.mapper.GetInput()

        # From _set get IDs of the selected items
        types = []
        faces = vtk.vtkCellArray()
        points = vtk.vtkPoints() # all points of the ugrid
        point_num = 0
        for cell_num, face_num in _set:
            cell = ugrid.GetCell(int(cell_num) - 1)
            face_num = int(face_num[1:]) # face number from .inp file: S1=1, S2=2 etc.

            if cell.GetCellType() in [12,25]: # hex
                face_num = [4,5,2,1,3,0][face_num - 1]
                face = cell.GetFace(face_num) # vtkCell made from face: quad or triangle
            elif cell.GetCellType() in [10,24]: # tet
                face_num = [3,0,1,2][face_num - 1]
                face = cell.GetFace(face_num) # vtkCell made from face: quad or triangle
            elif cell.GetCellType() == 13: # wedge
                face_num = face_num - 1
                face = cell.GetFace(face_num) # vtkCell made from face: quad or triangle
            elif cell.GetCellType() in [5,22,9,23]: # shells
                face_num = face_num - 1
                face = cell
            else:
                return

            _type = face.GetCellType() # 5, 9, 22, or 23

            # Linearize highlighted element
            nop = face.GetNumberOfPoints() # always even
            if not cell.IsLinear():
                nop = int(nop/2) # take only corners for second order elements
            if _type in [9, 23]: # quad
                facet = vtk.vtkQuad()
                types.append(9) # linear quad
            elif _type in [5, 22]: # triangle
                facet = vtk.vtkTriangle()
                types.append(5) # linear triangle
            else:
                return

            for i in range(nop):
                p = face.GetPoints().GetPoint(i)
                points.InsertNextPoint(p)
                facet.GetPointIds().SetId(i, point_num)
                point_num += 1
            faces.InsertNextCell(facet)

        # Construct new mesh from selected nodes/cells
        selected = vtk.vtkUnstructuredGrid()
        selected.SetPoints(points)
        selected.SetCells(types, faces)

        selectedMapper = vtk.vtkDataSetMapper()
        selectedMapper.SetInputData(selected)
        selectedActor = vtk.vtkActor()
        selectedActor.SetMapper(selectedMapper)
        selectedActor.GetProperty().SetColor(1.0, 0.6, 0.0) # color of selected nodes/cells
        selectedActor.GetProperty().SetPointSize(7) # size of selected nodes (if any)
        # selectedActor.GetProperty().EdgeVisibilityOn()

        # self.actor.SetVisibility(False)

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
