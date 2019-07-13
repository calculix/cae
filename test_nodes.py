import vtk

def main():
    pointSource = vtk.vtkPointSource()
    pointSource.SetNumberOfPoints(50)
    pointSource.Update()

    ids = vtk.vtkIdTypeArray()
    ids.SetNumberOfComponents(1)

    i = 10
    while i < 20:
        ids.InsertNextValue(i)
        i += 1

    selectionNode = vtk.vtkSelectionNode()
    selectionNode.SetFieldType(vtk.vtkSelectionNode.POINT)
    selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES)
    selectionNode.SetSelectionList(ids)

    selection = vtk.vtkSelection()
    selection.AddNode(selectionNode)

    extractSelection = vtk.vtkExtractSelection()
    extractSelection.SetInputConnection(0, pointSource.GetOutputPort())
    extractSelection.SetInputData(1, selection)
    extractSelection.Update()

    # In selection
    selected = vtk.vtkUnstructuredGrid()
    selected.ShallowCopy(extractSelection.GetOutput())

    selectionNode.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)
    extractSelection.Update()

    inputMapper = vtk.vtkDataSetMapper()
    inputMapper.SetInputConnection(pointSource.GetOutputPort())
    inputActor = vtk.vtkActor()
    inputActor.SetMapper(inputMapper)

    selectedMapper = vtk.vtkDataSetMapper()
    selectedMapper.SetInputData(selected)
    selectedActor = vtk.vtkActor()
    selectedActor.SetMapper(selectedMapper)
    selectedActor.GetProperty().SetColor(1.0, 0.6, 0.0)
    selectedActor.GetProperty().SetPointSize(7)

    window = vtk.vtkRenderWindow()
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    renderer = vtk.vtkRenderer()
    window.AddRenderer(renderer)

    renderer.AddActor(inputActor)
    renderer.AddActor(selectedActor)
    window.Render()
    interactor.Start()

if __name__ == '__main__':
    main()
