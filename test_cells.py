import vtk

def main():
    sphereSource = vtk.vtkSphereSource()
    sphereSource.Update()

    ids = vtk.vtkIdTypeArray()
    ids.SetNumberOfComponents(1)

    # Specify that we want to extract cells 10 through 19
    i = 10
    while i < 20:
        ids.InsertNextValue(i)
        i += 1

    selectionNode = vtk.vtkSelectionNode()
    selectionNode.SetFieldType(vtk.vtkSelectionNode.CELL)
    selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES)
    selectionNode.SetSelectionList(ids)

    selection = vtk.vtkSelection()
    selection.AddNode(selectionNode)

    extractSelection = vtk.vtkExtractSelection()
    extractSelection.SetInputConnection(0, sphereSource.GetOutputPort())
    extractSelection.SetInputData(1, selection)
    extractSelection.Update()

    # In selection
    selected = vtk.vtkUnstructuredGrid()
    selected.ShallowCopy(extractSelection.GetOutput())

    inputMapper = vtk.vtkDataSetMapper()
    inputMapper.SetInputConnection(sphereSource.GetOutputPort())
    inputActor = vtk.vtkActor()
    inputActor.SetMapper(inputMapper)

    selectedMapper = vtk.vtkDataSetMapper()
    selectedMapper.SetInputData(selected)

    selectedActor = vtk.vtkActor()
    selectedActor.SetMapper(selectedMapper)
    selectedActor.GetProperty().SetColor(1.0, 0.6, 0.0)


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