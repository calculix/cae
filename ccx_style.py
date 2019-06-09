# -*- coding: utf-8 -*-

import vtk, logging


class MouseInteractorHighLightActor(vtk.vtkInteractorStyleTrackballCamera):

    def __init__(self, parent=None):
        self.AddObserver('MiddleButtonPressEvent', self.middle_button_press_event)
        self.AddObserver('LeftButtonPressEvent', self.leftButtonPressEvent)
        self.AddObserver('RightButtonPressEvent', self.rightButtonPressEvent)
        self.LastPickedActor = None
        self.LastPickedProperty = vtk.vtkProperty()


    def leftButtonPressEvent(self, obj, event):
        clickPos = self.GetInteractor().GetEventPosition()

        picker = vtk.vtkPropPicker()
        picker.Pick(clickPos[0], clickPos[1], 0, self.GetDefaultRenderer())

        # Get the new
        self.NewPickedActor = picker.GetActor()

        # If something was selected
        if self.NewPickedActor:
            # If we picked something before, reset its property
            if self.LastPickedActor:
                self.LastPickedActor.GetProperty().DeepCopy(self.LastPickedProperty)

            # Save the property of the picked actor so that we can
            # restore it next time
            self.LastPickedProperty.DeepCopy(self.NewPickedActor.GetProperty())
            # Highlight the picked actor by changing its properties
            self.NewPickedActor.GetProperty().SetColor(1,0,0)
            self.NewPickedActor.GetProperty().SetDiffuse(1.0)
            self.NewPickedActor.GetProperty().SetSpecular(0.0)

            # TODO it selects all cells, not only one

            # Save the last picked actor
            self.LastPickedActor = self.NewPickedActor

        logging.info('Left button pressed')
        # self.log_widget.ensureCursorVisible() # scroll text to the end

        self.OnLeftButtonDown()
        return


    def middle_button_press_event(self, obj, event):
        logging.info('Middle button pressed')
        # self.log_widget.ensureCursorVisible() # scroll text to the end

        self.OnMiddleButtonDown()
        return


    def rightButtonPressEvent(self, obj, event):
        logging.info('Right button pressed')
        # self.log_widget.ensureCursorVisible() # scroll text to the end

        self.OnRightButtonDown()
        return