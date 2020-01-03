#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, December 2019
    Distributed under GNU General Public License v3.0

    Application's settings.
    Attributes values are maintained in config/Settings_linux.env.
    User dialog form is config/SettingsDialog.xml - use Qt Designer to edit.
"""


from path import Path
import os, sys, logging
from PyQt5 import QtWidgets, uic
import clean


# Session settings object used everywhere in the code
class Settings():


    # Read settings from file
    def __init__(self):
        self.p = Path() # calculate absolute paths

        # Try to read settings file
        try:
            with open(self.p.settings, 'r') as f:
                exec(f.read())
            if len(self.__dict__) < 2:
                raise Exception

        # Apply default values
        except:
            self.path_start_model = os.path.join(self.p.examples, 'default.inp')

            # Windows
            if os.name=='nt':
                self.path_ccx = os.path.join(self.p.bin, 'ccx_2.16_MT.exe')
                self.path_cgx = 'C:\\cgx.exe'
                self.path_paraview = 'C:\\Program Files\\ParaView\\bin\\paraview.exe'
                self.path_editor = 'C:\\Windows\\System32\\notepad.exe'

            # Linux
            else:
                self.path_ccx = os.path.join(self.p.bin, 'ccx_2.16_MT')
                self.path_cgx = '/usr/local/bin/cgx'
                self.path_paraview = '/opt/ParaView/bin/paraview'
                self.path_editor = '/snap/bin/code'

            self.logging_level = 'INFO'
            self.vtk_view = 'WithEdges'
            self.show_maximized = False
            self.show_empty_keywords = True
            self.expanded = True
            self.vtk_show_axes = True
            self.vtk_parallel_view = True
            self.show_help = True
            self.show_vtk = True


    # Open dialog window and pass settings
    def open(self):
        self.__init__() # re-read settings from file
        sd = SettingsDialog(settings=self)

        # Warning about Cygwin DLLs
        if os.name=='nt':
            logging.warning('In Windows ccx binary may not work if placed outside \'bin\' directory. It needs Cygwin DLLs!')

        # Get response from dialog window
        if sd.exec(): # == 1 if user pressed 'OK'
            sd.save()
            self.__init__() # read settings from file
            logging.warning('For some settings to take effect application\'s restart may be needed.')


    # Automatic save current settings during the workflow
    def save(self):
        # Pass values to dialog and save
        sd = SettingsDialog(settings=self)
        sd.save()



# User dialog window with all setting attributes: menu File->Settings
class SettingsDialog(QtWidgets.QDialog):


    # Create dialog window
    def __init__(self, settings=None):

        # Load UI form
        QtWidgets.QDialog.__init__(self)
        self.p = Path() # calculate absolute paths
        uic.loadUi(self.p.settings_xml, self) # load default settings from

        # Push settings values to the form
        if settings:
            for attr, value in settings.__dict__.items():
                try:
                    widget = self.findChild(QtWidgets.QCheckBox, attr)
                    widget.setChecked(value)
                except:
                    try:
                        widget = self.findChild(QtWidgets.QLineEdit, attr)
                        widget.setText(value)
                    except:
                        try:
                            widget = self.findChild(QtWidgets.QComboBox, attr)
                            widget.setCurrentText(value)
                        except:
                            pass


    # Save settings updated via or passed to dialog
    def save(self):
        with open(self.p.settings, 'w') as f:

            # Iterate over class attributes
            for attr, value in self.__dict__.items():
                class_name = value.__class__.__name__
                if class_name in ['QCheckBox', 'QLineEdit', 'QComboBox']:

                    # Write settings to file
                    if class_name == 'QCheckBox':
                        text = str(value.isChecked())
                    elif class_name == 'QLineEdit':
                        text = value.text()
                        text = '\'' + self.p.abspath(text) + '\'' # covert path to absolute
                        if '\\' in text: # reconstruct path for Windows
                            text = '\\\\'.join(text.split('\\'))
                        value = self.__dict__['label_' + attr]
                    elif class_name == 'QComboBox':
                        text = '\'' + value.currentText() + '\''
                        value = self.__dict__['label_' + attr]

                    line = 'self.{} = {}'.format(attr, text)
                    comment = value.text()
                    f.write('# ' + comment + '\n')
                    f.write(line + '\n\n')



# Test module
if __name__ == '__main__':

    # Create application
    app = QtWidgets.QApplication(sys.argv)

    # Create and open settings window
    settings = Settings()
    settings.open()

    # Clean cached files
    p = Path() # calculate absolute paths
    clean.cache(p.src)
