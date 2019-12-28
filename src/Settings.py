# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, December 2019
    Distributed under GNU General Public License v3.0

    Application's settings.
    Attributes values are maintained through the settings ENV file.
    User dialog form is maintained through the settings.XML file.
    You may edit settings.XML from Qt Designer.
"""


from Path import Path
import os, sys, logging
from gui.SettingsDialog import SettingsDialog


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
            print('Error reading ENV settings. Applying default values.')
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
            self.show_maximized = True
            self.show_empty_keywords = True
            self.expanded = True
            self.vtk_show_axes = True
            self.vtk_parallel_view = True
            self.show_help = True
            self.show_vtk = True


    # Open dialog window and pass settings
    def open(self):
        self.__init__() # re-read settings from file
        dialog = SettingsDialog(settings=self)

        # Warning about Cygwin DLLs
        if os.name=='nt':
            logging.warning('In Windows ccx binary may not work if placed outside \'bin\' directory. It needs Cygwin DLLs!')

        # Get response from dialog window
        if dialog.exec(): # == 1 if user pressed 'OK'
            dialog.save()
            self.__init__() # read settings from file
            logging.warning('For some settings to take effect application\'s restart may be needed.')


    # Automatic save current settings during the workflow
    def save(self):
        # Pass values to dialog and save
        settings = SettingsDialog(settings=self)
        settings.save()



# Test module
if __name__ == '__main__':

    from PyQt5 import QtWidgets
    from clean import cleanCache

    # Create application
    app = QtWidgets.QApplication(sys.argv)

    # Create and open settings window
    settings = Settings()
    settings.open()

    # Clean cached files
    p = Path() # calculate absolute paths
    cleanCache(p.src)
