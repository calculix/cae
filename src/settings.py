#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, June 2020
Distributed under GNU General Public License v3.0

Application's settings.
Attributes values are maintained in config/Settings_*.py.
User dialog form is config/SettingsDialog.xml - use Qt Designer to edit. """

# Standard modules
import os
import sys
import logging
import traceback

# External modules
from PyQt5 import QtWidgets, uic

# My modules
import path
import clean


# Session settings object used everywhere in the code
class Settings():

    # Read settings from file
    def __init__(self, p):
        self.p = p

        # Try to read settings file
        try:
            with open(self.p.settings, 'r') as f:
                exec(f.read())
            if len(self.__dict__) < 2:
                raise Exception

        # Apply default values
        except:

            # Windows
            if os.name=='nt':
                ext = '.exe'
                self.path_paraview = 'C:\\Program Files\\ParaView\\bin\\paraview.exe'
                self.path_editor = 'C:\\Windows\\System32\\notepad.exe'

            # Linux
            else:
                ext = ''
                self.path_paraview = '/opt/ParaView/bin/paraview'
                self.path_editor = '/usr/bin/gedit'

            self.path_ccx = os.path.join(self.p.bin, 'ccx' + ext)
            self.path_cgx = os.path.join(self.p.bin, 'cgx' + ext)
            self.start_model = os.path.join(self.p.examples, 'default.inp')
            self.logging_level = 'DEBUG'
            self.show_empty_keywords = True
            self.expanded = True
            self.start_cgx_by_default = True
            self.align_windows = True

    # Open dialog window and pass settings
    def open(self):
        self.__init__(self.p) # re-read settings from file
        sd = SettingsDialog(self.p, settings=self)

        # Warning about Cygwin DLLs
        if os.name=='nt':
            logging.warning('In Windows CCX/CGX binaries may not work if placed outside \'bin\' directory. They need Cygwin DLLs!')

        # Get response from dialog window
        if sd.exec(): # == 1 if user pressed 'OK'
            sd.save()
            self.__init__(self.p) # read settings from file
            logging.warning('For some settings to take effect application\'s restart may be needed.')

    # Automatic saving of current settings during the workflow
    def save(self):
        """ Pass values to dialog and save
        This method could produce redundant PyQt debug logging """
        sd = SettingsDialog(self.p, settings=self)
        sd.save()
    def save_bad(self):
        """ This method erases comments from the settings file """
        with open(self.p.settings, 'w') as f:
            for attr, value in self.__dict__.items():
                if attr == 'p':
                    continue
                if type(value) ==  str:
                    line = 'self.{} = \'{}\''.format(attr, value)
                else:
                    line = 'self.{} = {}'.format(attr, value)
                f.write(line + '\n\n')


# User dialog window with all setting attributes: menu File->Settings
class SettingsDialog(QtWidgets.QDialog):

    # Create dialog window
    def __init__(self, p, settings=None):
        self.p = p

        # Switch off logging
        hh = logging.getLogger().handlers
        logging.getLogger().handlers = []

        # Load UI form
        # Produces huge amount of redundant debug logs
        QtWidgets.QDialog.__init__(self)
        uic.loadUi(self.p.settings_xml, self) # load default settings

        # Switch on logging
        for h in hh:
            logging.getLogger().addHandler(h)

        # Push settings values to the form
        if settings:
            for attr, value in settings.__dict__.items():
                widget = self.findChild(QtWidgets.QCheckBox, attr)
                if widget is not None:
                    widget.setChecked(value)
                    continue

                widget = self.findChild(QtWidgets.QLineEdit, attr)
                if widget is not None:
                    widget.setText(value)
                    continue

                widget = self.findChild(QtWidgets.QComboBox, attr)
                if widget is not None:
                    widget.setCurrentText(value)
                    continue

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


# Run test
if __name__ == '__main__':
    clean.screen()
    logging.basicConfig(level=0, format='%(message)s')

    # Create application
    app = QtWidgets.QApplication(sys.argv)

    # Create and open settings window
    p = path.Path()
    s = Settings(p)
    s.save()
    s.open()

    # Clean cached files
    clean.cache(p.src)
