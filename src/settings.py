#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Application settings.
Attributes values are maintained in config/Settings_*.py.
User dialog form is config/SettingsDialog.xml - use Qt Designer to edit.
"""

# Standard modules
import os
import sys
import logging
import webbrowser

# External modules
from PyQt5 import QtWidgets, uic, QtCore

# My modules
import path
import tests


class Settings:
    """Session settings object is used everywhere in the code."""

    def __init__(self):
        """Read settings from file or apply defaults."""
        try:
            # Try to read settings file
            with open(path.p.settings, 'r') as f:
                exec(f.read())
            if len(self.__dict__) < 2:
                raise Exception

        except:
            # Apply default values
            if os.name=='nt':
                # Windows
                ext = '.exe'
                self.path_paraview = 'C:\\Program Files\\ParaView\\bin\\paraview.exe'
                self.path_editor = 'C:\\Windows\\System32\\notepad.exe'
            else:
                # Linux
                ext = ''
                self.path_paraview = '/usr/bin/paraview'
                self.path_editor = '/usr/bin/gedit'

            self.start_model = os.path.join(path.p.examples, 'default.inp')
            self.default_web_browser = 'internal'
            self.logging_level = 'DEBUG'
            self.show_empty_keywords = True
            self.expanded = True
            self.start_cgx_by_default = True
            self.align_windows = True
            self.show_help = False

    def open(self):
        """Open dialog window and pass settings."""
        self.__init__() # re-read settings from file
        sd = SettingsDialog(settings=self)

        # Get response from dialog window
        if sd.exec(): # == 1 if user pressed 'OK'
            sd.save()
            self.__init__() # read settings from file
            logging.warning('For some settings to take effect application\'s restart may be needed.')

    def save(self):
        """Automatic saving of current settings during the workflow."""
        sd = SettingsDialog(settings=self)
        sd.save()


class SettingsDialog(QtWidgets.QDialog):
    """User dialog window with all
    setting attributes: menu File->Settings.
    """

    def __init__(self, settings=None):
        """Create dialog window."""

        # Load UI form - produces huge amount of redundant debug logs
        logging.disable() # switch off logging
        super().__init__() # create dialog window
        uic.loadUi(path.p.settings_xml, self) # load default settings
        logging.disable(logging.NOTSET) # switch on logging

        self.add_widget_for_default_web_browser()

        # Actions
        self.path_paraview_button.clicked.connect(
            lambda: self.select_path(self.path_paraview))
        self.path_editor_button.clicked.connect(
            lambda: self.select_path(self.path_editor))
        self.start_model_button.clicked.connect(
            lambda: self.select_path(self.start_model))

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

    def add_widget_for_default_web_browser(self):
        """QComboBox for default web browser to open help."""
        webbrowser.get() # initialize web browsers
        exclude = ['xdg-open', 'gvfs-open', 'x-www-browser']
        wb_list = [x for x in webbrowser._tryorder if x not in exclude]
        wb_list = sorted(wb_list)
        self.default_web_browser.addItems(wb_list)

    def save(self):
        """Save settings updated via or passed to dialog."""
        with open(path.p.settings, 'w') as f:

            # Iterate over class attributes
            for attr, value in self.__dict__.items():
                class_name = value.__class__.__name__
                if class_name in ['QCheckBox', 'QLineEdit', 'QComboBox']:

                    # Write settings to file
                    if class_name == 'QCheckBox':
                        text = str(value.isChecked())
                    elif class_name == 'QLineEdit':
                        text = value.text()
                        text = '\'' + path.p.abspath(text) + '\'' # covert path to absolute
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

    def select_path(self, path_edit):
        """Open file dialog to select path."""
        file_name = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select path', '', '*')[0]
        if len(file_name):
            path_edit.setText(file_name)


@tests.test_wrapper()
def test():

    # Create application
    app = QtWidgets.QApplication(sys.argv)

    # Create and open settings window
    global s
    s.open()


s = Settings()


if __name__ == '__main__':
    test() # run test
