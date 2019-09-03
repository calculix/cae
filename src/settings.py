# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, September 2019
    Distributed under GNU General Public License v3.0

    Application's settings.
    Attributes values are maintained through the settings ENV file.
    User dialog form is maintained through the settings.XML file.
    You may edit settings.XML from Qt Designer.
"""


from path import Path
import os, sys, logging, re
from PyQt5 import QtWidgets, uic


# Session settings object used everywhere in the code
class Settings():


    # Read settings from file
    def __init__(self):
        self.p = Path() # calculate absolute pathes

        # Try to read settings file
        try:
            with open(self.p.settings, 'r') as f:
                exec(f.read())
            if len(self.__dict__) < 2:
                raise Exception

        # Apply default values from UI form
        except:
            print('Error reading ENV settings. Overwriting file with default values.')
            self = SettingsDialog()
            self.save()


    # Open dialog window and pass settings
    def open(self):
        dialog = SettingsDialog(push=True)

        # Warning about Cygwin DLLs
        if os.name=='nt':
            logging.warning('In Windows ccx binary may not work if placed outside \'bin\' directory. It needs Cygwin DLLs!')

        # Get response from dialog window
        if dialog.exec(): # == 1 if user pressed 'OK'
            dialog.save()
            logging.warning('For some settings to take effect application\'s restart may be needed.')
        

# User dialog window with all setting attributes: menu File->Settings
class SettingsDialog(QtWidgets.QDialog):


    # Create dialog window
    def __init__(self, push=False):

        # Load UI form
        QtWidgets.QDialog.__init__(self)
        self.p = Path() # calculate absolute pathes
        uic.loadUi(self.p.settings_xml, self) # load default settings from

        # Push settings values to the form
        if push:
            settings = Settings()
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


    # Save settings updated via dialog
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
                        if '\\' in text: # reconstruct path for Windows
                            text = '\\\\'.join(text.split('\\'))
                        text = '\'' + self.p.abspath(text) + '\'' # covert path to absolute
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
    from clean import cleanCache
    p = Path() # calculate absolute pathes
    cleanCache(p.src)
