# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, December 2019
    Distributed under GNU General Public License v3.0

    Application's settings.
    Attributes values are maintained through the settings ENV file.
    User dialog form is maintained through the settings XML file.
    You may edit settings.XML from Qt Designer.
"""


from Path import Path
from PyQt5 import QtWidgets, uic


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

