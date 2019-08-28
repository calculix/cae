# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    Application's settings.
    Attributes values are maintained through the ccx_settings.ENV file.
    User dialog form is maintained through the ccx_settings.XML file.
    You may edit ccx_settings.XML directly from Qt Designer.
"""


import os, sys, logging, re
from PyQt5 import QtWidgets, uic


# Session settings object used everywhere in the code
class Settings():


    # Read settings from file
    def __init__(self):
        path = os.path.dirname(sys.argv[0])
        if path.endswith('src'):
            path = path[:-3]
        self.file_name = os.path.join(path, 'src', 'ccx_settings.env') # full path
        f = open(self.file_name).read()
        self.lines = f.split('\n')
        exec(f)


    # Automatically save settings during the workflow
    def save(self):
        with open(self.file_name, 'w') as f:
            new_line = False
            for line in self.lines:
                match = re.search('^self\.(\S+)\s*=\s*(.+)', line)
                if match:
                    param = match.group(1)
                    value = match.group(2)
                    if value.startswith('\'') and value.endswith('\''):
                        value = '\'' + getattr(self, param) + '\''
                    else:
                        value = getattr(self, param)
                    line = 'self.{} = {}'.format(param, value)
                if new_line:
                    f.write('\n')
                f.write(line)
                new_line = True
        logging.info('Settings saved.')


    # Open dialog window and pass settings
    def open(self):
        dialog = Dialog()

        # Get response from dialog window
        if dialog.exec() == Dialog.Accepted: # if user pressed 'OK'
            dialog.onOk()
            logging.warning('For some settings to take effect application\'s restart may be needed.')


# User dialog window with all setting attributes: menu File->Settings
class Dialog(QtWidgets.QDialog):


    # Create dialog window
    def __init__(self):

        # Load UI form
        super(Dialog, self).__init__()
        path = os.path.dirname(sys.argv[0])
        if path.endswith('src'):
            path = path[:-3]
        ccx_settings_xml = os.path.join(path, 'src',  'ccx_settings.xml') # full path
        uic.loadUi(ccx_settings_xml, self)

        # Push settings values to the form
        self.settings = Settings() # read settings from file
        for attr, value in self.settings.__dict__.items():
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
    def onOk(self):
        with open(self.settings.file_name, 'w') as f:

            # Iterate over class attributes
            for attr, value in self.__dict__.items():
                if value.__class__.__name__ in ['QCheckBox', 'QLineEdit', 'QComboBox']:

                    # Update session settings
                    setattr(self.settings, attr, value)

                    # Write settings in file
                    if value.__class__.__name__ == 'QCheckBox':
                        line = 'self.{} = {}'.format(attr, value.isChecked())
                        comment = value.text()
                    if value.__class__.__name__ == 'QLineEdit':
                        line = 'self.{} = \'{}\''.format(attr, value.text())
                        value = self.__dict__['label_' + attr]
                        comment = value.text()
                    if value.__class__.__name__ == 'QComboBox':
                        line = 'self.{} = \'{}\''.format(attr, value.currentText())
                        value = self.__dict__['label_' + attr]
                        comment = value.text()
                    f.write('# ' + comment + '\n')
                    f.write(line + '\n\n')
                    logging.debug(line)


# Test module
if __name__ == '__main__':
    logging.getLogger().setLevel(logging.DEBUG)

    # Create application
    app = QtWidgets.QApplication(sys.argv)

    # Create and open settings window
    Settings().open()

    # Execute application
    a = app.exec_()

    # Clean cached files
    from .clean import cleanCache
    cleanCache()

    # Exit application
    sys.exit(a)
