# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, UJV Rez, June 2019.
    Distributed under GNU General Public License, version 2.

    Dialog window to edit CalculiX keyword
    Called from CAE via double click on keyword in the treeView
    Actually, here we define a keyword's implementation: its name and INP_code
"""


from PyQt5 import QtWidgets, uic, QtCore


class Dialog(QtWidgets.QDialog):

    def __init__(self, keyword): # here 'keyword' is ccx_dom.keyword object
        self.keyword = keyword # needed to pass to other functions

        # Create dialog window
        super(Dialog, self).__init__()

        # Load basic form
        uic.loadUi('ccx_dialog.ui', self)

        # For each keyword's argument create name and value widgets
        index = 0 # row number for vertical layout
        for argument in keyword.items:
            if argument.item_type != 'argument':
                continue

            # Argument's values
            if len(argument.values):
                # Predefined values to be chosen
                argument_values_widget = QtWidgets.QComboBox()
                argument_values_widget.addItems(argument.values)

                # Assign event to update textEdit widget
                argument_values_widget.currentIndexChanged.connect(self.onChange)

                # QComboBox doesn't expand by default
                sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
                sizePolicy.setHorizontalStretch(1) # expand horizontally
                argument_values_widget.setSizePolicy(sizePolicy)
            elif '=' in argument.name:
                # Values to be entered
                argument_values_widget = QtWidgets.QLineEdit()

                # Assign event to update textEdit widget
                argument_values_widget.textChanged.connect(self.onChange)
            else:
                # Flag to be checked
                argument_values_widget = QtWidgets.QCheckBox()
                argument.name += ' ' # shift checkbox a little bit for nice view

                # Assign event to update textEdit widget
                argument_values_widget.clicked.connect(self.onChange)

            # Mark required argument
            if argument.required:
                argument_required_widget = QtWidgets.QLabel()
                argument_required_widget.setText('Required:')
                argument_required_widget.setStyleSheet('color:Red;')
                self.vertical_layout.insertWidget(index, argument_required_widget)
                index += 1 # first time

            # Argument's name
            if '|' in argument.name:
                # Mutually exclusive arguments
                argument_name_widget = QtWidgets.QComboBox()
                argument_name_widget.addItems(argument.name.split('|'))

                # Assign event to update textEdit widget
                argument_name_widget.currentIndexChanged.connect(self.onChange)
            else:
                argument_name_widget = QtWidgets.QLabel()
                argument_name_widget.setText(argument.name)

            # Keep name and values in horizontal layout
            horizontal_layout = QtWidgets.QHBoxLayout()
            horizontal_layout.setContentsMargins(0, 0, 0, 20) # bottom margin
            horizontal_layout.addWidget(argument_name_widget)
            horizontal_layout.addWidget(argument_values_widget)
            horizontal_layout.setAlignment(QtCore.Qt.AlignLeft)

            self.vertical_layout.insertLayout(index, horizontal_layout) # add widgets to dialog window
            index += 1 # second time

        # Fill textEdit widget with default keyword's configuration
        self.onChange(None)

        # Actions
        self.buttonBox.accepted.connect(self.onOk)
        self.buttonBox.button(QtWidgets.QDialogButtonBox.Reset).clicked.connect(self.onReset)


    # Update piece of INP-code in the textEdit widget
    def onChange(self, event):
        arguments = {} # name:value
        string = self.keyword.name
        j = 0
        for i, widget in enumerate(self.children()):
            if i>3: # skip widgets created manually in ccx_dialog.ui

                # Get text from widget: argument's name and value
                text = '' # clear text from prev. iteration
                if widget.__class__.__name__ == 'QLabel':
                    text = widget.text()
                elif widget.__class__.__name__ == 'QLineEdit':
                    text = widget.text()
                elif widget.__class__.__name__ == 'QComboBox':
                    text = widget.currentText()
                elif widget.__class__.__name__ == 'QCheckBox':
                    if widget.isChecked():
                        text = 'QCheckBox'

                # print(i, j, widget.__class__.__name__, text)

                value = '' # clear value from prev. iteration
                if not 'Required' in text:
                    if (j % 2) == 0:
                        # Argument's name
                        name = text # name is always present
                    else:
                        # Argument's value
                        value = text

                        # flag goes without value, only flag name
                        if len(value.strip()):
                            if value == 'QCheckBox':
                                value = ''
                            arguments[name.strip()] = value
                    j += 1

        # Generate text for textEdit widget
        for name, value in arguments.items():
            if self.keyword.from_new_line:
                string += '\n' + name + value # argument goes from new line
            else:
                string += ', ' + name + value # argument goes inline

        self.textEdit.setText(string)
        # print(string)


    # Reset textEdit widget to initial state
    def onReset(self):
        self.textEdit.setText(self.keyword.name)


    # Return piece of created code for the .inp-file
    def onOk(self):
        super(Dialog, self).accept()
        return self.textEdit.toPlainText().strip().split('\n')
