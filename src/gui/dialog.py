#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Dialog window to create/edit keyword's implementation.
Called via double click on keyword in the treeView.
Here we define a keyword's implementation: its name and inp_code. """

# Standard modules
import os
import sys
import re
import math
import time
import logging
import webbrowser

# External modules
try:
    from PyQt5 import QtWidgets, uic, QtCore, QtGui
except:
    msg = 'Please, install PyQt5 with command:\n'\
        + 'pip3 install PyQt5'
    sys.exit(msg)

# My modules
try:
    from model.kom import ItemType
    import gui
except:
    sys_path = os.path.abspath(__file__)
    sys_path = os.path.dirname(sys_path)
    sys_path = os.path.join(sys_path, '..')
    sys_path = os.path.normpath(sys_path)
    sys_path = os.path.realpath(sys_path)
    sys.path.insert(0, sys_path)
    import clean
    import path
    import gui


class KeywordDialog(QtWidgets.QDialog):

    """
    p - Path
    s - Settings
    f - Window Factory
    """
    def __init__(self, p, s, f, KOM, item):
        self.p = p
        self.s = s
        self.f = f
        self.mw = f.mw
        self.item = item # needed to pass to other functions
        self.widgets = [] # list of created widgets

        # Switch off logging
        hh = logging.getLogger().handlers
        logging.getLogger().handlers = []

        # Create dialog, load form and align window
        super(KeywordDialog, self).__init__()
        uic.loadUi(self.p.dialog_xml, self)

        # Switch on logging
        for h in hh:
            logging.getLogger().addHandler(h)

        if self.s.align_windows:
            size = QtWidgets.QDesktopWidget().availableGeometry()
            # TODO check this bug in Windows 10
            if os.name=='nt': # just bug with window border padding
                width = math.floor(size.width() / 3) - 22
                height = size.height() - 55
            else:
                width = math.floor(size.width() / 3)
                height = size.height()
            self.setGeometry(0, 0, width, height)

        # Add window icon (different for each keyword)
        icon_name = item.name.replace('*', '') + '.png'
        icon_name = icon_name.replace(' ', '_')
        icon_name = icon_name.replace('-', '_')
        icon_path = os.path.join(self.p.img, 'icon_' + icon_name.lower())
        icon = QtGui.QIcon(icon_path)
        self.setWindowIcon(icon)

        # New implementation: draw full form for keyword's arguments
        if self.item.itype == ItemType.KEYWORD:
            self.setWindowTitle('New ' + self.item.name)

            # For each keyword's argument create name and value widgets
            row_number = 0 # row number for vertical layout
            for argument in self.item.items:
                if argument.itype != ItemType.ARGUMENT:
                    continue

                argument_values_items = argument.items
                for ag in argument.name.split('|'):
                    logging.debug('\nArgument ' + ag)

                    # Try to get existing implementations for argument.name
                    keyword = KOM.get_keyword_by_name('*' + ag)
                    if keyword is not None:
                        """
                            For example, add names of *AMPLITUDE implementations,
                                if argument.name is 'AMPLITUDE'
                        """
                        argument_values_items = ['']
                        # Example: ELSET argument in *ELSET keyword
                        if ag != self.item.name.upper()[1:]:
                            implementations = [item.name for item in keyword.get_implementations()]
                            logging.debug('\tKeyword ' + keyword.name)
                            logging.debug('\t\tImplementations ' + str(implementations))
                            logging.debug('\t\tArgument items ' + str(argument.items))
                            if len(implementations) and not len(argument.items):
                                argument.form = 'QComboBox'
                                if len(implementations) == 1:
                                    argument_values_items = implementations
                                if len(implementations) > 1:
                                    argument_values_items.extend(implementations)

                # Argument's values
                if argument.form == 'QComboBox':
                    argument_name_text = argument.name + '='

                    # Predefined values to be chosen
                    argument_values_widget = QtWidgets.QComboBox()
                    argument_values_widget.addItems(argument_values_items)

                    # Assign event to update textEdit widget
                    argument_values_widget.currentIndexChanged.connect(self.onChange)

                    # QComboBox doesn't expand by default
                    sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
                    sizePolicy.setHorizontalStretch(1) # expand horizontally
                    argument_values_widget.setSizePolicy(sizePolicy)

                elif argument.form == 'QLineEdit':
                    argument_name_text = argument.name + '='

                    # Values to be entered
                    argument_values_widget = QtWidgets.QLineEdit()

                    # Assign event to update textEdit widget
                    argument_values_widget.textChanged.connect(self.onChange)

                elif argument.form == 'QCheckBox':
                    argument_name_text = argument.name + ' ' # shift checkbox a little bit for nice view

                    # Flag to be checked
                    argument_values_widget = QtWidgets.QCheckBox()

                    # Assign event to update textEdit widget
                    argument_values_widget.clicked.connect(self.onChange)

                # Mark required argument
                if argument.required:
                    argument_required_widget = QtWidgets.QLabel()
                    argument_required_widget.setText('Required:')
                    argument_required_widget.setStyleSheet('color:Red;')
                    self.vertical_layout.insertWidget(row_number, argument_required_widget)
                    row_number += 1 # first time

                # Mutually exclusive arguments
                if '|' in argument.name:
                    argument_name_widget = QtWidgets.QComboBox()
                    if argument.form == 'QCheckBox':
                        arg_names = argument.name.split('|')
                    else:
                        arg_names = [n + '=' for n in argument.name.split('|')]
                    argument_name_widget.addItems(arg_names)

                    # Assign event to update textEdit widget
                    argument_name_widget.currentIndexChanged.connect(self.onChange)

                else:
                    argument_name_widget = QtWidgets.QLabel()
                    argument_name_widget.setText(argument_name_text)

                # Keep name and values in horizontal layout
                horizontal_layout = QtWidgets.QHBoxLayout()
                horizontal_layout.setContentsMargins(0, 0, 0, 20) # bottom margin
                horizontal_layout.addWidget(argument_name_widget)
                horizontal_layout.addWidget(argument_values_widget)
                horizontal_layout.setAlignment(QtCore.Qt.AlignLeft)

                # Save name and values for processing in self.onChange()
                self.widgets.append(argument_name_widget)
                self.widgets.append(argument_values_widget)

                # Add widgets to dialog window
                self.vertical_layout.insertLayout(row_number, horizontal_layout)
                row_number += 1 # second time

            # Fill textEdit widget with default keyword's configuration
            self.onChange(None)

        # Edit implementation: draw only textEdit
        if self.item.itype == ItemType.IMPLEMENTATION:
            self.setWindowTitle('Edit ' + self.item.name)
            for line in self.item.inp_code:
                self.textEdit.append(line)

        # Actions
        self.buttonBox.accepted.connect(self.onOk)
        self.buttonBox.button(QtWidgets.QDialogButtonBox.Reset).clicked.connect(self.onReset)
        self.buttonBox.helpRequested.connect(self.onHelp)

    # Update piece of INP-code in the textEdit widget
    def onChange(self, event):
        arguments = {} # name:value
        for i, widget in enumerate(self.widgets):

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

            # logging.debug('{} {} {}'.format(i, widget.__class__.__name__, text))

            value = '' # clear value from prev. iteration
            if not 'Required' in text:
                if (i % 2) == 0:
                    # Argument's name
                    name = text # name is always present
                else:
                    # Argument's value
                    value = text

                    # Flag goes without value, only flag name
                    if len(value.strip()):
                        if value == 'QCheckBox':
                            value = ''
                        arguments[name.strip()] = value

        # Generate text for textEdit widget
        if self.item.itype == ItemType.KEYWORD:
            string = self.item.name
            for name, value in arguments.items():
                if self.item.from_new_line:
                    string += '\n' + name + value # argument goes from new line
                else:
                    string += ', ' + name + value # argument goes inline
        if self.item.itype == ItemType.IMPLEMENTATION:
            string = self.item.parent.name

        self.textEdit.setText(string)

    # Reset textEdit widget to initial state
    def onReset(self):
        for i, widget in enumerate(self.widgets):
            if (i % 2) == 1: # iterate over values not labels
                if widget.__class__.__name__ == 'QLineEdit':
                    widget.setText('') # empty is default
                elif widget.__class__.__name__ == 'QComboBox':
                    widget.setCurrentIndex(0) # this row is default
                elif widget.__class__.__name__ == 'QCheckBox':
                    widget.setChecked(False) # uncheck is default
        self.onChange(None)

    # Return piece of created code for the .inp-file
    def onOk(self):
        super(KeywordDialog, self).accept()
        return self.textEdit.toPlainText().strip().split('\n')

    # Open HTML help page in a default web browser
    # TODO Works with Chrome and Firefox - check others
    def onHelp(self):
        url = self.get_url()
        if not os.path.isfile(url):
            msg = 'Help page does not exist:\n{}'.format(
                os.path.relpath(url, start=self.p.app_home_dir))
            logging.error(msg)
            return
        wb = webbrowser.get()
        cmd = wb.name + ' --new-window ' + url
        self.f.sw = gui.window.SlaveWindow(cmd, wb.name)
        self.f.sw.run()
        self.f.create_connection(2)

    # Get URL to the local help page
    def get_url(self):
        if self.item.itype == ItemType.KEYWORD:
            keyword_name = self.item.name[1:] # cut star
        if self.item.itype == ItemType.IMPLEMENTATION:
            keyword_name = self.item.parent.name[1:] # cut star

        # Avoid spaces in html page names
        html_page_name = re.sub(r'[ -]', '_', keyword_name)
        url = os.path.join(self.p.doc, html_page_name + '.html')
        return url


# TODO Get default web browser in Windows:
# https://stackoverflow.com/questions/19037216/how-to-get-a-name-of-default-browser-using-python
# from winreg import *
# with OpenKey(HKEY_CURRENT_USER, r"Software\\Microsoft\\Windows\\Shell\\Associations\\UrlAssociations\\http\\UserChoice") as key:
#     browser = QueryValueEx(key, 'Progid')[0]

# Get default web browser
# TODO Log web browser during startup tests
def test():
    print('Primary web browser:')
    wb = webbrowser.get()
    print(wb.name, wb.basename)

    print('\nAll:')
    for wb in webbrowser._tryorder:
        print(wb)

# Run test
if __name__ == '__main__':
    clean.screen()
    test()
    clean.cache()
