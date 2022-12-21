#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2023
Distributed under GNU General Public License v3.0

TODO Run internal help in a separated slave window. It will allow to
align dialog and slave help browser via connection.align_windows().

Dialog window to create/edit keyword implementation.
Called via double click on keyword in the treeView.
Keyword implementation is defined by its name and inp_code.
Dialog is created via Factory class, run_master_dialog() method.
So this dialog is a master, help webbrowser (if any) is a slave.
"""

# Standard modules
import os
import sys
import logging

# External modules
from PyQt5 import QtWidgets, uic, QtCore, QtGui, QtWebEngineWidgets
from PyQt5.QtWidgets import QMessageBox

# My modules
sys_path = os.path.abspath(__file__)
sys_path = os.path.dirname(sys_path)
sys_path = os.path.join(sys_path, '..')
sys_path = os.path.normpath(sys_path)
sys_path = os.path.realpath(sys_path)
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
from path import p
from settings import s
from model.kom import ItemType, KWL, KWT
import gui.window


class MyWidget(QtWidgets.QWidget):
    """My custom widget from which to construct the dialog.
    MyWidget is used to visualize Arguments.
    """

    def __init__(self, argument, widgets, pad):
        assert type(argument.name) is str, 'Wrong name type: {}'.format(type(argument.name))
        assert type(widgets) is list, 'Wrong widgets type: {}'.format(type(widgets))
        self.name = argument.name
        self.widgets = widgets
        self.required = argument.get_required()
        self.newline = argument.get_newline()
        super().__init__()

        self.vertical_layout = QtWidgets.QVBoxLayout()
        self.vertical_layout.setContentsMargins(pad, 0, 0, 0)
        if argument.comment:
            comment_label = QtWidgets.QLabel(argument.comment)
            comment_label.setStyleSheet('color:Gray;')
            self.vertical_layout.addWidget(comment_label)

        self.horizontal_layout = QtWidgets.QHBoxLayout()
        self.horizontal_layout.setContentsMargins(0, 0, 0, 10) # bottom margin
        self.horizontal_layout.setAlignment(QtCore.Qt.AlignLeft)

        self.label = None
        if '|' in self.name:
            """Mutually exclusive arguments
            name='FREQUENCY|TIME POINTS'
            """
            self.label = QtWidgets.QComboBox()
            self.label.addItems(self.name.split('|'))
            self.label.text = self.label.currentText
            self.label.my_signal = self.label.currentIndexChanged
        elif self.name:
            self.label = QtWidgets.QLabel(self.name)
            self.label.my_signal = self.label.linkHovered
        if self.label:
            widgets.insert(0, self.label)

        # Mark required argument
        if self.required:
            required_label = QtWidgets.QLabel()
            required_label.setText('*')
            required_label.setStyleSheet('color:Red;')
            widgets.insert(0, required_label)

        if pad:
            self.setEnabled(False)
        for w in widgets:
            self.horizontal_layout.addWidget(w)
        self.vertical_layout.addLayout(self.horizontal_layout)
        self.setLayout(self.vertical_layout)

    def setEnabled(self, status):
        for w in self.widgets:
            w.setEnabled(status)


class Combo(MyWidget):
    """Combo box with label."""

    def __init__(self, argument, pad):
        self.w = QtWidgets.QComboBox()
        self.w.addItems(argument.value.split('|'))

        # QComboBox doesn't expand by default
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(1) # expand horizontally
        self.w.setSizePolicy(sizePolicy)

        super().__init__(argument, [self.w], pad)
        self.my_signal = self.w.currentIndexChanged
        self.text = self.w.currentText

    def reset(self):
        self.w.setCurrentIndex(0)


class Line(MyWidget):
    """Line edit with label."""

    def __init__(self, argument, pad):
        self.w = QtWidgets.QLineEdit()
        self.w.setText(argument.value)
        super().__init__(argument, [self.w], pad)
        self.my_signal = self.w.textChanged
        self.text = self.w.text

    def reset(self):
        self.w.setText('')


class Check(MyWidget):
    """Checkbox with label."""

    def __init__(self, argument, pad):
        self.w = QtWidgets.QCheckBox()
        super().__init__(argument, [self.w], pad)
        self.my_signal = self.w.clicked

    def text(self):
        if self.w.isChecked():
            return self.__class__.__name__
        else:
            return ''

    def reset(self):
        self.w.setChecked(False)


class SelectFileWidget(MyWidget):
    """Custom widget to select files. With label."""

    def __init__(self, argument, pad):
        self.line_edit = QtWidgets.QLineEdit()
        self.push_button = QtWidgets.QPushButton('...', None)
        self.push_button.clicked.connect(self.get_file)
        self.push_button.setFixedSize(30, 30)
        super().__init__(argument, [self.line_edit, self.push_button], pad)
        self.my_signal = self.line_edit.textChanged
        self.text = self.line_edit.text

    def get_file(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Single File', '', '*.inp')[0]
        self.line_edit.setText(fname)

    def reset(self):
        self.line_edit.setText('')


class KeywordDialog(QtWidgets.QDialog):

    @gui.window.init_wrapper()
    def __init__(self, item):
        """Load form and show the dialog."""
        self.info = None # WindowInfo will be set in @init_wrapper
        self.item = item # the one was clicked in the treeView
        # self.widgets = [] # list of created widgets
        self.arguments = []

        # Load UI form - produces huge amount of redundant debug logs
        logging.disable() # switch off logging
        super().__init__() # create dialog window
        uic.loadUi(p.dialog_xml, self) # load empty dialog form
        logging.disable(logging.NOTSET) # switch on logging

        # Set window icon (different for each keyword)
        # TODO Test if it is Windows-specific
        icon_name = self.item.name.replace('*', '') + '.png'
        icon_name = icon_name.replace(' ', '_')
        icon_name = icon_name.replace('-', '_')
        icon_path = os.path.join(p.img, 'icon_' + icon_name.lower())
        icon = QtGui.QIcon(icon_path)
        self.setWindowIcon(icon)
        self.setWindowTitle(self.item.name)

        # Build widgets
        if item.itype == ItemType.IMPLEMENTATION:
            # Fill textEdit with implementation's inp_code
            for line in item.inp_code:
                self.textEdit.append(line)
        elif item.itype == ItemType.KEYWORD:
            # Create widgets for each keyword argument
            self.arguments = KWL.get_keyword_by_name(item.name).get_arguments()
            self.build_argument_widgets(self.arguments)
            self.change(None) # fill textEdit widget with default inp_code

        # Generate html help page from official manual
        self.doc = QtWebEngineWidgets.QWebEngineView()
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(2) # expand horizontally
        self.doc.setSizePolicy(sizePolicy)
        self.url = self.get_help_url()
        self.show_help = s.show_help

        self.show()

    def build_argument_widgets(self, arguments, pad=0):
        for argument in reversed(arguments):
            self.build_argument_widgets(argument.get_arguments(), pad=pad+20)
            if argument.use:
                # Try to get existing implementations
                # NOTE Keywords from the tree/list are different instances
                kwl = KWL.get_keyword_by_name(argument.use)
                kwt = KWT.keyword_dic[kwl.get_path2()]
                implementations = [impl.name for impl in kwt.get_implementations()]
                argument_widget = Combo(argument, pad)
                argument_widget.w.addItems(implementations)
            else:
                argument_widget = eval(argument.form)(argument, pad)
            self.vertical_layout.insertWidget(0, argument_widget)
            # self.widgets.insert(0, argument_widget)
            argument.widget = argument_widget

            # Connect signals to slots
            argument_widget.my_signal.connect(self.change)
            if argument_widget.label:
                argument_widget.label.my_signal.connect(self.change)

    def change(self, data, arguments=None, append=False):
        """Update piece of INP-code in the textEdit widget."""
        args = {} # name:value
        if not append:
            self.textEdit.setText(self.item.name)
        if arguments is None:
            arguments = self.arguments
        for a in arguments:
            w = a.widget

            name = ''
            if w.label:
                name = w.label.text()
            value = w.text() # argument value

            if value:
                # Checkbox goes without value, only flag name
                if value == Check.__name__:
                    value = ''
                elif name:
                    name += '='
                args[name] = value
                if w.newline:
                    txt = '\n' + name + value # argument goes from new line
                else:
                    txt = ', ' + name + value # argument goes inline
                self.textEdit.setText(self.textEdit.toPlainText() + txt)

            # Activate children argument widgets
            for arg in a.get_arguments():
                arg.widget.setEnabled(bool(w.text()))
            self.change(data, a.get_arguments(), append=True)

    def reset(self):
        """Reset textEdit widget to initial state."""
        if arguments is None:
            arguments = self.arguments
        for a in arguments:
            w = a.widget
        # for w in self.widgets:
            if hasattr(w, 'reset'):
                w.reset()

    def accept(self, arguments=None):
        """Check if all required fields are filled."""
        if arguments is None:
            arguments = self.arguments
        for a in arguments:
            w = a.widget
        # for w in self.widgets:
            if w.isEnabled() and w.required:
                # name = w.name.replace('=', '') # argument name
                name = w.name # argument name
                value = w.text() # argument value
                if not value:
                    msg = name + ' is required!'
                    QMessageBox.warning(self, 'Warning', msg)
                    return
        super().accept()

    def ok(self):
        """Return piece of created code for the .inp-file."""
        return self.textEdit.toPlainText().strip().split('\n')

    def get_help_url(self):
        """Get URL to the local doc page."""
        if self.item.itype == ItemType.KEYWORD:
            keyword_name = self.item.name[1:] # cut star
        if self.item.itype == ItemType.IMPLEMENTATION:
            keyword_name = self.item.parent.name[1:] # cut star

        # Avoid spaces and hyphens in html page names
        import re
        html_page_name = re.sub(r'[ -]', '_', keyword_name)
        url = os.path.join(p.doc, html_page_name + '.html')
        return url

    def show_hide_internal_help(self, click):
        """Show / Hide HTML help."""
        size = QtWidgets.QApplication.primaryScreen().availableSize()
        import math
        w = math.floor(size.width() / 3)
        h = self.geometry().height()
        if click:
            self.show_help = not self.show_help
        else:
            self.show_help = s.show_help

        # To show or not to show
        if self.show_help:
            self.doc.load(QtCore.QUrl.fromLocalFile(self.url)) # load help document
            self.setMaximumWidth(size.width())
            self.setMinimumWidth(size.width())
            self.resize(size.width(), h)
            self.horizontal_layout.addWidget(self.doc)
            self.buttonBox.button(QtWidgets.QDialogButtonBox.Help)\
                .setText('Hide help')
        else:
            self.doc.setParent(None) # remove widget
            self.buttonBox.button(QtWidgets.QDialogButtonBox.Help)\
                .setText('Help')
            self.setMaximumWidth(w)
            self.setMinimumWidth(500)
            self.resize(w, h)


def test_dialog():
    """Create keyword dialog"""
    app = QtWidgets.QApplication(sys.argv)
    item = KWL.get_keyword_by_name('*BOUNDARY')
    from gui.window import df
    df.run_master_dialog(item) # 0 = cancel, 1 = ok


if __name__ == '__main__':
    test_dialog()
