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
from model.kom import ItemType, KWL, KWT, Implementation
import gui.window


class SettingsDialog(QtWidgets.QDialog):
    """User dialog window with all
    setting attributes: menu File->Settings.
    """

    def __init__(self):
        """Create dialog window."""

        # Load UI form - produces huge amount of redundant debug logs
        logging.disable() # switch off logging
        super().__init__() # create dialog window
        uic.loadUi(p.settings_xml, self) # load default settings
        logging.disable(logging.NOTSET) # switch on logging

        self.add_widget_for_default_web_browser()

        # Actions
        self.path_paraview_button.clicked.connect(
            lambda: self.select_path(self.path_paraview))
        self.path_editor_button.clicked.connect(
            lambda: self.select_path(self.path_editor))
        self.start_model_button.clicked.connect(
            lambda: self.select_path(self.start_model))

        # Push values to the form
        for attr_name, attr_value in s.__dict__.items():
            widget = self.findChild(QtWidgets.QCheckBox, attr_name)
            if widget is not None:
                widget.setChecked(attr_value)
                continue

            widget = self.findChild(QtWidgets.QLineEdit, attr_name)
            if widget is not None:
                widget.setText(attr_value)
                continue

            widget = self.findChild(QtWidgets.QComboBox, attr_name)
            if widget is not None:
                widget.setCurrentText(attr_value)
                continue

    def add_widget_for_default_web_browser(self):
        """QComboBox for default web browser to open help.
        NOTE webbrowser does not automatically recognize
        all installed browsers in Windows.
        """
        import webbrowser
        webbrowser.get() # initialize web browsers
        exclude = ['xdg-open', 'gvfs-open', 'x-www-browser']
        wb_list = [x for x in webbrowser._tryorder if x not in exclude]
        wb_list = sorted(wb_list)
        self.default_web_browser.addItems(wb_list)

    def save(self):
        """Save settings updated via or passed to dialog.
        Update global settings variable 's'.
        """
        with open(p.settings, 'w') as f:

            # Iterate over class attributes
            for attr_name, attr_value in self.__dict__.items():
                class_name = attr_value.__class__.__name__
                if class_name in ['QCheckBox', 'QLineEdit', 'QComboBox']:

                    # Write settings to file
                    if class_name == 'QCheckBox':
                        setting_value = attr_value.isChecked()
                        text = str(setting_value)
                    elif class_name == 'QLineEdit':
                        text = attr_value.text()
                        text = '\'' + p.abspath(text) + '\'' # covert path to absolute
                        if '\\' in text: # reconstruct path for Windows
                            text = '\\\\'.join(text.split('\\'))
                        setting_value = text[1:-1] # cut quotes
                        attr_value = self.__dict__['label_' + attr_name]
                    elif class_name == 'QComboBox':
                        text = '\'' + attr_value.currentText() + '\''
                        setting_value = text
                        attr_value = self.__dict__['label_' + attr_name]

                    s.__setattr__(attr_name, setting_value)
                    line = 'self.{} = {}'.format(attr_name, text)
                    comment = attr_value.text()
                    f.write('# ' + comment + '\n')
                    f.write(line + '\n\n')

    def select_path(self, path_edit):
        """Open file dialog to select path."""
        file_name = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select path', '', '*')[0]
        if len(file_name):
            path_edit.setText(file_name)

    def open(self):
        """Open dialog window and pass settings."""
        if self.exec(): # == 1 if user pressed 'OK'
            self.save()
            # s.__init__() # read settings from file
            logging.warning('For some settings to take effect application\'s restart may be needed.')


def build_widgets(dialog, arguments, parent_layout):
    """Build widgets for direct children of the Keyword.
    Recursion for nested arguments/groups is implemented
    inside GroupWidget/ArgumentWidget - not here.
    """
    for a in arguments:
        txt = '{} has neither "form" nor "user" attributes'.format(a)
        assert hasattr(a, 'form') or hasattr(a, 'use'), txt

        if hasattr(a, 'form') and a.form:
            form = a.form
        if hasattr(a, 'use') and a.use:
            form = Combo.__name__
        a.widget = eval(form)(dialog, a) # argument's widget
        parent_layout.addWidget(a.widget)


def change(dialog, arguments=[], append=False):
    """Update piece of INP-code in the textEdit when
    a signal is emitted in any of argument's widgets.
    """
    if dialog.item.itype == ItemType.IMPLEMENTATION:
        dialog.textEdit.clear()
        return
    if not append:
        dialog.textEdit.setText(dialog.item.name)
    if not arguments and not append:
        arguments = dialog.item.get_arguments()
    for a in arguments:
        if a.widget is None:
            continue
        w = a.widget
        old_value = dialog.textEdit.toPlainText()
        new_value = w.text() if w.isEnabled() else '' # argument value
        if old_value.endswith('\n') and new_value.startswith(', '):
            new_value = new_value[2:]
        if w.__class__.__name__ != Empty.__name__:
            if old_value.endswith(', ') and new_value.startswith(', '):
                new_value = new_value[2:]
        dialog.textEdit.setText(old_value + new_value)

        # Recursively walk through the whole keyword arguments
        args = a.get_arguments()
        if args:
            change(dialog, args, append=True)


def reset(arguments):
    """Reset argument's widgets to initial state."""
    for a in arguments:
        w = a.widget
        if hasattr(w, 'reset'):
            w.reset()


class KeywordDialog(QtWidgets.QDialog):

    @gui.window.init_wrapper()
    def __init__(self, item):
        """Load form and show the dialog.
        'item' is one of: Keyword or Implementation
        """
        # Load UI form - produces huge amount of redundant debug logs
        logging.disable() # switch off logging
        super().__init__() # create dialog window
        uic.loadUi(p.dialog_xml, self) # load empty dialog form
        logging.disable(logging.NOTSET) # switch on logging

        self.item = item # the one was clicked in the treeView
        self.info = None # WindowInfo will be set in @init_wrapper
        self.arguments = []
        self.make_screenshot = False

        # Set window icon (different for each keyword)
        # TODO Test if it is Windows-specific
        icon_name = item.name.replace('*', '') + '.png'
        icon_name = icon_name.replace(' ', '_')
        icon_name = icon_name.replace('-', '_')
        icon_path = os.path.join(p.img, 'icon_' + icon_name.lower())
        icon = QtGui.QIcon(icon_path)
        self.setWindowIcon(icon)
        self.setWindowTitle(item.name)

        # Fill textEdit with implementation's inp_code
        if item.itype == ItemType.IMPLEMENTATION:
            for line in item.inp_code:
                self.textEdit.append(line)

        # Create widgets for each keyword argument
        elif item.itype == ItemType.KEYWORD:
            kw = KWL.get_keyword_by_name(item.name)
            self.arguments = kw.get_arguments()
            build_widgets(self, kw.get_arguments(), self.widgets_layout)
            change(self) # fill textEdit widget with default inp_code

        # Generate html help page from official manual
        self.doc = QtWebEngineWidgets.QWebEngineView()
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(2) # expand horizontally
        self.doc.setSizePolicy(sizePolicy)
        self.url = self.get_help_url()
        self.show_help = s.show_help

        self.show()

    def reset(self):
        """Reset all widgets to initial state."""
        reset(self.arguments)
        change(self) # update QTextEdit with INP-code

    def accept(self, arguments=None, depth=0):
        """Check if all required fields are filled.
        Also save dialog's screenshot if module is ran from tests.
        """
        # TODO Temporarily check is disabled. Checks are needed.
        super().accept()
        return

        ok = True
        if arguments is None:
            arguments = self.item.get_arguments()
        for a in arguments:
            w = a.widget
            if w is not None and w.isEnabled() and w.required:
                name = a.name # argument name
                value = w.text() # argument value
                if not value:
                    msg = 'Fill all required fields'
                    if name:
                        msg += ': ' + name
                    else:
                        msg += '!'
                    QMessageBox.warning(self, 'Warning', msg)
                    return False
            ok = ok and self.accept(a.get_arguments(), depth+1)
        if depth:
            return ok
        if depth==0 and ok:
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

    def screenshot(self):
        """Make screenshot of the dialog with all widgets and window title.
        Keyword name without an asterisk (*) is a PNG file name.
        """
        screen = QtWidgets.QApplication.primaryScreen()
        screenshot = screen.grabWindow(0, self.pos().x(), self.pos().y(), 
            self.size().width(), self.size().height())
        fname = os.path.join(p.img, 'KeywordDialog', self.item.name[1:] + '.png')
        screenshot.save(fname, 'png')


class ArgumentWidget(QtWidgets.QWidget):
    """ArgumentWidget is used to visualize Arguments."""

    def __init__(self, dialog, argument, widgets, show_label=True):
        assert type(argument.name) is str, 'Wrong name type: {}'.format(type(argument.name))
        assert type(widgets) is list, 'Wrong widgets type: {}'.format(type(widgets))
        self.name = argument.name
        self.widgets = widgets
        self.required = argument.get_required()
        self.newlines = argument.get_newlines()
        self.readonly = argument.get_readonly()
        super().__init__()

        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.setContentsMargins(0, 0, 0, 0)

        self.required_label = QtWidgets.QLabel()
        self.required_label.setText('*')
        self.required_label.setStyleSheet('color:Red;')

        if argument.comment:
            comment_layout = QtWidgets.QHBoxLayout()
            comment_layout.setContentsMargins(0, 0, 0, 0) # bottom margin
            comment_layout.setAlignment(QtCore.Qt.AlignLeft)
            comment_label = QtWidgets.QLabel(argument.comment)
            comment_label.setStyleSheet('color: Blue;')

            if self.required:
                comment_layout.addWidget(self.required_label)

            comment_layout.addWidget(comment_label)
            self.v_layout.addLayout(comment_layout)

        self.horizontal_layout = QtWidgets.QHBoxLayout()
        self.horizontal_layout.setContentsMargins(0, 0, 0, 10) # bottom margin
        self.horizontal_layout.setAlignment(QtCore.Qt.AlignLeft)
        self.v_layout.addLayout(self.horizontal_layout)

        self.label = None
        if '|' in self.name:
            """Mutually exclusive arguments
            name='FREQUENCY|TIME POINTS'
            """
            self.label = QtWidgets.QComboBox()
            self.label.addItems(self.name.split('|'))
            self.label.text = self.label.currentText
            self.label.currentIndexChanged.connect(lambda: change(dialog))
        elif self.name:
            self.label = QtWidgets.QLabel(self.name)
            self.label.linkHovered.connect(lambda: change(dialog))
        if show_label and self.label is not None:
            widgets.insert(0, self.label)

        # Mark required argument
        if self.required and not argument.comment:
            widgets.insert(0, self.required_label)

        for w in widgets:
            if hasattr(w, 'setReadOnly'):
                w.setReadOnly(self.readonly)
            self.horizontal_layout.addWidget(w)
        self.setLayout(self.v_layout)

        # Recursion for nested arguments/groups
        build_widgets(dialog, argument.get_arguments(), self.v_layout)

    def setEnabled(self, status):
        for w in self.widgets:
            w.setEnabled(status)
        if status:
            self.required_label.setStyleSheet('color:Red;')
        else:
            self.required_label.setStyleSheet('color:Gray;')

    def text(self):
        if not self.w.isEnabled():
            return ''
        if not self.isEnabled():
            return ''
        sep = ''
        if not self.newlines:
            sep = ', '
        ct = ''
        if hasattr(self.w, 'currentText'):
            ct = self.w.currentText()
        if hasattr(self.w, 'toPlainText'):
            ct = self.w.toPlainText()
        elif hasattr(self.w, 'text'):
            ct = self.w.text()
        if ct and not ct.startswith('Create'):
            if self.label is not None and self.label.text():
                return self.newlines + sep + self.label.text() + '=' + ct
            elif self.name:
                return self.newlines + sep + self.name + '=' + ct
            else:
                return self.newlines + sep + ct
        else:
            return self.newlines


class GroupWidget(QtWidgets.QWidget):
    """Custom widget container - a Group in kw_list.xml.
    Unite arguments and apply on them horizontal (HBox)
    o vertical (VBox) layout.
    """
    def __init__(self, group, layout):
        self.name = group.name
        self.required = group.get_required()
        self.newlines = group.get_newlines()
        for a in group.get_arguments():
            a.required = a.required or self.required
        super().__init__()

        # Add comment label
        v_layout = QtWidgets.QVBoxLayout()
        v_layout.setContentsMargins(0, 0, 0, 0)
        v_layout.insertLayout(0, layout)
        if group.comment:
            comment_label = QtWidgets.QLabel(group.comment)
            comment_label.setStyleSheet('color: Blue;')
            v_layout.insertWidget(0, comment_label)
        self.setLayout(v_layout)

    def text(self, *args):
        return self.newlines

    def reset(self):
        reset(self.arguments)


# TODO It is widget for argument - should inherit ArgumentWidget
class Table(QtWidgets.QWidget):
    """Custom QTableWidget."""

    def __init__(self, dialog, argument):
        self.dialog = dialog
        self.required = argument.get_required()
        self.arguments = argument.get_arguments()
        self.newlines = argument.get_newlines()
        super().__init__()
        layout = QtWidgets.QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)

        # Add comment label
        if argument.comment:
            comment_label = QtWidgets.QLabel(argument.comment)
            comment_label.setStyleSheet('color: Blue;')
            layout.insertWidget(0, comment_label)

        self.w = QtWidgets.QTableWidget()
        self.w.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.w.setAlternatingRowColors(True)
        self.w.setColumnCount(len(self.arguments))
        if self.required:
            headers = ['*' + a.comment for a in self.arguments]
            self.w.horizontalHeader().setStyleSheet('color:Red;')
        else:
            headers = [a.comment for a in self.arguments]
        self.w.setHorizontalHeaderLabels(headers)
        self.w.setRowCount(0)
        self.w.verticalHeader().hide()
        layout.addWidget(self.w)

        self.row_add()

        l = QtWidgets.QHBoxLayout()
        s = QtWidgets.QSpacerItem(0, 0, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        b1 = QtWidgets.QPushButton()
        b1.setText('+')
        b1.setFixedSize(30, 30)
        b2 = QtWidgets.QPushButton()
        b2.setText('-')
        b2.setFixedSize(30, 30)

        b1.clicked.connect(self.row_add)
        b2.clicked.connect(self.row_rem)
        l.addSpacerItem(s)
        l.addWidget(b1)
        l.addWidget(b2)

        layout.addLayout(l)
        self.setLayout(layout)

    def row_add(self):
        i = self.w.rowCount()
        self.w.insertRow(i)
        self.w.setRowHeight(i, 20)
        for j in range(self.w.columnCount()):
            self.w.itemChanged.connect(lambda: change(self.dialog))

    def row_rem(self):
        i = self.w.rowCount()
        if i > 1:
            self.w.removeRow(self.w.rowCount()-1)

    def text(self):
        txt = self.newlines
        for i in range(self.w.rowCount()):
            for j in range(self.w.columnCount()):
                try:
                    txt += self.w.item(i, j).text() + ', '
                except:
                    pass
            txt += '\n'
        return txt.rstrip()[:-1]

    def reset(self):
        self.w.clearContents()
        self.w.setRowCount(1)


# TODO It is widget for argument - should inherit ArgumentWidget
class Group(QtWidgets.QWidget):
    """QGroupBox with Argument widgets.
    Used for argument with arguments inside.
    *CLOAD
    """
    def __init__(self, dialog, argument, box_layout):
        self.required = argument.get_required()
        self.argument = argument
        self.arguments = argument.get_arguments()
        for a in self.arguments:
            a.required = self.required
        super().__init__()

        box_layout.setContentsMargins(20, 10, 8, 0)
        build_widgets(dialog, self.arguments, box_layout)

        self.gbox = QtWidgets.QGroupBox()
        self.gbox.setCheckable(True)
        self.gbox.setTitle(argument.name)
        self.gbox.clicked.connect(lambda: change(dialog))
        self.gbox.setLayout(box_layout)
        self.reset()

        v_layout = QtWidgets.QVBoxLayout()
        v_layout.setContentsMargins(0, 0, 0, 0)
        v_layout.addWidget(self.gbox)
        self.setLayout(v_layout)

    def text(self):
        txt = ''
        if self.gbox.isEnabled() and self.gbox.isChecked():
            txt = ', ' + self.argument.name
        return txt + self.argument.get_newlines() # TODO newlines in the begining?

    def reset(self):
        self.gbox.setChecked(self.required)
        reset(self.arguments)


class HGroup(Group):
    """QGroupBox with horizontal layout."""

    def __init__(self, dialog, argument):
        super().__init__(dialog, argument, QtWidgets.QHBoxLayout())


class VGroup(Group):
    """QGroupBox with vertical layout."""

    def __init__(self, dialog, argument):
        super().__init__(dialog, argument, QtWidgets.QVBoxLayout())


class Tabs(GroupWidget):
    """A Group drawn with QTabWidget.
    Widgets from all tabs will participate in the final INP code.
    """
    def __init__(self, dialog, group):
        self.arguments = group.get_arguments()
        layout = QtWidgets.QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        super().__init__(group, layout)

        self.w = QtWidgets.QTabWidget()
        for gr in self.arguments:
            tab = QtWidgets.QWidget()
            l = QtWidgets.QVBoxLayout()
            tab.setLayout(l)
            self.w.addTab(tab, gr.name)
            build_widgets(dialog, [gr], l) # build a.widget
            s = QtWidgets.QSpacerItem(0, 0, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
            l.addSpacerItem(s)

        layout.addWidget(self.w)


class Tabs2(Tabs):
    """Only widgets from active tab will participate in the final INP code."""

    def __init__(self, dialog, group):
        self.dialog = dialog
        super().__init__(dialog, group)
        for i,a in enumerate(group.get_arguments()):
            a.widget.setDisabled(i!=0)
        self.w.currentChanged.connect(self.tab_activated)

    def tab_activated(self, index):
        """Needed to deactivate tab widgets not to participate in dialog.accept()"""
        for i,a in enumerate(self.arguments):
            a.widget.setDisabled(i!=index)
        change(self.dialog)


class Box(GroupWidget):
    """GroupWidget."""

    def __init__(self, dialog, group, layout):
        self.arguments = group.get_arguments()
        layout.setContentsMargins(0, 0, 0, 0)
        super().__init__(group, layout)
        build_widgets(dialog, self.arguments, layout)


class HBox(Box):
    """GroupWidget with horizontal layout."""

    def __init__(self, dialog, group):
        super().__init__(dialog, group, QtWidgets.QHBoxLayout())


class VBox(Box):
    """GroupWidget with vertical layout."""

    def __init__(self, dialog, group):
        super().__init__(dialog, group, QtWidgets.QVBoxLayout())


class Or(GroupWidget):
    """GroupWidget with radiobuttons."""

    def __init__(self, dialog, group, layout):
        self.arguments = group.get_arguments()
        layout.setContentsMargins(0, 0, 0, 0)
        self.aw = {}
        for i, a in enumerate(self.arguments):
            hl = QtWidgets.QHBoxLayout()
            rb = QtWidgets.QRadioButton()
            rb.setChecked(not bool(i))
            hl.addWidget(rb)
            a.required = group.required
            build_widgets(dialog, [a], hl) # build a.widget
            self.aw[i] = (a, rb)
            rb.toggled.connect(self.toggle)
            rb.toggled.connect(lambda: change(dialog))
            layout.addLayout(hl)
        super().__init__(group, layout)
        self.toggle()

    def toggle(self):
        for (a, rb) in self.aw.values():
            a.widget.setEnabled(rb.isChecked())


class Grid(GroupWidget):
    """Grid of checkboxes."""

    def __init__(self, dialog, argument):
        self.newlines = argument.get_newlines()
        self.checkboxes = []
        l = QtWidgets.QGridLayout()
        row = 0
        col = 0
        for a in argument.value.split(','):
            cb = QtWidgets.QCheckBox(a.strip())
            self.checkboxes.append(cb)
            cb.toggled.connect(lambda: change(dialog))
            l.addWidget(cb, row, col)
            col += 1
            if col % 4 == 0:
                col = 0
                row += 1
        super().__init__(argument, l)

    def reset(self):
        for cb in self.checkboxes:
            cb.setChecked(False)

    def text(self):
        txt = ''
        for cb in self.checkboxes:
            if cb.isChecked():
                txt += ', ' + cb.text()
        if self.newlines:
            return self.newlines + txt[2:]
        else:
            return txt


class OrH(Or):
    """GroupWidget with radiobuttons. Horizontal layout."""

    def __init__(self, dialog, group):
        l = QtWidgets.QHBoxLayout()
        super().__init__(dialog, group, l)


class OrV(Or):
    """GroupWidget with radiobuttons. Vertical layout."""

    def __init__(self, dialog, group):
        l = QtWidgets.QVBoxLayout()
        super().__init__(dialog, group, l)


class Repl(ArgumentWidget):
    """ArgumentWidget with horizontal layout. Accepts string with
    asterisks as input. Allows user to replace asterisks with numbers.
    Contains labels and inputs.
    """
    def __init__(self, dialog, argument):
        self.argument = argument
        self.widgets = []
        self.text_widgets = []
        self.edits = []
        parts = argument.value.split('*')
        for i,s in enumerate(parts):
            l1 = QtWidgets.QLabel(s)
            self.widgets.append(l1)
            self.text_widgets.append(l1)
            if i < len(parts) - 1:
                l2 = QtWidgets.QLineEdit()
                l2.setFixedWidth(110)
                l2.setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
                l2.setValidator(QtGui.QIntValidator())
                l2.textChanged.connect(lambda: change(dialog))
                self.widgets.append(l2)
                self.text_widgets.append(l2)
                self.edits.append(l2)
        super().__init__(dialog, argument, self.widgets)

    def text(self):
        entered_text = ''
        for w in self.edits:
            entered_text += w.text()
        txt = self.argument.get_newlines()
        if self.argument.name:
            txt += self.argument.name + '='
        for w in self.text_widgets:
            txt += w.text()
        if entered_text:
            return ', ' + txt
        else:
            return ''

    def reset(self):
        for w in self.edits:
            w.clear()


class Combo(ArgumentWidget):
    """QComboBox widget with label."""

    def __init__(self, dialog, argument):
        self.dialog = dialog
        self.argument = argument
        self.w = QtWidgets.QComboBox()

        # Try to get existing implementations
        if hasattr(argument, 'use') and argument.use:
            implementations = [impl.name for impl in KWT.get_implementations(argument.use)]
            self.w.addItem('')
            if implementations:
                self.w.addItems(implementations)
            self.w.addItem('Create ' + argument.use)

        if argument.value:
            for v in argument.value.split('|'):
                self.w.addItem(v.strip())

        # QComboBox doesn't expand by default
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(1) # expand horizontally
        self.w.setSizePolicy(sizePolicy)

        self.w.currentIndexChanged.connect(self.index_changed)
        super().__init__(dialog, argument, [self.w])

    def reset(self):
        self.w.setCurrentIndex(0)

    def index_changed(self, *args):
        if hasattr(self.argument, 'use') and self.argument.use and self.w.currentText().startswith('Create'):
            i = args[0]
            kwt_item = KWT.get_top_keyword_by_name(KWT.root, self.argument.use)
            kwl_item = KWL.get_keyword_by_name(self.argument.use)

            # Exec dialog and recieve answer
            # Process response from dialog window if user pressed 'OK'
            d = KeywordDialog(kwl_item)
            if d.exec(): # 0 = cancel, 1 = ok
                # The generated piece of .inp code for the CalculiX input file
                inp_code = d.ok() # list of strings

                # Create implementation object for keyword
                impl = Implementation(kwt_item, inp_code) # create keyword implementation

                # Add impl to the tree
                from gui.tree import t
                t.put_implementation(impl)

                # Add and item to the drop-down list
                self.w.removeItem(i)
                self.w.addItem(impl.name)
                self.w.setCurrentIndex(i)

        change(self.dialog)


class Use2(ArgumentWidget):
    """A widget for multiple argument values, allow to choose Implementation type:
    <Argument form='Use2'>*NSET|*ELSET</Argument>
    ... and shows a combo with all implementations of the chosen type.

    Two combo boxes: left one allows to select implementation in the right one.
    """
    def __init__(self, dialog, argument):
        self.argument = argument
        self.use = QtWidgets.QComboBox()
        self.values = argument.value.split('|')
        for v in self.values:
            self.use.addItem(v)

        # Try to get existing implementations
        self.w = QtWidgets.QComboBox()
        self.draw_implementations()

        # QComboBox doesn't expand by default
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(1) # expand horizontally
        self.w.setSizePolicy(sizePolicy)

        self.w.currentIndexChanged.connect(lambda: change(dialog))
        self.use.currentIndexChanged.connect(self.draw_implementations)
        super().__init__(dialog, argument, [self.use, self.w])

    def draw_implementations(self, *args):
        self.w.clear()
        kw = self.values[self.use.currentIndex()]
        implementations = [impl.name for impl in KWT.get_implementations(kw)]
        self.w.addItem('')
        if implementations:
            self.w.addItems(implementations)
        self.w.addItem('Create ' + kw)

    def reset(self):
        self.use.setCurrentIndex(0)
        self.w.setCurrentIndex(0)


class Text(ArgumentWidget):
    """QTextEdit widget with label."""

    def __init__(self, dialog, argument):
        self.argument = argument
        self.w = QtWidgets.QTextEdit()
        self.w.setMinimumSize(200, 52)
        super().__init__(dialog, argument, [self.w])
        self.setText(argument.value)
        self.setEnabled = self.w.setEnabled
        self.w.textChanged.connect(lambda: change(dialog))

    def reset(self):
        self.setText(self.argument.value)

    def setText(self, text):
        self.w.clear()
        for line in text.split('\\n'):
            self.w.append(line)
        font = self.w.document().defaultFont()
        fontMetrics = QtGui.QFontMetrics(font)
        textSize = fontMetrics.size(0, self.w.toPlainText())
        textHeight = textSize.height() + 10 # Need to tweak
        self.w.setMaximumHeight(textHeight)


class Empty(ArgumentWidget):
    """No widget. Adds space into the INP code."""

    def __init__(self, dialog, argument):
        self.newlines = argument.get_newlines()
        super().__init__(dialog, argument, [])

    def text(self):
        return self.newlines + ', '


class Bool(ArgumentWidget):
    """Checkbox widget with label."""

    def __init__(self, dialog, argument):
        self.argument = argument
        self.status = '|' in argument.name
        if self.status:
            self.w = QtWidgets.QCheckBox()
        else:
            self.w = QtWidgets.QCheckBox(argument.name)
        self.w.clicked.connect(lambda: change(dialog))
        super().__init__(dialog, argument, [self.w], show_label=self.status)
        self.reset()

    def text(self):
        # TODO newlines?
        if self.w.isEnabled() and self.w.isChecked():
            if self.status:
                return ', ' + self.label.text()
            else:
                return ', ' + self.argument.name
        else:
            return ''

    def reset(self):
        self.w.setChecked(self.argument.get_required())


class SelectFile(ArgumentWidget):
    """A custom widget to select files. With label."""

    def __init__(self, dialog, argument):
        self.argument = argument
        self.w = QtWidgets.QLineEdit()
        self.push_button = QtWidgets.QPushButton('...', None)
        self.push_button.clicked.connect(self.get_file)
        self.push_button.setFixedSize(30, 30)
        self.w.textChanged.connect(lambda: change(dialog))
        super().__init__(dialog, argument, [self.w, self.push_button])

    def get_file(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Single File', '', '*.inp;;*.frd;;*.mtx')[0]
        self.w.setText(fname)

    def reset(self):
        self.w.setText(self.argument.value)


class Line(ArgumentWidget):
    """QLineEdit widget with label."""

    def __init__(self, dialog, argument):
        self.argument = argument
        self.w = QtWidgets.QLineEdit()
        self.w.setText(argument.value)
        self.w.textChanged.connect(lambda: change(dialog))
        self.setEnabled = self.w.setEnabled
        super().__init__(dialog, argument, [self.w])

    def reset(self):
        self.w.setText(self.argument.value)


class Int(Line):
    """Text widget accepting int number."""

    def __init__(self, dialog, argument):
        super().__init__(dialog, argument)
        self.w.setValidator(QtGui.QIntValidator())


class Float(Line):
    """Text widget accepting float number."""

    def __init__(self, dialog, argument):
        super().__init__(dialog, argument)
        self.w.setValidator(QtGui.QDoubleValidator())


def test_dialog():
    """Prepare logging."""
    logging.getLogger().setLevel(logging.NOTSET)
    fmt = logging.Formatter('%(levelname)s: %(message)s')
    for h in logging.getLogger().handlers:
        h.setFormatter(fmt)

    """Create keyword dialog."""
    app = QtWidgets.QApplication(sys.argv)
    item = KWL.get_keyword_by_name('*BOUNDARY')
    from gui.window import df
    df.run_master_dialog(item) # 0 = cancel, 1 = ok


if __name__ == '__main__':
    test_dialog()
