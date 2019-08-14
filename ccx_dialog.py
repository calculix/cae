# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    Dialog window to edit CalculiX keyword.
    Called via double click on keyword in the treeView.
    Here we define a keyword's implementation: its name and INP_code.
"""


from PyQt5 import QtWidgets, uic, QtCore, QtWebEngineWidgets
import sys, os, re, logging, ccx_dom


# Load HTML help into QWebEngineView
def saveHTML(item):
    USE_CACHED_HTML = True # if False cached html will NOT be used

    # Get keyword name
    if item.item_type == ccx_dom.item_type.KEYWORD:
        keyword_name = item.name[1:] # cut star
    if item.item_type == ccx_dom.item_type.IMPLEMENTATION:
        keyword_name = item.parent.name[1:] # cut star
    
    # Avoid spaces in html page names
    html_page_name = keyword_name.replace(' ', '_')

    folder = os.path.dirname(sys.argv[0])
    folder = os.path.abspath(folder)
    folder = os.path.join(folder, 'doc')
    url = os.path.join(folder, html_page_name + '.html')

    # Generate html file if it wasn't created previously
    if not os.path.isfile(url) or not USE_CACHED_HTML:

        # Open 'ccx.html' and find link to keyword's page
        href = 'ccx.html'
        with open(os.path.join(folder, href), 'r') as f:
            for line in f.readlines():
                match = re.search('node\d{3}\.html.{3}' + keyword_name, line) # regex to match href
                try:
                    href = match.group(0)[:12]
                    break
                except:
                    pass

        # Read html of the keyword's page
        html = '<html><head><link rel="stylesheet" type="text/css" href="style.css"/></head><body>'
        with open(os.path.join(folder, href), 'r') as f:
            append = False
            cut_breakline = True
            for line in f.readlines():
                if '<!--End of Navigation Panel-->' in line:
                    append = True
                    continue
                if '<HR>' in  line:
                    break
                if '<PRE>' in line:
                    cut_breakline = False
                if '</PRE>' in line:
                    cut_breakline = True
                if append:
                    if cut_breakline:
                        line = line[:-1] + ' ' # replace '\n' with space
                    html += line
        html += '</body></html>'
        html = re.sub('<A.+?\">', '', html) # '?' makes it not greedy
        html = html.replace('</A>', '')
        with open(url, 'w') as f:
            f.write(html)
        
    return url


class Dialog(QtWidgets.QDialog):


    def __init__(self, DOM, item):

        # Create dialog window
        super(Dialog, self).__init__()

        # Load basic form
        uic.loadUi('ccx_dialog.ui', self)

        self.widgets = [] # list of created widgets
        self.item = item # needed to pass to other functions

        # New implementation: draw full form for keyword's arguments
        if self.item.item_type == ccx_dom.item_type.KEYWORD:
            self.setWindowTitle('New ' + self.item.name)

            # For each keyword's argument create name and value widgets
            index = 0 # row number for vertical layout
            for argument in self.item.items:
                if argument.item_type != ccx_dom.item_type.ARGUMENT:
                    continue

                # Remove braces if any
                argument.name = re.sub('[\(\)]', '', argument.name)
                
                # Try to get existing implementations for argument.name
                keyword_name = '*' + argument.name[:-1] # cut '='
                keyword = DOM.getKeywordByName(keyword_name)
                if keyword:
                    """
                        For example, add names of *AMPLITUDE implementations,
                            if argument.name is 'AMPLITUDE'
                    """
                    implementations = [item.name for item in keyword.items \
                        if item.item_type == ccx_dom.item_type.IMPLEMENTATION]
                    if len(implementations):
                        logging.debug('{} {}'.format(keyword.name, implementations))
                        argument.items.extend(implementations)
                        argument.items.insert(0, '') # first row needed empty

                # Argument's values
                if len(argument.items):
                    # Predefined values to be chosen
                    argument_values_widget = QtWidgets.QComboBox()
                    argument_values_widget.addItems(argument.items)

                    # Assign event to update textEdit widget
                    argument_values_widget.currentIndexChanged.connect(self.onChange)

                    # QComboBox doesn't expand by default
                    sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
                    sizePolicy.setHorizontalStretch(1) # expand horizontally
                    argument_values_widget.setSizePolicy(sizePolicy)
                elif argument.name.endswith('='):
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

                # Mutually exclusive arguments
                if '|' in argument.name:

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

                # Save name and values for processing in self.onChange()
                self.widgets.append(argument_name_widget)
                self.widgets.append(argument_values_widget)

                # Add widgets to dialog window
                self.vertical_layout.insertLayout(index, horizontal_layout)
                index += 1 # second time

            # Fill textEdit widget with default keyword's configuration
            self.onChange(None)

        # Edit implementation: draw only textEdit
        if self.item.item_type == ccx_dom.item_type.IMPLEMENTATION:
            self.setWindowTitle('Edit ' + self.item.name)
            for line in self.item.INP_code:
                self.textEdit.append(line)

        # Generate html help page from official manual
        url = saveHTML(item)
        self.doc.load(QtCore.QUrl.fromLocalFile(url)) # load help document

        # Actions
        self.buttonBox.accepted.connect(self.onOk)
        self.buttonBox.button(QtWidgets.QDialogButtonBox.Reset).clicked.connect(self.onReset)


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
        if self.item.item_type == ccx_dom.item_type.KEYWORD:
            string = self.item.name
            for name, value in arguments.items():
                if self.item.from_new_line:
                    string += '\n' + name + value # argument goes from new line
                else:
                    string += ', ' + name + value # argument goes inline
        if self.item.item_type == ccx_dom.item_type.IMPLEMENTATION:
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
        super(Dialog, self).accept()
        return self.textEdit.toPlainText().strip().split('\n')
