# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, UJV Rez, June 2019.
    Distributed under GNU General Public License, version 2.

    Methods to work with main window's treeView widget.
"""


from PyQt5 import QtGui
import ccx_dialog, ccx_dom, ccx_log


class tree:

    def __init__(self, treeView, textEdit):

        # Configure logging
        self.logger = ccx_log.logger(textEdit)

        # Generate CalculiX DOM based on keywords hierarchy from ccx_dom.txt
        try:
            self.DOM = ccx_dom.DOM() # DOM is generated once per session
            self.logger.info('CalculiX object model generated.')
        except:
            self.logger.error('Can\'t generate keywords hierarchy!')

        # Needed for 'doubleClicked' method
        self.treeView = treeView

        # Now generate treeView items
        self.generateTreeView()


    # Recursively generate treeView widget items based on DOM
    def generateTreeView(self):
        self.model = QtGui.QStandardItemModel()
        parent = self.model.invisibleRootItem() # top element in QTreeView
        self.addToTree(parent, self.DOM.root.items) # pass root - group 'Model'
        self.treeView.setModel(self.model)
        self.treeView.expandAll() # expanded looks better


    # Used with generateTreeView() - implements recursion
    def addToTree(self, parent, children):
        """
            parent is QtGui.QStandardItem
            children are items of ccx_dom.DOM object
        """
        for item in children:
            if (item.item_type == 'keyword') or (item.item_type == 'group'):
                tree_element = QtGui.QStandardItem(item.name)
                tree_element.setData(item)
                parent.appendRow(tree_element)

                # Draw keyword's children for its implementation
                for i in range(len(item.implementations)):
                    impl = item.implementations[i] # keyword implementation object
                    e = QtGui.QStandardItem(impl.name)
                    e.setData(impl)
                    tree_element.setText(item.name + ' (' + str(len(item.implementations)) + ')')
                    tree_element.appendRow(e)
                    self.addToTree(e, item.items)

                # Do not draw keyword's children if it doesn't have implementations
                if item.item_type == 'group':
                    self.addToTree(tree_element, item.items)


    # Double click on treeView item: edit the keyword via dialog
    def doubleClicked(self, index):
        item = self.treeView.model().itemFromIndex(index) # treeView item, we obtain it from 'index'
        item = item.data() # now it is ccx_dom.group, ccx_dom.keyword or ccx_dom.implementation 

        # Only double clicking on ccx_dom.keyword creates dialog, not on ccx_dom.group
        if item.item_type == 'keyword':
            dialog = ccx_dialog.Dialog(item) # create dialog window and and pass ccx_dom.keyword object
            if dialog.exec_() == ccx_dialog.Dialog.Accepted: # if user pressed 'OK'

                # The generated piece of .inp code for the CalculiX input file
                INP_code = dialog.onOk()
                for line in INP_code.split('\n'):
                    self.logger.info(line) # show it

                # Create implementation object
                ccx_dom.implementation(item, INP_code)

                # Update treeView widget
                self.generateTreeView()

        # Click on keyword's implementation: show piece of INP_code
        elif item.item_type == 'implementation':
            for line in item.INP_code.split('\n'):
                self.logger.info(line) # show it
