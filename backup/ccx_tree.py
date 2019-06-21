# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, June 2019.
    Distributed under GNU General Public License, version 2.

    Methods to work with main window's treeView widget.
"""


from PyQt5 import QtWidgets, QtCore, QtGui
import ccx_dialog, ccx_dom, ccx_log


class tree:

    def __init__(self, treeView, textEdit):

        # Configure logging
        self.logger = ccx_log.logger(textEdit)

        # Generate CalculiX DOM based on keywords hierarchy from ccx_dom.txt
        self.DOM = ccx_dom.DOM(textEdit) # DOM is generated once per session during start

        # Make treeView available for other methods
        self.treeView = treeView

        # Expanded / collapsed flag - used in other methods
        self.expanded = False # start with collapsed tree items

        # Hide / show tree keywords without implementations 
        self.show_empty = False # start with show all tree items

        # Now generate treeView items
        self.generateTreeView()

        # Actions
        self.treeView.doubleClicked.connect(self.doubleClicked)
        self.treeView.customContextMenuRequested.connect(self.rightClicked)


    # Recursively generate treeView widget items based on DOM
    def generateTreeView(self):
        self.model = QtGui.QStandardItemModel()
        parent = self.model.invisibleRootItem() # top element in QTreeView
        self.addToTree(parent, self.DOM.root.items) # pass root - group 'Model'
        self.treeView.setModel(self.model)

        # Keep expanded/collapsed state
        if self.expanded:
            self.treeView.expandAll()
        else:
            self.treeView.collapseAll()

        # self.logger.info('TreeView generated.')


    # Used with generateTreeView() - implements recursion
    def addToTree(self, parent, children):
        """
            parent is QtGui.QStandardItem
            children are items of ccx_dom.DOM group or keyword
        """
        for item in children:

            if item.item_type == 'group':

                # Check if there are keywords with implementations
                is_empty = True
                for keyword in item.items:
                    if len(keyword.implementations):
                        is_empty = False
                        break

                if not is_empty:
                    tree_element = QtGui.QStandardItem(item.name)
                    tree_element.setData(item)

                    # Add icon to each group in tree
                    icon_name = item.name.replace('*', '') + '.png'
                    icon = QtGui.QIcon('./icons/' + icon_name.lower())
                    tree_element.setIcon(icon)

                    parent.appendRow(tree_element)
                    self.addToTree(tree_element, item.items)

            if item.item_type == 'keyword':
                if self.show_empty or len(item.implementations):
                    tree_element = QtGui.QStandardItem(item.name)
                    tree_element.setData(item)

                    # Add icon to each keyword in tree
                    icon_name = item.name.replace('*', '') + '.png'
                    icon = QtGui.QIcon('./icons/' + icon_name.lower())
                    tree_element.setIcon(icon)

                    parent.appendRow(tree_element)
                    if not len(item.implementations):
                        self.addToTree(tree_element, item.items)

                # Draw keyword's children for its implementation
                for i in range(len(item.implementations)):
                    impl = item.implementations[i] # keyword implementation object
                    e = QtGui.QStandardItem(impl.name)
                    e.setData(impl)

                    tree_element.setText(item.name + ' (' + str(len(item.implementations)) + ')')
                    tree_element.appendRow(e)
                    self.addToTree(e, item.items)


    # Double click on treeView item: edit the keyword via dialog
    def doubleClicked(self, index):
        item = self.treeView.model().itemFromIndex(index) # treeView item obtained from 'index'
        item = item.data() # now it is ccx_dom.group, ccx_dom.keyword or ccx_dom.implementation 

        # Double click on ccx_dom.group doesn't create dialog
        if item.item_type == 'keyword' \
            or item.item_type == 'implementation':

            # Create dialog window and pass item
            dialog = ccx_dialog.Dialog(item)

            # Get response from dialog window
            if dialog.exec_() == ccx_dialog.Dialog.Accepted: # if user pressed 'OK'

                # The generated piece of .inp code for the CalculiX input file
                INP_code = dialog.onOk() # list of strings
                for line in INP_code:
                    self.logger.info(line) # show it

                # Create implementation object for keyword
                if item.item_type == 'keyword':
                    ccx_dom.implementation(item, INP_code)

                # Replace implementation object with a new one
                elif item.item_type == 'implementation':
                    keyword = item.keyword
                    item.remove()
                    ccx_dom.implementation(keyword, INP_code)

                # Update treeView widget
                self.generateTreeView()


    # Context menu for right click
    def rightClicked(self):
        self.myMenu = QtWidgets.QMenu('Menu', self.treeView)

        action_show_hide = QtWidgets.QAction('Show/Hide empty containers', self.treeView)
        self.myMenu.addAction(action_show_hide)
        action_show_hide.triggered.connect(self.actionShowHide)

        action_expand_collapse = QtWidgets.QAction('Expand/Collapse all', self.treeView)
        self.myMenu.addAction(action_expand_collapse)
        action_expand_collapse.triggered.connect(self.actionExpandCollapse)

        try: # catch out of index
            index = self.treeView.selectedIndexes()[0] # selected item index
            item = self.treeView.model().itemFromIndex(index) # treeView item obtained from 'index'
            item = item.data() # now it is ccx_dom.group, ccx_dom.keyword or ccx_dom.implementation 
            if item.item_type == 'implementation':
                action_delete_implementation = QtWidgets.QAction('Delete', self.treeView)
                self.myMenu.addAction(action_delete_implementation)
                action_delete_implementation.triggered.connect(self.actionDeleteImplementation)
        except:
            pass

        self.myMenu.exec_(QtGui.QCursor.pos())


    # Show/Hide empty treeView items
    def actionShowHide(self):
        self.show_empty = not(self.show_empty)
        self.generateTreeView()


    # Expand or collapsea all treeView items
    def actionExpandCollapse(self):
        if self.expanded:
            self.treeView.collapseAll()
        else:
            self.treeView.expandAll()
        self.expanded = not(self.expanded)


    # Delete keyword's implementation from DOM
    def actionDeleteImplementation(self, item):
        index = self.treeView.selectedIndexes()[0] # selected item index
        item = self.treeView.model().itemFromIndex(index) # treeView item obtained from 'index'
        item = item.data() # now it is ccx_dom.group, ccx_dom.keyword or ccx_dom.implementation 

        if item.item_type == 'implementation':
            item.remove()

            # Update treeView items
            self.generateTreeView()
