# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, June 2019.
    Distributed under GNU General Public License, version 2.

    Methods to work with main window's treeView widget.
"""


from PyQt5 import QtWidgets, QtCore, QtGui
import ccx_dialog, ccx_log
from ccx_dom import *


class tree:

    def __init__(self, CAE):
        self.DOM = CAE.DOM
        self.treeView = CAE.treeView

        # Hide / show tree keywords without implementations 
        self.show_empty = False # start with show all tree items

        # Now generate treeView items
        self.model = QtGui.QStandardItemModel()
        self.generateTreeView()

        # Actions
        self.treeView.doubleClicked.connect(self.doubleClicked)
        self.treeView.customContextMenuRequested.connect(self.rightClicked)
        self.treeView.expanded.connect(self.treeViewExpanded)
        self.treeView.collapsed.connect(self.treeViewCollapsed)


    # Recursively generate treeView widget items based on DOM
    def generateTreeView(self):
        self.model.clear() # remove all items and data
        branch = self.model.invisibleRootItem() # top element in QTreeView
        self.addToTree(branch, self.DOM.root.items) # pass root groups
        self.treeView.setModel(self.model)


    # Used with generateTreeView() - implements recursion
    def addToTree(self, branch, items):

        for item in items:

            # Add to the tree only needed item_types
            if item.item_type not in \
                [item_type.GROUP, item_type.KEYWORD, item_type.IMPLEMENTATION]:
                continue

            # Check if there are keywords with implementations
            if self.show_empty \
                    or item.countImplementations() \
                    or item.item_type == item_type.IMPLEMENTATION:

                tree_element = QtGui.QStandardItem(item.name)
                tree_element.setData(item)
                branch.appendRow(tree_element)

                # Expand / collapse
                if item.expanded:
                    self.treeView.expand(tree_element.index())
                else:
                    self.treeView.collapse(tree_element.index())

                # Add icon to each keyword in tree
                icon_name = item.name.replace('*', '') + '.png'
                icon = QtGui.QIcon('./icons/' + icon_name.lower())
                tree_element.setIcon(icon)

                # Organize recursion
                impls = item.getImplementations()
                if len(impls):
                    self.addToTree(tree_element, impls)
                else:
                    self.addToTree(tree_element, item.items)


    # Double click on treeView item: edit the keyword via dialog
    def doubleClicked(self):
        index = self.treeView.selectedIndexes()[0] # selected item index
        tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
        item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION

        # Double click on GROUP doesn't create dialog
        if item.item_type == item_type.KEYWORD or item.item_type == item_type.IMPLEMENTATION:

            # Create dialog window and pass item
            dialog = ccx_dialog.Dialog(item)

            # Get response from dialog window
            if dialog.exec_() == ccx_dialog.Dialog.Accepted: # if user pressed 'OK'

                # The generated piece of .inp code for the CalculiX input file
                INP_code = dialog.onOk() # list of strings

                # Create implementation object for keyword
                if item.item_type == item_type.KEYWORD:
                    impl = implementation(item, INP_code)
                    impl_element = QtGui.QStandardItem(impl.name)
                    impl_element.setData(impl)
                    tree_element.appendRow(impl_element)

                # Replace implementation object with a new one
                elif item.item_type == item_type.IMPLEMENTATION:
                    # Remove old one
                    parent = tree_element.parent() # parent treeView item
                    keyword = parent.data() # parent keyword for implementation
                    keyword.items.remove(item) # remove implementation from keyword's items

                    # Add new one
                    impl = implementation(keyword, INP_code, name=item.name)
                    tree_element.setData(impl)


    # Context menu for right click
    def rightClicked(self):
        self.myMenu = QtWidgets.QMenu('Menu', self.treeView)

        try: # catch out of index
            index = self.treeView.selectedIndexes()[0] # selected item index
            item = self.model.itemFromIndex(index) # treeView item obtained from 'index'
            item = item.data() # now it is GROUP, KEYWORD or IMPLEMENTATION
            if item.item_type == item_type.IMPLEMENTATION:

                # 'Edit' action
                action_edit_implementation = QtWidgets.QAction('Edit', self.treeView)
                self.myMenu.addAction(action_edit_implementation)
                action_edit_implementation.triggered.connect(self.doubleClicked)

                # 'Delete' action
                action_delete_implementation = QtWidgets.QAction('Delete', self.treeView)
                self.myMenu.addAction(action_delete_implementation)
                action_delete_implementation.triggered.connect(self.actionDeleteImplementation)

            if item.item_type == item_type.KEYWORD:

                # 'Create' action
                action_create_implementation = QtWidgets.QAction('Create', self.treeView)
                self.myMenu.addAction(action_create_implementation)
                action_create_implementation.triggered.connect(self.doubleClicked)

            # Add splitter
            self.myMenu.addSeparator()

        except:
            pass

        if self.show_empty:
            title = 'Hide empty containers'
        else:
            title = 'Show empty containers'
        action_show_hide = QtWidgets.QAction(title, self.treeView)
        self.myMenu.addAction(action_show_hide)
        action_show_hide.triggered.connect(self.actionShowHide)

        action_expand_collapse = QtWidgets.QAction('Collapse all', self.treeView)
        self.myMenu.addAction(action_expand_collapse)
        action_expand_collapse.triggered.connect(self.actionCollapseAll)

        action_expand_collapse = QtWidgets.QAction('Expand all', self.treeView)
        self.myMenu.addAction(action_expand_collapse)
        action_expand_collapse.triggered.connect(self.actionExpandAll)

        self.myMenu.exec_(QtGui.QCursor.pos())


    # Show/Hide empty treeView items
    def actionShowHide(self):
        self.show_empty = not(self.show_empty)
        self.generateTreeView()


    # Expand or collapse all treeView items
    def actionCollapseAll(self):
        self.treeView.collapseAll()
    def actionExpandAll(self):
        self.treeView.expandAll()


    # Delete keyword's implementation from DOM
    def actionDeleteImplementation(self):
        index = self.treeView.selectedIndexes()[0] # selected item index
        tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
        item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION

        if item.item_type == item_type.IMPLEMENTATION:

            # Confirmation dialog to delete implementation
            answer = QtWidgets.QMessageBox.question(None,
                item.name, 'OK to delete ' + item.name + '?',
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                QtWidgets.QMessageBox.Yes)

            # If confirmed
            if answer == QtWidgets.QMessageBox.Yes:

                parent = tree_element.parent() # parent treeView item
                keyword = parent.data() # parent keyword for implementation
                keyword.items.remove(item) # remove implementation from keyword's items
                parent.removeRow(tree_element.row()) # remove row in treeView

                # Hide empty branch from tree
                def hideParent(branch):

                    # To hide current item/brunch it should be empty 'keyword' or 'group'
                    if not self.show_empty \
                        and not branch.hasChildren() \
                        and branch.data().item_type != item_type.IMPLEMENTATION:

                        # Hide current item/brunch from tree via calling parent.removeRow
                        parent = branch.parent()
                        if not parent:
                            parent = self.model.invisibleRootItem()
                        parent.removeRow(branch.row())

                        if parent != self.model.invisibleRootItem():
                            hideParent(parent)

                hideParent(parent)


    # Change DOM item's 'expanded' variable when user interacts with treeView
    def treeViewExpanded(self, index):
        tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
        item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION
        item.expanded = True
    def treeViewCollapsed(self, index):
        tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
        item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION
        item.expanded = False
