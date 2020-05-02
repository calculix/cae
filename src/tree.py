#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, December 2019
Distributed under GNU General Public License v3.0

Methods to work with main window's treeView widget. """

# Standard modules
import re
import os
import sys
import logging

# External modules
from PyQt5 import QtWidgets, QtCore, QtGui

# My modules
import gui
import model
from model.kom import item_type, implementation

class Tree:

    """
    p - Path
    s - Settings
    w - Window
    m - Model
    """
    def __init__(self, p, s, w, m):
        self.p = p
        self.s = s
        self.w = w
        self.m = m
        self.model = QtGui.QStandardItemModel()
        self.w.treeView.setModel(self.model)

    # Delete keyword's implementation in the treeView by pressing 'Delete' button
    def keyPressEvent(self, e):
        if e.key() == QtCore.Qt.Key_Delete:
            self.actionDeleteImplementation()

    # Recursively generate treeView widget items based on KOM
    def generateTreeView(self, m):
        self.model.clear() # remove all items and data from tree
        parent_element = self.model.invisibleRootItem() # top element in QTreeView
        self.addToTree(parent_element, self.m.KOM.root.items) # pass top level groups

    # Used with generateTreeView() - implements recursion
    def addToTree(self, parent_element, items):

        for item in items:

            # Add to the tree only needed item_types
            if item.item_type not \
                in [item_type.GROUP,
                    item_type.KEYWORD,
                    item_type.IMPLEMENTATION]:
                continue

            # Check if there are keywords with implementations
            if self.s.show_empty_keywords \
                    or item.name == 'Job' \
                    or item.countImplementations() \
                    or item.item_type == item_type.IMPLEMENTATION:

                # Create tree_element
                tree_element = QtGui.QStandardItem(item.name)
                tree_element.setData(item)
                parent_element.appendRow(tree_element)

                # Set text color for tree_element
                brush = QtGui.QBrush()
                brush.setColor(QtCore.Qt.gray)
                if item.isActive():
                    brush.setColor(QtCore.Qt.black)
                    # Bold font for implementations
                    if item.item_type == item_type.IMPLEMENTATION:
                        font = QtGui.QFont()
                        font.setBold(True)
                        tree_element.setFont(font)
                tree_element.setForeground(brush)

                # Expand / collapse
                if item.expanded:
                    self.w.treeView.expand(tree_element.index())
                else:
                    self.w.treeView.collapse(tree_element.index())

                # Add icon to each keyword in tree
                icon_name = item.name.replace('*', '') + '.png'
                icon_name = icon_name.replace(' ', '_')
                icon_name = icon_name.replace('-', '_')
                icon_path = os.path.join(self.p.img, 'icon_' + icon_name.lower())
                icon = QtGui.QIcon(icon_path)
                tree_element.setIcon(icon)

                # # Append job name
                # if item.name == 'Job':
                #     self.job_element = tree_element
                #     self.appendJobName()

                # Organize recursion
                impls = item.getImplementations()
                if len(impls):
                    self.addToTree(tree_element, impls)
                else:
                    self.addToTree(tree_element, item.items)

    # Double click on treeView item: edit the keyword via dialog
    def doubleClicked(self):
        index = self.w.treeView.selectedIndexes()[0] # selected item index
        tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
        item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION

        # Double click on GROUP doesn't create dialog
        if item and item.item_type in \
            [item_type.KEYWORD, item_type.IMPLEMENTATION]:

            if item.active:

                # Create dialog window and pass item
                dialog = gui.keyword_dialog.KeywordDialog(self.m.KOM, item)

                # Get response from dialog window
                if dialog.exec() == gui.keyword_dialog.KeywordDialog.Accepted: # if user pressed 'OK'

                    # The generated piece of .inp code for the CalculiX input file
                    INP_code = dialog.onOk() # list of strings

                    # Create implementation object for keyword
                    if item.item_type == item_type.KEYWORD:
                        impl = implementation(item, INP_code) # create keyword's implementation

                        # Regenerate tree_element's children
                        tree_element.removeRows(0, tree_element.rowCount()) # remove all children
                        self.addToTree(tree_element, item.getImplementations()) # add only implementations

                        # Reparse mesh or constraints
                        # self.m.Mesh.reparse(INP_code)
                        reparsed = model.parsers.mesh.Mesh(INP_code=INP_code, old=self.m.Mesh)
                        self.m.Mesh.updateWith(reparsed)
                        self.clicked() # rehighlight

                    # Replace implementation object with a new one
                    elif item.item_type == item_type.IMPLEMENTATION:
                        # Remove old one
                        parent = tree_element.parent() # parent treeView item
                        keyword = parent.data() # parent keyword for implementation
                        keyword.items.remove(item) # remove implementation from keyword's items

                        # Add new one
                        impl = implementation(keyword, INP_code, name=item.name)
                        tree_element.setData(impl)

                        # Reparse mesh or constraints
                        reparsed = model.parsers.mesh.Mesh(INP_code=INP_code, old=self.m.Mesh)
                        self.m.Mesh.updateWith(reparsed)
                        self.clicked() # rehighlight

            else:
                logging.warning('Please, create ' + item.getParentKeywordName() + ' first.')

    # Highlight node sets, element sets or surfaces
    def clicked(self):

        # Debug for Ctrl+Click
        if not len(self.w.treeView.selectedIndexes()):
            return

        index = self.w.treeView.selectedIndexes()[0] # selected item index
        tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
        item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION

        # Highlight entities
        if item and item.item_type == item_type.IMPLEMENTATION:
            ipn_up = item.parent.name.upper()
            lead_line = item.INP_code[0]
            _set = []

            # Highlight mesh entities
            if ipn_up == '*NSET' or ipn_up == '*NODE':
                match = re.search('NSET\s*=\s*([\w\-]*)', lead_line.upper())
                if match: # if there is NSET attribute
                    name = lead_line[match.start(1):match.end(1)] # node set name
                    if name in self.m.Mesh.nsets:
                        self.w.post('plot n ' + name)

            elif ipn_up == '*ELSET' or ipn_up == '*ELEMENT':
                match = re.search('ELSET\s*=\s*([\w\-]*)', lead_line.upper())
                if match: # if there is ELSET attribute
                    name = lead_line[match.start(1):match.end(1)] # element set name
                    if name in self.m.Mesh.elsets:
                        self.w.post('plot e ' + name)

            elif ipn_up == '*SURFACE':

                # Surface type - optional attribute
                stype = 'ELEMENT' # 'ELEMENT' or 'NODE'
                match = re.search('\*SURFACE\s*,.*TYPE\s*=\s*(\w*)', lead_line.upper())
                if match:
                    stype = lead_line[match.start(1):match.end(1)]

                match = re.search('NAME\s*=\s*([\w\-]*)', lead_line.upper())
                name = lead_line[match.start(1):match.end(1)] # surface name
                if stype == 'ELEMENT':
                    self.w.post('plot f ' + name)
                elif stype=='NODE':
                    self.w.post('plot f ' + name)

            # Highlight Loads & BC
            elif ipn_up in ['*BOUNDARY', '*CLOAD', '*CFLUX']:
                for line in item.INP_code[1:]:
                    line = line.strip()
                    n = line.replace(',', ' ').split()[0]
                    try:
                        # Single node number
                        _set.append(int(n))
                    except ValueError as err:
                        # Nodes in node set
                        _set.extend([n.num for n in self.m.Mesh.nsets[n].items])
                        pass
                # self.w.VTK.highlight(set(_set), 1) # 1 = vtk.vtkSelectionNode.POINT

        # else:
        #     self.w.deselect_cgx_sets()

    # Context menu for right click
    def rightClicked(self):
        self.myMenu = QtWidgets.QMenu('Menu', self.w.treeView)

        try:
            index = self.w.treeView.selectedIndexes()[0] # selected item index
            tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
            item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION

            # Context menu for any keyword and implementations
            if item:
                if item.item_type == item_type.IMPLEMENTATION:

                    # 'Edit' action
                    action_edit_implementation = QtWidgets.QAction('Edit', self.w.treeView)
                    self.myMenu.addAction(action_edit_implementation)
                    action_edit_implementation.triggered.connect(self.doubleClicked)

                    # 'Delete' action
                    action_delete_implementation = QtWidgets.QAction('Delete', self.w.treeView)
                    self.myMenu.addAction(action_delete_implementation)
                    action_delete_implementation.triggered.connect(self.actionDeleteImplementation)

                if item.item_type == item_type.KEYWORD:

                    # 'Create' action
                    action_create_implementation = QtWidgets.QAction('Create', self.w.treeView)
                    self.myMenu.addAction(action_create_implementation)
                    action_create_implementation.triggered.connect(self.doubleClicked)

            # Add splitter
            self.myMenu.addSeparator()

        except IndexError:
            pass

        # Context menu elements which always present
        if self.s.show_empty_keywords:
            action_show_hide = QtWidgets.QAction('Hide empty containers', self.w.treeView)
        else:
            action_show_hide = QtWidgets.QAction('Show empty containers', self.w.treeView)
        self.myMenu.addAction(action_show_hide)
        action_show_hide.triggered.connect(self.actionShowHide)

        action_expand_collapse = QtWidgets.QAction('Collapse all', self.w.treeView)
        self.myMenu.addAction(action_expand_collapse)
        action_expand_collapse.triggered.connect(self.actionCollapseAll)

        action_expand_collapse = QtWidgets.QAction('Expand all', self.w.treeView)
        self.myMenu.addAction(action_expand_collapse)
        action_expand_collapse.triggered.connect(self.actionExpandAll)

        self.myMenu.exec_(QtGui.QCursor.pos())

    # Show/Hide empty treeView items
    def actionShowHide(self):
        self.s.show_empty_keywords = not(self.s.show_empty_keywords)
        self.s.save() # save 'show_empty_keywords' value in settings
        self.generateTreeView(self.m)

    # Expand or collapse all treeView items
    def actionCollapseAll(self):
        self.w.treeView.collapseAll()
        self.s.expanded = False
        self.s.save()
    def actionExpandAll(self):
        self.w.treeView.expandAll()
        self.s.expanded = True
        self.s.save()

    # Delete keyword's implementation from KOM
    def actionDeleteImplementation(self):
        index = self.w.treeView.selectedIndexes()[0] # selected item index
        tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
        item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION

        if item and item.item_type == item_type.IMPLEMENTATION:

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

                # Regenerate parent's children
                parent.removeRows(0, parent.rowCount()) # remove all children
                self.addToTree(parent, keyword.items)

                if not self.s.show_empty_keywords:
                    self.hideParent(parent)
    def hideParent(self, tree_element):

        # To hide current item/brunch it should be empty 'keyword' or 'group'
        if not self.s.show_empty_keywords \
            and not tree_element.hasChildren() \
            and tree_element.data().item_type != item_type.IMPLEMENTATION:

            # Hide current item/brunch from tree via calling parent.removeRow
            parent = tree_element.parent()
            if not parent:
                parent = self.model.invisibleRootItem()
            parent.removeRow(tree_element.row())

            if parent != self.model.invisibleRootItem():
                self.hideParent(parent)

    # Change KOM item's 'expanded' variable when user interacts with treeView
    def treeViewExpanded(self, index):
        tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
        item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION
        if item:
            item.expanded = True
    def treeViewCollapsed(self, index):
        tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
        item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION
        if item:
            item.expanded = False
