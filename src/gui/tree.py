#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Methods to work with main window treeView widget.
"""

# Standard modules
import re
import os
import sys
import logging

# External modules
from PyQt5 import QtWidgets, QtCore, QtGui

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
from model.parsers.mesh import Mesh
from model.kom import ItemType, Implementation
from model.kom import KOM
from model import m
from gui.window import factory
import tests


class Tree:

    def __init__(self):
        self.model = QtGui.QStandardItemModel()
        if hasattr(factory.mw, 'treeView'):
            factory.mw.treeView.setModel(self.model)

    def keyPressEvent(self, e):
        """Delete keyword implementation in the
        treeView by pressing 'Delete' button.
        """
        if e.key() == QtCore.Qt.Key_Delete:
            self.actionDeleteImplementation()

    def generateTreeView(self):
        """Recursively generate treeView widget items based on KOM."""
        self.model.clear() # remove all items and data from tree
        parent_element = self.model.invisibleRootItem() # top element in QTreeView
        self.addToTree(parent_element, KOM.root.items) # pass top level groups

    def addToTree(self, parent_element, items):
        """Used with generateTreeView() - implements recursion."""
        for item in items:

            # Add to the tree only needed item_types
            allowed = [ItemType.GROUP, ItemType.KEYWORD,
                ItemType.IMPLEMENTATION]
            if item.itype not in allowed:
                continue

            # Check if there are keywords with implementations
            if not s.show_empty_keywords \
                and not item.count_implementations() \
                and item.itype != ItemType.IMPLEMENTATION:
                continue

            # Create tree_element
            tree_element = QtGui.QStandardItem(item.name)
            tree_element.setData(item)
            parent_element.appendRow(tree_element)

            # Set text color for tree_element
            brush = QtGui.QBrush()
            brush.setColor(QtCore.Qt.gray)
            if item.is_active():
                brush.setColor(QtCore.Qt.black)
                # Bold font for implementations
                if item.itype == ItemType.IMPLEMENTATION:
                    font = QtGui.QFont()
                    font.setBold(True)
                    tree_element.setFont(font)
            tree_element.setForeground(brush)

            # Expand / collapse
            if item.expanded:
                factory.mw.treeView.expand(tree_element.index())
            else:
                factory.mw.treeView.collapse(tree_element.index())

            # Add icon to each keyword in tree
            icon_name = item.name.replace('*', '') + '.png'
            icon_name = icon_name.replace(' ', '_')
            icon_name = icon_name.replace('-', '_')
            icon_path = os.path.join(p.img, 'icon_' + icon_name.lower())
            icon = QtGui.QIcon(icon_path)
            tree_element.setIcon(icon)

            # Organize recursion
            impls = item.get_implementations()
            if len(impls):
                self.addToTree(tree_element, impls)
            else:
                self.addToTree(tree_element, item.items)

    def doubleClicked(self):
        """Double click on treeView item: edit the keyword via dialog."""
        index = factory.mw.treeView.selectedIndexes()[0] # selected item index
        tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
        item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION

        # Double click on GROUP doesn't create dialog
        allowed_types = [ItemType.KEYWORD, ItemType.IMPLEMENTATION]
        if not item or item.itype not in allowed_types:
            return

        if not item.active:
            kw_name = item.get_parent_keyword_name()
            msg = 'Please, create {} first.'.format(kw_name)
            logging.warning(msg)
            return

        # Exec dialog and recieve answer
        # Process response from dialog window if user pressed 'OK'
        if factory.run_master_dialog(item): # 0 = cancel, 1 = ok

            # The generated piece of .inp code for the CalculiX input file
            inp_code = factory.mw.ok() # list of strings

            # Create implementation object for keyword
            if item.itype == ItemType.KEYWORD:
                impl = Implementation(item, inp_code) # create keyword implementation

                # Regenerate tree_element children
                tree_element.removeRows(0, tree_element.rowCount()) # remove all children
                self.addToTree(tree_element, item.get_implementations()) # add only implementations

                # Reparse mesh or constraints
                # TODO Use Mesh without m
                # m.Mesh.reparse(inp_code)
                reparsed = Mesh(icode=inp_code, old=m.Mesh)
                if m.Mesh is not None:
                    m.Mesh.updateWith(reparsed)
                self.clicked() # rehighlight

            # Replace implementation object with a new one
            elif item.itype == ItemType.IMPLEMENTATION:
                # Remove old one
                parent = tree_element.parent() # parent treeView item
                keyword = parent.data() # parent keyword for implementation
                keyword.items.remove(item) # remove implementation from keyword items

                # Add new one
                impl = Implementation(keyword, inp_code, name=item.name)
                tree_element.setData(impl)

                # Reparse mesh or constraints
                reparsed = Mesh(icode=inp_code, old=m.Mesh)
                m.Mesh.updateWith(reparsed)
                self.clicked() # rehighlight

    def clicked(self):
        """Highlight node sets, element sets or surfaces."""
        
        if not hasattr(factory.mw, 'treeView'):
            return

        # Debug for Ctrl+Click
        if not len(factory.mw.treeView.selectedIndexes()):
            return

        # Highlight only when INP is opened
        if factory.sw is None:
            return
        if not (p.path_cgx + ' -c ') in factory.sw.cmd:
            return

        index = factory.mw.treeView.selectedIndexes()[0] # selected item index
        tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
        item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION

        # Highlight entities
        if item and item.itype == ItemType.IMPLEMENTATION:
            ipn_up = item.parent.name.upper()
            _set = []

            i = 0
            while item.inp_code[i].startswith('**'):
                i += 1
            lead_line = item.inp_code[i]

            # Highlight mesh entities
            if ipn_up == '*NSET' or ipn_up == '*NODE':
                match = re.search('NSET\s*=\s*([\w\!\#\%\$\&\"\'\(\)\*\=\+\-\.\/\:\;\<\>\?\@\[\]\^\_\`\{\\\|\}\~]*)', lead_line.upper())
                if match: # if there is NSET attribute
                    name = lead_line[match.start(1):match.end(1)] # node set name
                    if name in m.Mesh.nsets:
                        factory.connection.post('plot n ' + name)

            elif ipn_up == '*ELSET' or ipn_up == '*ELEMENT':
                match = re.search('ELSET\s*=\s*([\w\!\#\%\$\&\"\'\(\)\*\=\+\-\.\/\:\;\<\>\?\@\[\]\^\_\`\{\\\|\}\~]*)', lead_line.upper())
                if match: # if there is ELSET attribute
                    name = lead_line[match.start(1):match.end(1)] # element set name
                    if name in m.Mesh.elsets:
                        factory.connection.post('plot e ' + name)

            elif ipn_up == '*SURFACE':

                # Surface type - optional attribute
                stype = 'ELEMENT' # 'ELEMENT' or 'NODE'
                match = re.search('\*SURFACE\s*,.*TYPE\s*=\s*(\w*)', lead_line.upper())
                if match:
                    stype = lead_line[match.start(1):match.end(1)]

                match = re.search('NAME\s*=\s*([\w\!\#\%\$\&\"\'\(\)\*\=\+\-\.\/\:\;\<\>\?\@\[\]\^\_\`\{\\\|\}\~]*)', lead_line.upper())
                name = lead_line[match.start(1):match.end(1)] # surface name
                if stype == 'ELEMENT':
                    factory.connection.post('plot f ' + name)
                elif stype=='NODE':
                    factory.connection.post('plot f ' + name)

            # Highlight Loads & BC
            elif ipn_up in ['*BOUNDARY', '*CLOAD', '*CFLUX']:
                for line in item.inp_code[1:]:
                    line = line.strip()
                    n = line.replace(',', ' ').split()[0]
                    try:
                        # Single node number
                        _set.append(int(n))
                    except ValueError as err:
                        # Nodes in node set
                        _set.extend([n.num for n in m.Mesh.nsets[n].items])
                        pass
                # factory.mw.VTK.highlight(set(_set), 1) # 1 = vtk.vtkSelectionNode.POINT

        # else:
        #     factory.mw.deselect_cgx_sets()

    def rightClicked(self):
        """Context menu for right click."""
        self.myMenu = QtWidgets.QMenu('Menu', factory.mw.treeView)

        try:
            index = factory.mw.treeView.selectedIndexes()[0] # selected item index
            tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
            item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION

            # Context menu for any keyword and implementations
            if item:
                if item.itype == ItemType.IMPLEMENTATION:

                    # 'Edit' action
                    action_edit_implementation = QtWidgets.QAction('Edit', factory.mw.treeView)
                    self.myMenu.addAction(action_edit_implementation)
                    action_edit_implementation.triggered.connect(self.doubleClicked)

                    # 'Delete' action
                    action_delete_implementation = QtWidgets.QAction('Delete', factory.mw.treeView)
                    self.myMenu.addAction(action_delete_implementation)
                    action_delete_implementation.triggered.connect(self.actionDeleteImplementation)

                if item.itype == ItemType.KEYWORD:

                    # 'Create' action
                    action_create_implementation = QtWidgets.QAction('Create', factory.mw.treeView)
                    self.myMenu.addAction(action_create_implementation)
                    action_create_implementation.triggered.connect(self.doubleClicked)

            # Add splitter
            self.myMenu.addSeparator()

        except IndexError:
            pass

        # Context menu elements which always present
        if s.show_empty_keywords:
            action_show_hide = QtWidgets.QAction('Hide empty containers', factory.mw.treeView)
        else:
            action_show_hide = QtWidgets.QAction('Show empty containers', factory.mw.treeView)
        self.myMenu.addAction(action_show_hide)
        action_show_hide.triggered.connect(self.actionShowHide)

        action_collapse_expand = QtWidgets.QAction('Collapse/expand all', factory.mw.treeView)
        self.myMenu.addAction(action_collapse_expand)
        action_collapse_expand.triggered.connect(self.action_collapse_expand_all)

        self.myMenu.exec_(QtGui.QCursor.pos())

    def actionShowHide(self):
        """Show/Hide empty treeView items."""
        s.show_empty_keywords = not(s.show_empty_keywords)
        s.save() # save 'show_empty_keywords' value in settings
        self.generateTreeView()

    def action_collapse_expand_all(self):
        """Expand or collapse all treeView items."""
        if s.expanded:
            factory.mw.treeView.collapseAll()
        else:
            factory.mw.treeView.expandAll()
        s.expanded = not s.expanded
        s.save()

    def actionDeleteImplementation(self):
        """Delete keyword implementation from KOM."""

        def hide_parent(tree_element):
            # To hide current item/brunch it should be empty 'keyword' or 'group'
            if not s.show_empty_keywords \
                and not tree_element.hasChildren() \
                and tree_element.data().itype != ItemType.IMPLEMENTATION:

                # Hide current item/brunch from tree via calling parent.removeRow
                parent = tree_element.parent()
                if not parent:
                    parent = self.model.invisibleRootItem()
                parent.removeRow(tree_element.row())

                if parent != self.model.invisibleRootItem():
                    hide_parent(parent)

        index = factory.mw.treeView.selectedIndexes()[0] # selected item index
        tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
        item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION

        if item and item.itype == ItemType.IMPLEMENTATION:

            # Confirmation dialog to delete implementation
            answer = QtWidgets.QMessageBox.question(None,
                item.name, 'OK to delete ' + item.name + '?',
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                QtWidgets.QMessageBox.Yes)

            # If confirmed
            if answer == QtWidgets.QMessageBox.Yes:

                parent = tree_element.parent() # parent treeView item
                keyword = parent.data() # parent keyword for implementation
                keyword.items.remove(item) # remove implementation from keyword items

                # Regenerate parent children
                parent.removeRows(0, parent.rowCount()) # remove all children
                self.addToTree(parent, keyword.items)

                if not s.show_empty_keywords:
                    hide_parent(parent)

    def expanded_or_collapsed(self, index):
        """Change KOM item 'expanded' variable when user interacts with treeView."""
        tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
        item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION
        if item:
            item.expanded = not item.expanded


# Create treeView items based on Keyword Object Model
t = Tree()


@tests.test_wrapper()
def test():
    """Start keyword editor and generate the tree."""
    import log
    log.stop_logging()

    # Create application
    app = QtWidgets.QApplication(sys.argv)

    # Show main window
    factory.run_master()

    # Create treeView items based on KOM
    t = Tree()

    # Actions
    if hasattr(factory.mw, 'treeView'):
        factory.mw.treeView.doubleClicked.connect(t.doubleClicked)
        factory.mw.treeView.clicked.connect(t.clicked)
        factory.mw.treeView.customContextMenuRequested.connect(t.rightClicked)
        factory.mw.treeView.expanded.connect(t.expanded_or_collapsed)
        factory.mw.treeView.collapsed.connect(t.expanded_or_collapsed)

    t.generateTreeView()

    # Execute application
    app.exec()


if __name__ == '__main__':
    test() # run test
