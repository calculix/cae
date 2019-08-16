# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    Methods to work with main window's treeView widget.
    Depends on ccx_cae.CAE.
"""


from PyQt5 import QtWidgets, QtCore, QtGui
import ccx_dialog, ccx_dom, re, os, logging


class tree:


    def __init__(self, CAE):
        self.CAE = CAE

        # Read application's global settings
        self.settings = CAE.settings

        # Now generate treeView items
        self.model = QtGui.QStandardItemModel()
        self.CAE.treeView.setModel(self.model)
        self.generateTreeView()

        # Actions
        self.CAE.treeView.doubleClicked.connect(self.doubleClicked)
        self.CAE.treeView.clicked.connect(self.clicked)
        self.CAE.treeView.customContextMenuRequested.connect(self.rightClicked)
        self.CAE.treeView.expanded.connect(self.treeViewExpanded)
        self.CAE.treeView.collapsed.connect(self.treeViewCollapsed)


    # Recursively generate treeView widget items based on DOM
    def generateTreeView(self):
        self.model.clear() # remove all items and data from tree
        parent_element = self.model.invisibleRootItem() # top element in QTreeView
        self.addToTree(parent_element, self.CAE.DOM.root.items) # pass top level groups


    # Used with generateTreeView() - implements recursion
    def addToTree(self, parent_element, items):

        for item in items:

            # Add to the tree only needed item_types
            if item.item_type not \
                in [ccx_dom.item_type.GROUP, 
                    ccx_dom.item_type.KEYWORD, 
                    ccx_dom.item_type.IMPLEMENTATION]:
                continue

            # Check if there are keywords with implementations
            if self.settings.show_empty_keywords \
                    or item.name == 'Job' \
                    or item.countImplementations() \
                    or item.item_type == ccx_dom.item_type.IMPLEMENTATION:

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
                    if item.item_type == ccx_dom.item_type.IMPLEMENTATION:
                        font = QtGui.QFont()
                        font.setBold(True)
                        tree_element.setFont(font)
                tree_element.setForeground(brush)

                # Expand / collapse
                if item.expanded:
                    self.CAE.treeView.expand(tree_element.index())
                else:
                    self.CAE.treeView.collapse(tree_element.index())

                # Add icon to each keyword in tree
                icon_name = item.name.replace('*', '') + '.png'
                icon = QtGui.QIcon('./icons/' + icon_name.lower())
                tree_element.setIcon(icon)

                # Append job name
                if item.name == 'Job':
                    job_element = QtGui.QStandardItem(self.CAE.job.name)
                    tree_element.appendRow(job_element)

                # Organize recursion
                impls = item.getImplementations()
                if len(impls):
                    self.addToTree(tree_element, impls)
                else:
                    self.addToTree(tree_element, item.items)


    # Double click on treeView item: edit the keyword via dialog
    def doubleClicked(self):
        index = self.CAE.treeView.selectedIndexes()[0] # selected item index
        tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
        item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION

        # Double click on GROUP doesn't create dialog
        if item and item.item_type in \
            [ccx_dom.item_type.KEYWORD, ccx_dom.item_type.IMPLEMENTATION]:

            if item.active:

                # Create dialog window and pass item
                dialog = ccx_dialog.Dialog(self.CAE.DOM, item)

                # Get response from dialog window
                if dialog.exec() == ccx_dialog.Dialog.Accepted: # if user pressed 'OK'

                    # The generated piece of .inp code for the CalculiX input file
                    INP_code = dialog.onOk() # list of strings

                    # Create implementation object for keyword
                    if item.item_type == ccx_dom.item_type.KEYWORD:
                        impl = ccx_dom.implementation(item, INP_code) # create keyword's implementation

                        # Regenerate tree_element's children
                        tree_element.removeRows(0, tree_element.rowCount()) # remove all children
                        self.addToTree(tree_element, item.getImplementations()) # add only implementations

                    # Replace implementation object with a new one
                    elif item.item_type == ccx_dom.item_type.IMPLEMENTATION:
                        # Remove old one
                        parent = tree_element.parent() # parent treeView item
                        keyword = parent.data() # parent keyword for implementation
                        keyword.items.remove(item) # remove implementation from keyword's items

                        # Add new one
                        impl = ccx_dom.implementation(keyword, INP_code, name=item.name)
                        tree_element.setData(impl)

            else:
                logging.warning('Please, create ' + item.parent.name + ' first.')


    # Highlight node sets, element sets or surfaces
    def clicked(self):
        self.CAE.VTK.actionSelectionClear() # clear selection
        index = self.CAE.treeView.selectedIndexes()[0] # selected item index
        tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
        item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION

        if item:
            ipn = item.parent.name.upper()

            if item.item_type == ccx_dom.item_type.IMPLEMENTATION:
                lead_line = item.INP_code[0].upper()
                _set = []

                # Hightlight mesh entities
                if ipn == '*NSET' or ipn == '*NODE':
                    match = re.search('NSET\s*=\s*(\w*)', lead_line)
                    if match: # if there if NSET attribute
                        name = match.group(1) # node set name
                        _set = [self.CAE.VTK.node2point[n.num] \
                                for n in self.CAE.mesh.nsets[name].nodes]
                        self.CAE.VTK.highlight(_set, 1) # 1 = vtk.vtkSelectionNode.POINT
                elif ipn == '*ELSET' or ipn == '*ELEMENT':
                    match = re.search('ELSET\s*=\s*(\w*)', lead_line)
                    if match: # if there if ELSET attribute
                        name = match.group(1) # element set name
                        _set = [self.CAE.VTK.element2cell[e.num] \
                                for e in self.CAE.mesh.elsets[name].elements]
                        self.CAE.VTK.highlight(_set, 0) # 0 = vtk.vtkSelectionNode.CELL
                elif ipn == '*SURFACE':

                    # Surface type - optional attribute
                    stype = 'ELEMENT' # 'ELEMENT' or 'NODE'
                    match = re.search('\*SURFACE\s*,.*TYPE\s*=\s*(\w*)', lead_line)
                    if match:
                        stype = match.group(1)

                    name = re.search('NAME\s*=\s*(\w*)', lead_line).group(1) # surface name
                    if stype == 'ELEMENT':
                        _set = self.CAE.mesh.surfaces[name + stype].set
                        self.CAE.VTK.highlightSURFACE(_set)
                    elif stype=='NODE':
                        _set = [self.CAE.VTK.node2point[n.num] \
                                for n in self.CAE.mesh.surfaces[name + stype].set]
                        self.CAE.VTK.highlight(_set, 1) # 1 = vtk.vtkSelectionNode.POINT

                # Hightlight Loads & BC
                elif ipn in ['*BOUNDARY', '*CLOAD', '*CFLUX']:
                    for line in item.INP_code[1:]:
                        line = line.strip().upper()
                        n = line.replace(',', ' ').split()[0]
                        try:
                            # Single node number
                            _set.append(int(n))
                        except ValueError as err:
                            # Nodes in node set
                            _set.extend([self.CAE.VTK.node2point[n.num] \
                                for n in self.CAE.mesh.nsets[n].nodes])
                    self.CAE.VTK.highlight(set(_set), 1) # 1 = vtk.vtkSelectionNode.POINT
                elif ipn in ['*BOUNDARYF', '*MASS FLOW', '*DLOAD', '*DFLUX',
                            '*RADIATE', '*FILM', '*MODAL DAMPING']:
                    pass


    # Context menu for right click
    def rightClicked(self):
        self.myMenu = QtWidgets.QMenu('Menu', self.CAE.treeView)

        try:
            index = self.CAE.treeView.selectedIndexes()[0] # selected item index
            tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
            item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION

            # Context menu for any keyword and implementations
            if item:
                if item.item_type == ccx_dom.item_type.IMPLEMENTATION:

                    # 'Edit' action
                    action_edit_implementation = QtWidgets.QAction('Edit', self.CAE.treeView)
                    self.myMenu.addAction(action_edit_implementation)
                    action_edit_implementation.triggered.connect(self.doubleClicked)

                    # 'Delete' action
                    action_delete_implementation = QtWidgets.QAction('Delete', self.CAE.treeView)
                    self.myMenu.addAction(action_delete_implementation)
                    action_delete_implementation.triggered.connect(self.actionDeleteImplementation)

                if item.item_type == ccx_dom.item_type.KEYWORD:

                    # 'Create' action
                    action_create_implementation = QtWidgets.QAction('Create', self.CAE.treeView)
                    self.myMenu.addAction(action_create_implementation)
                    action_create_implementation.triggered.connect(self.doubleClicked)

            # Context menu for Job
            elif tree_element.text() == self.CAE.job.name:

                # Write input file
                action = QtWidgets.QAction('Write input', self.CAE.treeView)
                self.myMenu.addAction(action)
                action.triggered.connect(lambda: self.CAE.IE.exportINP())

                # Submit job
                action = QtWidgets.QAction('Submit', self.CAE.treeView)
                if os.path.isfile(self.CAE.job.name):
                    action.setDisabled(False)
                else:
                    action.setDisabled(True)
                self.myMenu.addAction(action)
                action.triggered.connect(self.CAE.job.submit)

            # Add splitter
            self.myMenu.addSeparator()

        except IndexError:
            pass

        # Context menu elements which always present
        if self.settings.show_empty_keywords:
            action_show_hide = QtWidgets.QAction('Hide empty containers', self.CAE.treeView)
        else:
            action_show_hide = QtWidgets.QAction('Show empty containers', self.CAE.treeView)
        self.myMenu.addAction(action_show_hide)
        action_show_hide.triggered.connect(self.actionShowHide)

        action_expand_collapse = QtWidgets.QAction('Collapse all', self.CAE.treeView)
        self.myMenu.addAction(action_expand_collapse)
        action_expand_collapse.triggered.connect(self.actionCollapseAll)

        action_expand_collapse = QtWidgets.QAction('Expand all', self.CAE.treeView)
        self.myMenu.addAction(action_expand_collapse)
        action_expand_collapse.triggered.connect(self.actionExpandAll)

        self.myMenu.exec_(QtGui.QCursor.pos())


    # Show/Hide empty treeView items
    def actionShowHide(self):
        self.settings.show_empty_keywords = not(self.settings.show_empty_keywords)
        self.settings.save() # save 'show_empty_keywords' value in settings
        self.generateTreeView()


    # Expand or collapse all treeView items
    def actionCollapseAll(self):
        self.CAE.treeView.collapseAll()
        self.settings.expanded = False
        self.settings.save()
    def actionExpandAll(self):
        self.CAE.treeView.expandAll()
        self.settings.expanded = True
        self.settings.save()


    # Delete keyword's implementation from DOM
    def actionDeleteImplementation(self):
        index = self.CAE.treeView.selectedIndexes()[0] # selected item index
        tree_element = self.model.itemFromIndex(index) # treeView item obtained from 'index'
        item = tree_element.data() # now it is GROUP, KEYWORD or IMPLEMENTATION

        if item.item_type == ccx_dom.item_type.IMPLEMENTATION:

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

                if not self.settings.show_empty_keywords:
                    self.hideParent(parent)
    def hideParent(self, tree_element):

        # To hide current item/brunch it should be empty 'keyword' or 'group'
        if not self.settings.show_empty_keywords \
            and not tree_element.hasChildren() \
            and tree_element.data().item_type != ccx_dom.item_type.IMPLEMENTATION:

            # Hide current item/brunch from tree via calling parent.removeRow
            parent = tree_element.parent()
            if not parent:
                parent = self.model.invisibleRootItem()
            parent.removeRow(tree_element.row())

            if parent != self.model.invisibleRootItem():
                self.hideParent(parent)


    # Change DOM item's 'expanded' variable when user interacts with treeView
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
