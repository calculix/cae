# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, October 2019
    Distributed under GNU General Public License v3.0

    Main window class.
"""


from PyQt5 import QtWidgets, uic, QtCore, QtGui
import logging, os
from tree import tree
from VTK import VTK
from log import myLoggingHandler
from ie import importFile


# Main window
class MainWindow(QtWidgets.QMainWindow):


    # Create main window
    def __init__(self, p, settings):
        QtWidgets.QMainWindow.__init__(self) # create main window
        uic.loadUi(p.cae_xml, self) # load form

        # Configure logs to be shown in window
        logging.getLogger().addHandler(myLoggingHandler(self.textEdit))
        logging.getLogger().setLevel(settings.logging_level)

        # When logger is ready - check if settings read correctly
        if hasattr(settings, 'error_path'):
            logging.error('Error path in settings file: ' + settings.error_path + '. Loading default values.')

        # Create VTK widget
        if settings.show_vtk:
            self.VTK = VTK() # create everything for model visualization
            self.h_splitter.addWidget(self.VTK.widget)
            self.setMinimumSize(1280, 600)
            self.resize(1280, 720)
        else:
            self.toolBar.setParent(None) # hide toolbar

        # MainWindow actions
        if True:
            # self.treeView.keyPressEvent = self.keyPressEvent

            # File actions
            self.action_file_settings.triggered.connect(settings.open)
            self.action_file_exit.triggered.connect(QtWidgets.qApp.quit)

            # Help actions
            self.action_help_readme.triggered.connect(lambda:
                    self.help('https://github.com/imirzov/ccx_cae#calculix-cae'))
            self.action_help_yahoo.triggered.connect(lambda:
                    self.help('https://groups.yahoo.com/neo/groups/CALCULIX/conversations/topics/15616'))
            self.action_help_issues.triggered.connect(lambda:
                    self.help('https://github.com/imirzov/ccx_cae/issues'))

            # VTK actions
            if settings.show_vtk:
                # self.actionSelectionNodes.triggered.connect(self.VTK.actionSelectionNodes)
                # self.actionSelectionElements.triggered.connect(self.VTK.actionSelectionElements)
                # self.actionSelectionClear.triggered.connect(self.VTK.actionSelectionClear)
                self.actionViewParallel.triggered.connect(self.VTK.actionViewParallel)
                self.actionViewPerspective.triggered.connect(self.VTK.actionViewPerspective)
                self.actionViewFront.triggered.connect(self.VTK.actionViewFront)
                self.actionViewBack.triggered.connect(self.VTK.actionViewBack)
                self.actionViewTop.triggered.connect(self.VTK.actionViewTop)
                self.actionViewBottom.triggered.connect(self.VTK.actionViewBottom)
                self.actionViewLeft.triggered.connect(self.VTK.actionViewLeft)
                self.actionViewRight.triggered.connect(self.VTK.actionViewRight)
                self.actionViewIso.triggered.connect(self.VTK.actionViewIso)
                self.actionViewFit.triggered.connect(self.VTK.actionViewFit)
                self.actionViewWireframe.triggered.connect(self.VTK.actionViewWireframe)
                self.actionViewSurface.triggered.connect(self.VTK.actionViewSurface)
                self.actionViewSurfaceWithEdges.triggered.connect(self.VTK.actionViewSurfaceWithEdges)


    # Delete keyword's implementation in the treeView by pressing 'Delete' button
    # def keyPressEvent(self, e):
    #     if e.key() == QtCore.Qt.Key_Delete:
    #         self.tree.actionDeleteImplementation()


    # Open links from the Help menu
    def help(self, link):
        url = QtCore.QUrl(link)
        logging.info('Going to ' + link)
        if not QtGui.QDesktopServices.openUrl(url):
            logging.warning('Can\'t open url: ' + link)
