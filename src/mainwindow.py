# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, December 2019
    Distributed under GNU General Public License v3.0

    Main window class.
"""


from PyQt5 import QtWidgets, uic, QtCore, QtGui
import logging
from VTK import VTK
from log import myLoggingHandler


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


    # TODO Delete keyword's implementation in the treeView by pressing 'Delete' button
    # def keyPressEvent(self, e):
    #     if e.key() == QtCore.Qt.Key_Delete:
    #         self.tree.actionDeleteImplementation()


    # Open links from the Help menu
    def help(self, link):
        url = QtCore.QUrl(link)
        logging.info('Going to ' + link)
        if not QtGui.QDesktopServices.openUrl(url):
            logging.warning('Can\'t open url: ' + link)
