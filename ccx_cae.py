# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    CalculiX CAE
    Main window application
    INP parser and writer
"""


import sys, os, argparse, vtk
from PyQt5 import QtWidgets, uic, QtCore, QtGui, Qt
import ccx_mesh, ccx_tree, ccx_log, ccx_vtk, ccx_dom, ccx_inp


class CAE(QtWidgets.QMainWindow):


    # Create main window
    def __init__(self):

        # Create main window
        QtWidgets.QMainWindow.__init__(self)

        # Load form
        uic.loadUi('ccx_cae.ui', self)

        # Configure logging
        self.logger = ccx_log.logger(self)

        # Create VTK widget
        self.VTK = ccx_vtk.VTK(self) # create everything for model visualization
        self.vl.addWidget(self.VTK.widget) # add vtk_widget to the form
        self.frame.setLayout(self.vl) # apply layout: it will expand vtk_widget to the frame size

        # Generate CalculiX DOM based on keywords hierarchy from ccx_dom.txt
        self.DOM = ccx_dom.DOM(self) # empty DOM w/o implementations

        # Generate empty CalculiX model and treeView items
        self.tree = ccx_tree.tree(self)

        # Default start model could be chosen with command line parameter
        parser = argparse.ArgumentParser()
        parser.add_argument("--mesh", "-mesh",
                            help="Mesh .inp file",
                            type=str, default='./examples/acou1.inp')
        args = parser.parse_args()

        # Import default ugrid
        self.inp = ccx_inp.inp(self)
        self.inp.importINP(args.mesh)

        # Actions
        self.actionImportINP.triggered.connect(self.inp.importINP)
        self.actionExportINP.triggered.connect(self.inp.exportINP)
        self.treeView.keyPressEvent = self.keyPressEvent


    # Delete keyword's implementation in the treeView by pressing 'Delete' button
    def keyPressEvent(self, e):
        if e.key() == QtCore.Qt.Key_Delete:
            self.tree.actionDeleteImplementation()


if __name__ == '__main__':

    # Clean cached files before start
    os.system('py3clean .')

    app = QtWidgets.QApplication(sys.argv)
    window = CAE()
    window.show() # window.showMaximized() or window.show()
    sys.exit(app.exec_())
