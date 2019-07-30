# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    CalculiX CAE
    Main window application
"""


import sys, os, argparse
from PyQt5 import QtWidgets, uic, QtCore, QtGui
import ccx_mesh, ccx_cae_tree, ccx_cae_log, ccx_vtk, ccx_dom, ccx_cae_ie


class CAE(QtWidgets.QMainWindow):


    # Create main window
    def __init__(self):

        # Create main window
        QtWidgets.QMainWindow.__init__(self)

        # Load form
        uic.loadUi('ccx_cae.ui', self)

        # Configure logging
        self.logger = ccx_cae_log.logger(self)

        # Create VTK widget
        self.VTK = ccx_vtk.VTK() # create everything for model visualization
        self.vl.addWidget(self.VTK.widget) # add vtk_widget to the form
        self.frame.setLayout(self.vl) # apply layout: it will expand vtk_widget to the frame size

        self.mesh = None # mesh from .inp-file - will be parsed in ccx_cae_ie.py
        self.IE = ccx_cae_ie.IE(self) # import/export of .inp-file
        self.DOM = ccx_dom.DOM() # empty DOM w/o implementations
        self.logger.messages(self.DOM.msg_list) # process list of log messages

        # Create/regenerate treeView items: empty model or with implementations
        self.tree = ccx_cae_tree.tree(self)

        # Default start model could be chosen with command line parameter
        parser = argparse.ArgumentParser()
        parser.add_argument("--mesh", "-mesh",
                            help="Mesh .inp file", type=str,
                            # default='')
                            # default='ccx_mesh.inp')
                            default='./examples/achtel2.inp')
        args = parser.parse_args()
        if len(args.mesh): # import default ugrid
            msgs = self.IE.importINP(args.mesh)
            self.logger.messages(msgs)

        # Actions
        self.actions()


    # Actions
    def actions(self):
        self.treeView.keyPressEvent = self.keyPressEvent

        # VTK actions
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
