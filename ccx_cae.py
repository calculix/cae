# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    CalculiX CAE - main window.
    How to run:
        python3 ccx_cae.py -inp ccx_mesh.inp
"""


import sys, os, argparse, logging, shutil, subprocess
os.environ['PATH'] += os.path.dirname(sys.executable) # Pyinstaller bug in Windows
from PyQt5 import QtWidgets, uic, QtCore, QtGui
import ccx_cae_tree, ccx_vtk, ccx_dom, ccx_cae_ie, ccx_settings, ccx_job, ccx_log


# Main window
class CAE(QtWidgets.QMainWindow):


    # Create main window
    def __init__(self, settings, default_start_model):
        QtWidgets.QMainWindow.__init__(self) # create main window
        uic.loadUi('ccx_cae.ui', self) # load form

        # Configure logs to be shown in window
        logging.getLogger().addHandler(ccx_log.myLoggingHandler(self))
        logging.getLogger().setLevel(settings.logging_level)

        # Create VTK widget
        self.VTK = ccx_vtk.VTK() # create everything for model visualization
        self.vl.addWidget(self.VTK.widget) # add vtk_widget to the form

        self.mesh = None # mesh from .inp-file - will be parsed in ccx_cae_ie.py
        self.IE = ccx_cae_ie.IE(self, settings) # import/export of .inp-file
        self.DOM = ccx_dom.DOM() # empty DOM w/o implementations
        self.job = ccx_job.Job(settings, default_start_model) # create job object
        self.tree = ccx_cae_tree.tree(self) # create treeView items based on DOM
        self.IE.importFile(default_start_model) # import default ugrid

        # Actions
        if True:
            self.treeView.keyPressEvent = self.keyPressEvent

            # VTK actions
            self.actionSelectionNodes.triggered.connect(self.VTK.actionSelectionNodes)
            self.actionSelectionElements.triggered.connect(self.VTK.actionSelectionElements)
            self.actionSelectionClear.triggered.connect(self.VTK.actionSelectionClear)
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


# Here application starts
if __name__ == '__main__':

    # Create application
    app = QtWidgets.QApplication(sys.argv)

    # Read application's global settings
    settings = ccx_settings.Settings()

    # Default start model could be chosen with command line parameter
    parser = argparse.ArgumentParser()
    parser.add_argument('-inp', type=str, help='your .inp file',
                        default=settings.default_start_model)
    args = parser.parse_args()

    # Create and show main window
    window = CAE(settings, args.inp)
    if settings.showMaximized:
        window.showMaximized()
    else:
        window.show()

    # Execute application
    a = app.exec_()

    # Clean cached files
    if os.path.isdir('__pycache__'):
        shutil.rmtree('__pycache__') # works in Linux as in Windows

    # Exit application
    sys.exit(a)
