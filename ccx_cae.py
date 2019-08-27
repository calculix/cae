# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    CalculiX CAE - main window.
    How to run:
        python3 ccx_cae.py -inp ccx_mesh.inp
"""

import sys, os

# Update enviroment variable PATH: pyinstaller bug in Windows
home_dir = os.path.dirname(sys.argv[0]) # app. home directory
if home_dir not in os.environ['PATH']:
    if not os.environ['PATH'].endswith(os.pathsep):
        os.environ['PATH'] += os.pathsep
    os.environ['PATH'] += home_dir

import argparse, logging, shutil, subprocess, clean
from PyQt5 import QtWidgets, uic, QtCore, QtGui
import ccx_cae_tree, ccx_vtk, ccx_kom, ccx_cae_ie, ccx_settings, ccx_job, ccx_log


# Main window
class CAE(QtWidgets.QMainWindow):


    # Create main window
    def __init__(self, settings, path_start_model):
        QtWidgets.QMainWindow.__init__(self) # create main window
        ui = os.path.join(os.path.dirname(sys.argv[0]), 'ccx_cae.xml') # full path
        uic.loadUi(ui, self) # load form

        # Configure logs to be shown in window
        logging.getLogger().addHandler(ccx_log.myLoggingHandler(self))
        logging.getLogger().setLevel(settings.logging_level)

        # Abs. path to the path_start_model
        if len(path_start_model):
            path_start_model = os.path.abspath(path_start_model)
            if not os.path.isfile(path_start_model):
                path_start_model = os.path.join(os.path.dirname(sys.argv[0]),
                    os.path.basename(path_start_model))

        # Create VTK widget
        if settings.show_vtk:
            self.VTK = ccx_vtk.VTK() # create everything for model visualization
            self.h_splitter.addWidget(self.VTK.widget)
            self.setMinimumSize(1280, 600)
            self.resize(1280, 720)
        else:
            self.toolBar.setParent(None) # hide toolbar

        self.mesh = None # mesh from .inp-file - will be parsed in ccx_cae_ie.py
        self.IE = ccx_cae_ie.IE(self, settings) # import/export of .inp-file
        self.KOM = ccx_kom.KOM() # empty KOM w/o implementations
        self.job = ccx_job.Job(settings, path_start_model) # create job object
        self.tree = ccx_cae_tree.tree(self) # create treeView items based on KOM
        if len(path_start_model):
            self.IE.importFile(path_start_model) # import default start model

        # Actions
        if True:
            self.treeView.keyPressEvent = self.keyPressEvent
            self.action_file_settings.triggered.connect(settings.open)

            # Job actions
            self.action_job_write_input.triggered.connect(self.IE.writeInput)
            self.action_job_edit_inp.triggered.connect(self.job.editINP)
            self.action_job_open_subroutine.triggered.connect(self.job.openSubroutine)
            self.action_job_rebuild_ccx.triggered.connect(self.job.rebuildCCX)
            self.action_job_submit.triggered.connect(self.job.submit)
            self.action_job_view_log.triggered.connect(self.job.viewLog)
            self.action_job_open_cgx.triggered.connect(self.job.openCGX)
            self.action_job_export_vtu.triggered.connect(self.job.exportVTU)
            self.action_job_open_paraview.triggered.connect(self.job.openParaview)

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
    def keyPressEvent(self, e):
        if e.key() == QtCore.Qt.Key_Delete:
            self.tree.actionDeleteImplementation()


if __name__ == '__main__':

    # Create application
    app = QtWidgets.QApplication(sys.argv)

    # Read application's global settings
    settings = ccx_settings.Settings()

    # Default start model could be chosen with command line parameter
    parser = argparse.ArgumentParser()
    parser.add_argument('-inp', type=str, help='your .inp file',
                        default=settings.path_start_model)
    args = parser.parse_args()

    # Create and show main window
    window = CAE(settings, args.inp)
    if settings.show_maximized:
        window.showMaximized()
    else:
        window.show()

    # Execute application
    a = app.exec_()

    # Clean cached files
    clean.cache()

    # Exit application
    sys.exit(a)
