# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, September 2019
    Distributed under GNU General Public License v3.0

    CalculiX CAE - main window.
    How to run:
        python3 cae.py -inp model.inp
"""


# Pyinstaller bug in Windows: append 'app_home_dir' and 'src' directories to PATH
from path import Path
p = Path() # calculate absolute pathes
p.append_to_PATH([p.app_home_dir, p.src])

# Main imports
import os, sys, argparse, logging, shutil, subprocess
from PyQt5 import QtWidgets, uic, QtCore, QtGui
from cae_tree import tree
from VTK import VTK
from kom import KOM
from cae_ie import IE
from settings import Settings 
from job import Job
from log import myLoggingHandler
from clean import cleanCache


# Main window
class CAE(QtWidgets.QMainWindow):


    # Create main window
    def __init__(self, settings, path_start_model):
        QtWidgets.QMainWindow.__init__(self) # create main window
        ui = os.path.join(p.config, 'cae.xml') # full path
        uic.loadUi(ui, self) # load form

        # Configure logs to be shown in window
        logging.getLogger().addHandler(myLoggingHandler(self.textEdit))
        logging.getLogger().setLevel(settings.logging_level)

        # When logger is ready - check if settings read correctly
        if hasattr(settings, 'error'):
            logging.error('Error reading ENV settings file. Default values used.')

        # Abs. path to the path_start_model
        if len(path_start_model):
            path_start_model = os.path.join(p.app_home_dir, path_start_model)

        # Create VTK widget
        if settings.show_vtk:
            self.VTK = VTK() # create everything for model visualization
            self.h_splitter.addWidget(self.VTK.widget)
            self.setMinimumSize(1280, 600)
            self.resize(1280, 720)
        else:
            self.toolBar.setParent(None) # hide toolbar

        self.mesh = None # mesh from .inp-file - will be parsed in cae_ie.py
        self.IE = IE(self, settings) # import/export of .inp-file
        self.KOM = KOM() # empty KOM w/o implementations
        # TODO try to reorder to omit double job calling/renaming
        self.job = Job(settings, path_start_model) # create job object
        self.tree = tree(self, settings) # create treeView items based on KOM
        if len(path_start_model):
            self.IE.importFile(path_start_model) # import default start model

        # Actions
        if True:
            self.treeView.keyPressEvent = self.keyPressEvent

            # File actions
            self.action_file_import.triggered.connect(self.IE.importFile)
            self.action_file_settings.triggered.connect(settings.open)
            self.action_file_exit.triggered.connect(QtWidgets.qApp.quit)

            # Job actions
            self.action_job_write_input.triggered.connect(self.IE.writeInput)
            self.action_job_edit_inp.triggered.connect(self.job.editINP)
            self.action_job_open_subroutine.triggered.connect(self.job.openSubroutine)
            self.action_job_rebuild_ccx.triggered.connect(self.job.rebuildCCX)
            self.action_job_submit.triggered.connect(self.job.submit)
            self.action_job_view_log.triggered.connect(self.job.viewLog)
            self.action_job_open_cgx.triggered.connect(self.job.openCGX)
            self.action_job_export_vtu.triggered.connect(self.job.exportVTU)
            self.action_job_open_paraview.triggered.connect(self.job.openParaView)

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
    def keyPressEvent(self, e):
        if e.key() == QtCore.Qt.Key_Delete:
            self.tree.actionDeleteImplementation()


    # Open README.md on the GitHub
    def help(self, link):
        url = QtCore.QUrl(link)
        logging.info('Going to ' + link)
        if not QtGui.QDesktopServices.openUrl(url):
            logging.warning('Can\'t open url: ' + link)


if __name__ == '__main__':

    # Create application
    app = QtWidgets.QApplication(sys.argv)

    # Read application's global settings
    settings = Settings()

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
    cleanCache(p.src)

    # Exit application
    sys.exit(a)
