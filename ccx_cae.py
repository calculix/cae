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
import ccx_cae_tree, ccx_vtk, ccx_dom, ccx_cae_ie, ccx_settings, ccx_job


# Logging handler
class myLoggingHandler(logging.Handler):


    # Initialization
    def __init__(self, CAE):
        super().__init__() # create handler
        self.textEdit = CAE.textEdit
        self.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))


    # Sends log messages to CAE's textEdit widget
    def emit(self, LogRecord):
        msg_text = self.format(LogRecord)

        # Message color depending on logging level
        color = {
                'DEBUG':'Gray',
                'INFO':'Black',
                'WARNING':'Blue',
                'ERROR':'Red',
            }[LogRecord.levelname]

        self.textEdit.append('<p style=\'color:{0}; margin:0px;\'>{1}</p>'.format(color, msg_text))
        self.textEdit.moveCursor(QtGui.QTextCursor.End) # scroll text to the end


# Main window
class CAE(QtWidgets.QMainWindow):


    # Create main window
    def __init__(self, settings):

        # Create main window
        QtWidgets.QMainWindow.__init__(self)

        # Load form
        uic.loadUi('ccx_cae.ui', self)

        # # Read application's global settings
        # settings = ccx_settings.Settings()

        # Configure logging
        logging.getLogger().addHandler(myLoggingHandler(self))
        logging.getLogger().setLevel(settings.logging_level)

        # Create VTK widget
        self.VTK = ccx_vtk.VTK() # create everything for model visualization
        self.vl.addWidget(self.VTK.widget) # add vtk_widget to the form

        self.mesh = None # mesh from .inp-file - will be parsed in ccx_cae_ie.py
        self.IE = ccx_cae_ie.IE(self) # import/export of .inp-file
        self.DOM = ccx_dom.DOM() # empty DOM w/o implementations

        # Default start model could be chosen with command line parameter
        parser = argparse.ArgumentParser()
        parser.add_argument('-inp', type=str, help='your .inp file',
                            default=settings.default_start_model)
        args = parser.parse_args()

        # Create job object
        self.job = ccx_job.Job(settings, args.inp)

        # Create/regenerate treeView items: empty model or with implementations
        self.tree = ccx_cae_tree.tree(self)

        # Import default ugrid
        self.IE.importFile(args.inp)

        # Actions
        self.actions()


    # Actions
    def actions(self):
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

    # Create main window
    window = CAE(settings)
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
