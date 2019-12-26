# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, October 2019
    Distributed under GNU General Public License v3.0

    CalculiX CAE - main function.
    How to run:
        python3 cae.py -inp model.inp


    Что почитать:
    https://habr.com/ru/post/276593/
    https://fktpm.ru/file/84-sovershennyj-kod.pdf
    Google:
        model view controller python
        шаблоны проектирования python
"""


# Pyinstaller bug in Windows: append 'app_home_dir' and 'src' directories to PATH
from path import Path
p = Path() # calculate absolute paths
p.append_to_PATH([p.app_home_dir, p.src])

# Main imports
import os, sys, argparse, shutil, subprocess
from PyQt5 import QtWidgets
from settings import Settings
from mainwindow import MainWindow
from model import Model
from tree import tree
from ie import importFile, writeInput
from clean import cleanCache


if __name__ == '__main__':
    # TODO: have a look at https://github.com/davidfraser/pyan
    # from pycallgraph import PyCallGraph
    # from pycallgraph import Config
    # from pycallgraph import GlobbingFilter
    # from pycallgraph.output import GraphvizOutput
    # modules = [m[:-3]+'*' for m in os.listdir(p.src) if m.endswith('.py')] + ['MainWindow*']
    # config = Config()
    # config.trace_filter = GlobbingFilter(
    #     include=modules, exclude=['logging*', '*FileFinder'])
    # graphviz = GraphvizOutput(output_file='architecture.png')
    # with PyCallGraph(output=graphviz, config=config):

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
        mw = MainWindow(p, settings)
        if settings.show_maximized:
            mw.showMaximized()
        else:
            mw.show()

        # Generate FEM model
        m = Model(settings, args.inp)

        # Create treeView items based on KOM
        t = tree(p, settings, mw.treeView, mw.VTK, m.KOM)

        # Abs. path to INP file
        if len(args.inp):
            path_to_inp = os.path.join(p.app_home_dir, args.inp)
            importFile(p, settings, m, mw, t, path_to_inp) # import default start model

        # MainWindow actions
        if True:

            # File actions
            mw.action_file_import.triggered.connect(
                lambda: importFile(p, settings, m, mw, t))

            # Job actions
            mw.action_job_write_input.triggered.connect(writeInput)
            mw.action_job_edit_inp.triggered.connect(m.job.editINP)
            mw.action_job_open_subroutine.triggered.connect(m.job.openSubroutine)
            mw.action_job_rebuild_ccx.triggered.connect(m.job.rebuildCCX)
            mw.action_job_submit.triggered.connect(m.job.submit)
            mw.action_job_view_log.triggered.connect(m.job.viewLog)
            mw.action_job_open_cgx.triggered.connect(m.job.openCGX)
            mw.action_job_export_vtu.triggered.connect(m.job.exportVTU)
            mw.action_job_open_paraview.triggered.connect(m.job.openParaView)



        # Execute application
        a = app.exec()

        # Clean cached files
        cleanCache(p.src)

        # Exit application
        sys.exit(a)
