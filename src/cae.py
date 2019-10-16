# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, October 2019
    Distributed under GNU General Public License v3.0

    CalculiX CAE - main function.
    How to run:
        python3 cae.py -inp model.inp
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
from ie import importFile
from clean import cleanCache


if __name__ == '__main__':
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
        m = Model(p, settings, args.inp)

        # Create treeView items based on KOM
        t = tree(p, settings, mw.treeView, mw.VTK, m.KOM)



        # Execute application
        a = app.exec()

        # Clean cached files
        cleanCache(p.src)

        # Exit application
        sys.exit(a)
