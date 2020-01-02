# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, January 2020
    Distributed under GNU General Public License v3.0

    CalculiX CAE - main module.
    Creates main objects and starts the application.

    How to run:
        python3 cae.py
        python3 cae.py -inp model.inp
"""


# Pyinstaller bug in Windows: append 'app_home_dir' and 'src' directories to PATH
from Path import Path
p = Path() # calculate absolute paths
p.append_to_PATH([p.app_home_dir, p.src])

# Main imports
import os, sys, argparse
from PyQt5 import QtWidgets
from Settings import Settings
from actions import actions
from gui.MainWindow import MainWindow
from model.Model import Model
from model.Job import Job
from ie import importFile
from Tree import Tree
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
    s = Settings()

    # Default start model could be chosen with command line parameter
    parser = argparse.ArgumentParser()
    parser.add_argument('-inp', type=str, help='your .inp file',
                        default=s.path_start_model)
    args = parser.parse_args()



    # Create and show main window
    w = MainWindow(p, s)
    if s.show_maximized:
        w.showMaximized()
    else:
        w.show()

    m = Model() # generate FEM model
    j = Job(p, s, args.inp) # create job object
    t = Tree(p, s, w, m) # create treeView items based on KOM

    # Abs. path to INP file
    if len(args.inp):
        path_to_inp = os.path.join(p.app_home_dir, args.inp)
        importFile(s, w, m, t, j, path_to_inp) # import default start model


    # MainWindow actions
    actions(s, w, m, t, j)


    # Execute application
    a = app.exec()

    # Recursively clean cached files in all subfolders
    cleanCache(p.src)

    # Exit application
    sys.exit(a)
