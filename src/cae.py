# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, December 2019
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
from layer.ie import importFile
from layer.Tree import Tree
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
    t = Tree(p, settings, mw, m)
    # Abs. path to INP file
    if len(args.inp):
        path_to_inp = os.path.join(p.app_home_dir, args.inp)
        importFile(settings, mw, m, t, path_to_inp) # import default start model

    # MainWindow actions
    actions(settings, mw, m, t)


    # Execute application
    a = app.exec()

    # Recursively clean cached files in all subfolders
    cleanCache(p.src)

    # Exit application
    sys.exit(a)
