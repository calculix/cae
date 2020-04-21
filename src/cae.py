#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, January 2020
Distributed under GNU General Public License v3.0

CalculiX CAE - main module.
Creates main objects and starts the application.

How to run:
python3 cae.py
python3 cae.py -inp model.inp """


# System imports
import os
import sys
import time
import argparse
import logging

# Requirements
try:
    from PyQt5 import QtWidgets
except:
    msg = 'Please, install PyQt5 with command:\n'\
        + 'pip3 install PyQt5'
    sys.exit(msg)

# My modules
import settings
import actions
import gui
import model
import ie
import tree
import clean
import path


# Pyinstaller bug in Windows:
# append 'app_home_dir' and 'src' directories to PATH
p = path.Path() # calculate absolute paths
p.append_to_PATH([p.app_home_dir, p.src])


if __name__ == '__main__':
    clean.screen()

    # Exit if OS is not Linux or Windows
    if os.name not in ['nt', 'posix']:
        msg = 'SORRY, {} OS is not supported.'
        sys.exit(msg.format(os.name))

    # # Draw apps architecture
    # from pycallgraph import PyCallGraph
    # from pycallgraph import Config
    # from pycallgraph import GlobbingFilter
    # from pycallgraph.output import GraphvizOutput
    # modules = [m[:-3]+'*' for m in os.listdir(p.src) if m.endswith('.py')] + ['Window*']
    # config = Config()
    # config.trace_filter = GlobbingFilter(
    #     include=modules, exclude=['logging*', '*FileFinder'])
    # graphviz = GraphvizOutput(output_file='architecture.png')
    # with PyCallGraph(output=graphviz, config=config):

    # Create application
    app = QtWidgets.QApplication(sys.argv)

    # Read application's global settings
    s = settings.Settings()

    # Configure global logging level
    logging.getLogger().setLevel(s.logging_level)

    # Default start model (INP file)
    # could be chosen with command line parameter
    parser = argparse.ArgumentParser()
    parser.add_argument('-inp',
        type=str, help='your .inp file',
        default=s.start_model)
    args = parser.parse_args()
    start_model = ''
    if len(args.inp):
        start_model = os.path.join(p.app_home_dir, args.inp)

    # Show CAE and get window ID
    if os.name=='nt':
        w = gui.window.Windows_window(p, s)
    else:
        w = gui.window.Linux_window(p, s)
    w.show()
    w.initialize()
    w.wid1 = w.get_wid(w.windowTitle())
    if w.wid1 is None:
        msg = 'ERROR! Can\'t get {} window ID.'
        sys.exit(msg.format(w.windowTitle()))

    # Main block
    m = model.model.Model() # generate FEM model
    j = model.job.Job(p, s, args.inp) # create job object
    t = tree.Tree(p, s, w, m) # create treeView items based on KOM
    ie.importFile(s, w, m, t, j, start_model) # import default model
    actions.actions(s, w, m, t, j) # window actions

    # Execute application
    app.exec()

    # Kill CGX after CAE exit
    gui.cgx.kill()

    # Recursively clean cached files in all subfolders
    clean.cache(p.src)
