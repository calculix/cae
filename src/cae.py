#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, July 2020
Distributed under GNU General Public License v3.0

CalculiX CAE - main module.
Creates main objects and starts the application.

How to run:
Ctrl+F5 from VSCode
or:
python3 ./src/cae.py
python3 ./src/cae.py -inp yourmodel.inp """

# Standard modules
import os
import sys
import time
import argparse
import logging

# External modules
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
import tree
import clean
import path
import importer

# Pyinstaller bug in Windows:
# append 'app_home_dir' and 'src' directories to PATH
p = path.Path() # calculate absolute paths
p.append_to_PATH([p.app_home_dir, p.src])

if __name__ == '__main__':
    start_time = time.perf_counter()
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
    s = settings.Settings(p)

    # Configure global logging level
    logging.getLogger().setLevel(s.logging_level)

    # Default start model (INP file)
    # could be chosen with command line parameter
    parser = argparse.ArgumentParser()
    parser.add_argument('-inp',
        type=str, help='your .inp file',
        default=s.start_model)
    args = parser.parse_args()

    # Show CAE window and get window ID
    # A new logger's handler is created here
    if os.name=='nt':
        w = gui.window.Windows_window(p, s)
    else:
        w = gui.window.Linux_window(p, s)
    w.show()
    w.wid1 = w.get_wid('CalculiX CAE')
    if s.align_windows:
        w.align()

    # Main block
    m = model.Model() # generate FEM model
    t = tree.Tree(p, s, w, m) # create treeView items based on KOM
    j = model.job.Job(p, s, w) # create job object
    actions.actions(p, s, w, m, t, j) # window actions

    # Import default model
    if len(args.inp):
        start_model = os.path.join(p.app_home_dir, args.inp)
        importer.import_file(p, s, w, m, t, j, start_model)

    # Or start empty
    else:
        logging.warning('No default start model specified.')
        m.KOM = model.kom.KOM(p, s)
        t.generateTreeView(m)

    logging.info('Started in {:.1f} seconds.\n'
        .format(time.perf_counter() - start_time))

    # Execute application
    app.exec()

    # Kill CGX after CAE exit
    gui.cgx.kill(w)

    # Recursively clean cached files in all subfolders
    clean.cache(p.src)
