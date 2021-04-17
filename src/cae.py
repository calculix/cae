#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Main module.
Creates main objects and starts the application.

How to run:
Ctrl+F5 from VSCode
or:
python3 ./src/cae.py
python3 ./src/cae.py -inp yourmodel.inp """

# Run checks before start
import checks
import logging
ch = checks.Check()
ch.check_os()
ch.check_python()
ch.check_requirements()
logging.getLogger().handlers = [] # stop logging into checks.log

# Standard modules
import os
import sys
import time
import argparse

# External modules
from PyQt5 import QtWidgets

# My modules
import path
import clean
import settings
import gui
import model
import tree
import importer
import actions

# Pyinstaller bug in Windows:
# append 'app_home_dir' and 'src' directories to PATH
p = path.Path() # calculate absolute paths
p.append_to_PATH([p.app_home_dir, p.src])

start_time = time.perf_counter()
clean.screen()

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
if os.name == 'nt':
    w = gui.window.MainWindowWindows(p, s)
else:
    w = gui.window.MainWindowLinux(p, s)
w.show()

# Main block
m = model.Model() # generate FEM model
t = tree.Tree(p, s, w, m) # create treeView items based on KOM
j = model.job.Job(p, s, w, m) # create job object
i = importer.Importer(p, s, w, m, t, j)
actions.actions(p, s, w, m, t, j, i) # window actions

# Import default model
if len(args.inp):
    start_model = os.path.join(p.app_home_dir, args.inp)
    i.import_file(start_model)

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
