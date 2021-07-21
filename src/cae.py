#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Main module. Creates main objects and starts the application.
NOTE QtWebEngineWidgets has to be imported before the app start.

How to run:
Ctrl+F5 from VSCode
or:
python3 ./src/cae.py
python3 ./src/cae.py -inp yourmodel.inp
"""

import clean
clean.screen()

# First time logging configure
import logging
logging.basicConfig(level=logging.NOTSET)

# Run some important checks before start
import checks
checks.run_startup_checks()
clean.screen()

# Standard modules
import os
import sys
import time
import argparse

start_time = time.perf_counter()

# External modules
from PyQt5 import QtWidgets, QtWebEngineWidgets

# My modules
import path
import settings
import gui.job

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

# Default start model (INP file)
# could be chosen with command line parameter
parser = argparse.ArgumentParser()
parser.add_argument('-inp', type=str,
    help='your .inp file',
    default=settings.s.start_model)
args = parser.parse_args()

# Show main window with text logging handler
# TODO Avoid terminal logs
from gui.window import factory
factory.run_master(path.p.main_xml)

# Main block
# TODO Do not create inctances here - do it in appropriate modules
import importer
i = importer.Importer() # prepare to import model
import actions
actions.actions(i) # window actions

# Import default model
if len(args.inp):
    start_model = os.path.join(path.p.app_home_dir, args.inp)
    i.import_file(start_model)

# Or start empty
else:
    logging.warning('No default start model specified.')
    from tree import t
    t.generateTreeView()

logging.info('Started in {:.1f} seconds.\n'
    .format(time.perf_counter() - start_time))

# Execute application
app.exec()

# Kill CGX after CAE exit
factory.kill_slave()

# Recursively clean cached files in all subfolders
clean.cache(path.p.src)
