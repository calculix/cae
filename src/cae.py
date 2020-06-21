#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, June 2020
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
import file_tools
import clean
import path

# Pyinstaller bug in Windows:
# append 'app_home_dir' and 'src' directories to PATH
p = path.Path() # calculate absolute paths
p.append_to_PATH([p.app_home_dir, p.src])

""" INP importer:
- Enrich KOM with implementations from parsed file.
- Generate new tree with keyword implementations.
- Parse mesh.
- Open model in CGX.

We can get here via menu File -> Import
or directly after application start.

s - Settings
w - Window
m - Model
t - Tree
j - Job """
def import_file(s, w, m, t, j, file_name=None):
    if file_name is None or len(file_name)==0:
        file_name = QtWidgets.QFileDialog.getOpenFileName(w, \
            'Import INP/UNV file', j.dir, \
            'INP (*.inp);;UNV (*.unv)')[0]

    if file_name is not None and len(file_name):

        # Rename job before tree regeneration
        j.initialize(file_name[:-4] + '.inp')

        # Convert UNV to INP
        if file_name.lower().endswith('.unv'):
            j.convertUNV()
            if not os.path.isfile(j.inp):
                logging.error('Can not convert\n' + j.unv)
                return

        # Open model in CGX
        if s.run_cgx_on_start:
            j.cgx_inp(s, w)

        # Show model name in window's title
        w.setWindowTitle('CalculiX CAE - ' + j.name)

        # Generate new KOM without implementations
        m.KOM = model.kom.KOM()

        # Parse INP and enrich KOM with parsed objects
        logging.info('Loading model\n{}'.format(j.inp))
        lines = file_tools.read_lines(j.inp)
        def importer(INP_doc, KOM):
            keyword_chain = []
            impl_counter = {}
            for i in range(len(INP_doc)):
                line = INP_doc[i]

                # Parse keyword
                if line.startswith('*'):
                    logging.debug(line)

                    # Distinguish 'NODE' and 'NODE PRINT'
                    if ',' in line:
                        keyword_name = line.split(',')[0]
                    else:
                        keyword_name = line

                    # Find KOM keyword path corresponding to keyword_chain
                    keyword_chain.append(keyword_name)
                    path = KOM.getPath(keyword_chain)
                    if path:
                        logging.debug('path found: ' + str([item.name for item in path]))

                        # Read INP_code for the current keyword
                        INP_code = [line] # line is stripped in mesh.py
                        while i+1 < len(INP_doc) and \
                            not INP_doc[i+1].startswith('*'): # here will be no comments - they are removed in mesh.py
                            INP_code.append(INP_doc[i+1])
                            i += 1

                        # Create keyword implementations
                        impl = None
                        path_as_string = '' # string representation of 'path' accounting for implementations
                        for j in range(len(path)):
                            # Choose where to create implementation
                            if impl:
                                # Implementation will be created inside another implementation
                                item = impl.getItemByName(path[j].name)
                            else:
                                # Implementation will be created inside keyword or group
                                item = path[j]

                            path_as_string += '/' + item.name
                            if j == len(path) - 1: # last item is always keyword
                                # Create implementation (for example, MATERIAL-1)
                                impl = model.kom.implementation(item, INP_code)
                                # logging.debug('1')
                            elif item.item_type == model.kom.item_type.KEYWORD:
                                # If for this keyword implementation was created previously
                                counter = impl_counter[path_as_string] - 1
                                impl = item.items[counter] # first implementation, for example, STEP-1
                                path_as_string += '/' + impl.name
                                # logging.debug('2')
                            else:
                                impl = item
                                # logging.debug('3')

                        # Count implementation
                        if path_as_string in impl_counter:
                            # If current keyword already has implementations
                            impl_counter[path_as_string] += 1
                        else:
                            # If first implementation was created for current keyword
                            impl_counter[path_as_string] = 1

                    else:
                        logging.warning('Wrong keyword {}.'.format(keyword_name))
        importer(lines, m.KOM) # pass whole INP-file to the parser

        # Add parsed implementations to the tree
        t.generateTreeView(m)

        # Parse mesh
        m.Mesh = model.parsers.mesh.Mesh(INP_file=j.inp)

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
    if s.align_windows:
        w.align()

    # How the application is executed:
    # from source code or as binary
    if sys.argv[0].endswith('.py'):
        logging.info('Running from sources.')
    else:
        logging.info('Running from binaries.')

    # Main block
    m = model.Model() # generate FEM model
    t = tree.Tree(p, s, w, m) # create treeView items based on KOM
    j = model.job.Job(p) # create job object
    actions.actions(p, s, w, m, t, j) # window actions
    import_file(s, w, m, t, j, start_model) # import default model

    if s.run_cgx_on_start:
        w.post(s.model_view)
    logging.info('Started in {:.1f} seconds.'
        .format(time.perf_counter() - start_time))

    # Execute application
    app.exec()

    # Kill CGX after CAE exit
    gui.cgx.kill()

    # Recursively clean cached files in all subfolders
    clean.cache(p.src)
