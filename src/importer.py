#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, July 2020
Distributed under GNU General Public License v3.0

CalculiX CAE file importer:
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

# Standard modules
import logging

# External modules
try:
    from PyQt5 import QtWidgets
except:
    msg = 'Please, install PyQt5 with command:\n'\
        + 'pip3 install PyQt5'
    sys.exit(msg)

# My modules
import gui
import model
import file_tools

def import_inp(s, INP_doc, KOM):
    keyword_chain = []
    impl_counter = {}
    for i in range(len(INP_doc)):
        line = INP_doc[i]

        # Parse keyword
        if line.startswith('*'):
            logging.debug('\n' + line)

            # Distinguish 'NODE' and 'NODE PRINT'
            if ',' in line:
                keyword_name = line.split(',')[0]
            else:
                keyword_name = line

            # Find KOM keyword path corresponding to keyword_chain
            keyword_chain.append(keyword_name)
            path, msg = KOM.get_path(keyword_chain)
            logging.debug(msg)
            if path is not None:
                logging.debug('path found: ' + ', '.join([item.name for item in path]))

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
                        impl = model.kom.implementation(s, item, INP_code)
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

def import_file(p, s, w, m, t, j, file_name=''):
    if len(file_name) == 0:
        file_name = QtWidgets.QFileDialog.getOpenFileName(w, \
            'Import INP/UNV file', j.dir, \
            'INP (*.inp);;UNV (*.unv)')[0]

    if file_name is not None and len(file_name):
        gui.cgx.kill(w) # close old CGX
        w.textEdit.clear() # clear logs in the textEdit
        gui.log.stop_stdout_readers(w)

        # Rename job before tree regeneration
        # A new logger's handler is created here
        j.initialize(file_name[:-4] + '.inp')

        # Convert UNV to INP
        if file_name.lower().endswith('.unv'):
            j.convert_unv()
            if not os.path.isfile(j.inp):
                logging.error('Can not convert\n' + j.unv)
                return

        # Show model name in window's title
        w.setWindowTitle('CalculiX CAE - ' + j.name)

        # Generate new KOM without implementations
        m.KOM = model.kom.KOM(p, s)

        # Parse INP and enrich KOM with parsed objects
        logging.info('Loading model\n{}'.format(j.inp))
        lines = file_tools.read_lines(j.inp)
        import_inp(s, lines, m.KOM) # pass whole INP-file to the parser

        # Add parsed implementations to the tree
        t.generateTreeView(m)

        # Parse mesh
        m.Mesh = model.parsers.mesh.Mesh(INP_file=j.inp)

        # Open a new non-empty model in CGX
        if not s.start_cgx_by_default:
            logging.warning('"Settings -> Start CGX by default" is unchecked.')
            return
        if not len(m.Mesh.nodes):
            logging.warning('Empty mesh, CGX will not start!')
            return
        j.cgx_inp(m)