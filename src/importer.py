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
import traceback
import os

# External modules
try:
    from PyQt5 import QtWidgets
except:
    msg = 'Please, install PyQt5 with command:\n'\
        + 'pip3 install PyQt5'
    sys.exit(msg)

# My modules
import model
import gui
import file_tools
import clean
import tests

# Split inp_doc on blocks
def split_on_blocks(inp_doc, KOM):
    keyword_blocks = []
    i = 0
    while i < len(inp_doc):
        inp_code = [inp_doc[i], ]
        while i+1 < len(inp_doc):
            if inp_doc[i+1].upper().startswith(KOM.keyword_names):
                break
            i += 1
            inp_code.append(inp_doc[i])
        if len(inp_code):
            keyword_blocks.append(inp_code)
        i += 1

    # Omit adding comments to the end of block
    for i in range(len(keyword_blocks) - 1):
        inp_code = keyword_blocks[i]
        while len(inp_code) and inp_code[-1].startswith('**'):
            line = inp_code.pop()
            keyword_blocks[i+1].insert(0, line)
        if not len(inp_code):
            keyword_blocks.pop(i)

    # print_blocks(keyword_blocks)
    return keyword_blocks

def print_blocks(keyword_blocks):
    for inp_code in keyword_blocks:
        print_block(inp_code)

def print_block(inp_code):
    for i in range(len(inp_code)):
        print(inp_code[i])
        if i == 5:
            break
    print()
    print(len(inp_code), '---')
    print()

def import_inp(s, inp_doc, KOM):
    keyword_chain = []
    impl_counter = {}

    for inp_code in split_on_blocks(inp_doc, KOM):
        keyword_name = None
        msg = 'Error parsing INP code:\n'
        for line in inp_code:
            msg += line + '\n'
            if line.upper().startswith(KOM.keyword_names):

                # Distinguish 'NODE' and 'NODE PRINT'
                if ',' in line:
                    keyword_name = line.split(',')[0]
                else:
                    keyword_name = line

                logging.debug('\n' + line)
                break
        if keyword_name is None:
            logging.error(msg)
            continue

        # Find KOM keyword path corresponding to keyword_chain
        keyword_chain.append(keyword_name)
        path, msg = KOM.get_path(keyword_chain)
        logging.debug(msg)
        if path is None:
            logging.warning('Wrong keyword {}.'.format(keyword_name))
        else:
            logging.debug('path found: ' + ', '.join([item.name for item in path]))

            # Create keyword implementations
            impl = None
            path_as_string = '' # string representation of 'path' accounting for implementations
            for j in range(len(path)):
                # Choose where to create implementation
                if impl:
                    # Implementation will be created inside another implementation
                    item = impl.get_item_by_name(path[j].name)
                else:
                    # Implementation will be created inside keyword or group
                    item = path[j]

                path_as_string += '/' + item.name
                if j == len(path) - 1: # last item is always keyword
                    # Create implementation (for example, MATERIAL-1)
                    impl = model.kom.Implementation(s, item, inp_code)
                    # logging.debug('1')
                elif item.itype == model.kom.ItemType.KEYWORD:
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

    return KOM

def import_file(p, s, w, m, t, j, file_name=''):
    if len(file_name) == 0:
        file_name = QtWidgets.QFileDialog.getOpenFileName(w, \
            'Import INP/UNV file', j.dir, \
            'INP (*.inp);;UNV (*.unv)')[0]

    if file_name is not None and len(file_name):
        w.textEdit.clear()

        # Rename job before tree regeneration
        # A new logger's handler is created here
        j.__init__(p, s, w, m, file_name[:-4] + '.inp')

        w.stop_stdout_readers()

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

        # gui.cgx.kill(w)
        has_nodes = len(m.Mesh.nodes)
        gui.cgx.open_inp(w, j.inp, has_nodes)

# Run test
if __name__ == '__main__':
    clean.screen()
    os.chdir(os.path.dirname(__file__))
    print = tests.print

    # Prepare logging
    log_file = __file__[:-3] + '.log'
    logging.getLogger().addHandler(tests.myHandler(log_file))
    logging.getLogger().setLevel(logging.WARNING)

    limit = 3000 # how many files to process
    counter = 1
    kom_xml = '../config/kom.xml'
    # examples_dir = '../../examples/abaqus/eif'
    examples_dir = '../../examples/ccx_2.16.test'

    print(log_file, 'IMPORTER (KEYWORDS PARSER) TEST\n\n')
    examples = tests.scan_all_files_in(examples_dir, '.inp', limit)
    for file_name in examples:
        print(log_file, '\n{} {}'.format(counter, file_name))
        inp_doc = file_tools.read_lines(file_name)

        # Build new clean/empty keyword object model
        k = model.kom.KOM(None, None, kom_xml)

        try:
            # Parse inp_doc end enrich existing KOM
            k = import_inp(None, inp_doc, k)
        except:
            logging.error(traceback.format_exc())
        counter += 1

    clean.cache()
