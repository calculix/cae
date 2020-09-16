#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2020
Distributed under GNU General Public License v3.0

File importer:
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

# TODO Fix import:
# *Material, name=MATERIAL

# Standard modules
import os
import io
import re
import sys
import time
import logging
import traceback

# External modules
try:
    from PyQt5 import QtWidgets
except:
    msg = 'Please, install PyQt5 with command:\n'\
        + 'pip3 install PyQt5'
    print(msg)
    raise SystemExit # the best way to exit

# My modules
import model
import gui
import file_tools
import clean
import tests


class KeywordBlock:

    def __init__(self, keyword_name, inp_code):
        self.inp_code = inp_code
        self.keyword_name = keyword_name

    def print_debug_info(self):
        i = 0
        for j in range(len(self.inp_code)):
            if not self.inp_code[j].startswith('**'):
                break
            else:
                i += 1
        if i > 0:
            comment = '\n'.join(self.inp_code[:i])
            logging.debug(comment)
        self.inp_code = '\n'.join(self.inp_code[i:i+1])
        logging.debug(self.inp_code + '\n')


# Split inp_doc on blocks
def split_on_blocks(inp_doc):
    keyword_blocks = []
    i = 0
    regex = r'^\*[\w\s-]+'
    while i < len(inp_doc):
        match = re.match(regex, inp_doc[i])
        if match is not None:
            keyword_name = match.group(0).strip()

            # Comments before the block
            j = 0
            while i > j and inp_doc[i-j-1].startswith('**'):
                j += 1
            start = i - j # index of block start

            i += 1
            while i < len(inp_doc):
                match = re.match(regex, inp_doc[i])
                if match is not None:
                    i -= 1
                    break
                i += 1

            # Comments after the block (if not EOF)
            if i < len(inp_doc) - 1:
                while inp_doc[i].startswith('**'):
                    i -= 1
            end = i # index of block end

            inp_code = inp_doc[start:end+1]
            kwb = KeywordBlock(keyword_name, inp_code)
            keyword_blocks.append(kwb)

        i += 1

    return keyword_blocks

def import_inp(s, inp_doc, KOM):
    parent = KOM.root
    impl_counter = {}

    for kwb in split_on_blocks(inp_doc):

        # Create implementations (for example, MATERIAL-1)
        kw = None
        while kw is None and parent is not None: # root has None parent
            kw = KOM.get_top_keyword_by_name(parent, kwb.keyword_name)
            if kw is not None:
                parent = model.kom.Implementation(s, kw, kwb.inp_code)
            else:
                parent = parent.parent
        if kw is None:
            parent = KOM.root
            logging.warning('Misplaced or wrong keyword {}.'\
                .format(kwb.keyword_name))

    # KOM.test()
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
        w.setWindowTitle('CalculiX Advanced Environment - ' + j.name)

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

# Test importer on all CalculiX examples
if __name__ == '__main__':
    clean.screen()
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    start_time = time.perf_counter()
    print = tests.print

    # Prepare logging
    log_file = __file__[:-3] + '.log'
    h = tests.myHandler(log_file) # remove old log file
    log_capture_string = io.StringIO()
    ch = logging.StreamHandler(log_capture_string)
    ch.setLevel(logging.DEBUG)
    fmt = logging.Formatter('%(levelname)s: %(message)s')
    ch.setFormatter(fmt)
    logging.getLogger().addHandler(ch)

    limit = 50000 # how many files to process
    # examples_dir = '../../examples/ccx/test'
    # examples_dir = '../../examples/abaqus/eif'
    # examples_dir = '../../examples/yahoo'
    examples_dir = '../../examples'
    counter = 0
    kom_xml = '../config/kom.xml'

    print(log_file, 'IMPORTER (KEYWORDS PARSER) TEST\n\n')
    examples = tests.scan_all_files_in(examples_dir, '.inp', limit)

    lines_count = 0
    for file_name in examples:
        counter += 1
        relpath = os.path.relpath(file_name, start=os.getcwd())
        inp_doc = file_tools.read_lines(file_name)

        # Build new clean/empty keyword object model
        k = model.kom.KOM(None, None, kom_xml)

        try:
            # Parse inp_doc end enrich existing KOM
            k = import_inp(None, inp_doc, k)
        except:
            logging.error(traceback.format_exc())

        # Log only problematic INP files
        log_contents = log_capture_string.getvalue()
        if len(log_contents) != lines_count:
            log_contents = log_contents[lines_count:]
            lines_count = log_capture_string.tell()
            print(log_file, '\n{} {}'.format(counter, relpath))
            print(log_file, log_contents)
    log_capture_string.close()

    msg = '\n{} INP files. Total time {}.'
    delta = tests.time_delta(start_time)
    print(log_file, msg.format(len(examples), delta))
    clean.cache()
