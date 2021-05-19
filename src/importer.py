#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2021
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

# Standard modules
import os
import io
import re
import sys
import time
import logging

# External modules
from PyQt5 import QtWidgets

# My modules
import model
import gui
import tests


# Keyword block
class Block:

    def __init__(self, keyword_name, comments, lead_line, data_lines):
        self.keyword_name = keyword_name # string
        self.comments = comments # list
        self.lead_line = lead_line # string
        self.data_lines = data_lines # list
    
    def get_inp_code(self):
        inp_code = []
        inp_code.extend(self.comments)
        inp_code.append(self.lead_line)
        inp_code.extend(self.data_lines)
        return inp_code

    def print_debug_info(self):

        sys.stdout.write('\nCOMMENTS:\n')
        for line in self.comments:
            sys.stdout.write(line + '\n')

        sys.stdout.write('LEAD:\n')
        sys.stdout.write(self.lead_line + '\n')

        sys.stdout.write('DATA:\n')
        for line in self.data_lines:
            sys.stdout.write(line + '\n')


class Importer:

    def __init__(self, s, f, m, t, j):
        self.p = s.getp() # path
        self.s = s # settings
        self.f = f # window factory
        if f is not None:
            self.w = f.mw # master window
        else:
            self.w = None
        self.m = m # model
        self.t = t # tree
        self.j = j # job
        self.keyword_blocks = []

    # Split inp_doc on blocks
    def split_on_blocks(self, inp_doc):
        self.keyword_blocks = []
        i = 0
        regex = r'^\*[\w\s-]+'
        while i < len(inp_doc):
            match = re.match(regex, inp_doc[i])
            if match is not None:
                keyword_name = match.group(0).strip()

                # Comments before the block
                comments = []
                j = 0 # amount of comment lines
                while i > j and inp_doc[i-j-1].startswith('**'):
                    j += 1
                    comments.insert(0, inp_doc[i-j])

                # Lead line - a line(s) with keyword
                lead_line = inp_doc[i].rstrip()
                while lead_line.endswith(','):
                    i += 1
                    lead_line = lead_line + ' ' + inp_doc[i].rstrip()

                i += 1
                start = i # start of data lines
                while i < len(inp_doc):
                    match = re.match(regex, inp_doc[i])
                    if match is not None:
                        i -= 1
                        break # reached next keyword
                    i += 1

                # Comments after the block (if not EOF)
                if i < len(inp_doc) - 1:
                    while inp_doc[i].startswith('**'):
                        i -= 1
                end = i # index of block end

                data_lines = inp_doc[start:end+1]
                b = Block(keyword_name, comments, lead_line, data_lines)
                # b.print_debug_info()
                self.keyword_blocks.append(b)

            i += 1

    # Create keyword implementations
    def import_inp(self):
        parent = self.m.KOM.root
        impl_counter = {}
        messages = []

        for kwb in self.keyword_blocks:

            # Create implementations (for example, MATERIAL-1)
            kw = None
            while kw is None and parent is not None: # root has None parent
                kw = self.m.KOM.get_top_keyword_by_name(parent, kwb.keyword_name)
                if kw is not None:
                    parent = model.kom.Implementation(self.s, kw, kwb.get_inp_code())
                else:
                    parent = parent.parent
            if kw is None:
                parent = self.m.KOM.root
                msg = 'Misplaced or wrong keyword {}.'\
                    .format(kwb.keyword_name)
                if msg not in messages:
                    messages.append(msg)
                    logging.warning(msg)

    # Main method in the class
    def import_file(self, file_name):
        if file_name is None and \
            (self.w is None or self.j is None):
            raise SystemExit

        if file_name is None:
            file_name = QtWidgets.QFileDialog.getOpenFileName(self.w, \
                'Import INP/UNV file', self.j.dir, \
                'INP (*.inp);;UNV (*.unv)')[0]

        if file_name is not None and len(file_name):
            self.w.textEdit.clear()

            # Rename job before tree regeneration
            # A new logger's handler is created here
            self.j.__init__(self.s, self.f,
                self.m, file_name[:-4] + '.inp')

            self.f.stop_stdout_readers()

            # Convert UNV to INP
            if file_name.lower().endswith('.unv'):
                self.j.convert_unv()
                if not os.path.isfile(self.j.inp):
                    logging.error('Can not convert\n' + self.j.unv)
                    return

            # Show model name in window's title
            self.w.setWindowTitle('CalculiX Advanced Environment - ' \
                + self.j.name)

            # Generate new KOM without implementations
            self.m.KOM = model.kom.KOM(self.s)

            # Get INP code and split it on blocks
            logging.info('Loading model\n{}'.format(self.j.inp))
            inp_doc = read_lines(self.j.inp)
            self.split_on_blocks(inp_doc) # fill keyword_blocks

            # Parse INP and enrich KOM with parsed objects
            self.import_inp()

            # Add parsed implementations to the tree
            self.t.generateTreeView(self.m)

            # Parse mesh
            self.m.Mesh = model.parsers.mesh.Mesh(blocks=self.keyword_blocks)

            # Open a new non-empty model in CGX
            if not self.s.start_cgx_by_default:
                msg = '"Settings -> Start CGX by default" is unchecked.'
                logging.warning(msg)
                return
            if not len(self.m.Mesh.nodes):
                msg = 'Empty mesh, CGX will not start!'
                logging.warning(msg)
                return

            has_nodes = len(self.m.Mesh.nodes)
            gui.cgx.open_inp(self.p, self.f, self.j.inp, has_nodes)


# Recurcively reads all the file lines and its includes.
# Does not omit comments and empty lines.
def read_lines(INP_file):
    INP_file = os.path.abspath(INP_file)
    if not os.path.isfile(INP_file):
        msg_text = 'File not found: ' + INP_file
        logging.error(msg_text)
        return []

    lines = []
    with open(INP_file, 'r', errors='ignore') as f:
        for line in f.readlines():
            line = line.strip()
            lines.append(line)

            # Append lines from include file
            if line.upper().startswith('*INCLUDE'):
                inc_file = line.split('=')[1].strip()
                inc_file = os.path.normpath(
                    os.path.join(os.path.dirname(INP_file), inc_file))
                lines.extend(read_lines(inc_file))

    return lines

# Test importer on all CalculiX examples
@tests.test_wrapper()
def test():
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    start_time = time.perf_counter()
    print = tests.print
    m = model.Model() # generate FEM model

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

    print(log_file, 'IMPORTER (KEYWORDS PARSER) TEST\n')
    examples = tests.scan_all_files_in(examples_dir, '.inp', limit)

    lines_count = 0
    for file_name in examples:
        counter += 1

        # Build new clean/empty keyword object model
        m.KOM = model.kom.KOM(None, None, kom_xml='../config/kom.xml')

        # Parse inp_doc end enrich existing KOM
        i = Importer(None, None, None, m, None, None)
        inp_doc = read_lines(file_name)
        i.split_on_blocks(inp_doc) # fill keyword_blocks
        i.import_inp()

        # Log only problematic INP files
        log_contents = log_capture_string.getvalue()
        if len(log_contents) != lines_count:
            log_contents = log_contents[lines_count:]
            lines_count = log_capture_string.tell()
            relpath = os.path.relpath(file_name, start=os.getcwd())
            print(log_file, '{} {}'.format(counter, relpath))
            print(log_file, log_contents)
    log_capture_string.close()

    msg = '\n{} INP files.'
    print(log_file, msg.format(len(examples)))

# Run test
if __name__ == '__main__':
    test()
