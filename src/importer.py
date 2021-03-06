#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

File importer:
- Enrich KOM with implementations from parsed file.
- Generate new tree with keyword implementations.
- Parse mesh.
- Open model in CGX.

We can get here via menu File -> Import
or directly after application start.

TODO
Only file list is logged. No importer messages.
"""

# Standard modules
import os
import re
import sys
import logging

# External modules
from PyQt5 import QtWidgets

# My modules
from model import m
from model.kom import KOM, KeywordObjectModel, Implementation
import log
import tests


class Block:
    """Keyword block."""

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

    def __init__(self):
        from gui.window import factory
        self.w = factory.mw # master window
        self.keyword_blocks = []

    def split_on_blocks(self, inp_doc):
        """Split inp_doc on blocks."""
        self.keyword_blocks = []
        i = 0
        regex = r'^\*[\w\s-]+'
        while i < len(inp_doc):
            match = re.match(regex, inp_doc[i])
            if match is not None:
                keyword_name = match.group(0).strip()

                # Comments before the block
                comments = []
                counter = 0 # amount of comment lines
                while i > counter and inp_doc[i-counter-1].startswith('**'):
                    counter += 1
                    comments.insert(0, inp_doc[i-counter])

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

    def import_inp(self):
        """Create keyword implementations."""
        parent = KOM.root
        messages = []

        for kwb in self.keyword_blocks:

            # Create implementations (for example, MATERIAL-1)
            kw = None
            while kw is None and parent is not None: # root has None parent
                kw = KOM.get_top_keyword_by_name(parent, kwb.keyword_name)
                if kw is not None:
                    parent = Implementation(kw, kwb.get_inp_code())
                else:
                    parent = parent.parent
            if kw is None:
                parent = KOM.root
                msg = 'Misplaced or wrong keyword {}.'\
                    .format(kwb.keyword_name)
                if msg not in messages:
                    messages.append(msg)
                    logging.warning(msg)

    def import_file(self, file_name):
        """Main method in the class"""
        from gui.job import j
        if file_name is None and self.w is None:
            msg = 'file_name and self.w are None.'
            raise SystemExit(msg)

        if file_name is None:
            file_name = QtWidgets.QFileDialog.getOpenFileName(self.w, \
                'Import INP/UNV file', j.dir, \
                'INP (*.inp);;UNV (*.unv)')[0]

        if file_name is not None and len(file_name):
            self.w.textEdit.clear()

            # Rename job before tree regeneration
            # A new logger handler is created here
            # TODO Do not call job __init__ twice
            j.__init__(file_name[:-4] + '.inp')

            from gui import stdout
            stdout.stop_readers()

            # Convert UNV to INP
            if file_name.lower().endswith('.unv'):
                j.convert_unv()
                if not os.path.isfile(j.inp):
                    logging.error('Can not convert\n' + j.unv)
                    return

            # Show model name in window title
            title = 'CalculiX Advanced Environment - ' + j.name
            self.w.setWindowTitle(title)

            # Generate new KOM without implementations
            KOM = KeywordObjectModel()

            # Get INP code and split it on blocks
            logging.info('Loading model\n{}'.format(j.inp))
            inp_doc = read_lines(j.inp)
            self.split_on_blocks(inp_doc) # fill keyword_blocks

            # Parse INP and enrich KOM with parsed objects
            self.import_inp()

            # Add parsed implementations to the tree
            from gui.tree import t
            t.generateTreeView()

            # Parse mesh
            # TODO Use Mesh instead of m.Mesh
            from model.parsers import mesh
            m.Mesh = mesh.Mesh(blocks=self.keyword_blocks)

            # Open a new non-empty model in CGX
            from settings import s
            if not s.start_cgx_by_default:
                msg = '"Settings -> Start CGX by default" is unchecked.'
                logging.warning(msg)
                return
            if not len(m.Mesh.nodes):
                msg = 'Empty mesh, CGX will not start!'
                logging.warning(msg)
                return

            has_nodes = len(m.Mesh.nodes)
            from gui import cgx
            cgx.open_inp(j.inp, has_nodes)


# Prepare to import model
i = Importer()


def read_lines(INP_file):
    """Recurcively reads all the INP file lines and its includes.
    Does not omit comments and empty lines.
    """
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


@tests.test_wrapper()
def test():
    """Test importer on all CalculiX examples."""

    # Prepare logging
    log_file = __file__[:-3] + '.log'
    log.stop_logging()
    log.add_my_handler(logging.WARNING)
    log.print_to_file(log_file, 'IMPORTER (KEYWORDS PARSER) TEST')

    limit = 50000 # how many files to process
    examples_dir = '../../examples'
    # examples_dir = '../../examples/ccx/test'
    # examples_dir = '../../examples/abaqus/eif'
    # examples_dir = '../../examples/yahoo'
    examples = tests.scan_all_files_in(examples_dir, '.inp', limit)

    counter = 0
    for file_name in examples:
        counter += 1
        relpath = os.path.relpath(file_name, start=os.getcwd())
        log.print_to_file(log_file, '\n{} {}'.format(counter, relpath))

        # Build new clean/empty keyword object model
        KOM = KeywordObjectModel()

        # Parse inp_doc end enrich existing KOM
        i = Importer()
        inp_doc = read_lines(file_name)
        i.split_on_blocks(inp_doc) # fill keyword_blocks
        i.import_inp()

    msg = '\n{} INP files.'
    log.print_to_file(log_file, msg.format(len(examples)))


if __name__ == '__main__':
    test() # run test
