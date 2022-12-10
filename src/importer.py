#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2023
Distributed under GNU General Public License v3.0

File importer:
- Enrich KWT with implementations from parsed file.
- Generate new tree with keyword implementations.
- Parse mesh.
- Open model in CGX.

We can get here via menu File -> Import
or directly after application start.

TODO Only file list is logged. No importer messages.
"""

# Standard modules
import os
import re
import sys
import time
import logging
import pathlib

# External modules
from PyQt5 import QtWidgets

# My modules
from path import p
from model import m
from model.kom import KWT, KeywordTree, Implementation
import log
from utils import tests


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
        from gui.window import wf
        self.w = wf.mw # master window
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

    def parse_blocks(self):
        """Create keyword implementations in self.keyword_blocks."""
        parent = KWT.root
        messages = []

        for kwb in self.keyword_blocks:

            # Create implementations (for example, MATERIAL-1)
            kw = None
            while kw is None and parent is not None: # root has None parent
                kw = KWT.get_top_keyword_by_name(parent, kwb.keyword_name)
                if kw is not None:
                    parent = Implementation(kw, kwb.get_inp_code())
                else:
                    parent = parent.parent
            if kw is None:
                parent = KWT.root
                msg = 'Misplaced or wrong keyword {}.'\
                    .format(kwb.keyword_name)
                if msg not in messages:
                    messages.append(msg)
                    logging.warning(msg)

    def import_file(self, file_name):
        """Main method in the class.
        NOTE file_name could be None
        """
        from gui.job import j
        if file_name is None and self.w is None:
            msg = 'file_name and self.w are None.'
            logging.warning(msg)
            raise SystemExit(msg)

        if file_name is None:
            file_name = QtWidgets.QFileDialog.getOpenFileName(self.w, \
                'Import INP/FBD/FBL/UNV file', p.examples, \
                'INP (*.inp);;FBD/FBL (*.fbd *.fbl);;UNV (*.unv)')[0]

        if file_name is not None and len(file_name):
            # self.w.textEdit.clear()

            # Update job instance before the tree regeneration
            # A new logger handler is created here
            j.generate(file_name[:-4] + '.inp')

            # Show model name in window title
            title = 'CalculiX Advanced Environment - ' + j.name
            self.w.setWindowTitle(title)

            from gui import stdout
            stdout.stop_readers()

            # Generate new KWT without implementations
            KWT.__init__()

            inp_files = []
            if file_name.lower().endswith('.inp'):
                inp_files.append(j.inp)

            # Convert UNV to INP
            if file_name.lower().endswith('.unv'):
                j.convert_unv()
                inp_files.append(j.inp)
                if not os.path.isfile(j.inp):
                    logging.error('Can not convert\n' + j.unv)
                    return

            # Pass FBD to CGX
            from gui import cgx
            if file_name.lower().endswith(('.fbd', '.fbl')):
                """Get list of newly created or updated files."""
                ext = ('.fbd', '.fbl', '.inp', '.nam', '.msh', '.sur', '.con', '.dlo', '.bou')
                flist_before = {}
                for f in os.listdir(os.path.dirname(file_name)):
                    fname = pathlib.Path(f)
                    flist_before[f] = fname.stat().st_ctime
                cgx.restart_and_read_fbd(file_name)
                # cgx.restart_empty()
                # time.sleep(1)
                # cgx.read_fbd_file(file_name)
                time.sleep(1)
                flist_after = os.listdir(os.path.dirname(file_name))
                for f in flist_after:
                    if not f.endswith(ext):
                        continue
                    if not f in flist_before:
                        inp_files.append(f)
                    else:
                        modified_before = flist_before[f]
                        fname = pathlib.Path(f)
                        modified_after = fname.stat().st_ctime
                        if modified_before != modified_after:
                            inp_files.append(f)

            # Get INP code from all model files
            logging.info('Loading model files...')
            inp_doc = []
            for f in inp_files:
                logging.debug(f)
                inp_doc.extend(read_lines(f))

            # Split INP code on blocks - fill self.keyword_blocks
            self.split_on_blocks(inp_doc)

            # Parse keyword_blocks and enrich KWT with parsed objects
            self.parse_blocks()

            # Write generated INP code as file and reopen it normally
            if file_name.lower().endswith(('.fbd', '.fbl')):
                j.write_input(KWT.get_inp_code_as_lines(), file_name=j.inp)

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

            if os.path.isfile(j.inp):
                has_nodes = len(m.Mesh.nodes)
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
        KWT = KeywordTree()

        # Parse inp_doc end enrich existing KWT
        i = Importer()
        inp_doc = read_lines(file_name)
        i.split_on_blocks(inp_doc) # fill keyword_blocks
        i.parse_blocks()

    msg = '\n{} INP files.'
    log.print_to_file(log_file, msg.format(len(examples)))


if __name__ == '__main__':
    test() # run test
