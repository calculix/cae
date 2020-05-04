#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, May 2020
Distributed under GNU General Public License v3.0

Test mesh parser on all CalculiX examples.
Ctrl + F5 to Run. """

# Standard modules
import os
import sys
import time
import logging
import PyQt5

# My modules
from model.parsers import mesh
import clean


# How many files to process
limit = 50000000


log_file = os.path.join(os.path.dirname(__file__), 'tests.log')


# Configure logging to emit messages via 'print' method
class myHandler(logging.Handler):

    def __init__(self):
        super().__init__()
        fmt = logging.Formatter('%(levelname)s: %(message)s')
        self.setFormatter(fmt)

        # Remove old log file
        if os.path.isfile(log_file):
            os.remove(log_file)

    def emit(self, LogRecord):
        print(self.format(LogRecord))


# Redefine print method to write logs to file
def print(*args):
    line = ' '.join([str(arg) for arg in args])
    line = line.rstrip() + '\n'
    with open(log_file, 'a') as f:
        f.write(line)
    sys.stdout.write(line)


# List all .ext-files here and in all subdirectories
def scan_all_files_in(start_folder, ext):
    all_files = []
    for f in os.scandir(start_folder):
        if f.is_dir():
            for ff in scan_all_files_in(f.path, ext):
                all_files.append(ff)
        elif f.is_file() and f.name.endswith(ext):
            all_files.append(f.path)
    return sorted(all_files)[:limit]


# Test
if __name__ == '__main__':
    start_time = time.perf_counter()

    # Prepare logging
    logging.getLogger().addHandler(myHandler())
    logging.getLogger().setLevel(logging.INFO)

    print('MESH PARSER TEST\n\n')
    counter = 0
    for file_name in scan_all_files_in('examples', '.inp'):

        # Skip some files
        if 'materials.inp' in file_name:
            continue
        if 'default.inp' in file_name:
            continue

        counter += 1
        relpath = os.path.relpath(file_name, start=os.getcwd())
        print('\n{}\n{}: {}'.format('='*50, counter, relpath))

        # Parse mesh and plot it in VTK
        app = PyQt5.QtWidgets.QApplication([])
        m = mesh.Mesh(INP_file=file_name) # parse mesh

    print('\nTotal {:.1f} seconds.'
        .format(time.perf_counter() - start_time))
    clean.cache()
