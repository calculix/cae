#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, January 2020
    Distributed under GNU General Public License v3.0

    Test for all CalculiX examples.
    Ctrl + F5 to Run.
"""


from PyQt5 import QtWidgets
import time, logging, os

from gui import vtk_widget
from model.parsers import mesh
import clean
from log import myHandler
from log import print


# How many files to process
limit = 50000000


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

    counter = 0
    for file_name in scan_all_files_in('examples', '.inp'):

        # Skip some files
        if 'materials.inp' in file_name:
            continue
        if 'default.inp' in file_name:
            continue

        counter += 1
        print('\n{}\n{}: {}'.format('='*50, counter, file_name))

        # Parse mesh and plot it in VTK
        app = QtWidgets.QApplication([])
        m = mesh.Mesh(INP_file=file_name) # parse mesh
        # vtk_widget.VTK().plotMesh(m)

    print('\nTotal {:.1f} seconds.'
        .format(time.perf_counter() - start_time))
    clean.cache()
