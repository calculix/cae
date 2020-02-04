#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, January 2020
    Distributed under GNU General Public License v3.0

    Test mesh parser on all CalculiX examples.
    Ctrl + F5 to Run.
"""


import os
import time
import logging
import PyQt5

from gui import vtk_widget
from model.parsers import mesh
from log import myHandler
from log import print
import clean


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

        # TODO Doesn't work!
        # VTK = vtk_widget.VTK()
        # VTK.plotMesh(m)

    print('\nTotal {:.1f} seconds.'
        .format(time.perf_counter() - start_time))
    clean.cache()
