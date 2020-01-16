#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, January 2020
    Distributed under GNU General Public License v3.0

    Test for all CalculiX examples.
    Ctrl + F5 to Run.
"""


from PyQt5 import QtWidgets
import time, logging, glob, os

from gui import vtk_widget
from model.parsers import mesh
import clean


# Configure logging to emit messages via 'print' method
class myHandler(logging.Handler):

    def __init__(self):
        super().__init__()
        self.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))

    def emit(self, LogRecord):
        msg_text = self.format(LogRecord)
        print(msg_text)


# Parse mesh and plot it in VTK
def test(file_name):
    app = QtWidgets.QApplication([])
    m = mesh.Mesh(INP_file=file_name) # parse mesh
    vtk_widget.VTK().plotMesh(m)


# Redefine print method to write logs to file
def print(*args):
    line = ' '.join([str(arg) for arg in args])
    with open('src/tests.log', 'a') as f:
        f.write(str(line) + '\n')


# Test
if __name__ == '__main__':
    try: os.remove('src/tests.log')
    except: pass
    start_time = time.perf_counter()
    logging.getLogger().addHandler(myHandler())
    logging.getLogger().setLevel(logging.DEBUG)
    file_list = glob.glob('examples/**/*.inp', recursive=True)
    file_list = sorted(file_list)

    i = 0
    for file_name in file_list:

        # Skip some files
        if 'mkraska' in file_name:
            continue
        if 'materials.inp' in file_name:
            continue

        print('\n' + '='*50 + '\n{0}: {1}'.format(i+1, file_name))
        test(file_name)
        i += 1
        # if i==10: break # 10 files only

    print('\nTotal {:.1f} seconds.'
        .format(time.perf_counter() - start_time))
    clean.cache()
