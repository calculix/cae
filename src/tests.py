# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License, version 2.

    Test for all CalculiX examples.
    Run with command:
        python3 tests.py > tests.log
"""

import os, sys, time, logging, shutil
import VTK
from mesh import Parse
from PyQt5 import QtWidgets


class myHandler(logging.Handler):

    def __init__(self):
        super().__init__()
        self.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))

    def emit(self, LogRecord):
        msg_text = self.format(LogRecord)
        print(msg_text)


class Tester:

    def __init__(self):

        # Configure logging
        logging.getLogger().addHandler(myHandler())
        logging.getLogger().setLevel(logging.DEBUG) # control the logging level

        DIRECTORY = os.path.abspath('examples')
        start = time.perf_counter() # start time

        file_list = [file_name for file_name in os.listdir(DIRECTORY) if file_name.lower().endswith('.inp')]
        for i, file_name in enumerate(file_list):
            print('\n' + '='*50 + '\n{0}: {1}'.format(i+1, file_name))

            self.test(os.path.join(DIRECTORY, file_name))

            # break # one file only
            # if i==50: break # 10 files only

        print('\nTotal {:.1f} seconds'.format(time.perf_counter()-start)) # end time

    def test(self, file_name):
        app = QtWidgets.QApplication(sys.argv)

        # Create VTK widget
        VTK = VTK.VTK()

        # Parse mesh and convert it to ugrid
        mesh = Parse(file_name) # parse mesh
        ugrid = VTK.mesh2ugrid(mesh)

        # Plot ugrid in VTK
        if ugrid:
            VTK.mapper.SetInputData(ugrid) # ugrid is our mesh data
            VTK.actionViewIso() # iso view after import


if __name__ == '__main__':

    Tester()

    # Delete cached files
    if os.path.isdir('__pycache__'):
        shutil.rmtree('__pycache__') # works in Linux as in Windows
