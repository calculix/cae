# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    Different tests for all example input files.
    Run with command:
        python3 test.py > test.log
"""

import os, sys, time, logging
import ccx_vtk, ccx_mesh
from PyQt5 import QtWidgets


class myHandler(logging.Handler):

    def __init__(self, parent):
        super().__init__()

    def emit(self, record):
        msg = self.format(record)
        print(msg)


class Tester:


    # Test mesh parser for all example input files
    def test_mesh_parser(self, file_name):
        print('-'*17 + '\ntest_mesh_parser:\n' + '-'*17)
        mesh = ccx_mesh.Parse(file_name) # parse mesh


    # Test vtk camera for all example files
    def test_vtk_camera(self, file_name):
        # print('-'*16 + '\ntest_vtk_camera:\n' + '-'*16)
        app = QtWidgets.QApplication(sys.argv)

        # Create VTK widget
        VTK = ccx_vtk.VTK()

        # Parse mesh and convert it to ugrid
        mesh = ccx_mesh.Parse(file_name) # parse mesh
        ugrid = VTK.mesh2ugrid(mesh)

        # Plot ugrid in VTK
        if ugrid:
            VTK.mapper.SetInputData(ugrid) # ugrid is our mesh data
            VTK.actionViewIso() # iso view after import


    def __init__(self):

        # Configure logging
        formatter = logging.Formatter('%(levelname)s: %(message)s')
        handler = myHandler(self)
        handler.setFormatter(formatter)
        logging.getLogger().addHandler(handler)
        logging.getLogger().setLevel(logging.INFO) # control the logging level

        DIRECTORY = './examples/'
        start = time.perf_counter() # start time

        file_list = [file_name for file_name in os.listdir(DIRECTORY) if file_name.endswith('.inp')]
        for i, file_name in enumerate(file_list):
            print('\n' + '='*50 + '\n{0}: {1}'.format(i+1, file_name))

            # Test mesh parser
            # self.test_mesh_parser(DIRECTORY + file_name)

            # Test vtk camera
            self.test_vtk_camera(DIRECTORY + file_name)

            # break # one file only
            # if i==50: break # 10 files only

        print('\nTotal {:.1f} seconds'.format(time.perf_counter()-start)) # end time


if __name__ == '__main__':

    # Clean cached files before start
    os.system('py3clean .')

    Tester()
