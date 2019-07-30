# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    Different tests for all example input files.
    Run with command:
        python3 test.py > test.log
"""

import os, sys, vtk
import ccx_vtk, ccx_dom, ccx_cae_ie, ccx_mesh, ccx_cae_log
from PyQt5 import QtWidgets, uic, QtCore, QtGui


# Test mesh parser for all example input files
def test_mesh_parser(logger, file_name):
    print('-'*17 + '\ntest_mesh_parser:\n' + '-'*17)
    mesh = ccx_mesh.Parse(file_name) # parse mesh
    logger.messages(mesh.msg_list) # process list of log messages


# Test vtk camera for all example files
def test_vtk_camera(logger, file_name):
    print('-'*16 + '\ntest_vtk_camera:\n' + '-'*16)
    app = QtWidgets.QApplication(sys.argv)

    # Create VTK widget
    VTK = ccx_vtk.VTK()

    # Parse mesh and convert it to ugrid
    mesh = ccx_mesh.Parse(file_name) # parse mesh
    msgs, ugrid = VTK.mesh2ugrid(mesh)

    # Plot ugrid in VTK
    if ugrid:
        VTK.mapper.SetInputData(ugrid) # ugrid is our mesh data
        VTK.actionViewIso() # iso view after import
        logger.messages(VTK.msg_list) # process list of log messages


if __name__ == '__main__':

    # Clean cached files before start
    os.system('py3clean .')

    logger = ccx_cae_log.logger()
    DIRECTORY = './examples/'

    file_list = [file_name for file_name in os.listdir(DIRECTORY) if file_name.endswith('.inp')]
    for i, file_name in enumerate(file_list):
        print('\n' + '='*50 + '\n{0}: {1}'.format(i+1, file_name))

        # Test mesh parser
        test_mesh_parser(logger, DIRECTORY + file_name)

        # Test vtk camera
        test_vtk_camera(logger, DIRECTORY + file_name)

        # break # one file only
        # if i==50: break # 10 files only
