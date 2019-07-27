# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    Test mesh parser for all example input files.
    Run with command:
        python3 test.py > test.txt
"""

import os, ccx_mesh, ccx_log


class TestMesh:

    def __init__(self):

        DIRECTORY = './examples/'
        for file_name in os.listdir(DIRECTORY):
            print('='*50)
            print(file_name)
            
            # Parse mesh
            mesh = ccx_mesh.Parse(DIRECTORY + file_name)
            
            # Process list of INFO and ERROR messages
            ccx_log.logger().messages(mesh.msg_list)

            # break # one file only


if __name__ == '__main__':

    # Clean cached files before start
    os.system('py3clean .')

    # Run test
    TestMesh()
