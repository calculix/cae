#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2023
Distributed under GNU General Public License v3.0

Methods for cleaning up temporary/unused files/folders etc.
"""

import os
import sys
import shutil


def screen():
    """Clean screen."""
    os.system('cls' if os.name=='nt' else 'clear')


def cache(folder=None):
    """Recursively delete cached files in all subfolders."""
    if not folder:
        folder = os.getcwd()
    pycache = os.path.join(folder, '__pycache__')
    if os.path.isdir(pycache):
        shutil.rmtree(pycache) # works in Linux as in Windows

    # Recursively clear cache in child folders
    for f in os.scandir(folder):
        if f.is_dir():
            cache(f.path)


def files(startFolder=None):
    """Cleaup trash files in startFolder and all subfolders."""
    if not startFolder:
        startFolder = os.getcwd()
    extensions = (  '.12d', '.cvg', '.dat', '.vwf', '.out', '.nam', '.inp1', '.inp2',
                    '.sta', '.equ', '.eig', '.stm', '.mtx', '.net', '.inp0', '.rin',
                    '.fcv', 'dummy' )
    for f in os.scandir(startFolder):
        if f.is_dir(): # if folder
            files(f.path)
        elif f.is_file() and f.name.endswith(extensions):
            try:
                os.remove(f.path)
                sys.__stdout__.write('Delelted: ' + f.path + '\n')
            except:
                sys.__stdout__.write(f.path + ': ' + sys.exc_info()[1][1] + '\n')


def results():
    """Cleaup old result files."""
    extensions = ('.frd', '.vtk', '.vtu')
    for f in os.scandir('.'):
        if f.name.endswith(extensions):
            try:
                os.remove(f.path)
                sys.__stdout__.write('Delelted: ' + f.path + '\n')
            except:
                sys.__stdout__.write(f.path + ': ' + sys.exc_info()[1][1] + '\n')
