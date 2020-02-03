#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, January 2020
Distributed under GNU General Public License v3.0
Methods for cleaning up temporary/unused files/folders """

import os
import sys
import shutil


# Clean screen
def screen():
    os.system('cls' if os.name=='nt' else 'clear')


# Recursively delete cached files in all subfolders
def cache(folder=None):
    if not folder:
        folder = os.getcwd()
    pycache = os.path.join(folder, '__pycache__')
    if os.path.isdir(pycache):
        shutil.rmtree(pycache) # works in Linux as in Windows

    # Recursively clear cache in child folders
    for f in os.scandir(folder):
        if f.is_dir():
            cache(f.path)


# Cleaup trash files in startFolder and all subfolders
def files(startFolder=None):
    extensions = (  '.12d', '.cvg', '.dat', '.vwf', '.out', '.nam', '.inp1', '.inp2',
                    '.sta', '.equ', '.eig', '.stm', '.mtx', '.net', '.inp0', '.rin',
                    '.fcv', 'dummy' )
    if not startFolder:
        startFolder = os.getcwd()
    for f in os.scandir(startFolder):
        if f.is_dir(): # if folder
            files(f.path)
        elif f.is_file() and f.name.endswith(extensions):
            try:
                os.remove(f.path)
                sys.__stdout__.write('Delelted: ' + f.path + '\n')
            except:
                sys.__stdout__.write(f.path + ': ' + sys.exc_info()[1][1] + '\n')


# Cleaup old result files
def results():
    extensions = ('.frd', '.vtk', '.vtu')
    for f in os.scandir('.'):
        if f.name.endswith(extensions):
            try:
                os.remove(f.path)
                sys.__stdout__.write('Delelted: ' + f.path + '\n')
            except:
                sys.__stdout__.write(f.path + ': ' + sys.exc_info()[1][1] + '\n')
