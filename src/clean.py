#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    Methods for cleaning up temporary/unused files/folders.
"""


import os, sys, shutil


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
    for f in os.listdir(folder):
        f = os.path.join(folder, f)
        if os.path.isdir(f):
            cache(f)


# Cleaup trash files in startFolder and all subfolders
def files(startFolder=None):
    extensions = (  '.12d', '.cvg', '.dat', '.vwf', '.out', '.nam', '.inp1', '.inp2',
                    '.sta', '.log', '.equ', '.eig', '.stm', '.mtx', '.net', '.inp0',
                    '.rin', '.fcv'  )
    if not startFolder:
        startFolder = os.getcwd()
    for f in os.listdir(startFolder):
        f = os.path.join(startFolder, f)
        if os.path.isdir(f): # if folder
            files(f)
        elif f.endswith(extensions):
            try:
                os.remove(f)
                sys.__stdout__.write('Delelted: ' + f + '\n')
            except:
                sys.__stdout__.write(f + ': ' + sys.exc_info()[1][1] + '\n')


# Cleaup old result files
def results():
    extensions = ('.frd', '.vtk', '.vtu')
    for f in os.listdir('.'):
        if f.endswith(extensions):
            try:
                os.remove(f)
                sys.__stdout__.write('Delelted: ' + f + '\n')
            except:
                sys.__stdout__.write(f + ': ' + sys.exc_info()[1][1] + '\n')
