# -*- coding: utf-8 -*-

""" © Ihor Mirzov, September 2019
Distributed under GNU General Public License v3.0

Prepare binaries for publishing:
- python3 make_release.py
or 'Ctrl+F5' from VSCode """

import os
import shutil
import datetime
import PyInstaller.__main__

def copy(src, dst, skip):
    for f in os.listdir(src):
        if f!='dist' and not f.endswith(skip):
            src_path = os.path.join(src, f)
            dst_path = os.path.join('dist', src, f)

            if os.path.isdir(src_path):
                if not os.path.isdir(dst_path):
                    os.mkdir(dst_path)
                copy(src_path, dst_path, skip)

            if os.path.isfile(src_path):
                shutil.copy2(src_path, dst_path)

if __name__ == '__main__':
    if os.name=='nt':
        op_sys = '_windows'
        skip = ('_linux', '_linux.env', '.sh', '.desktop',
            'ccx', 'cgx', 'unv2ccx', 'ccx2paraview')
        extension = '.exe' # binary extension in OS
        TEMP = 'C:\\Windows\\Temp\\'
    else:
        op_sys = '_linux'
        skip = ('_windows', '_windows.env', '.bat', '.exe', '.dll')
        extension = '' # binary extension in OS
        TEMP = '/tmp/'

    PROJECT_NAME = os.path.split(os.getcwd())[-1] # name of project's folder
    DATE = '_' + datetime.datetime.now().strftime('%Y%m%d')
    ARCH = os.path.join('./releases', PROJECT_NAME + DATE + op_sys)

    # Remove prev. trash
    if os.path.isdir('./dist'):
        shutil.rmtree('./dist')

    # Run pyinstaller to create binaries
    args = [
        './src/cae.py',
        '--workpath=' + TEMP,   # temp dir
        '-w',                   # no console during app run
        ]
    PyInstaller.__main__.run(args)

    # Delete cached files
    if os.path.isdir('./src/__pycache__'):
        shutil.rmtree('./src/__pycache__') # works in Linux as in Windows

    # Delete .spec file
    if os.path.isfile('cae.spec'):
        os.remove('cae.spec')

    # Rename ./dist/cae to ./dist/src
    shutil.move('./dist/cae', './dist/src')

    # Prepare skip list
    with open('.gitignore', 'r') as f:
        lines = f.readlines()
    for i in range(len(lines)):
        skip += (lines[i].rstrip().lstrip('*'), )
    skip += ('.git', '.gitignore', '.py',
        'tests.log', 'requirements.txt',
        'dist', 'gui', 'model')

    # Copy files and folders from sources to 'dist'
    copy('.', 'dist', skip)

    # Make archive
    if os.path.isfile(ARCH + '.zip'):
        os.remove(ARCH + '.zip') # delete old

    # Complress whole directory
    shutil.make_archive(ARCH, 'zip', 'dist')

    # Remove unneeded files and folders
    shutil.rmtree(TEMP + 'cae')
    shutil.rmtree(os.path.abspath('dist'))
