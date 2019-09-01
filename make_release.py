# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    Prepare binaries for publishing:
        python3 release_binaries.py
"""


import os, shutil, datetime, subprocess
import PyInstaller.__main__


if __name__ == '__main__':

    if os.name=='nt':
        op_sys = '_windows'
        skip_files = ('_linux', '.sh', '.desktop')
        extension = '.exe' # binary extension in OS
        TEMP = 'C:\\Windows\\Temp\\'
    else:
        op_sys = '_linux'
        skip_files = ('_windows', '.bat', '.exe', '.dll')
        extension = '' # binary extension in OS
        TEMP = '/tmp/'

    PROJECT_NAME = os.path.split(os.getcwd())[-1] # name of project's folder
    DIRECTORY = os.path.join(os.path.abspath('dist'), 'cae')
    DATE = '_' + datetime.datetime.now().strftime('%Y%m%d')
    ARCH = os.path.join('..', PROJECT_NAME + DATE + op_sys)

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

    # Copy some files and folders from sources
    skip_files += ('dist', 'src', '.py', '.git', '.vscode',
            '.gitignore')
    for f in os.listdir():
        # All dirs and files except Python sources
        if not f.endswith(skip_files):
            if os.path.isdir(f):
                shutil.copytree(f, os.path.join('dist', f),
                    ignore=shutil.ignore_patterns('*.py'))
            else:
                shutil.copy2(f, os.path.join('dist', f))

    # Make archive
    if os.path.isfile(ARCH + '.zip'):
        os.remove(ARCH + '.zip') # delete old

    # Complress whole directory
    shutil.make_archive(ARCH, 'zip', 'dist')

    # Remove unneeded files and folders
    shutil.rmtree(TEMP + 'cae')
    shutil.rmtree(os.path.abspath('dist'))
