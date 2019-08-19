# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    Prepare CCX_CAE binaries for publishing.
"""


import os, shutil, datetime, subprocess
import PyInstaller.__main__


if __name__ == '__main__':
    PROJECT_NAME = 'ccx_cae'
    DIRECTORY = './dist/' + PROJECT_NAME + '/'
    DATE = '_' + datetime.datetime.now().strftime('%Y%m%d')
    op_sys = '_windows' if os.name=='nt' else '_ubuntu'
    ARCH = '../' + PROJECT_NAME + op_sys + DATE
    TEMP = 'C:\\Windows\\Temp\\' if os.name=='nt' else '/tmp/'

    # Run pyinstaller to create binaries
    PyInstaller.__main__.run([PROJECT_NAME + '.py', '--workpath=' + TEMP])

    # Delete cached files
    if os.path.isdir('__pycache__'):
        shutil.rmtree('__pycache__') # works in Linux as in Windows

    # Delete .spec file
    if os.path.isfile(PROJECT_NAME + '.spec'):
        subprocess.run('rm ' + PROJECT_NAME + '.spec', shell=True)

    # Copy some files and folders from sources
    copy_files = []
    skip_files = ['dist', '.git']
    for file_name in os.listdir():
        # Append all dirs and files except Python sources
        if not file_name in skip_files and \
           not file_name.endswith('.py'):
            copy_files.append(file_name)
    for f in copy_files:
        subprocess.run('cp -r ' + f + ' ' + DIRECTORY + f, shell=True)

    # Make archive
    if os.path.isfile(ARCH + '.zip'):
        subprocess.run('rm ' + ARCH + '.zip', shell=True) # delete old
    shutil.make_archive(ARCH, 'zip', DIRECTORY) # create new

    # Remove unneeded files and folders
    subprocess.run('rm -r ' + TEMP + PROJECT_NAME + ' ./dist', shell=True)
