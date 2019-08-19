# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    Prepare binaries for publishing.
"""


import os, shutil, datetime
import PyInstaller.__main__


if __name__ == '__main__':

    DIRECTORY = './dist/ccx_cae/'
    DATE = datetime.datetime.now().strftime('%Y%m%d')
    op_sys = os.system('windows' if os.name=='nt' else 'ubuntu')
    ARCH = 'ccx_cae_' + op_sys + '_' + DATE
    TEMP = os.system('C:\Windows\Temp\' if os.name=='nt' else '/tmp/')

    # Run pyinstaller to create binaries
    PyInstaller.__main__.run(['ccx_cae.py', '--workpath=' + TEMP])

    # Copy some files and folders from sources
    copy_files = ['doc', 'icons', 'examples',
                'ccx_cae.ui', 'ccx_dialog.ui',
                'ccx_dom.inp', 'ccx_mesh.inp', 'Materials.inp',
                'Settings.env']
    for f in copy_files:
        os.system('cp -r ' + f + ' ' + DIRECTORY + f)

    # Make archive
    os.system('rm ' + ARCH + '.zip') # delete old
    shutil.make_archive(ARCH, 'zip', DIRECTORY) # create new

    # Remove unneeded files and folders
    os.system('rm -r ' + TEMP + 'ccx_cae ./dist')
    os.system('rm ccx_cae.spec')

    # Delete cached files
    if os.path.isdir('__pycache__'):
        shutil.rmtree('__pycache__') # works in Linux as in Windows
