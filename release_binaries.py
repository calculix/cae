# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    Prepare binaries for publising.
"""


import os, shutil, datetime
import PyInstaller.__main__


if __name__ == '__main__':

    DIRECTORY = './dist/ccx_cae/'
    DATE = datetime.datetime.now().strftime('%Y%m%d')
    ARCH = 'ccx_cae_binaries_' + DATE

    # Run pyinstaller to create binaries
    PyInstaller.__main__.run([
        'ccx_cae.py',
        '--workpath=/tmp',
    ])

    # Copy some files and folders from sources
    copy_files = ['doc', 'icons', 'ccx_cae.ui', 'ccx_dialog.ui', 'ccx_dom.inp', 'ccx_mesh.inp']
    for f in copy_files:
        os.system('cp -r ' + f + ' ' + DIRECTORY + f)

    # Make archive
    os.system('rm ' + ARCH + '.zip') # delete old
    shutil.make_archive(ARCH, 'zip', DIRECTORY) # create new

    # Remove unneeded files and folders
    os.system('rm -r /tmp/ccx_cae ./dist')
    os.system('rm ccx_cae.spec')
    os.system('py3clean .') # clear cache
