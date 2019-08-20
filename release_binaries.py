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

    one_file = False
    PROJECT_NAME = os.path.split(os.getcwd())[-1] # name of project's folder

    DIRECTORY = os.path.join(os.path.abspath('dist'), PROJECT_NAME)
    DATE = '_' + datetime.datetime.now().strftime('%Y%m%d')
    op_sys = '_windows' if os.name=='nt' else '_ubuntu'
    ARCH = os.path.join('..', PROJECT_NAME + op_sys + DATE)
    TEMP = 'C:\\Windows\\Temp\\' if os.name=='nt' else '/tmp/'
    extension = ('.exe' if os.name=='nt' else '') # file extension in OS

    # Run pyinstaller to create binaries
    args = [PROJECT_NAME + '.py', '--workpath=' + TEMP]
    if one_file:
        args.append('--onefile')
    PyInstaller.__main__.run(args)

    # Delete cached files
    if os.path.isdir('__pycache__'):
        shutil.rmtree('__pycache__') # works in Linux as in Windows

    # Delete .spec file
    if os.path.isfile(PROJECT_NAME + '.spec'):
        os.remove(PROJECT_NAME + '.spec')

    # Copy some files and folders from sources
    if not one_file:
        skip_files = ['tests', 'dist', '.git']
        for f in os.listdir():
            # All dirs and files except Python sources
            if not f in skip_files and not f.lower().endswith('.py'):
                if os.path.isdir(f):
                    shutil.copytree(f, os.path.join(DIRECTORY, f))
                else:
                    shutil.copy2(f, os.path.join(DIRECTORY, f))

    # Make archive
    if os.path.isfile(ARCH + '.zip'):
        os.remove(ARCH + '.zip') # delete old
    if not one_file:
        # Complress whole directory
        shutil.make_archive(ARCH, 'zip', DIRECTORY)
    else:
        # Complress single file
        import zipfile
        with zipfile.ZipFile(ARCH + '.zip', 'w') as archive:
            archive.write(DIRECTORY + extension,
                arcname=os.path.basename(DIRECTORY) + extension)

    # Remove unneeded files and folders
    shutil.rmtree(TEMP + PROJECT_NAME)
    shutil.rmtree(os.path.abspath('dist'))
