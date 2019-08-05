# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    Prepare source code for publising.
"""


import os, datetime


if __name__ == '__main__':

    DIRECTORY = os.getcwd()
    SKIP = ['examples', '.git', 'test.py', 'test.txt',
            'release_binaries.py', os.path.basename(sys.executable)]
    DATE = datetime.datetime.now().strftime('%Y%m%d')
    ARCH = '../ccx_cae_sources_' + DATE + '.zip'

    files = ' '
    for f in os.listdir(DIRECTORY): # iterate over files in DIRECTORY
        if not f in SKIP:
            files = files + ' ' + f

    os.system('rm ' + ARCH) # delete old source archive
    os.system('python -m zipfile -c ' + ARCH + files) # compress command
