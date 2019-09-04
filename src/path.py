# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, September 2019
    Distributed under GNU General Public License v3.0

    Utility to calculate absolute pathes to the application's main folders.
"""


import os, sys


class Path:


    def __init__(self):
        op_sys = 'windows' if os.name=='nt' else 'linux' # OS name

        # Application's home directory - the one with README.md and LICENSE
        self.app_home_dir = os.path.abspath(
                os.path.join(os.path.dirname(sys.argv[0]), '..'))

        self.config = os.path.join(self.app_home_dir, 'config')
        self.cae_xml = os.path.join(self.config, 'cae.xml')
        self.kom_xml = os.path.join(self.config, 'kom.xml')
        self.settings_xml = os.path.join(self.config, 'settings.xml')
        self.dialog_xml = os.path.join(self.config, 'dialog.xml')

        self.bin = os.path.join(self.app_home_dir, 'bin')
        self.ccx = os.path.join(self.app_home_dir, 'ccx_' + op_sys, 'ccx_free_form_fortran')
        self.settings = os.path.join(self.config, 'settings_' + op_sys + '.env')
        self.doc = os.path.join(self.app_home_dir, 'doc')
        self.examples = os.path.join(self.app_home_dir, 'examples')
        self.img = os.path.join(self.app_home_dir, 'img')
        self.src = os.path.join(self.app_home_dir, 'src')


    # Pyinstaller bug in Windows: append 'app_home_dir' and 'src' directories to PATH
    def append_to_PATH(self, pathes):
        if not os.environ['PATH'].endswith(os.pathsep):
            os.environ['PATH'] += os.pathsep
        for path in pathes:
            if path not in os.environ['PATH']:
                os.environ['PATH'] += path
                os.environ['PATH'] += os.pathsep


    # Convert relative path to absolute and check
    def abspath(self, rel):

        # We do not know if rel is really relative path
        if os.path.isfile(os.path.join(self.app_home_dir, rel)):
            return os.path.join(self.app_home_dir, rel)

        # If rel is absolute path - return as is
        else:
            return rel