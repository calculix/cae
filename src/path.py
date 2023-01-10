#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2023
Distributed under GNU General Public License v3.0

Absolute paths to the application main folders.
Static class whos fields shouldn't change on the run.
It is initialized on the app start - and that's all.
"""

import os
import platform


class Path:

    def __init__(self):
        self.op_sys = platform.system().lower() # OS name

        # Windows
        if os.name == 'nt':
            self.extension = '.exe' # file extension in OS

        # Linux
        else:
            self.extension = '' # file extension in OS

        # Application home directory - the one with README.md and LICENSE
        self.app_home_dir = os.path.normpath(
            os.path.join(os.path.dirname(
                os.path.realpath(__file__)), '..'))

        self.config = os.path.join(self.app_home_dir, 'config')
        self.main_xml = os.path.join(self.config, 'Window.xml')
        self.kw_list_xml = os.path.join(self.config, 'kw_list.xml')
        self.kw_tree_xml = os.path.join(self.config, 'kw_tree.xml')
        self.settings_xml = os.path.join(self.config, 'SettingsDialog.xml')
        self.dialog_xml = os.path.join(self.config, 'KeywordDialog.xml')
        self.version_txt = os.path.join(self.config, 'version.txt')

        self.requirements = os.path.join(self.config, 'requirements_' + self.op_sys + '.txt')
        self.bin = os.path.join(self.app_home_dir, 'bin')
        self.ccx = os.path.join(self.app_home_dir, 'ccx_' + self.op_sys, 'src')
        self.settings = os.path.join(self.config, 'Settings_' + self.op_sys + '.py')
        self.doc = os.path.join(self.app_home_dir, 'doc')
        self.examples = os.path.join(self.app_home_dir, 'examples')
        self.img = os.path.join(self.app_home_dir, 'img')
        self.src = os.path.join(self.app_home_dir, 'src')
        self.path_ccx = os.path.join(self.bin, 'ccx' + self.extension)
        self.path_cgx = os.path.join(self.bin, 'cgx' + self.extension)
        self.log = os.path.join(self.app_home_dir, 'cae.log')

    def append_to_PATH(self, paths):
        """Pyinstaller bug in Windows:
        append 'app_home_dir' and 'src' directories to PATH.
        """
        if not os.environ['PATH'].endswith(os.pathsep):
            os.environ['PATH'] += os.pathsep
        for path in paths:
            if path not in os.environ['PATH']:
                os.environ['PATH'] += path
                os.environ['PATH'] += os.pathsep

    def abspath(self, rel):
        """Convert relative path to absolute and check."""

        # We do not know if rel is really relative path
        if os.path.isfile(os.path.join(self.app_home_dir, rel)):
            return os.path.join(self.app_home_dir, rel)

        # If rel is absolute path - return as is
        else:
            return rel


p = Path()
