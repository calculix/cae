#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2023
Distributed under GNU General Public License v3.0

Application settings.
Attributes values are maintained in config/Settings_*.py.
User dialog form is config/SettingsDialog.xml - use Qt Designer to edit.
"""

# Standard modules
import os

# My modules
from path import p


class Settings:
    """Session settings object is used everywhere in the code."""

    def __init__(self):
        """Read settings from file or apply defaults."""
        try:
            # Try to read settings file
            with open(p.settings, 'r') as f:
                exec(f.read())
            if len(self.__dict__) < 2:
                raise Exception

        except:
            # Apply default values
            if os.name=='nt':
                # Windows
                self.path_paraview = 'C:\\Program Files\\ParaView\\bin\\paraview.exe'
                self.path_editor = 'C:\\Windows\\System32\\notepad.exe'
            else:
                # Linux
                self.path_paraview = '/usr/bin/paraview'
                self.path_editor = '/usr/bin/gedit'

            self.start_model = os.path.join(p.examples, 'default.inp')
            self.default_web_browser = 'internal'
            self.logging_level = 'DEBUG'
            self.show_empty_keywords = True
            self.expanded = True
            self.start_cgx_by_default = True
            self.align_windows = False
            self.show_help = False
            self.perform_startup_checks = True


s = Settings()
