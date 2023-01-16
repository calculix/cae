#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2023
Distributed under GNU General Public License v3.0

Test for src/path.py
"""

# Standard modules
import os
import sys
import unittest

# My modules
sys_path = os.path.abspath(__file__)
sys_path = os.path.dirname(sys_path)
sys_path = os.path.join(sys_path, '..')
sys_path = os.path.normpath(sys_path)
sys_path = os.path.realpath(sys_path)
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
from path import p


class TestPath(unittest.TestCase):

    def test_path(self):
        """Ensure that all app's folders and files are present"""
        skip = ('op_sys', 'extension', '__')
        for attr in dir(p):
            if attr.startswith(skip):
                continue
            a = getattr(p, attr)
            if type(a) is not str:
                continue
            if 'Settings_' in a:
                continue
            status = os.path.isfile(a) or os.path.isdir(a)
            if not status:
                print('\nNot found:', a)
            self.assertTrue(status)


if __name__ == '__main__':
    unittest.main()
