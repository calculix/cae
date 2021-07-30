#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Test for src/settings.py
"""

# Standard modules
import os
import sys
import unittest

# External modules
from PyQt5 import QtWidgets

# My modules
sys_path = os.path.abspath(__file__)
sys_path = os.path.dirname(sys_path)
sys_path = os.path.join(sys_path, '..')
sys_path = os.path.normpath(sys_path)
sys_path = os.path.realpath(sys_path)
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
from settings import s


class Test(unittest.TestCase):

    # Create and open settings window
    def test_settings(self):
        app = QtWidgets.QApplication(sys.argv)
        s.open()
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
