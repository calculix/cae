#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2023
Distributed under GNU General Public License v3.0

Test for src/gui/dialog.py
"""

# Standard modules
import os
import sys
import unittest

# External modules
from PyQt5 import QtWidgets, QtWebEngineWidgets

# My modules
sys_path = os.path.abspath(__file__)
sys_path = os.path.dirname(sys_path)
sys_path = os.path.join(sys_path, '..')
sys_path = os.path.normpath(sys_path)
sys_path = os.path.realpath(sys_path)
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
from model.kom import KOM
from gui.dialog import KeywordDialog


class TestDialog(unittest.TestCase):

    def test_dialog(self):
        """Create keyword dialog"""
        app = QtWidgets.QApplication(sys.argv)
        i = KOM.get_keyword_by_name('*NODE')
        try:
            d = KeywordDialog(i)
        except:
            d = None
        self.assertTrue(d is not None)


if __name__ == '__main__':
    unittest.main()
