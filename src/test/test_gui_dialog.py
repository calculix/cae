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
from model.kom import KWL
from gui.dialog import KeywordDialog


class TestDialog(unittest.TestCase):

    def test_dialog(self):
        """Create keyword dialog"""
        app = QtWidgets.QApplication(sys.argv)
        try:
            item = KWL.get_keyword_by_name('*AMPLITUDE')
            d = KeywordDialog(item)
            # from gui.window import df
            # df.run_master_dialog(item) # 0 = cancel, 1 = ok
        except:
            d = None
        self.assertTrue(d is not None)


def cycle_keyword_dialogs():
    from gui.window import df
    app = QtWidgets.QApplication(sys.argv)

    from model.kom import KWL
    for i, k in enumerate(KWL.keywords):
        print(i, k.name)

        # Create keyword dialog
        if not df.run_master_dialog(k): # 0 = cancel, 1 = ok
            break

if __name__ == '__main__':
    # unittest.main()
    cycle_keyword_dialogs()
