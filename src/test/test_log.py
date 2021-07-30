#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Test for src/log.py
"""

# Standard modules
import os
import sys
import logging
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
from gui.window import MasterWindow


app = None
window = None


class qAppTest(unittest.TestCase):
    """Helper class to provide QApplication instances."""

    def setUp(self):
        """Creates the QApplication instance."""
        global app
        if app is None:
            app = QtWidgets.QApplication([])

        global window
        if window is None:
            window = MasterWindow()

        # # Configure logging level
        # for h in logging.getLogger().handlers:
        #     h.setLevel(logging.NOTSET)
        #     print(h.name, logging.getLevelName(h.level))


class Test(qAppTest):
    """Test if logging works."""

    def test0_window(self):
        global window
        self.assertTrue(hasattr(window, 'textEdit'))

    def test1_debug(self):
        global window
        logging.debug('debug')
        txt = window.textEdit.toPlainText()
        self.assertTrue('debug' in txt)

    def test2_info(self):
        global window
        logging.info('info')
        txt = window.textEdit.toPlainText()
        self.assertTrue('info' in txt)

    def test3_warning(self):
        global window
        logging.warning('warning')
        txt = window.textEdit.toPlainText()
        self.assertTrue('warning' in txt)

    def test4_error(self):
        global window
        logging.error('error')
        txt = window.textEdit.toPlainText()
        self.assertTrue('error' in txt)


if __name__ == '__main__':
    unittest.main()
