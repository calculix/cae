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
from log import stop_logging
from gui.window import wf


# Create application
app = QtWidgets.QApplication(sys.argv)

stop_logging()

# Show main window and add_text_handler()
wf.run_master()
for h in logging.getLogger().handlers:
    print(h.name)


class Test(unittest.TestCase):

    def test_logging_debug(self):
        'Test if logging.debug works'
        logging.debug('zxc')
        txt = wf.mw.textEdit.toPlainText()
        self.assertTrue('zxc' in txt)

    def test_logging_info(self):
        'Test if logging.info works'
        logging.info('qwe')
        txt = wf.mw.textEdit.toPlainText()
        self.assertTrue('qwe' in txt)

    def test_logging_warning(self):
        'Test if logging.warning works'
        logging.warning('rty')
        txt = wf.mw.textEdit.toPlainText()
        self.assertTrue('rty' in txt)

    def test_logging_error(self):
        'Test if logging.error works'
        logging.error('asd')
        txt = wf.mw.textEdit.toPlainText()
        self.assertTrue('asd' in txt)


if __name__ == '__main__':
    unittest.main()
