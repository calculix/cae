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
# import log
from gui.window import MasterWindow


# Show main window
# log.stop_logging()
app = QtWidgets.QApplication([])
mw = MasterWindow()

# # Configure logging level
# for h in logging.getLogger().handlers:
#     h.setLevel(logging.NOTSET)
#     print(h.name, logging.getLevelName(h.level))

logging.debug('debug')
logging.info('info')
logging.warning('warning')
logging.error('error')

txt = mw.textEdit.toPlainText()
print(txt)
# app.exec()


class Test(unittest.TestCase):
    """Test if logging works."""

    def test1_debug(self):
        global txt
        self.assertTrue('debug' in txt)

    def test2_info(self):
        global txt
        self.assertTrue('info' in txt)

    def test3_warning(self):
        global txt
        self.assertTrue('warning' in txt)

    def test4_error(self):
        global txt
        self.assertTrue('error' in txt)


if __name__ == '__main__':
    unittest.main()
