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
from PyQt5 import QtWidgets, QtWebEngineWidgets, QtCore

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
    """Open KeywordDialog, click widgets to test them,
    repeat the same for the next keyword (if clicked OK),
    or quit (if clicked Cancel).
    """
    app = QtWidgets.QApplication(sys.argv)
    for i, k in enumerate(KWL.keywords):
        print(i, k.name)
        d = KeywordDialog(k)
        d.make_screenshot = True
        if not d.exec(): # 0 = cancel, 1 = ok
            break


class DialogPngCreator():
    """Export KeywordDialog renders.
    It allows to check what will look like dialogs for all the keywords.
    """

    def __init__(self):
        app = QtWidgets.QApplication(sys.argv)
        self.counter = 0
        self.start_index = 0
        self.end_index = 0 # len(KWL.keywords) - 1
        self.stack = []
        self.timer1 = QtCore.QTimer()
        self.timer2 = QtCore.QTimer()
        self.timer1.timeout.connect(self.open_dialogs)
        self.timer2.timeout.connect(self.save_dialogs)
        self.timer1.start(300)
        app.exec_()

    def open_dialogs(self):
        k = KWL.keywords[self.start_index + self.counter]
        d = KeywordDialog(k)
        print(self.start_index + self.counter, d.item.name)
        self.stack.append(d)
        self.counter += 1
        if self.start_index + self.counter > self.end_index:
            QtCore.QTimer.singleShot(1000, self.start_timer2)
            self.timer1.stop()
            print('timer1 stopped')

    def start_timer2(self):
        print('wait')
        self.timer2.start(1000)

    def save_dialogs(self):
        d = self.stack.pop()
        self.counter -= 1
        print(self.start_index + self.counter, d.item.name)
        d.screenshot()
        d.close()
        if self.counter <= 0:
            self.timer2.stop()
            print('timer2 stopped')


if __name__ == '__main__':
    # unittest.main()
    # cycle_keyword_dialogs()
    # DialogPngCreator()
    pass
