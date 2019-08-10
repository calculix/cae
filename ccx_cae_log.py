# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, August 2019.
    Distributed under GNU General Public License, version 2.

    Logging handler.
    Sends log messages to CAE's textEdit widget.
    Each logging level has own message color.
"""


import logging
from PyQt5 import QtGui


class myHandler(logging.Handler):

    def __init__(self, CAE):
        super().__init__()
        self.textEdit = CAE.textEdit
        self.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))

    def emit(self, LogRecord):
        msg_text = self.format(LogRecord)

        # Message color depending on logging level
        color = {
                'DEBUG':'Gray',
                'INFO':'Black',
                'WARNING':'Blue',
                'ERROR':'Red',
            }[LogRecord.levelname]

        self.textEdit.append('<p style=\'color:{0}; margin:0px;\'>{1}</p>'.format(color, msg_text))
        self.textEdit.moveCursor(QtGui.QTextCursor.End) # scroll text to the end
