# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    Logging methods.
    Depends on ccx_cae.CAE.
"""


from PyQt5 import QtGui
from enum import Enum


# Enums for 'msg_type' variable
class msgType(Enum):
    INFO = 0
    WARNING = 1
    ERROR = 2


# Message for logger
class msg:
    
    def __init__(self, msg_type, msg_text):
        self.msg_type = msg_type
        self.msg_text = msg_text


class logger:

    def __init__(self, CAE=None):
        self.CAE = CAE


    # Short syntax for all messages
    def info(self, msg_text):
        self.message(msg(msgType.INFO, msg_text))
    def warning(self, msg_text):
        self.message(msg(msgType.WARNING, msg_text))
    def error(self, msg_text):
        self.message(msg(msgType.ERROR, msg_text))


    # Process one message
    def message(self, msg):

        # Message color depending on status
        color = {msgType.INFO:'Black', msgType.WARNING:'Blue', msgType.ERROR:'Red'}[msg.msg_type]

        # Message marker depending on status
        status = {msgType.INFO:'', msgType.WARNING:'WARNING: ', msgType.ERROR:'ERROR: '}[msg.msg_type]

        print('{0}{1}'.format(status, msg.msg_text))
        if self.CAE: # could be None in tests
            self.CAE.textEdit.append('<p style=\'color:{0}; margin:0px;\'>{1}{2}</p>'.format(color, status, msg.msg_text))
            self.CAE.textEdit.moveCursor(QtGui.QTextCursor.End) # scroll text to the end


    # Process list of messages
    def messages(self, msg_list):
        for msg in msg_list:
            self.message(msg)
