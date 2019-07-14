# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    Logging methods
"""


from PyQt5 import QtGui
from enum import Enum


# Enums for 'msg_type' variable
class msgType(Enum):
    INFO = 0
    ERROR = 1


# Message for logger
class msg:
    
    def __init__(self, msg_type, msg_text):
        self.msg_type = msg_type
        self.msg_text = msg_text


class logger:

    def __init__(self, CAE=None):
        self.CAE = CAE
    

    # Process one message
    def message(self, msg):
        if msg.msg_type == msgType.INFO:
            self.info(msg.msg_text)
        if msg.msg_type == msgType.ERROR:
            self.error(msg.msg_text)


    # Process list of messages
    def messages(self, msg_list):
        for msg in msg_list:
            self.message(msg)


    # Info log with Black font color
    def info(self, msg):
        print(msg)
        if self.CAE: # could be None in tests
            self.CAE.textEdit.append('<p style=\'color:Black; margin:0px;\'>' + msg + '</p>')
            self.CAE.textEdit.moveCursor(QtGui.QTextCursor.End) # scroll text to the end


    # Error log with Red font color
    def error(self, msg):
        print('ERROR!', msg)
        if self.CAE: # could be None in tests
            self.CAE.textEdit.append('<p style=\'color:Red; margin:0px;\'>ERROR! ' + msg + '</p>')
            self.CAE.textEdit.moveCursor(QtGui.QTextCursor.End) # scroll text to the end
