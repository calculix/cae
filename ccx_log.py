# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, June 2019.
    Distributed under GNU General Public License, version 2.

    Logging methods
"""


from PyQt5 import QtGui


class logger:


    def __init__(self, CAE):
        # Here we shall output logs
        self.textEdit = CAE.textEdit


    def info(self, msg):
        self.textEdit.append('<p style=\'color:Black; margin:0px;\'>' + msg + '</p>')
        self.textEdit.moveCursor(QtGui.QTextCursor.End) # scroll text to the end


    def error(self, msg):
        self.textEdit.append('<p style=\'color:Red; margin:0px;\'>ERROR! ' + msg + '</p>')
        self.textEdit.moveCursor(QtGui.QTextCursor.End) # scroll text to the end

