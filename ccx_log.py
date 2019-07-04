# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    Logging methods
"""


from PyQt5 import QtGui


class logger:


    def __init__(self, CAE):
        self.CAE = CAE


    # Info log with Black font color
    def info(self, msg):
        self.CAE.textEdit.append('<p style=\'color:Black; margin:0px;\'>' + msg + '</p>')
        self.CAE.textEdit.moveCursor(QtGui.QTextCursor.End) # scroll text to the end


    # Error log with Red font color
    def error(self, msg):
        self.CAE.textEdit.append('<p style=\'color:Red; margin:0px;\'>ERROR! ' + msg + '</p>')
        self.CAE.textEdit.moveCursor(QtGui.QTextCursor.End) # scroll text to the end

