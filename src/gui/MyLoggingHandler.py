# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, September 2019
    Distributed under GNU General Public License v3.0

    Custon logging handler + log utilities
"""


import logging, re
from PyQt5 import QtGui


# Logging handler
class MyLoggingHandler(logging.Handler):


    # Initialization
    def __init__(self, textEdit):
        super().__init__() # create handler
        self.textEdit = textEdit
        fmt = logging.Formatter('%(levelname)s, %(module)s: %(message)s')
        self.setFormatter(fmt)


    # Sends log messages to MainWindow's textEdit widget
    def emit(self, LogRecord):

        # Message color depending on logging level
        color = {
                'NOTSET':   'Green',
                'DEBUG':    'Gray',
                'INFO':     'Black',
                'WARNING':  'Blue',
                'ERROR':    'Red',
                'CRITICAL': 'Pink',
            }

        if LogRecord.levelname in color:
            msg_text = self.format(LogRecord)
            color = color[LogRecord.levelname]
        else:
            msg_text = LogRecord.getMessage()
            color = 'Brown'

        self.textEdit.append('<p style=\'color:{0}; margin:0px;\'>{1}</p>'.format(color, msg_text))
        self.textEdit.moveCursor(QtGui.QTextCursor.End) # scroll text to the end


# Process one stdout message
def logLine(line):
    logging_level = {
            'NOTSET': 0,
            'DEBUG': 10,
            'INFO': 20,
            'WARNING': 30,
            'ERROR': 40,
            'CRITICAL': 50,
        }
    if line.startswith(tuple(logging_level.keys())):
        match = re.search('(^\w+): (.+)', line) # levelname and message
        if match: # skip logging of empty strings
            level = match.group(1)
            msg_text = match.group(2)
            logging.log(logging_level[level], msg_text)
    else:
        logging.log(25, line)
