#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, September 2019
Distributed under GNU General Public License v3.0

Custon logging handler.
The Text handler is created once per session.
The File handler is created any time you open a model. """

import os
import re
import logging
from PyQt5 import QtGui


class MyLoggingHandler(logging.Handler):

    # Target could be textEdit or file
    def __init__(self, target):
        super().__init__() # create handler
        self.target = target
        fmt = logging.Formatter('%(levelname)s, %(module)s: %(message)s')
        self.setFormatter(fmt)


# Handler to show logs in the CAE's textEdit
class MyTextLoggingHandler(MyLoggingHandler):

    # Sends log messages to Window's textEdit widget
    def emit(self, LogRecord):

        # Message color depending on logging level
        color = {
            'NOTSET':   'Black',
            'DEBUG':    'Gray',
            'INFO':     'Green',
            'WARNING':  'Blue',
            'ERROR':    'Red',
            'CRITICAL': 'Pink', }
        if LogRecord.levelname in color:
            color = color[LogRecord.levelname]
        else:
            color = 'Brown'

        # Move newlines before the levelname
        msg_text = LogRecord.getMessage()
        while msg_text.startswith('\n'):
            self.log_msg('<p></p>')
            msg_text = msg_text[1:]

        # Keep all newlines in message
        msg_text = msg_text.replace('\n', '<br/>')

        msg_text = '<p style=\'margin:0px;\'><span style=\'color:{};\'>{}</span>: {}</p>'\
            .format(color, LogRecord.levelname, msg_text)
        self.log_msg(msg_text)
    
    def log_msg(self, msg_text):
        self.target.append(msg_text)
        self.target.moveCursor(QtGui.QTextCursor.End) # scroll text to the end


# Handler to write the job's log file
class MyFileLoggingHandler(MyLoggingHandler):

    def emit(self, LogRecord):

        # Move newlines before the levelname
        msg_text = LogRecord.getMessage()
        while msg_text.startswith('\n'):
            self.log_msg('')
            msg_text = msg_text[1:]

        msg_text = '{}: {}'.format(LogRecord.levelname, msg_text)        
        self.log_msg(msg_text)

    def log_msg(self, msg_text):
        with open(self.target, 'a') as f:
            f.write(msg_text + '\n')


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


def add_text_handler(textEdit):
    logging.getLogger().handlers = []
    h = MyTextLoggingHandler(textEdit)
    logging.getLogger().addHandler(h)


def add_file_handler(log_file):
    if len(logging.getLogger().handlers) > 1:
        logging.getLogger().handlers.pop()
    h = MyFileLoggingHandler(log_file)
    logging.getLogger().addHandler(h)
    if os.path.exists(log_file):
        os.remove(log_file)
