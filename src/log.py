#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Custom logging handlers.
The Text handler is created once per session and prints into CAE textEdit.
The File handler is created any time you open a model.
Log file has the same name as a model.
"""

# Standard modules
import os
import logging

# External modules
from PyQt5 import QtGui, QtWidgets

# My modules
from settings import s
from utils import tests


mh = 'MyHandler'
mtlh = 'MyTextLoggingHandler'
mflh = 'MyFileLoggingHandler'
mslh = 'MyStreamLoggingHandler'
fmt = logging.Formatter('%(levelname)s: %(message)s')


# def initialize_logging(level=s.logging_level):
#     global fmt
#     logging.basicConfig(level=level, format=fmt._fmt)


def get_logging_info():
    hh = logging.getLogger().handlers
    if not len(hh):
        msg = 'No logging handlers.'
    else:
        msg = 'Total {} logging handlers:'.format(len(hh))
        for h in hh:
            msg += '\n' + str(h.get_name())
    print(msg)


def print_to_file(log_file, *args):
    """Redefine print method to write logs to file and stdout."""
    line = ' '.join([str(arg) for arg in args])
    if log_file is not None:
        with open(log_file, 'a') as f:
            f.write(line + '\n')
    print(line)


class myHandler(logging.Handler):
    """Handler to emit messages via print_to_file()."""

    def __init__(self, log_file):
        super().__init__()
        global fmt
        self.setFormatter(fmt)
        self.log_file = log_file

        # Remove old log file
        if os.path.isfile(log_file):
            os.remove(log_file)

    def emit(self, LogRecord):
        print_to_file(self.log_file, self.format(LogRecord))


class MyLoggingHandler(logging.Handler):
    """Target could be textEdit or file."""

    def __init__(self, target):
        super().__init__() # create handler
        self.target = target
        global fmt
        self.setFormatter(fmt)


class MyTextLoggingHandler(MyLoggingHandler):
    """Called only once on startup.
    Handler to show logs in master window textEdit.
    """

    def emit(self, LogRecord):
        """Sends log messages to Window textEdit widget."""

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
            color = '#2F4F4F' # DARKSLATEGRAY

        msg = LogRecord.getMessage()

        # Move leading and trailing newlines outside the message
        leading_newlines = ''
        trailing_newlines = ''
        while msg.startswith('\n'):
            leading_newlines += '<br/>'
            msg = msg[1:]
        while msg.endswith('\n'):
            trailing_newlines += '<br/>'
            msg = msg[:-1]

        # Keep all newlines inside the message
        msg = msg.replace('\n', '<br/>')

        if LogRecord.levelname == 'Level 25':
            msg = '<p style=\"margin:0px; color:{};\">{}{}{}</p>'\
                .format(color, leading_newlines, msg, trailing_newlines)
        else:
            msg = '{}<span style=\"color:{};\">{}</span>: {}{}'\
                .format(leading_newlines, color, LogRecord.levelname, msg, trailing_newlines)
            msg = '<p style=\"margin:0px;\">{}</p>'.format(msg)

        # Print message and scroll it to the end
        self.target.append(msg)
        self.target.moveCursor(QtGui.QTextCursor.End)


class MyFileLoggingHandler(MyLoggingHandler):
    """Called at each file import.
    Handler to write logs into file.
    """

    def emit(self, LogRecord):

        # Move newlines before the levelname
        msg = LogRecord.getMessage()
        # print(msg + '\n')
        while msg.startswith('\n'):
            with open(self.target, 'a') as f:
                f.write('\n')
            msg = msg[1:]

        msg = '{}: {}'.format(LogRecord.levelname, msg)
        if msg.startswith('Level 25: '):
            msg = msg[10:]


        # Append message to the log file
        with open(self.target, 'a') as f:
            f.write(msg + '\n')


def stop_logging():
    """Stop all logging handlers, even default stdout."""
    for h in logging.getLogger().handlers:
        h.close()
    logging.getLogger().handlers = []
    logging.disable() # switch off logging


def remove_handler_by_name(name):
    hh = logging.getLogger().handlers
    for h in hh:
        if h.name == name:
            h.close()
            hh.remove(h)
            break


def add_my_handler(level=s.logging_level):
    """Handler to emit messages via print_to_file().
    Combination of file and stream handlers.
    """
    import inspect
    caller = os.path.realpath(inspect.stack()[1][1])
    log_file = caller[:-3] + '.log'
    logging.basicConfig(level=level)
    h = myHandler(log_file)
    global mh
    h.set_name(mh)
    h.setLevel(level)
    logging.getLogger().addHandler(h)


def remove_my_handler():
    """Remove all MyHandler handlers."""
    global mh
    remove_handler_by_name(mh)


def add_text_handler(textEdit, level=s.logging_level):
    """Handler to show logs in master window textEdit."""
    logging.basicConfig(level=level)
    h = MyTextLoggingHandler(textEdit)
    global mtlh
    h.set_name(mtlh)
    h.setLevel(level)
    logging.getLogger().addHandler(h)


def remove_text_handler():
    """Remove all textEdit handlers."""
    global mtlh
    remove_handler_by_name(mtlh)


def add_file_handler(log_file, level=s.logging_level):
    """Handler to write logs into file."""
    logging.basicConfig(level=level)
    h = MyFileLoggingHandler(log_file)
    global mflh
    h.set_name(mflh)
    h.setLevel(level)
    logging.getLogger().addHandler(h)


def remove_file_handler():
    """Remove all file handlers."""
    global mflh
    remove_handler_by_name(mflh)


def add_stream_handler(level=s.logging_level):
    """Handler to write logs into stream."""
    logging.basicConfig(level=level)
    h = logging.StreamHandler()
    global mslh
    h.set_name(mslh)
    h.setLevel(level)
    global fmt
    h.setFormatter(fmt)
    logging.getLogger().addHandler(h)


def remove_stream_handler():
    """Remove all stream handlers."""
    global mslh
    remove_handler_by_name(mslh)


"""
    NOTE Doesn't write to file system messages and uncatched errors

    h = myHandler(log_file) # remove old log file
    log_capture_string = io.StringIO()
    ch = logging.StreamHandler(log_capture_string)
    ch.setLevel(logging.DEBUG)
    global fmt
    ch.setFormatter(fmt)
    logging.getLogger().addHandler(ch)

    log_contents = log_capture_string.getvalue()
    if len(log_contents) != lines_count:
        log_contents = log_contents[lines_count:]
        lines_count = log_capture_string.tell()
        relpath = os.path.relpath(file_name, start=os.getcwd())
        print_to_file(log_file, '{} {}'.format(counter, relpath))
        print_to_file(log_file, log_contents)
    log_capture_string.close()
"""
