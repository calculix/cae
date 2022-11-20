#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2023
Distributed under GNU General Public License v3.0

Custom logging handlers.
The File handler is created any time you open a model.
Log file has the same name as a model.
"""

# Standard modules
import os
import sys
import logging

# My modules
from settings import s
from path import p


mh = 'MyHandler'
mflh = 'MyFileLoggingHandler'
mslh = 'MyStreamLoggingHandler'
fmt = logging.Formatter('%(levelname)s: %(message)s')


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


class MyFileLoggingHandler(MyLoggingHandler):
    """Called at each file import.
    Handler to write logs into file.
    """

    def emit(self, LogRecord):

        # Move newlines before the levelname
        msg = LogRecord.getMessage()
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


def colored(r, g, b, text):
    return "\033[38;2;{};{};{}m{}\033[38;2;255;255;255m".format(r, g, b, text)


class MyStreamLoggingHandler(MyLoggingHandler):
    """Handler to show logs in terminal."""

    def emit(self, LogRecord):

        # Message color depending on logging level
        color = {
            'NOTSET':   (255, 255, 255), # White
            'DEBUG':    (128, 128, 255), # Blue
            'INFO':     (128, 255, 128), # Green
            'WARNING':  (255, 255,   0), # Yellow
            'ERROR':    (255,   0,   0), # Red
            'CRITICAL': (255, 192, 203), # Pink
            }
        if LogRecord.levelname in color:
            color = color[LogRecord.levelname]
        else:
            color = (47, 79, 79) # DARKSLATEGRAY

        msg = LogRecord.getMessage()

        # Move leading and trailing newlines outside the message
        leading_newlines = ''
        trailing_newlines = ''
        while msg.startswith('\n'):
            leading_newlines += '\n'
            msg = msg[1:]
        while msg.endswith('\n'):
            trailing_newlines += '\n'
            msg = msg[:-1]

        if LogRecord.levelname == 'Level 25':
            msg = leading_newlines + msg + trailing_newlines
        else:
            msg = leading_newlines + colored(*color, LogRecord.levelname) \
                + ': ' + msg + trailing_newlines

        print(msg)


def stop_logging():
    """Stop all logging handlers, even default stdout."""
    l = logging.getLogger()
    for h in l.handlers:
        h.close()
    l.handlers = []
    l.addHandler(logging.NullHandler())
    l.propagate = False


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
    Is used in test() functions.
    """
    import inspect
    caller = os.path.realpath(inspect.stack()[1][1])
    log_file = caller[:-3] + '.log'
    global fmt
    logging.basicConfig(level=level, format=fmt._fmt)
    h = myHandler(log_file)
    global mh
    h.set_name(mh)
    h.setLevel(level)
    logging.getLogger().addHandler(h)


def add_file_handler(log_file, level=s.logging_level):
    """Handler to write logs into file."""
    h = MyFileLoggingHandler(log_file)
    global mflh
    h.set_name(mflh)
    h.setLevel(level)
    h.setFormatter(fmt)
    logging.getLogger().addHandler(h)


def remove_file_handler():
    """Remove all file handlers."""
    global mflh
    remove_handler_by_name(mflh)
    if os.path.exists(p.log):
        os.remove(p.log)


def add_stream_handler(level=s.logging_level):
    """Handler to write logs into stream."""
    global fmt
    h = MyStreamLoggingHandler(sys.stdout)
    global mslh
    h.set_name(mslh)
    h.setLevel(level)
    h.setFormatter(fmt)
    logging.getLogger().addHandler(h)
