#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Custon logging handlers for catching CCX and CGX output.
The Text handler is created once per session
and prints into CAE's textEdit.
The File handler is created any time you open a model.
Log file has the same name as a model. """

# TODO There is StreamHandler in import.py

# TODO Logging can't be switched off completely - 
# - see test() in tree.py for example.

# Standard modules
import os
import sys
import re
import time
import logging
import threading

# External modules
from PyQt5 import QtGui, QtWidgets

# My modules
sys_path = os.path.abspath(__file__)
sys_path = os.path.dirname(sys_path)
sys_path = os.path.join(sys_path, '..')
sys_path = os.path.normpath(sys_path)
sys_path = os.path.realpath(sys_path)
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
import tests
import gui
import settings
import path


# Configure logging to emit messages via 'print' method
class myHandler(logging.Handler):

    def __init__(self, log_file):
        super().__init__()
        fmt = logging.Formatter('%(levelname)s: %(message)s')
        self.setFormatter(fmt)
        self.log_file = log_file

        # Remove old log file
        if os.path.isfile(log_file):
            os.remove(log_file)

    def emit(self, LogRecord):
        print(self.log_file, self.format(LogRecord))


# Redefine print method to write logs to file
def print(log_file, *args):
    line = ' '.join([str(arg) for arg in args]) + '\n'
    if log_file is not None:
        with open(log_file, 'a') as f:
            f.write(line)
    sys.stdout.write(line)


class MyLoggingHandler(logging.Handler):

    # Target could be textEdit or file
    def __init__(self, target):
        super().__init__() # create handler
        self.target = target
        fmt = logging.Formatter('%(levelname)s, %(module)s: %(message)s')
        self.setFormatter(fmt)


# Called only once on startup
# Handler to show logs in master window's textEdit
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
            color = '#2F4F4F' # DARKSLATEGRAY

        # Move newlines before the levelname
        msg = LogRecord.getMessage()
        while msg.startswith('\n'):
            self.target.append('<p></p>')
            msg = msg[1:]

        # Keep all newlines in message
        msg = msg.replace('\n', '<br/>')

        if LogRecord.levelname == 'Level 25':
            msg = '<p style=\"margin:0px; color:{};\">{}</p>'\
                .format(color, msg)
        else:
            msg = '<span style=\"color:{};\">{}</span>: {}'\
                .format(color, LogRecord.levelname, msg)
            msg = '<p style=\"margin:0px;\">{}</p>'.format(msg)

        # Print message and scroll it to the end
        self.target.append(msg)
        self.target.moveCursor(QtGui.QTextCursor.End)


# Called at each file import
# Handler to write logs into file
class MyFileLoggingHandler(MyLoggingHandler):

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

        # Write message to the log file
        with open(self.target, 'a') as f:
            f.write(msg + '\n')


# Only 2 logging hadlers are allowed
allowed_handlers = 2

# Switch off logging
def switch_off_logging():
    hh = logging.getLogger().handlers
    msg = 'Switching off {} logging handlers'.format(len(hh))
    logging.debug(msg)
    logging.getLogger().handlers = []
    return hh

# Switch on logging
def switch_on_logging(hh):
    for h in hh:
        logging.getLogger().addHandler(h)
    msg = 'Switching on {} logging handlers'.format(len(hh))
    logging.debug(msg)

"""
# Emit error log if there are too many handlers
def check_amount_of_handlers():
    amount = len(logging.getLogger().handlers)
    if amount > allowed_handlers:
        msg = 'Exceeded amount of allowed logging handlers (2).'
        logging.error(msg)
        switch_off_logging()
    return amount
"""

# Handler to show logs in master window's textEdit
def add_text_handler(textEdit):
    remove_text_handler()
    h = MyTextLoggingHandler(textEdit)
    h.set_name = 'MyTextLoggingHandler'
    logging.getLogger().addHandler(h)

# Remove all textEdit handlers
def remove_text_handler():
    hh = logging.getLogger().handlers
    for h in hh:
        if h.name == 'MyTextLoggingHandler':
            hh.remove(h)

# Handler to write logs into file
def add_file_handler(log_file):
    remove_file_handler()
    h = MyFileLoggingHandler(log_file)
    h.set_name = 'MyFileLoggingHandler'
    logging.getLogger().addHandler(h)
    if os.path.exists(log_file):
        os.remove(log_file)

# Remove all file handlers
def remove_file_handler():
    hh = logging.getLogger().handlers
    for h in hh:
        if h.name == 'MyFileLoggingHandler':
            hh.remove(h)


# Reads CCX and CGX outputs
class StdoutReader:

    def __init__(self, stdout, prefix, read_output=True, connection=None):
        self.stdout = stdout
        self.prefix = prefix
        self.read_output = read_output
        self.connection = connection
        self.active = True
        self.name = None
    
    def stop(self):
        self.active = False

    # Process one non-empty stdout message
    def log_line(self, line):
        if self.read_output:
            logging_levels = {
                'NOTSET': 0,
                'DEBUG': 10,
                'INFO': 20,
                'WARNING': 30,
                'ERROR': 40,
                'CRITICAL': 50}
            if line.strip().startswith(tuple(logging_levels.keys())):
                match = re.search('^\s*(\w+):*(.+)', line) # levelname and message
                if match: # skip logging of empty strings
                    level = match.group(1)
                    line = match.group(2)
                    logging.log(logging_levels[level], line)
            else:
                logging.log(25, line)

    def filter_backspaces(self, line):
        return line

    # Infininte cycle to read process'es stdout
    def read_and_log(self):

        # Save from 100% CPU bug
        time.sleep(3)
        
        # Flush CGX buffer
        if self.connection is not None:
            self.connection.post(' ')

        # Read and log output
        while self.active:
            line = self.stdout.readline()
            line = line.replace(b'\r', b'') # for Windows
            if line == b' \n':
                continue
            if b'key: from string' in line:
                continue
            if line != b'':
                line = self.filter_backspaces(line)
                line = line.decode().rstrip()
                self.log_line(line)
                time.sleep(0.03) # CAE dies during fast logging
            else:
                # Here we get if CGX is closed
                logging.debug('StdoutReader stopped')
                break

        # Exit from function crashes app on textEdit scrolling!
        while self.active:
            time.sleep(1)

    # Read and log CGX stdout
    # Attention: it's infinite stdout reader!
    # An old reader quits if a new one is started
    def start(self):
        self.name = 'thread_{}_{}_{}'\
            .format(threading.active_count(),
                self.prefix, int(time.time()))
        t = threading.Thread(target=self.read_and_log,
            args=(), name=self.name, daemon=True)
        t.start()

        # List currently running threads
        t_names = sorted([t.name for t in threading.enumerate() \
            if t.name != threading.main_thread().name])
        msg = '\nLogging threads:\n' + '\n'.join(t_names) + '\n'
        logging.debug(msg)


class CgxStdoutReader(StdoutReader):

    def log_line(self, line):
        logging.log(25, line)

    # Posting to CGX window outputs to std line with backspaces:
    # pplploplotplot plot eplot e plot e Oplot e OUplot e OUTplot e OUT-plot e OUT-Splot e OUT-S3
    def filter_backspaces(self, line):
        if 8 in line:
            l = bytearray()
            i = 1
            b = line[-i]
            while b != 8:
                l.insert(0, b)
                i += 1
                b = line[-i]
            return l
        else:
            return line


# Run app's main window and test loggers
@tests.test_wrapper()
def test():

    # Create application
    app = QtWidgets.QApplication(sys.argv)

    # Calculate absolute paths
    p = path.Path()

    # Read application's global settings
    s = settings.Settings(p)

    # Show main window
    f = gui.window.Factory(s)
    f.run_master(p.main_xml)

    # Execute application
    app.exec()

# Run test
if __name__ == '__main__':
    test()
