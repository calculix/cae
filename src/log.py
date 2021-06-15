#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Custon logging handlers for catching CCX and CGX output.
The Text handler is created once per session
and prints into CAE textEdit.
The File handler is created any time you open a model.
Log file has the same name as a model. """

# Standard modules
import os
import sys
import re
import time
import logging
import inspect
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
import path
import tests

mh = 'MyHandler'
mtlh = 'MyTextLoggingHandler'
mflh = 'MyFileLoggingHandler'
mslh = 'MyStreamLoggingHandler'

# Redefine print method to write logs to file and stdout
def print(log_file, *args):
    line = ' '.join([str(arg) for arg in args]) + '\n'
    if log_file is not None:
        with open(log_file, 'a') as f:
            f.write(line)
    sys.stdout.write(line)


# Handler to emit messages via 'print' method
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


class MyLoggingHandler(logging.Handler):

    # Target could be textEdit or file
    def __init__(self, target):
        super().__init__() # create handler
        self.target = target
        fmt = logging.Formatter('%(levelname)s, %(module)s: %(message)s')
        self.setFormatter(fmt)


# Called only once on startup
# Handler to show logs in master window textEdit
class MyTextLoggingHandler(MyLoggingHandler):

    # Sends log messages to Window textEdit widget
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

        # Append message to the log file
        with open(self.target, 'a') as f:
            f.write(msg + '\n')


# Stop all logging handlers, even default stdout
def stop_logging():
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

# Handler to emit messages via 'print' method
# Combination of file and stream handlers
def add_my_handler(level=None):
    if level is None:
        import settings
        level = settings.s.logging_level
    caller = os.path.realpath(inspect.stack()[1][1])
    log_file = caller[:-3] + '.log'
    h = myHandler(log_file)
    h.set_name(mh)
    h.setLevel(level)
    logging.getLogger().addHandler(h)

# Remove all MyHandler handlers
def remove_my_handler():
    remove_handler_by_name(mh)

# Handler to show logs in master window textEdit
def add_text_handler(textEdit, level=None):
    if level is None:
        import settings
        level = settings.s.logging_level
    h = MyTextLoggingHandler(textEdit)
    h.set_name(mtlh)
    h.setLevel(level)
    logging.getLogger().addHandler(h)

# Remove all textEdit handlers
def remove_text_handler():
    remove_handler_by_name(mtlh)

# Handler to write logs into file
def add_file_handler(log_file, level=None):
    if level is None:
        import settings
        level = settings.s.logging_level
    h = MyFileLoggingHandler(log_file)
    h.set_name(mflh)
    h.setLevel(level)
    logging.getLogger().addHandler(h)
    """
    TODO Doesn't write to file system messages and uncatched errors

    h = log.myHandler(log_file) # remove old log file
    log_capture_string = io.StringIO()
    ch = logging.StreamHandler(log_capture_string)
    ch.setLevel(logging.DEBUG)
    fmt = logging.Formatter('%(levelname)s: %(message)s')
    ch.setFormatter(fmt)
    logging.getLogger().addHandler(ch)

    log_contents = log_capture_string.getvalue()
    if len(log_contents) != lines_count:
        log_contents = log_contents[lines_count:]
        lines_count = log_capture_string.tell()
        relpath = os.path.relpath(file_name, start=os.getcwd())
        print(log_file, '{} {}'.format(counter, relpath))
        print(log_file, log_contents)
    log_capture_string.close()
    """

# Remove all file handlers
def remove_file_handler():
    remove_handler_by_name(mflh)

# Handler to write logs into stream
def add_stream_handler(level=None):
    if level is None:
        import settings
        level = settings.s.logging_level
    h = logging.StreamHandler()
    h.set_name(mslh)
    h.setLevel(level)
    fmt = logging.Formatter('%(levelname)s: %(message)s')
    h.setFormatter(fmt)
    logging.getLogger().addHandler(h)

# Remove all stream handlers
def remove_stream_handler():
    remove_handler_by_name(mslh)


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


# Run main window and test loggers
@tests.test_wrapper()
def test():

    # Create application
    app = QtWidgets.QApplication(sys.argv)

    # Show main window
    import gui.window
    f = gui.window.Factory()
    f.run_master(path.p.main_xml)

    # Execute application
    app.exec()

# Run test
if __name__ == '__main__':
    test()
