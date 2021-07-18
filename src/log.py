#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Custom logging handlers.
The Text handler is created once per session and prints into CAE textEdit.
The File handler is created any time you open a model.
Log file has the same name as a model.

Classes for catching CCX and CGX output.
Difference between StdoutReaderLogger and CgxStdoutReaderLogger is that
first one deals with app which runs and quits, while CGX continues to run.
"""

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
import settings

mh = 'MyHandler'
mtlh = 'MyTextLoggingHandler'
mflh = 'MyFileLoggingHandler'
mslh = 'MyStreamLoggingHandler'
fmt = logging.Formatter('%(levelname)s: %(message)s')


# def initialize_logging(level=settings.s.logging_level):
#     global fmt
#     logging.basicConfig(level=level, format=fmt._fmt)


def get_logging_info():
    hh = logging.getLogger().handlers
    if not len(hh):
        msg = 'No logging handlers.'
    else:
        msg = 'Total {} logging handlers:'.format(len(hh))
        for h in hh:
            msg += '\n' + h.get_name()
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


def add_my_handler(level=settings.s.logging_level):
    """Handler to emit messages via print_to_file().
    Combination of file and stream handlers.
    """
    caller = os.path.realpath(inspect.stack()[1][1])
    log_file = caller[:-3] + '.log'
    h = myHandler(log_file)
    h.set_name(mh)
    h.setLevel(level)
    logging.getLogger().addHandler(h)


def remove_my_handler():
    """Remove all MyHandler handlers."""
    remove_handler_by_name(mh)


def add_text_handler(textEdit, level=settings.s.logging_level):
    """Handler to show logs in master window textEdit."""
    h = MyTextLoggingHandler(textEdit)
    h.set_name(mtlh)
    h.setLevel(level)
    logging.getLogger().addHandler(h)


def remove_text_handler():
    """Remove all textEdit handlers."""
    remove_handler_by_name(mtlh)


def add_file_handler(log_file, level=settings.s.logging_level):
    """Handler to write logs into file."""
    h = MyFileLoggingHandler(log_file)
    h.set_name(mflh)
    h.setLevel(level)
    logging.getLogger().addHandler(h)


def remove_file_handler():
    """Remove all file handlers."""
    remove_handler_by_name(mflh)


def add_stream_handler(level=settings.s.logging_level):
    """Handler to write logs into stream."""
    h = logging.StreamHandler()
    h.set_name(mslh)
    h.setLevel(level)
    global fmt
    h.setFormatter(fmt)
    logging.getLogger().addHandler(h)


def remove_stream_handler():
    """Remove all stream handlers."""
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



stdout_readers = []


def list_threads():
    """List currently running threads."""
    t_names = sorted([t.name for t in threading.enumerate() \
        if t.name != threading.main_thread().name])
    msg = '\nRunning threads:\n' + '\n'.join(t_names) + '\n'
    logging.debug(msg)


class StdoutReaderLogger:
    """Read and log output of a console app.
    Used in job.submit() and job.rebuild_ccx()."""

    def __init__(self, stdout, prefix, read_output=True, connection=None):
        self.stdout = stdout
        self.prefix = prefix
        self.read_output = read_output
        self.connection = connection
        self.active = True
        self.name = None

    def stop(self):
        self.active = False

    def log_line(self, line):
        """Process one non-empty stdout message."""
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

    def read_and_log(self):
        """Infininte cycle to read process'es stdout.
        CAE dies during fast logging.
        So we collect lines during 1 second - then log.
        """
        log_msg = ''
        start_time = time.perf_counter()
        # time.sleep(3) # Save from 100% CPU bug

        # Read and log output
        while self.active:
            line = self.stdout.readline() # freezes untile read a line
            line = line.replace(b'\r', b'') # for Windows
            if line == b' \n':
                continue
            if line != b'':
                line = line.decode().rstrip()
                log_msg += line + '\n'
                delta = time.perf_counter() - start_time
                if delta > 1:
                    start_time = time.perf_counter()
                    self.log_line(log_msg)
                    log_msg = ''
            else:
                # Here we get if there is no lines left in the stdout
                break

        if len(log_msg):
            self.log_line(log_msg)
        logging.debug('{} STOPPED.'.format(self.name))

    def start(self):
        """Read and log CGX stdout.
        Attention: it's infinite stdout reader!
        An old reader quits if a new one is started.
        """
        self.name = 'thread_{}_{}_{}'\
            .format(threading.active_count(),
                self.prefix, int(time.time()))
        t = threading.Thread(target=self.read_and_log,
            args=(), name=self.name, daemon=True)
        t.start()
        msg = '{} STARTED.'.format(self.name)
        logging.debug(msg)


class CgxStdoutReaderLogger(StdoutReaderLogger):
    """Read and log CGX output."""

    def log_line(self, line):
        logging.log(25, line)

    def filter_backspaces(self, line):
        """Posting to CGX window outputs to std line with backspaces:
        pplploplotplot plot eplot e plot e Oplot e OUplot e OUTplot e OUT-plot e OUT-Splot e OUT-S3
        """
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

    def read_and_log(self):
        """Infininte cycle to read CGX stdout.
        CAE dies during fast logging.
        So we may use time.sleep() after each log line.
        """

        # Flush CGX buffer
        if self.connection is not None:
            self.connection.post(' ')

        # Read and log output
        # time.sleep(3) # Save from 100% CPU bug
        while self.active:
            line = self.stdout.readline() # freezes untile read a line
            line = line.replace(b'\r', b'') # for Windows
            if line == b' \n':
                continue
            if b'key: from string' in line:
                continue
            if line != b'':
                line = self.filter_backspaces(line)
                line = line.decode().rstrip()
                self.log_line(line)
                time.sleep(0.03)
            else:
                # Here we get if there is no lines left in the stdout
                break

        logging.debug('{} STOPPED.'.format(self.name))


def start_stdout_reader(stdout, prefix, read_output):
    """Start stdout reading and logging thread."""
    sr = StdoutReaderLogger(stdout, 'read_stdout', read_output)
    global stdout_readers
    stdout_readers.append(sr)
    sr.start()


def start_cgx_stdout_reader(stdout, prefix, read_output, connection):
    """Start CGX stdout reading and logging thread."""
    sr = CgxStdoutReaderLogger(stdout, prefix, read_output, connection)
    global stdout_readers
    stdout_readers.append(sr)
    sr.start()


def stop_stdout_readers():
    """Quit all active logging threads."""
    global stdout_readers
    readers = [sr for sr in stdout_readers if sr.active]
    if len(readers):
        msg = 'Stopping threads:\n'
        for sr in readers:
            msg += sr.name + '\n'
        logging.debug(msg)
        for sr in readers:
            sr.stop()
        time.sleep(1)


@tests.test_wrapper()
def test():
    """Run main window and test loggers."""

    # Create application
    app = QtWidgets.QApplication(sys.argv)

    # Show main window
    stop_logging()
    import gui.window
    f = gui.window.Factory()
    f.run_master(path.p.main_xml) # has add_text_handler()

    logging.info('\nqwe\n')
    logging.warning('rty')

    # Execute application
    app.exec()


if __name__ == '__main__':
    test() # run test
