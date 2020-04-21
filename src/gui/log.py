#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, April 2020
Distributed under GNU General Public License v3.0

Custon logging handlers for catching CCX and CGX output.
The Text handler is created once per session
and prints into CAE's textEdit.
The File handler is created any time you open a model.
Log file has the same name as model. """


# Standard modules
import os
import sys
import re
import time
import logging
import queue
import threading
import multiprocessing

# External modules
from PyQt5 import QtGui, QtWidgets


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
        msg = LogRecord.getMessage()
        while msg.startswith('\n'):
            self.log_msg('<p></p>')
            msg = msg[1:]

        # Keep all newlines in message
        msg = msg.replace('\n', '<br/>')

        if LogRecord.levelname != 'Level 25':
            msg = '<span style=\'color:{};\'>{}</span>: {}'\
                .format(color, LogRecord.levelname, msg)
        msg = '<p style=\'margin:0px;\'>{}</p>'.format(msg)

        # Print message and scroll it to the end
        self.target.append(msg)
        self.target.moveCursor(QtGui.QTextCursor.End)
        self.flush()


# Handler to write the job's log file
class MyFileLoggingHandler(MyLoggingHandler):

    def emit(self, LogRecord):

        # Move newlines before the levelname
        msg = LogRecord.getMessage()
        while msg.startswith('\n'):
            self.log_msg('')
            msg = msg[1:]

        msg = '{}: {}'.format(LogRecord.levelname, msg)
        if msg.startswith('Level 25: '):
            msg = msg[10:]

        # Write message to the log file
        with open(self.target, 'a') as f:
            f.write(msg + '\n')
        self.flush()


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


# Process one non-empty stdout message
def logLine(line):
    if not len(line.strip()):
        return
    logging_levels = {
        'NOTSET': 0,
        'DEBUG': 10,
        'INFO': 20,
        'WARNING': 30,
        'ERROR': 40,
        'CRITICAL': 50}
    if line.startswith(tuple(logging_levels.keys())):
        match = re.search('(^\w+): (.+)', line) # levelname and message
        if match: # skip logging of empty strings
            level = match.group(1)
            line = match.group(2)
            logging.log(logging_levels[level], line)
    else:
        logging.log(25, line)


# Read and log pipe lines
# TODO Time spent is printed before logging is finished
# TODO CGX prints logs after close or after reopen
# TODO pplploplotplot plot eplot e plot e Oplot e OUplot e OUTplot e OUT-plot e OUT-Splot e OUT-S3
# TODO join old threads
def read_output(pipe, name):

    # Put pipe lines into queue while process is running
    def enqueue_output_1(pipe):
        while True:
            line = pipe.readline()
            if line != b'':
                line = line.decode().strip()
                logLine(line)
                time.sleep(0.04)
                pipe.flush()
                sys.stdout.flush()
            else:
                time.sleep(0.1)
        pipe.close()

    def enqueue_output_2(pipe, q):
        line = ''
        prev = ''
        for char in iter(pipe.read, b''):
            char = char.decode()
            if prev + char != '\n':
                line += char
                # print(char, end='')
            else:
                # print()
                q.put_nowait(line[:-1])
                line = ''
            prev = char
        q.put_nowait('END') # a mark to break while loop
        pipe.close()
        # print('END')

    def log_output(pipe, q):
        while True:
            try:
                # time.sleep(0.05)
                line = q.get_nowait()
                if line == 'END':
                    break
                logLine(line)
            except queue.Empty:
                QtWidgets.qApp.processEvents() # do not block GUI
                time.sleep(0.1) # reduce CPU usage

    # for t in threading.enumerate():
    #     if t.name != 'MainThread':
    #         t.join(timeout=5)
    t = threading.Thread(target=enqueue_output_1,
        args=(pipe, ), name=name, daemon=True)
    t.start()
    # list_threads()


# def list_threads():
#     n = threading.active_count()
#     t = threading.current_thread()
#     print('n={}, current thread is {}'.format(n, t.name))
#     for t in threading.enumerate():
#         print(t.name)
#     print()
