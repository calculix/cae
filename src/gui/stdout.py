#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Classes for catching CCX and CGX output.
Difference between StdoutReaderLogger and CgxStdoutReaderLogger is that
first one deals with app which runs and quits, while CGX continues to run.
"""

# Standard modules
import re
import time
import logging
import threading


stdout_readers = []


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
        time.sleep(3) # Save from crush and 100% CPU bug

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
        time.sleep(3) # Save from crush and 100% CPU bug

        # Flush CGX buffer
        if self.connection is not None:
            self.connection.post(' ')

        # Read and log output
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
