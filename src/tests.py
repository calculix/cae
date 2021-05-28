#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Utilities for testing.

Class Check tests user system configuration.
Run test: Ctrl+F5 from VSCode """

# TODO Use unittest in all module's tests

# Standard modules
import os
import sys
import time
import logging

# My modules
import clean


# TODO Move it to log.py
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

# Return spent time delta in format hh:mm:ss.s
def get_time_delta(start):
    delta = time.perf_counter() - start
    return '{:02d}:{:02d}:{:04.1f}'\
        .format(int(delta/3600), int(delta%3600/60), delta%3600%60)

# Log spent time delta in format hh:mm:ss.s
def log_time_delta(start, log_file=None):
    msg = '\nTotal ' + get_time_delta(start)
    print(log_file, msg)


# List all .ext-files here and in all subdirectories
def scan_all_files_in(start_folder, ext, limit=1000000):
    all_files = []
    for f in os.scandir(start_folder):
        if f.is_dir():
            for ff in scan_all_files_in(f.path, ext):
                all_files.append(ff)
        elif f.is_file() and f.name.endswith(ext):
            ff = os.path.normpath(f.path)
            all_files.append(ff)
    return sorted(all_files)[:limit]

# Common wrapper for test() method in different modules
def test_wrapper():
    def wrap(method):
        def fcn():
            start = time.perf_counter()
            fmt = '%(levelname)s: %(message)s'
            logging.basicConfig(level=logging.NOTSET, format=fmt)
            clean.screen()
            method()
            clean.cache()
            log_time_delta(start)
        return fcn
    return wrap

# TODO Invent some test
@test_wrapper()
def test():
    pass

# Run test
if __name__ == '__main__':
    test()
