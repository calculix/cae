#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, July 2020
Distributed under GNU General Public License v3.0

Utilities for testing """


# Standard modules
import os
import sys
import logging


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
    line = ' '.join([str(arg) for arg in args])
    line = line.rstrip() + '\n'
    with open(log_file, 'a') as f:
        f.write(line)
    sys.stdout.write(line)


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
