#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Utilities for testing. Some modules has test() method decorated
with test_wrapper(). The decorator is described here.
"""

# Standard modules
import os
import time
import logging
import inspect

# My modules
from utils import clean


def scan_all_files_in(start_folder, ext, limit=1000000):
    """List all .ext-files here and in all subdirectories."""
    all_files = []
    for f in os.scandir(start_folder):
        if f.is_dir():
            for ff in scan_all_files_in(f.path, ext):
                all_files.append(ff)
        elif f.is_file() and f.name.endswith(ext):
            ff = os.path.normpath(f.path)
            all_files.append(ff)
    return sorted(all_files)[:limit]


def get_time_delta(start):
    """Return spent time delta in format hh:mm:ss.s."""
    delta = time.perf_counter() - start
    return '{:02d}:{:02d}:{:04.1f}'\
        .format(int(delta/3600), int(delta%3600/60), delta%3600%60)


def log_time_delta(start, log_file=None):
    """Log spent time delta in format hh:mm:ss.s"""
    msg = '\nTotal ' + get_time_delta(start)
    if log_file is None:
        print(msg)
    else:
        import log
        log.print_to_file(log_file, msg)
        print(msg)


def test_wrapper():
    """Common wrapper for test() method in different modules..."""
    def wrap(method):
        def fcn():
            start = time.perf_counter()

            caller = os.path.realpath(inspect.stack()[1][1])
            os.chdir(os.path.dirname(caller))

            # Remove old log file
            log_file = caller[:-3] + '.log'
            if os.path.isfile(log_file):
                os.remove(log_file)

            # First time logging configure for tests
            fmt = logging.Formatter('%(levelname)s: %(message)s')
            logging.basicConfig(level=logging.NOTSET, format=fmt._fmt)

            clean.screen()
            method()
            clean.cache()

            log_time_delta(start, log_file)
        return fcn
    return wrap
