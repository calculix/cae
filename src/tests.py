#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Utilities for testing.

Class Check tests user system configuration.
Run test: Ctrl+F5 from VSCode """

# TODO Use unittest

# Standard modules
import os
import sys
import time
import logging
import importlib
import subprocess


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


# Convert seconds to format hh:mm:ss
def time_delta(start_time):
    delta = time.perf_counter() - start_time
    return time.strftime('%M:%S', time.gmtime(delta))


# Tests user system configuration
class Check:

    # Prepare logging
    def __init__(self):
        log_file = __file__[:-3] + '.log'
        h = myHandler(log_file)
        logging.getLogger().addHandler(h)
        logging.getLogger().setLevel(logging.DEBUG)
        print(log_file, 'STARTUP TESTS\n')

    # Exit if OS is not Linux or Windows
    def check_os(self):
        if os.name not in ['nt', 'posix']:
            msg = 'SORRY, {} OS is not supported.'.format(os.name)
            logging.warning(msg)
            raise SystemExit # the best way to exit

    # Check python version
    def check_python(self):
        """ sys.version_info:
        (3, 5, 2, 'final', 0) """
        v1 = sys.version_info[:3]
        v2 = (3, 4)
        if sys.version_info < v2:
            msg = 'You have Python {}.\n'.format(v1) \
                + 'Required version is {} or newer.'.format(v2)
            logging.error(msg)
            raise SystemExit # the best way to exit
        else:
            msg = 'Python version is {}.'.format(v1)
            logging.info(msg)

    # Check each package from requirements.txt
    def check_requirements(self):
        path = os.path.normpath(os.path.join(
            os.path.dirname(__file__),
            '..', 'requirements.txt'))
        with open(path) as f:
            requirements = f.readlines()
            for name in requirements:
                name = name.strip()
                if len(name):
                    self.check_package(name)

    # Check if package is installed
    def check_package(self, name):
        try:
            importlib.import_module(name)
            logging.info(name + ' OK')
        except ImportError:
            msg = 'Required package \'' + name \
                + '\' is not installed. Trying to fix...'
            logging.warning(msg)
            self.install_package(name)

    # Install package from code
    def install_package(self, name):
        try:
            cmd = [sys.executable, '-m', 'pip', 'install', name]
            subprocess.check_call(cmd)
        except:
            msg = 'Can not install required package \'' \
                + name + '\'. ' \
                + 'Please, install it manually.'
            logging.error(msg)
            raise SystemExit # the best way to exit


if __name__ == '__main__':
    ch = Check()
    ch.check_os()
    ch.check_python()
    ch.check_requirements()
    ch.check_package('qwe')
