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
import webbrowser

# My modules
import path
import clean


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

# Return spent time delta in format hh:mm:ss.s
def get_time_delta(start):
    delta = time.perf_counter() - start
    return '{:02d}:{:02d}:{:04.1f}'\
        .format(int(delta/3600), int(delta%3600/60), delta%3600%60)

# Log spent time delta in format hh:mm:ss.s
def log_time_delta(start, log_file=None):
    msg = '\nTotal ' + get_time_delta(start)
    print(log_file, msg)


# Tests user system configuration etc.
class Check:

    # Prepare logging
    def __init__(self):
        self.log_file = __file__[:-3] + '.log'

    # Initialize logging
    def start_logging(self):
        h = myHandler(self.log_file)
        logging.getLogger().addHandler(h)
        logging.getLogger().setLevel(logging.NOTSET) # 0
        print(self.log_file, 'STARTUP TESTS\n')

    # After all tests stop logging into tests.log
    def stop_logging(self):
        logging.getLogger().handlers = []

    # Exit if OS is not Linux or Windows
    def check_os(self):
        if os.name not in ['nt', 'posix']:
            msg = 'Sorry, {} OS is not supported.'.format(os.name)
            raise SystemExit(msg) # the best way to exit

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

    # Get default web browser
    def check_default_web_browser(self):
        wb = webbrowser.get()
        msg = 'Default web browser is {}.'.format(wb.name)
        logging.info(msg)

    # Check each package from requirements.txt
    def check_requirements(self):
        path = os.path.normpath(os.path.join(
            os.path.dirname(__file__),
            '..', 'requirements.txt'))
        requirements = []
        with open(path) as f:
            requirements = f.readlines()
        for name in requirements:
            name = name.strip()
            if len(name) and not self.check_required_package(name):
                raise SystemExit # the best way to exit

    # Check if package is installed
    def check_package(self, name):
        try:
            importlib.import_module(name)
            logging.info(name + ' OK')
            return True
        except ImportError:
            msg = 'Package \'' + name \
                + '\' is not installed.'
            logging.warning(msg)
            return False

    # Check if required package is installed
    def check_required_package(self, name):
        if self.check_package(name):
            return True
        else:
            logging.warning('Trying to fix...')
            return self.install_package(name)

    # Automatically install package
    def install_package(self, name):
        try:
            cmd = [sys.executable, '-m', 'pip', 'install', name]
            subprocess.check_call(cmd)
            return True
        except:
            msg = 'Can not install required package \'' \
                + name + '\'. ' \
                + 'Please, install it manually.'
            logging.error(msg)
            return False

    # Run all checks
    def check_all(self):
        self.check_os()
        self.check_python()
        self.check_default_web_browser()
        self.check_requirements() # containts SystemExit


def run():
    ch = Check()
    ch.start_logging()
    ch.check_all()
    ch.stop_logging()

if __name__ == '__main__':
    start = time.perf_counter()
    clean.screen()
    ch = Check()
    ch.start_logging()
    ch.check_all()
    ch.check_package('qwe')
    ok = ch.check_required_package('rty')
    clean.cache()
    log_time_delta(start)
