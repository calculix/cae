#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2021
Distributed under GNU General Public License v3.0

Checks user system configuration.
Run test: Ctrl+F5 from VSCode """

# TODO Use unittest

# Standard modules
import os
import sys
import logging
import importlib
import subprocess

# My modules
import tests

class Check:

    # Prepare logging
    def __init__(self):
        log_file = __file__[:-3] + '.log'
        h = tests.myHandler(log_file)
        logging.getLogger().addHandler(h)
        logging.getLogger().setLevel(logging.DEBUG)
        tests.print(log_file, 'STARTUP TESTS\n')

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
