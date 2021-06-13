#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

The class tests user system configuration.
Is called on the application startup.
Run test: Ctrl+F5 from VSCode """

# TODO Use unittest

# Standard modules
import os
import sys
import logging
import importlib
import subprocess
import webbrowser

# My modules
import tests
import log


# Tests user system configuration etc.
class Checks:

    # Exit if OS is not Linux or Windows
    @staticmethod
    def check_os():
        if os.name not in ['nt', 'posix']:
            msg = 'Sorry, {} OS is not supported.'.format(os.name)
            logging.error(msg)
            raise SystemExit # the best way to exit

    # Check python version
    @staticmethod
    def check_python():
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
    @staticmethod
    def check_default_web_browser():
        wb = webbrowser.get()
        msg = 'Default web browser is {}.'.format(wb.name)
        logging.info(msg)

    # Check each package from requirements.txt
    @staticmethod
    def check_requirements():
        path = os.path.normpath(os.path.join(
            os.path.dirname(__file__),
            '..', 'requirements.txt'))
        requirements = []
        with open(path) as f:
            requirements = f.readlines()
        for name in requirements:
            name = name.strip()
            if len(name) and not Checks.check_required_package(name):
                raise SystemExit # the best way to exit

    # Check if package is installed
    @staticmethod
    def check_package(name):
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
    @staticmethod
    def check_required_package(name):
        if Checks.check_package(name):
            return True
        else:
            logging.warning('Trying to fix...')
            return Checks.install_package(name)

    # Automatically install package
    @staticmethod
    def install_package(name):
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
    @staticmethod
    def check_all():
        Checks.check_os()
        Checks.check_python()
        Checks.check_default_web_browser()
        Checks.check_requirements() # containts SystemExit


# Tests running before the app start
def run_startup_checks():
    logging.basicConfig(level=logging.NOTSET)
    log.stop_logging()
    log_file = __file__[:-3] + '.log'
    log.add_my_handler()
    log.print(log_file, 'STARTUP TESTS\n')
    Checks.check_all()
    log.remove_my_handler()

# Run some checks
@tests.test_wrapper()
def test():
    log.stop_logging()
    log.add_my_handler()
    Checks.check_all()
    Checks.check_package('qwe')
    Checks.check_required_package('rty')

# Run test
if __name__ == '__main__':
    test()
