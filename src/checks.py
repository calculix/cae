#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Utilities to test user's system configuration.
Is called on the application startup.

Run test: Ctrl+F5 from VSCode.
"""


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


class Checks:
    """Tests user system configuration etc."""
    cmd = [sys.executable, '-m', 'pip', 'list']
    modules = subprocess.check_output(cmd)
    modules = modules.decode().rstrip()

    @staticmethod
    def check_os():
        """Exit if OS is not Linux or Windows."""
        if os.name not in ['nt', 'posix']:
            msg = 'Sorry, {} OS is not supported.'.format(os.name)
            logging.error(msg)
            raise SystemExit # the best way to exit

    @staticmethod
    def check_python():
        """Check python version.
        sys.version_info:
        (3, 5, 2, 'final', 0)
        """
        v1 = sys.version_info[:3]
        v2 = (3, 4)
        if sys.version_info < v2:
            msg = 'You have Python {}.{}.{}.\n'.format(*v1) \
                + 'Required version is {}.{} or newer.'.format(*v2)
            logging.error(msg)
            raise SystemExit # the best way to exit
        else:
            msg = 'Python version is {}.{}.{}.'.format(*v1)
            logging.info(msg)

    @staticmethod
    def check_default_web_browser():
        """Get default web browser."""
        wb = webbrowser.get()
        if 'nt' in os.name:
            wbname = wb.__class__.__name__
        else:
            wbname = wb.name
        msg = 'Default web browser is {}.'.format(wbname)
        logging.info(msg)

    @staticmethod
    def check_requirements():
        """Check each package from requirements.txt."""
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

    @classmethod
    def check_package(self, name):
        """Check if package is installed."""
        # try:
        #     importlib.import_module(name)
        if name.lower() in self.modules.lower():
            logging.info(name + ' OK')
            return True
        else:
        # except ImportError:
            msg = 'Package \'' + name \
                + '\' is not installed.'
            logging.warning(msg)
            return False

    @staticmethod
    def check_required_package(name):
        """Check if required package is installed."""
        if Checks.check_package(name):
            return True
        else:
            return Checks.install_package(name)

    @staticmethod
    def install_package(name, prefix=''):
        """Automatically install package."""
        try:
            logging.info('Installing ' + name + '...')
            cmd = [sys.executable, '-m', 'pip', prefix + 'install', name]
            # subprocess.check_call(cmd)
            msg = subprocess.check_output(cmd)
            msg = msg.decode().rstrip()
            logging.info(msg)
            return True
        except:
            msg = 'Can not {}install required package \''.format(prefix) \
                + name + '\'. ' \
                + 'Please, {}install it manually.'.format(prefix)
            logging.error(msg)
            return False

    @classmethod
    def uninstall_package(cls, name):
        """Automatically uninstall package."""
        cls.install_package(name, prefix='un')

    @staticmethod
    def check_all():
        """Run all checks."""
        log_file = __file__[:-3] + '.log'
        log.print_to_file(log_file, 'STARTUP TESTS\n')
        Checks.check_os()
        Checks.check_python()
        Checks.check_default_web_browser()
        Checks.check_requirements() # containts SystemExit


def run_startup_checks():
    """Tests running before the app start."""
    log.stop_logging()
    logging.disable(logging.NOTSET) # switch on logging
    log.add_my_handler()
    Checks.check_all()
    log.remove_my_handler()


@tests.test_wrapper()
def test():
    """Run some checks."""
    log.stop_logging()
    logging.disable(logging.NOTSET) # switch on logging
    log.add_my_handler()
    Checks.uninstall_package('unv2ccx')
    Checks.check_all() # install back 'unv2ccx'
    Checks.check_package('qwe')
    Checks.check_required_package('rty')

if __name__ == '__main__':
    test() # run test
