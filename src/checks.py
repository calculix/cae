#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2021
Distributed under GNU General Public License v3.0

Checks user system configuration.
Run test: Ctrl+F5 from VSCode """

# TODO Replace print with logging

# Standard modules
import os
import sys
import importlib
import subprocess

# My modules
import clean

# Exit if OS is not Linux or Windows
def check_os():
    if os.name not in ['nt', 'posix']:
        msg = 'SORRY, {} OS is not supported.'.format(os.name)
        print(msg)
        raise SystemExit # the best way to exit

# Check python version
def check_python():
    """ sys.version_info:
    (3, 5, 2, 'final', 0) """
    v1 = sys.version_info[:3]
    v2 = (3, 4)
    msg = 'ERROR! You have Python {}.\n'.format(v1) \
        + 'Required version is {} or newer.'.format(v2)
    if sys.version_info < v2:
        print(msg)
        raise SystemExit # the best way to exit
    else:
        msg = 'Python version is {}.'.format(v1)
        print(msg)

# Check each package from requirements.txt
def check_requirements():
    path = os.path.join(os.path.dirname(__file__),
        '..', 'requirements.txt')
    with open(path) as f:
        requirements = f.readlines()
        for name in requirements:
            name = name.strip()
            if len(name):
                check_package(name)

# Check if package is installed
def check_package(name):
    try:
        importlib.import_module(name)
        print(name, 'OK')
    except ImportError:
        msg = 'Required package \'' + name \
            + '\' is not installed. Trying to fix...'
        print(msg)
        install_package(name)

# Install package from code
def install_package(name):
    try:
        cmd = [sys.executable, '-m', 'pip', 'install', name]
        subprocess.check_call(cmd)
    except:
        msg = 'Can not install required package \'' \
            + name + '\'.\n' \
            + 'Please, install it manually.'
        print(msg)
        raise SystemExit # the best way to exit

if __name__ == '__main__':
    check_os()
    clean.screen()
    check_python()
    check_requirements()
    check_package('qwe')
