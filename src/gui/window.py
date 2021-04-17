#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

MainWindow class. Here also we keep links to another
app windows like CGX and web browser. Those links (WIDs)
are needed to align application windows.

For documentation on used functions refer to:
http://python-xlib.sourceforge.net/doc/html/python-xlib_21.html
https://github.com/python-xlib/python-xlib
https://github.com/asweigart/pyautogui """


# Standard modules
import os
import time
import logging
import subprocess
import inspect
import webbrowser
if 'nt' in os.name:
    import ctypes
    from ctypes import wintypes

# External modules
from PyQt5 import QtWidgets, uic

# My modules
import gui
import gui.connection


# Main window
class MainWindow(QtWidgets.QMainWindow):

    # Create main window
    """
    p - Path
    s - Settings
    """
    def __init__(self, p, s):
        self.p = p
        self.s = s
        self.stdout_readers = []
        self.slave_title = None
        self.slave_process = None # running process to send commands to
        self.wc = None

        QtWidgets.QMainWindow.__init__(self) # create main window
        uic.loadUi(p.main_xml, self) # load form
        self.size = QtWidgets.QDesktopWidget().availableGeometry()

        # Handler to show logs in the CAE's textEdit
        gui.log.add_text_handler(self.textEdit)

        # Window ID - to pass keycodes to
        # self.wid1 = None # cae window
        # self.wid2 = None # cgx window
        # self.wid3 = None # keyword dialog
        # self.wid4 = None # help (web browser)
        self.keyboardMapping = None

        # INP | FRD - corresponds to opened file
        self.mode = None

    # Run slave process
    # Close opened CGX (if any) and open a new one
    # Get window ID, align windows and post to CGX
    def run_slave(self, params, slave_title):
        if not os.path.isfile(self.p.path_cgx):
            logging.error('CGX not found:\n' \
                + self.p.path_cgx)
            return

        # Kill previous slave window
        self.kill_slave()

        cmd = self.p.path_cgx + ' ' + params

        # Run command - open new slave window - text editor or CGX
        # Open CGX without terminal/cmd window
        if self.p.op_sys == 'windows':
            self.slave_process = subprocess.Popen(cmd.split(),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True)
        if self.p.op_sys == 'linux':
            self.slave_process = subprocess.Popen(cmd.split(),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT)
        msg = '{} PID={}'.format(slave_title, self.slave_process.pid)
        logging.debug(msg)

        # Create connection between master and slave windows
        if os.name == 'nt':
            self.wc = gui.connection.WindowConnectionWindows(self, slave_title)
        else:
            self.wc = gui.connection.WindowConnectionLinux(self, slave_title)
        self.wc.connect()
        if self.s.align_windows:
            self.wc.align()

        # Start stdout reading and logging thread
        sr = gui.log.CgxStdoutReader(
            self.slave_process.stdout, 'read_cgx_stdout', self)
        self.stdout_readers.append(sr)
        sr.start()

        # Read config to align model to iso view
        file_name = os.path.join(self.p.config, 'iso.fbd')
        if os.path.isfile(file_name):
            self.wc.post('read ' + file_name)
        else:
            logging.error('No config file iso.fbd')

        # Read config to register additional colors
        # Those colors are needed to paint sets and surfaces
        file_name = os.path.join(self.p.config, 'colors.fbd')
        if os.path.isfile(file_name):
            self.wc.post('read ' + file_name)
        else:
            logging.error('No config file colors.fbd')

        # Caller fuction name: cgx_inp | cgx_frd
        self.mode = inspect.stack()[1].function

    # Kill all slave processes
    def kill_slave(self):
        if self.slave_process is not None:
            self.wc.wid2 = None
            count = 0
            while self.slave_process.poll() is None:
                try:
                    if os.name == 'nt':
                        os.system('TASKKILL /F /PID {} /T'.format(self.slave_process.pid))
                    else:
                        self.slave_process.kill()
                except:
                    logging.error(traceback.format_exc())
                time.sleep(0.1)
                count += 1
                if count >= 10:
                    break
            if self.slave_process.poll() is None:
                msg = 'Can not kill {}, PID={}.'\
                    .format(self.slave_title, self.slave_process.pid)
                logging.warning(msg)
            else:
                msg = 'Killed {}, PID={}.'\
                    .format(self.slave_title, self.slave_process.pid)
                logging.debug(msg)
                self.slave_process = None

    def stop_stdout_readers(self):
        readers = [sr for sr in self.stdout_readers if sr.active]
        if len(readers):
            msg = 'Stopping threads:\n'
            for sr in readers:
                msg += sr.name + '\n'
                sr.stop()
            logging.debug(msg)
            time.sleep(1)

    # Open links from the Help menu
    def help(self, url):
        logging.info('Going to\n' + url)
        if not webbrowser.open(url, new=2):
            logging.warning('Can\'t open url: ' + url)


class MainWindowLinux(MainWindow):

    # Initialize variables and map methods to GUI buttons
    def __init__(self, p, s):
        super(MainWindowLinux, self).__init__(p, s)


class MainWindowWindows(MainWindow):

    # Initialize variables and map methods to GUI buttons
    def __init__(self, p, s):
        super(MainWindowWindows, self).__init__(p, s)
