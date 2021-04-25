#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

MasterWindow class. Here also we keep links to another
app windows like CGX and web browser. Those links (WIDs)
are needed to align application windows.

For documentation on used functions refer to:
http://python-xlib.sourceforge.net/doc/html/python-xlib_21.html
https://github.com/python-xlib/python-xlib
https://github.com/asweigart/pyautogui """


# Standard modules
import os
import sys
import time
import logging
import subprocess
import inspect
import webbrowser
from shutil import which
if os.name == 'nt':
    import ctypes
    from ctypes import wintypes

# External modules
from PyQt5 import QtWidgets, uic

# My modules
try:
    # Normal run
    import gui
    import gui.connection
except:
    # Test run
    sys_path = os.path.abspath(__file__)
    sys_path = os.path.dirname(sys_path)
    sys_path = os.path.join(sys_path, '..')
    sys_path = os.path.normpath(sys_path)
    sys_path = os.path.realpath(sys_path)
    sys.path.insert(0, sys_path)
    import clean
    import gui
    import gui.connection
    import path
    import settings
    import window


# Main window
class MasterWindow(QtWidgets.QMainWindow):

    # Create main window
    """
    p - Path
    s - Settings
    """
    def __init__(self, p, s):
        self.p = p # global paths
        self.s = s # global settings
        self.stdout_readers = [] # for slave window
        self.slave_title = None
        self.slave_process = None # running process to send commands to
        self.connections = {} # 1:wc1, 2:wc2...

        QtWidgets.QMainWindow.__init__(self) # create main window
        uic.loadUi(p.main_xml, self) # load form

        # TODO Wrond architecture - its not window property
        self.desktopSize = QtWidgets.QDesktopWidget().availableGeometry()

        # Handler to show logs in the CAE's textEdit
        if hasattr(self, 'textEdit'): # skip for test_sendkeys
            gui.log.add_text_handler(self.textEdit)

        # Caller fuction name:
        # cgx.open_inp | cgx.open_frd | other
        self.mode = None

    # Run slave process
    # Close opened CGX (if any) and open a new one
    # Get window ID, align windows and post to CGX
    # TODO Move it to the SlaveWindow
    def run_slave(self, cmd, slave_title):
        self.slave_title = slave_title

        # Kill previous slave window
        self.kill_slave()

        # Run command - open new slave window - text editor or CGX
        # Open CGX without terminal/cmd window
        if self.p.op_sys == 'windows':
            shell=True
        if self.p.op_sys == 'linux':
            shell=False
        self.slave_process = subprocess.Popen(cmd.split(),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=shell)
        msg = '{} PID={}'.format(slave_title, self.slave_process.pid)
        logging.debug(msg)

        # Create connection between master and slave windows
        self.create_connection(1, slave_title)

        # Start stdout reading and logging thread
        # if hasattr(self, 'textEdit'): # skip for test_sendkeys
        sr = gui.log.CgxStdoutReader(
            self.slave_process.stdout, 'read_cgx_stdout', self)
        self.stdout_readers.append(sr)
        sr.start()

        # Caller fuction name:
        # cgx.open_inp | cgx.open_frd | other
        self.mode = inspect.stack()[1].function

    # Kill all slave processes
    def kill_slave(self):
        if self.slave_process is not None:
            self.connections[1].wid2 = None
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

    # Create connection between master and slave windows
    def create_connection(self, _id, slave_title):
        if os.name == 'nt':
            wc = gui.connection.WindowConnectionWindows(self, slave_title)
        else:
            wc = gui.connection.WindowConnectionLinux(self, slave_title)
        self.connections[_id] = wc
        wc.connect()
        if self.s.align_windows:
            wc.align()

    # Kill logging threads
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


# TODO
class SlaveWindow:
    
    def __init__(self, slave_title):
        self.slave_title = slave_title


# Keycodes sending to text editor and CGX
def test_sendkeys():

    # Create application
    app = QtWidgets.QApplication([])
    p = path.Path()
    p.main_xml = os.path.join(p.config, 'sendkeys.xml')
    s = settings.Settings(p)

    # Create master window and assign methods for buttons
    w = MasterWindow(p, s)
    w.show()

    # Configure global logging level
    logging.getLogger().setLevel(s.logging_level)
    fmt = '%(levelname)s: %(message)s'
    logging.basicConfig(format=fmt)

    # Map methods to GUI buttons
    if os.name == 'nt':
        cmd1 = 'notepad.exe'
        title1 = 'Untitled - Notepad'
    else:
        cmd1 = 'gedit'
        title1 = 'Untitled Document 1 - gedit'
        if which(cmd1) is None:
            cmd1 = 'kate'
            title1 = 'Untitled Document 1 - Kate'
            if which(cmd1) is None:
                cmd1 = ''
                msg = 'Neither gedit nor kate is available!'
                logging.error(msg)
    if len(cmd1):
        w.b0.clicked.connect(lambda: w.run_slave(cmd1, title1))
    else:
        w.b0.setDisabled(True)
    w.b1.clicked.connect(lambda: gui.cgx.open_inp(w, s.start_model, 1))
    w.b2.clicked.connect(lambda: w.connections[1].post('plot n all'))
    w.b3.clicked.connect(lambda: w.connections[1].post('plot e all'))
    w.customSend.clicked.connect(
        lambda: w.connections[1].post(w.customEdit.text()))
    w.customSend.clicked.connect(
        lambda: w.customEdit.clear())
    all_symbols = '!"#$%&\'()*+,-./1234567890:;<=>?@' \
        + 'ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`' \
        + 'abcdefghijklmnopqrstuvwxyz{|}~'
    w.testAllSymbols.clicked.connect(lambda: w.connections[1].post(all_symbols))

    # Execute application + exit
    a = app.exec()
    w.kill_slave()

# Run test
if __name__ == '__main__':
    clean.screen()
    test_sendkeys()
    clean.cache()