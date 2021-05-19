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
import traceback
from shutil import which
if os.name == 'nt':
    import ctypes
    from ctypes import wintypes

# External modules
from PyQt5 import QtWidgets, uic

# My modules
sys_path = os.path.abspath(__file__)
sys_path = os.path.dirname(sys_path)
sys_path = os.path.join(sys_path, '..')
sys_path = os.path.normpath(sys_path)
sys_path = os.path.realpath(sys_path)
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
import tests
import gui
import gui.connection
import path
import settings

# Common wrapper for MasterWindow/SlaveWindow __init__() method
def init_wrapper():
    def wrap(method):
        def fcn(self, arg):
            before = get_opened_windows()

            method(self, arg)

            # Wait for window to open ang get its info
            new_windows_infos = []
            start = time.perf_counter() # start time
            while not len(new_windows_infos):
                after = get_opened_windows()
                time.sleep(0.3)
                new_windows_infos = get_new_windows_infos(before, after)
                if time.perf_counter() - start > 3:
                    msg = 'New slave window starts too slowly.'
                    logging.error(msg)
                    break

            if len(new_windows_infos):
                wi = new_windows_infos[0]
                try:
                    wi.pid = self.process.pid
                except:
                    pass # MasterWindow has no process variable
                msg = wi.to_string()
                logging.debug(msg)
                self.info = wi

        return fcn
    return wrap


""" Window factory class - used to create master and slave windows
and bind/connect them together. Logging facility as killing methods
for slave window are also maintained here. """
class Factory:

    def __init__(self, s):
        self.s = s # global settings
        self.p = s.getp() # Path
        self.mw = None # master window
        self.sw = None # slave window
        self.stdout_readers = [] # for slave window
        args = [None, None]
        if os.name == 'nt':
            wc1 = gui.connection.WindowConnectionWindows(*args)
            wc2 = gui.connection.WindowConnectionWindows(*args)
        else:
            wc1 = gui.connection.WindowConnectionLinux(*args)
            wc2 = gui.connection.WindowConnectionLinux(*args)
        self.connection = None

    # Draw window form described in xml UI file
    def run_master(self, xml):
        self.mw = MasterWindow(xml)

    # Generate dialog window and show it
    def run_master_dialog(self, s, k, i):
        args = [s, k, i]
        self.mw = gui.dialog.KeywordDialog(args) # one argument for init_wrapper()
        
        # Actions
        self.mw.buttonBox.accepted.connect(self.mw.ok)
        self.mw.buttonBox.button(QtWidgets.QDialogButtonBox.Reset).clicked.connect(self.mw.reset)
        self.mw.buttonBox.helpRequested.connect(self.open_web_browser)

        d = self.mw.exec()
        self.kill_slave()
        return d

    # Open HTML help page in a default web browser
    # TODO Opens Chrome and Firefox - check others
    def open_web_browser(self):
        url = self.mw.get_help_url()
        if not os.path.isfile(url):
            msg = 'Help page does not exist:\n{}'.format(
                os.path.relpath(url, start=self.p.app_home_dir))
            logging.error(msg)
            return
        wb = webbrowser.get()
        cmd = wb.name + ' --new-window ' + url
        logging.info('Going to\n' + url)
        self.run_slave(cmd)
        self.create_connection()

    # Close opened slave window and open a new one
    def run_slave(self, cmd):
        self.kill_slave()
        self.sw = SlaveWindow(cmd)

        # Start stdout reading and logging thread
        if self.p.path_cgx in cmd:
            self.start_stdout_reader('read_cgx_stdout')

    # Kill all slave processes
    # TODO Firefox process is being killed correctly, but window remains
    def kill_slave(self):
        if self.sw is None:
            return
        if self.sw.process is None:
            return

        if self.p.path_cgx in self.sw.cmd:
            self.stop_stdout_readers()

        count = 0
        while self.sw.process.poll() is None:
            try:
                if os.name == 'nt':
                    os.system('TASKKILL /F /PID {} /T'.format(self.sw.process.pid))
                else:
                    self.sw.process.kill()
            except:
                logging.error(traceback.format_exc())
            time.sleep(0.1)
            count += 1
            if count >= 10:
                break
        if self.sw.process.poll() is None:
            msg = 'Can not kill {}, PID={}.'\
                .format(self.sw.info.wname, self.sw.process.pid)
            logging.warning(msg)
        else:
            msg = 'Killed {}, PID={}.'\
                .format(self.sw.info.wname, self.sw.process.pid)
            logging.debug(msg)
            self.sw.process = None
            self.sw.info = None

    # Start stdout reading and logging thread
    def start_stdout_reader(self, aim):
        sr = gui.log.CgxStdoutReader(self.sw.process.stdout, aim, self)
        self.stdout_readers.append(sr)
        sr.start()

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

    # Connect master and slave windows and align them
    def create_connection(self):
        if self.mw is None:
            logging.error('No master window.')
            return
        if self.sw is None:
            logging.error('No slave window.')
            return

        args = [self.mw, self.sw]
        if os.name == 'nt':
            self.connection = gui.connection.WindowConnectionWindows(*args)
        else:
            self.connection = gui.connection.WindowConnectionLinux(*args)

        if self.s.align_windows:
            self.connection.align_windows()


# Left window - main window of the app
class MasterWindow(QtWidgets.QMainWindow):

    # Load form and show the window
    @init_wrapper()
    def __init__(self, xml):
        self.info = None # WindowInfo will be set in @init_wrapper

        super().__init__() # create main window
        uic.loadUi(xml, self) # load form

        # Handler to show logs in the CAE's textEdit
        if hasattr(self, 'textEdit'): # skip for test_sendkeys
            gui.log.add_text_handler(self.textEdit)

        self.show()

    # Open links from the Help menu
    def help(self, url):
        if webbrowser.open(url, new=2):
            logging.info('Going to\n' + url)
        else:
            logging.warning('Can not open url\n' + url)


""" Right application's window - CGX, text editor,
web browser. Starts as a separate process """
class SlaveWindow:

    # Run slave process without terminal window, set window info
    @init_wrapper()
    def __init__(self, cmd):
        self.info = None # WindowInfo will be set in @init_wrapper
        self.cmd = cmd
        shell = (os.name == 'nt') # True or False
        self.process = subprocess.Popen(cmd.split(),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=shell)


# Snapshot window list
def get_opened_windows():
    args = [None, None]
    if os.name == 'posix':
        wc = gui.connection.WindowConnectionLinux(*args)
    elif os.name == 'nt':
        wc = gui.connection.WindowConnectionWindows(*args)
    else:
        logging.error('Unsupported OS.')
        raise SystemExit
    return wc.get_opened_windows()

# Gets two dictionaries with opened windows,
# compares and returns info for newly opened window.
# Only one new window is allowed.
def get_new_windows_infos(opened_windows_before, opened_windows_after):
    new_windows_infos = []
    for wi in opened_windows_after:
        if wi.wid not in [i.wid for i in opened_windows_before]:
            new_windows_infos.append(wi)
    if len(new_windows_infos) > 1:
        msg = 'Can\'t define slave WID: there is more than one newly opened window.'
        logging.error(msg)
        for wi in new_windows_infos:
            logging.debug(wi.to_string())
        raise SystemExit
    return new_windows_infos

# Keycodes sending to text editor and CGX
@tests.test_wrapper()
def test_sendkeys():
    p = path.Path()
    s = settings.Settings(p)
    xml = os.path.join(p.config, 'sendkeys.xml')

    # Create application
    app = QtWidgets.QApplication([])

    # Create master window
    f = Factory(s)
    f.run_master(xml)

    # Map methods to GUI buttons
    if os.name == 'nt':
        cmd1 = 'notepad.exe'
        # title1 = 'Untitled - Notepad'
    else:
        cmd1 = 'gedit'
        # title1 = 'Untitled Document 1 - gedit'
        if which(cmd1) is None:
            cmd1 = 'kate'
            # title1 = 'Untitled Document 1 - Kate'
            if which(cmd1) is None:
                cmd1 = ''
                msg = 'Neither gedit nor kate is available!'
                logging.error(msg)
                f.mw.b0.setDisabled(True)
    f.mw.b0.clicked.connect(lambda: f.run_slave(cmd1))
    f.mw.b0.clicked.connect(lambda: f.create_connection())
    f.mw.b1.clicked.connect(lambda: gui.cgx.open_inp(p, f, s.start_model, 1))
    f.mw.b2.clicked.connect(lambda: f.connection.post('plot n all'))
    f.mw.b3.clicked.connect(lambda: f.connection.post('plot e all'))
    f.mw.b4.clicked.connect(lambda: f.connection.post(f.mw.customEdit.text()))
    f.mw.b4.clicked.connect(lambda: f.mw.customEdit.clear())
    all_symbols = '!"#$%&\'()*+,-./1234567890:;<=>?@' \
        + 'ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`' \
        + 'abcdefghijklmnopqrstuvwxyz{|}~'
    f.mw.b5.clicked.connect(lambda: f.connection.post(all_symbols))

    # Execute application + exit
    a = app.exec()
    f.kill_slave()

# Run test
if __name__ == '__main__':
    test_sendkeys()
