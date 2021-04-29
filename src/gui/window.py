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


""" Window factory class - used to create master and slave windows
and bind/connect therm together. Logging facility as killing methods
for slave window are also maintained here. """
class Factory:

    def __init__(self, p, s):
        self.p = p # global paths
        self.s = s # global settings
        self.mw = MasterWindow(self.p.main_xml)
        self.sw = None # slave window
        self.stdout_readers = [] # for slave window
        wc1 = gui.connection.WindowConnection(None)
        wc2 = gui.connection.WindowConnection(None)
        self.connections = {1:wc1, 2:wc2}
        self.w = QtWidgets.QDesktopWidget().availableGeometry().width()
        self.h = QtWidgets.QDesktopWidget().availableGeometry().height()

        # Caller fuction name: cgx.open_inp | cgx.open_frd | other
        self.mode = None # TODO test it

    # Run slave process
    # Close opened CGX (if any) and open a new one
    # Get window ID, align windows and post to CGX
    def run_slave(self, cmd, title):
        self.kill_slave()
        self.sw = SlaveWindow(cmd, title)
        self.sw.run()

        # Create connection between master and slave windows
        self.create_connection(1)

        # Start stdout reading and logging thread
        self.start_stdout_reader('read_cgx_stdout')

        # Caller fuction name:
        # cgx.open_inp | cgx.open_frd | other
        self.mode = inspect.stack()[1].function

    # Kill all slave processes
    def kill_slave(self):
        if self.sw is None:
            return
        if self.sw.process is None:
            return

        self.stop_stdout_readers()

        self.connections[1].wid2 = None
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
                .format(self.sw.title, self.sw.process.pid)
            logging.warning(msg)
        else:
            msg = 'Killed {}, PID={}.'\
                .format(self.sw.title, self.sw.process.pid)
            logging.debug(msg)
            self.sw.process = None

    # Start stdout reading and logging thread
    # TODO BUG It logs prev command, not the last one
    def start_stdout_reader(self, aim):
        sr = gui.log.CgxStdoutReader(self.sw.process.stdout,
            aim, self)
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

    # Create connection between master and slave windows
    def create_connection(self, _id):
        if os.name == 'nt':
            wc = gui.connection.WindowConnectionWindows(self)
        elif os.name == 'posix':
            wc = gui.connection.WindowConnectionLinux(self)
        else:
            msg = 'Sorry, {} OS is not supported.'.format(os.name)
            raise SystemExit(msg) # the best way to exit
        self.connections[_id] = wc
        wc.connect()


# Left application's window - CAE
class MasterWindow(QtWidgets.QMainWindow):

    def __init__(self, xml):
        self.xml = xml

    # Load for and show the window
    def run(self):
        QtWidgets.QMainWindow.__init__(self) # create main window
        uic.loadUi(self.xml, self) # load form

        # Handler to show logs in the CAE's textEdit
        if hasattr(self, 'textEdit'): # skip for test_sendkeys
            gui.log.add_text_handler(self.textEdit)
        
        self.show()

    # Open links from the Help menu
    def help(self, url):
        logging.info('Going to\n' + url)
        if not webbrowser.open(url, new=2):
            logging.warning('Can\'t open url: ' + url)


# Right application's window - CGX, text editor, web browser
class SlaveWindow:

    # Run command - open new slave window - text editor or CGX
    def __init__(self, cmd, title):
        self.cmd = cmd
        self.title = title

        if os.name == 'nt':
            self.shell = True
        elif os.name == 'posix':
            self.shell = False
        else:
            msg = 'Sorry, {} OS is not supported.'.format(os.name)
            raise SystemExit(msg) # the best way to exit

    # Run slave process without terminal/cmd window
    def run(self):
        self.process = subprocess.Popen(self.cmd.split(),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=self.shell)
        msg = '{} PID={}'.format(self.title, self.process.pid)
        logging.debug(msg)


# Keycodes sending to text editor and CGX
def test_sendkeys():

    # Create application
    app = QtWidgets.QApplication([])
    p = path.Path()
    p.main_xml = os.path.join(p.config, 'sendkeys.xml')
    s = settings.Settings(p)

    # Create master window
    f = Factory(p, s)
    f.mw.run()

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
                f.mw.b0.setDisabled(True)
    f.mw.b0.clicked.connect(lambda: f.run_slave(cmd1, title1))
    f.mw.b1.clicked.connect(lambda: gui.cgx.open_inp(f, s.start_model, 1))
    f.mw.b2.clicked.connect(lambda: f.connections[1].post('plot n all'))
    f.mw.b3.clicked.connect(lambda: f.connections[1].post('plot e all'))
    f.mw.b4.clicked.connect(lambda: f.connections[1].post(f.mw.customEdit.text()))
    f.mw.b4.clicked.connect(lambda: f.mw.customEdit.clear())
    all_symbols = '!"#$%&\'()*+,-./1234567890:;<=>?@' \
        + 'ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`' \
        + 'abcdefghijklmnopqrstuvwxyz{|}~'
    f.mw.b5.clicked.connect(lambda: f.connections[1].post(all_symbols))

    # Execute application + exit
    a = app.exec()
    f.kill_slave()

# Run test
if __name__ == '__main__':
    clean.screen()
    test_sendkeys()
    clean.cache()