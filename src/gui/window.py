#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Class to work with master/slave windows. Keyword editor is always
master, CGX - slave. Keyword dialog is master, web browser - slave.
Window IDs (WIDs) are needed to align application windows.

For documentation on used functions refer to:
http://python-xlib.sourceforge.net/doc/html/python-xlib_21.html
https://github.com/python-xlib/python-xlib
https://github.com/asweigart/pyautogui

TODO
File "/media/ihor/WORK/Calculix/cae/src/gui/cgx.py", line 33, in read_fbd_file
factory.connection.post('read ' + file_name)
AttributeError: 'NoneType' object has no attribute 'post'
"""

# Standard modules
import os
import sys
import time
import logging
import subprocess
import webbrowser
import traceback
from shutil import which

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
import path
import settings
import tests
import gui.connection
import log
import gui.stdout


def init_wrapper():
    """Common wrapper for MasterWindow/SlaveWindow __init__() method."""
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
                msg = 'Starting new window:\n' + wi.to_string()
                logging.debug(msg)
                self.info = wi

        return fcn
    return wrap


class Factory:
    """Window factory class - used to create master and slave windows
    and bind/connect them together. Logging facility as killing methods
    for slave window are also maintained here.
    """

    def __init__(self):
        self.mw = None # master window
        self.sw = None # slave window
        args = [None, None]
        if os.name == 'nt':
            wc1 = gui.connection.WindowConnectionWindows(*args)
            wc2 = gui.connection.WindowConnectionWindows(*args)
        else:
            wc1 = gui.connection.WindowConnectionLinux(*args)
            wc2 = gui.connection.WindowConnectionLinux(*args)
        self.connection = None

    def run_master(self, xml=path.p.main_xml):
        """Draw window form described in xml UI file."""
        self.mw = MasterWindow(xml)

    def run_master_dialog(self, k, i):
        """Generate dialog window and show it."""
        args = [k, i]
        import gui.dialog
        self.mw = gui.dialog.KeywordDialog(args) # one argument for init_wrapper()

        # Actions
        self.mw.buttonBox.accepted.connect(self.mw.ok)
        self.mw.buttonBox.button(QtWidgets.QDialogButtonBox.Reset).clicked.connect(self.mw.reset)
        self.mw.buttonBox.helpRequested.connect(self.open_help)

        if settings.s.align_windows:
            gui.connection.align_master(self.mw)

        if settings.s.show_help:
            self.open_help(False)
        d = self.mw.exec()
        self.kill_slave()
        return d

    def open_help(self, click=True):
        """Open HTML help page in a default web browser."""
        if settings.s.default_web_browser == 'internal':
            self.mw.show_hide_internal_help(click)
        else:
            if not os.path.isfile(self.mw.url):
                msg = 'Help page does not exist:\n{}'.format(
                    os.path.relpath(self.mw.url, start=path.p.app_home_dir))
                logging.error(msg)
                return
            logging.info('Going to\n' + self.mw.url)
            wb = webbrowser.get()
            if 'nt' in os.name:
                """NOTE In Windows webbrowser can not
                return path to the browser executable."""
                wb.open_new(self.mw.url)
            else:
                cmd = wb.name + ' --new-window ' + self.mw.url
                self.run_slave(cmd)

    def run_slave(self, cmd):
        """Close opened slave window and open a new one."""
        self.kill_slave()
        self.sw = SlaveWindow(cmd)
        self.create_connection()

        # Start stdout reading and logging thread
        if path.p.path_cgx in cmd and hasattr(self.mw, 'textEdit'):
            html = self.mw.textEdit.toHtml()
            self.mw.textEdit.clear()
            self.mw.textEdit.setHtml(html)

            # Start stdout reading and logging thread
            if self.connection is not None:
                args = [self.sw.process.stdout,
                        'read_cgx_stdout',\
                        True, self.connection]
                gui.stdout.start_cgx_reader(*args)
            else:
                msg = 'Window connection not established.'
                logging.error(msg)

    def kill_slave(self):
        """Kill all slave processes."""
        if self.sw is None:
            return
        if self.sw.process is None:
            return
        if self.connection is None:
            return
        if path.p.path_cgx in self.sw.cmd:
            gui.stdout.stop_readers()

        # First try to close window
        # TODO Test Alt+F4 in Linux
        # TODO Firefox process is being killed correctly, but window remains
        # NOTE Dangerous method - closes everything
        if 'nt' in os.name:
            import ctypes
            ctypes.windll.user32.SetForegroundWindow(self.sw.info.wid)
            self.connection.send_hotkey('Alt_L', 'F4')
        # else:
        #     from Xlib import display, X
        #     self.d = display.Display()
        #     win = self.d.create_resource_object('window', self.sw.info.wid)
        #     win.set_input_focus(X.RevertToNone, X.CurrentTime)
        #     self.d.sync()
        #     self.connection.send_hotkey('Alt_L', 'F4')

        # Then kill its process
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


    def create_connection(self):
        """Connect master and slave windows and align them."""
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

        if settings.s.align_windows:
            self.connection.align_windows()


factory = Factory()


class MasterWindow(QtWidgets.QMainWindow):
    """Left window - main window of the app."""

    @init_wrapper()
    def __init__(self, xml=path.p.main_xml):
        """Load form and show the window."""
        self.info = None # WindowInfo will be set in @init_wrapper

        # Load UI form - produces huge amount of redundant debug logs
        logging.disable() # switch off logging
        super().__init__() # create main window
        uic.loadUi(xml, self) # load form
        logging.disable(logging.NOTSET) # switch on logging

        # Handler to show logs in the CAE textEdit
        if hasattr(self, 'textEdit'): # skip for test()
            log.add_text_handler(self.textEdit)

        self.show()

    def help(self, url):
        """Open links from the Help menu."""
        if webbrowser.open(url, new=2):
            logging.info('Going to\n' + url)
        else:
            logging.warning('Can not open url\n' + url)


class SlaveWindow:
    """Right application window - CGX, text editor,
    web browser. Starts as a separate process.
    """

    @init_wrapper()
    def __init__(self, cmd):
        """Run slave process without terminal window, set window info."""
        self.info = None # WindowInfo will be set in @init_wrapper
        self.cmd = cmd
        shell = (os.name == 'nt') # True or False
        self.process = subprocess.Popen(cmd.split(),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=shell)


def get_opened_windows():
    """Snapshot window list."""
    args = [None, None]
    if os.name == 'posix':
        wc = gui.connection.WindowConnectionLinux(*args)
    elif os.name == 'nt':
        wc = gui.connection.WindowConnectionWindows(*args)
    else:
        logging.error('Unsupported OS.')
        raise SystemExit
    return wc.get_opened_windows()


def get_new_windows_infos(opened_windows_before, opened_windows_after):
    """Gets two dictionaries with opened windows,
    compares and returns info for newly opened window.
    Only one new window is allowed.
    """
    new_windows_infos = []
    for wi in opened_windows_after:
        if wi.wid not in [i.wid for i in opened_windows_before]:
            new_windows_infos.append(wi)
    if len(new_windows_infos) > 1:
        msg = 'Can\'t connect to the slave.' \
            + '\nThere is more than one newly opened window.' \
            + '\nPlease, reopen it manually.'
        logging.error(msg)
        QtWidgets.QMessageBox.critical(None, 'Error', msg)
        for wi in new_windows_infos:
            logging.debug(wi.to_string())
    return new_windows_infos


@tests.test_wrapper()
def test():
    """Test sendkeys:
    keycodes sending to text editor and CGX.
    """
    import gui.cgx

    # Create application
    app = QtWidgets.QApplication([])

    # Create master window
    global factory
    xml = os.path.join(path.p.config, 'sendkeys.xml')
    factory.run_master(xml)

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
                factory.mw.b0.setDisabled(True)
    factory.mw.b0.clicked.connect(lambda: factory.run_slave(cmd1))
    factory.mw.b1.clicked.connect(lambda: gui.cgx.open_inp(settings.s.start_model, 1))
    factory.mw.b2.clicked.connect(lambda: factory.connection.post('plot n all'))
    factory.mw.b3.clicked.connect(lambda: factory.connection.post('plot e all'))
    factory.mw.b4.clicked.connect(lambda: factory.connection.post(factory.mw.customEdit.text()))
    factory.mw.b4.clicked.connect(lambda: factory.mw.customEdit.clear())
    all_symbols = '!"#$%&\'()*+,-./1234567890:;<=>?@' \
        + 'ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`' \
        + 'abcdefghijklmnopqrstuvwxyz{|}~'
    factory.mw.b5.clicked.connect(lambda: factory.connection.post(all_symbols))

    # Execute application + exit
    a = app.exec()
    factory.kill_slave()


if __name__ == '__main__':
    test()
