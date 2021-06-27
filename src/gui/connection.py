#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

TODO Invent some test for the module.

Classes for keycodes sending from master window to slave.
It is a layer between system libraries Xlib/ctypes and GUI.

mw = master window
sw = slave window
"""

# Standard modules
import os
import sys
import time
import math
import logging
try:
    import ctypes
    from ctypes import wintypes
except:
    from Xlib import display, protocol, X, XK
    from Xlib.ext.xtest import fake_input

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


def align_master(w):
    """Align master window - w could be QMainWindow or QDialog.
    frameSize() is window size, geometry() is the central widget size,
    width and height - desirable window size.
    """
    ag = QtWidgets.QApplication.primaryScreen().availableGeometry()
    w.move(ag.left(), ag.top())
    size = QtWidgets.QApplication.primaryScreen().availableSize()
    width = math.floor(size.width() / 3)
    height = size.height()
    dw = w.frameSize().width() - w.geometry().width() # delta width
    dh = w.frameSize().height() - w.geometry().height() # delta height
    w.resize(width-dw, height-dh)


class WindowInfo:

    def __init__(self, wid, pid, wname):
        self.wid = wid # window identificator
        self.hex_wid = '0x{}'.format(hex(wid)[2:].zfill(8))
        self.pid = pid # process ID
        self.wname = wname # window title

    def to_string(self):
        return '{}  {: 8d}  {}'\
            .format(self.hex_wid, self.pid, self.wname)


def post_wrapper(wc):
    """Common wrapper for WindowConnection.post()."""
    def wrap(method):
        def fcn(wc, cmd):
            if wc.get_slave_process() is not None \
                and wc.get_slave_process().poll() is None\
                and wc.wid2 is not None:

                if len(cmd):
                    # Drop previously pressed keys
                    if 'CalculiX GraphiX' in wc.get_slave_title():
                        method(wc, '\n')

                    # Post command
                    method(wc, cmd + '\n')

                    # Flush CGX buffer to get continuous output
                    if 'CalculiX GraphiX' in wc.get_slave_title():
                        method(wc, ' ')

                    return
                else:
                    logging.warning('Empty command.')

            else:
                logging.warning('No slave window.')
        return fcn
    return wrap


class WindowConnection:
    """Pair master and slave windows together
    to align and enable keycodes sending.
    """

    def __init__(self, mw, sw):
        self.opened_windows = [] # list of WindowInfo objects
        screen = QtWidgets.QApplication.primaryScreen()
        self.w = screen.availableSize().width()
        self.h = screen.availableSize().height()

        if mw is None or sw is None:
            return

        self.mw = mw # master window
        self.sw = sw # slave window
        # TODO Do we need self.wid?
        self.wid1 = self.get_master_wid() # master window id
        self.wid2 = self.get_slave_wid() # slave window id

    def log_opened_windows(self):
        msg = 'Window list:'
        if not len(self.opened_windows):
            self.get_opened_windows()
        for wi in self.opened_windows:
            msg += wi.to_string()
        logging.debug(msg)

    def get_slave_title(self):
        try:
            return self.sw.info.wname
        except:
            logging.error('Can not get slave title.')
            return ""

    def get_slave_process(self):
        try:
            return self.sw.process
        except:
            logging.warning('Can not get slave process.')
            return None

    def get_slave_wid(self):
        try:
            return self.sw.info.wid
        except:
            logging.error('Can not get slave wid.')
            return None

    def get_master_wid(self):
        try:
            return self.mw.info.wid
        except:
            logging.error('Can not get master wid.')
            return None


class WindowConnectionLinux(WindowConnection):

    def __init__(self, mw, sw):
        super().__init__(mw, sw)

        self.d = display.Display()
        self.screen = self.d.screen()
        self.root = self.screen.root

        # 0:lowercase, 1:shifted
        self.keyboardMapping = {
            '\t':('Tab', 0),
            '\n':('Return', 0),
            '\r':('Return', 0),
            '\b':('BackSpace', 0),
            ' ': ('space', 0),
            '!': ('exclam', 1),
            '#': ('numbersign', 1),
            '%': ('percent', 1),
            '$': ('dollar', 1),
            '&': ('ampersand', 1),
            '"': ('quotedbl', 1),
            "'": ('apostrophe', 0),
            '(': ('parenleft', 1),
            ')': ('parenright', 1),
            '*': ('asterisk', 1),
            '=': ('equal', 0),
            '+': ('plus', 1),
            ',': ('comma', 0),
            '-': ('minus', 0),
            '.': ('period', 0),
            '/': ('slash', 0),
            ':': ('colon', 1),
            ';': ('semicolon', 0),
            '<': ('less', 0),
            '>': ('greater', 1),
            '?': ('question', 1),
            '@': ('at', 1),
            '[': ('bracketleft', 0),
            ']': ('bracketright', 0),
            '^': ('asciicircum', 1),
            '_': ('underscore', 1),
            '`': ('grave', 0),
            '{': ('braceleft', 1),
            '\\':('backslash', 0),
            '|': ('bar', 1),
            '}': ('braceright', 1),
            '~': ('asciitilde', 1) }
        KEY_NAMES = [
            'Shift_L', 'Control_L', 'Alt_L', 'Pause', 'Caps_Lock',
            'Escape', 'Page_Up', 'Page_Down', 'End', 'Home',
            'Left', 'Up', 'Right', 'Down', 'Print',
            'Insert', 'Delete', 'Help', 'Super_L', 'Super_R',
            'KP_0', 'KP_1', 'KP_2', 'KP_3', 'KP_4', 'KP_5', 'KP_6',
            'KP_7', 'KP_8', 'KP_9', 'KP_Multiply', 'KP_Add',
            'KP_Separator', 'KP_Subtract', 'KP_Decimal', 'KP_Divide',
            'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12',
            'Num_Lock', 'Scroll_Lock', 'Shift_R', 'Control_R', 'Alt_R' ]
        for c in KEY_NAMES:
            self.keyboardMapping[c] = (c, 0)
        for c in '1234567890abcdefghijklmnopqrstuvwxyz':
            self.keyboardMapping[c] = (c, 0)
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            self.keyboardMapping[c] = (c, 1)

    @post_wrapper(WindowConnection)
    def post(self, cmd):
        """Activate window, send message to slave, deactivate."""

        def event(win, e, keycode, case):
            """Create X event."""
            return e(
                time=X.CurrentTime,
                root=self.root,
                window=win,
                same_screen=0, child=X.NONE,
                root_x=0, root_y=0, event_x=0, event_y=0,
                state=case, # 0=lowercase, 1=uppercase
                detail=keycode)

        def sendkey(win, symbol):
            """Key press + release."""
            if symbol in self.keyboardMapping:
                c = self.keyboardMapping[symbol][0]
                keysym = XK.string_to_keysym(c)
                code = self.d.keysym_to_keycode(keysym)
                case = self.keyboardMapping[symbol][1]

                e = event(win, protocol.event.KeyPress, code, case)
                win.send_event(e, propagate=True)
                e = event(win, protocol.event.KeyRelease, code, case)
                win.send_event(e, propagate=True)
            else:
                msg = 'Symbol {} is not supported.'.format(symbol)
                logging.warning(msg)

        win = self.d.create_resource_object('window', self.wid2)
        for symbol in cmd:
            sendkey(win, symbol)
        self.d.sync()

    def get_opened_windows(self):
        """Get dictionary with windows id, pid and wname."""
        self.opened_windows.clear()
        ids = self.root.get_full_property(
            self.d.intern_atom('_NET_CLIENT_LIST'),
            X.AnyPropertyType).value
        for _id_ in ids:
            w = self.d.create_resource_object('window', _id_)
            try:
                pid = w.get_full_property(
                    self.d.intern_atom('_NET_WM_PID'), X.AnyPropertyType).value[0]
                wname = w.get_full_property(
                    self.d.intern_atom('_NET_WM_NAME'), X.AnyPropertyType).value.decode()
            except:
                pid = 0
                wname = w.get_wm_name()
            wi = WindowInfo(_id_, pid, wname)
            self.opened_windows.append(wi)
        return self.opened_windows

    def send_hotkey(self, *keys):
        """Key press + release."""
        for key in keys:
            if key not in self.keyboardMapping:
                logging.warning('WARNING! Key {} is not supported.'.format(key))
                return
        for key in keys:
            c = self.keyboardMapping[key][0]
            keysym = XK.string_to_keysym(c)
            code = self.d.keysym_to_keycode(keysym)
            fake_input(self.d, X.KeyPress, code)
        for key in reversed(keys):
            c = self.keyboardMapping[key][0]
            keysym = XK.string_to_keysym(c)
            code = self.d.keysym_to_keycode(keysym)
            fake_input(self.d, X.KeyRelease, code)
        time.sleep(0.3)
        self.d.sync()

    def align_windows(self):
        """Align master and slave windows."""

        # Align master
        self.mw.setGeometry(0, 0, math.floor(self.w/3), self.h)

        # Align slave
        if self.wid2 is not None:
            win = self.d.create_resource_object('window', self.wid2)
            win.set_input_focus(X.RevertToNone, X.CurrentTime)
            self.send_hotkey('Super_L', 'Down') # restore window
            time.sleep(0.3)
            win.configure(x=math.ceil(self.w/3), y=0,
                width=math.floor(self.w*2/3), height=self.h)
        else:
            logging.error('Nothing to align - slave WID is None.')

        self.d.sync()


class WindowConnectionWindows(WindowConnection):

    def __init__(self, mw, sw):
        super().__init__(mw, sw)

        # 0:lowercase, 1:shifted
        self.keyboardMapping = {
            '\t':(9, 0),
            '\n':(13, 0),
            ' ': (32, 0),
            '!': (49, 1),
            '@': (50, 1),
            '#': (51, 1),
            '$': (52, 1),
            '%': (53, 1),
            '^': (54, 1),
            '&': (55, 1),
            '*': (56, 1),
            '(': (57, 1),
            ')': (48, 1),
            ';': (0xBA, 0),
            ':': (0xBA, 1),
            '=': (0xBB, 0),
            '+': (0xBB, 1),
            ',': (0xBC, 0),
            '<': (0xBC, 1),
            '-': (0xBD, 0),
            '_': (0xBD, 1),
            '.': (0xBE, 0),
            '>': (0xBE, 1),
            '/': (0xBF, 0),
            '?': (0xBF, 1),
            '`': (0xC0, 0),
            '~': (0xC0, 1),
            '[': (0xDB, 0),
            '{': (0xDB, 1),
            '\\':(0xDC, 0),
            '|': (0xDC, 1),
            ']': (0xDD, 0),
            '}': (0xDD, 1),
            '\'':(0xDE, 0),
            '\"':(0xDE, 1),
            'F4':(0x73, 0),
            'Alt_L': (0x12, 0),
            # 'VK_LWIN': (0x5B, 0),
            # 'VK_LEFT': (0x25, 0),
            # 'VK_UP':   (0x26, 0),
            # 'VK_RIGHT':(0x27, 0),
            # 'VK_DOWN': (0x28, 0),
            }
        for c in '0123456789':
            self.keyboardMapping[c] = (ord(c), 0)
        for c in 'abcdefghijklmnopqrstuvwxyz':
            self.keyboardMapping[c] = (ord(c)-32, 0)
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            self.keyboardMapping[c] = (ord(c), 1)

    @post_wrapper(WindowConnection)
    def post(self, cmd):
        """Activate window, send message to slave, deactivate."""

        def sendkey(symbol):
            """Key press (0) + release (2)."""
            if symbol not in self.keyboardMapping:
                msg = 'Symbol {} is not supported.'.format(symbol)
                logging.warning(msg)
                return
            code = self.keyboardMapping[symbol][0]
            case = self.keyboardMapping[symbol][1]
            if case == 1: # press Shift
                ctypes.windll.user32.keybd_event(0x10, 0, 0, 0)
            ctypes.windll.user32.keybd_event(code, 0, 0, 0) # key press
            ctypes.windll.user32.keybd_event(code, 0, 2, 0) # key release
            if case == 1: # release Shift
                ctypes.windll.user32.keybd_event(0x10, 0, 2, 0)

        # TODO Check if it's needed:
        ctypes.windll.user32.SetForegroundWindow(self.wid2)
        for symbol in cmd:
            time.sleep(0.001) # BUG doesn't work otherwise
            sendkey(symbol)
        # TODO Check if it's needed:
        ctypes.windll.user32.SetForegroundWindow(self.wid1)

    def get_opened_windows(self):
        self.opened_windows.clear()
        WNDENUMPROC = ctypes.WINFUNCTYPE(
            wintypes.BOOL,
            wintypes.HWND,    # _In_ hwnd
            wintypes.LPARAM,) # _In_ lParam

        @WNDENUMPROC
        def enum_proc(hwnd, lParam):
            if ctypes.windll.user32.IsWindowVisible(hwnd):
                pid = wintypes.DWORD()
                length = ctypes.windll.user32.GetWindowTextLengthW(hwnd) + 1
                buff = ctypes.create_unicode_buffer(length)
                ctypes.windll.user32.GetWindowTextW(hwnd, buff, length)
                wi = WindowInfo(hwnd, int(str(pid)[8:-1]), buff.value)
                self.opened_windows.append(wi)
            return True

        ctypes.windll.user32.EnumWindows(enum_proc, 0)
        return self.opened_windows

    # TODO Test Alt+F4 in Windows
    def send_hotkey(self, *keys):
        """Key press + release."""
        for key in keys:
            if key not in self.keyboardMapping:
                logging.warning('WARNING! Key {} is not supported.'.format(key))
                return
        for key in keys:
            code = self.keyboardMapping[key][0]
            ctypes.windll.user32.keybd_event(code, 0, 0, 0) # key press
        for key in reversed(keys):
            code = self.keyboardMapping[key][0]
            ctypes.windll.user32.keybd_event(code, 0, 2, 0) # key release
        time.sleep(0.3)

    def align_windows(self):
        """Align master and slave windows.
        NOTE CGX 'wpos' and 'wsize' works worse than 'ctypes'.
        TODO Split on align_master() and align_slave()
        TODO Move align_master() to the WindowConnection
        """
        # Align master
        align_master(self.mw)

        # Align slave
        ag = QtWidgets.QApplication.primaryScreen().availableGeometry()
        size = QtWidgets.QApplication.primaryScreen().availableSize()
        width = math.floor(size.width() / 3)
        height = size.height()
        if self.wid2 is not None:
            ctypes.windll.user32.ShowWindow(self.wid2, 9) # restore window
            time.sleep(0.3)
            ctypes.windll.user32.MoveWindow(self.wid2,
                ag.left() + width + 1, ag.top(),
                size.width() - width, height, True)
        else:
            logging.error('Nothing to align - slave WID is None.')


@tests.test_wrapper()
def test():
    pass


if __name__ == '__main__':
    test() # run test
