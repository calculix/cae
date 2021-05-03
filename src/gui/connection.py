#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Classes for keycodes sending from master window to slave.
    mw = master window
    sw = slave window
"""

# Standard modules
import os
import sys
import re
import time
import math
import logging
import webbrowser
if os.name == 'nt':
    import ctypes
    from ctypes import wintypes

# External modules
from PyQt5 import QtWidgets
if os.name == 'posix':
    import Xlib
    from Xlib import display, protocol, X, XK
    from Xlib.ext.xtest import fake_input

# My modules
import window
sys_path = os.path.abspath(__file__)
sys_path = os.path.dirname(sys_path)
sys_path = os.path.join(sys_path, '..')
sys_path = os.path.normpath(sys_path)
sys_path = os.path.realpath(sys_path)
sys.path.insert(0, sys_path)
import clean
import path
import settings
# if os.name == 'nt':
#     from gui.forcefocus import forceFocus

# Common wrapper for WindowConnection.get_wid()
def wid_wrapper(wc):
    def wrap(method):
        def fcn(wc, title):
            start = time.perf_counter() # start time
            wid = None
            while wid is None:
                time.sleep(0.1) # wait for window to start
                wid = method(wc, title)
                if time.perf_counter() - start > 5:
                    msg = 'Can\'t get \'{}\' window.'.format(title)
                    logging.error(msg)
                    wc.log_opened_windows()
                    msg = 'Communication with {} will not work.'
                    logging.warning(msg.format(title))
                    # TODO Make thread and check CGX availability periodically
                    return None
            msg = '{} WID=0x{}'.format(title, hex(wid)[2:].zfill(8))
            logging.debug(msg)
            return wid
        return fcn
    return wrap

# Common wrapper for WindowConnection.post()
def post_wrapper(wc):
    def wrap(method):
        def fcn(wc, cmd):
            if wc.sw.process is not None \
                and wc.sw.process.poll() is None\
                and wc.wid2 is not None:

                if len(cmd):
                    # Drop previously pressed keys
                    if 'cgx' in wc.sw.title:
                        method(wc, '\n')

                    # Post command
                    method(wc, cmd + '\n')

                    # Flush CGX buffer to get continuous output
                    if 'cgx' in wc.sw.title:
                        method(wc, ' ')

                    return
                else:
                    logging.warning('Empty command.')

            else:
                logging.warning('No slave window.')
        return fcn
    return wrap


# Pair master and slave windows together
# to align and enable keycodes sending
# TODO make possible to use QtWidgets.QDialog as master
class WindowConnection:

    def __init__(self, f):
        self.wid1 = None # master window id
        self.wid2 = None # slave window id
        self.opened_windows = {} # {wid:(pid, wname)}

        if f is None:
            return

        # TODO Pass all these arguemtns, not Factory
        self.mw = f.mw # master window
        self.sw = f.sw # slave window
        self.w = f.w
        self.h = f.h
        self.s = f.s # global settings

    def connect(self):
        self.wid1 = self.get_wid(self.mw.windowTitle())
        self.wid2 = self.get_wid(self.sw.title)
        if self.s.align_windows:
            self.align()

    def disconnect(self):
        self.wid1 = None
        self.wid2 = None

    def log_opened_windows(self):
        msg = 'Window list:'
        if not len(self.opened_windows):
            self.get_opened_windows()
        for _id_, (pid, wname) in self.opened_windows.items():
            msg += '\n0x{} {: 6d} {}'\
                .format(hex(_id_)[2:].zfill(8), pid, wname)
        logging.debug(msg)


class WindowConnectionLinux(WindowConnection):

    def __init__(self, f):
        super(WindowConnectionLinux, self).__init__(f)

        self.d = Xlib.display.Display()
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

    # Get window by title and return its ID
    @wid_wrapper(WindowConnection)
    def get_wid(self, title):
        ids = self.root.get_full_property(
            self.d.intern_atom('_NET_CLIENT_LIST'),
            X.AnyPropertyType).value
        wid = None
        for _id_ in ids:
            w = self.d.create_resource_object('window', _id_)
            try:
                wname = w.get_full_property(
                    self.d.intern_atom('_NET_WM_NAME'), X.AnyPropertyType).value.decode()
            except:
                wname = w.get_wm_name()
            regex = '(\S+ - )*' + re.escape(title.lower()) + '( - \S+)*'
            if re.fullmatch(regex, wname.lower()) is not None:
                wid = _id_
                break

            # Firefox
            # if 'new ' + title.lower() in wname.lower():
            #     wid = _id_
            #     break

        return wid

    # Activate window, send message to slave, deactivate
    @post_wrapper(WindowConnection)
    def post(self, cmd):

        # Create X event
        def event(win, e, keycode, case):
            return e(
                time=Xlib.X.CurrentTime,
                root=self.root,
                window=win,
                same_screen=0, child=Xlib.X.NONE,
                root_x=0, root_y=0, event_x=0, event_y=0,
                state=case, # 0=lowercase, 1=uppercase
                detail=keycode)

        # Key press + release
        def sendkey(win, symbol):
            if symbol in self.keyboardMapping:
                c = self.keyboardMapping[symbol][0]
                keysym = XK.string_to_keysym(c)
                code = self.d.keysym_to_keycode(keysym)
                case = self.keyboardMapping[symbol][1]

                e = event(win, Xlib.protocol.event.KeyPress, code, case)
                win.send_event(e, propagate=True)
                e = event(win, Xlib.protocol.event.KeyRelease, code, case)
                win.send_event(e, propagate=True)
            else:
                msg = 'Symbol {} is not supported.'.format(symbol)
                logging.warning(msg)

        win = self.d.create_resource_object('window', self.wid2)
        for symbol in cmd:
            sendkey(win, symbol)
        self.d.sync()

    # Get dictionary with windows id, pid and wname
    def get_opened_windows(self):
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
            self.opened_windows[_id_] = (pid, wname)
        return self.opened_windows
    
    # Key press + release
    def send_hotkey(self, *keys):
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

    # Align CGX and browser windows
    # CAE window is already aligned in __init__()
    def align(self):
        # Align master
        self.mw.setGeometry(0, 0, math.floor(self.w/3), self.h)
        # if self.wid1 is not None:
        #     win = self.d.create_resource_object('window', self.wid1)
        #     win.set_input_focus(Xlib.X.RevertToNone, Xlib.X.CurrentTime)
        #     self.send_hotkey('Super_L', 'Down') # restore window
        #     win.configure(x=0, y=0,
        #         width=math.floor(self.w/3), height=self.h)
        # else:
        #     logging.error('Master WID is None.')

        # if self.wid2 is not None:
        #     self.post('wpos {} {}'\
        #         .format(math.ceil(self.w/3), 0))
        #     self.post('wsize {} {}'\
        #         .format(math.floor(self.w*2/3), self.h))

        # Align slave
        if self.wid2 is not None:
            win = self.d.create_resource_object('window', self.wid2)
            win.set_input_focus(Xlib.X.RevertToNone, Xlib.X.CurrentTime)
            self.send_hotkey('Super_L', 'Down') # restore window
            win.configure(x=math.ceil(self.w/3), y=0,
                width=math.floor(self.w*2/3), height=self.h)
        else:
            logging.error('Slave WID is None.')

        self.d.sync()


class WindowConnectionWindows(WindowConnection):

    def __init__(self, f):
        super(WindowConnectionWindows, self).__init__(f)

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
            'VK_LWIN':  (0x5B, 0),
            'VK_LEFT':  (0x25, 0),
            'VK_UP':    (0x26, 0),
            'VK_RIGHT': (0x27, 0),
            'VK_DOWN':  (0x28, 0),
            }
        for c in '0123456789':
            self.keyboardMapping[c] = (ord(c), 0)
        for c in 'abcdefghijklmnopqrstuvwxyz':
            self.keyboardMapping[c] = (ord(c)-32, 0)
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            self.keyboardMapping[c] = (ord(c), 1)

    # If window found, its ID is returned
    @wid_wrapper(WindowConnection)
    def get_wid(self, title):
        self.wid = None
        WNDENUMPROC = ctypes.WINFUNCTYPE(
            wintypes.BOOL,
            wintypes.HWND,    # _In_ hwnd
            wintypes.LPARAM,) # _In_ lParam

        @WNDENUMPROC
        def enum_proc(hwnd, lParam):
            if ctypes.windll.user32.IsWindowVisible(hwnd):
                length = ctypes.windll.user32.GetWindowTextLengthW(hwnd) + 1
                buff = ctypes.create_unicode_buffer(length)
                ctypes.windll.user32.GetWindowTextW(hwnd, buff, length)
                regex = '(\S+ - )*' + re.escape(title.lower()) + '( - \S+)*'
                if re.fullmatch(regex, buff.value.lower()) is not None:
                    self.wid = hwnd
            return True

        ctypes.windll.user32.EnumWindows(enum_proc, 0)
        return self.wid

    # Activate window, send message to slave, deactivate
    @post_wrapper(WindowConnection)
    def post(self, cmd):

        # Key press (0) + release (2)
        def sendkey(symbol):
            if symbol in self.keyboardMapping:
                code = self.keyboardMapping[symbol][0]
                case = self.keyboardMapping[symbol][1]
                if case == 1: # press Shift
                    ctypes.windll.user32.keybd_event(0x10, 0, 0, 0)
                ctypes.windll.user32.keybd_event(code, 0, 0, 0)
                ctypes.windll.user32.keybd_event(code, 0, 2, 0)
                if case == 1: # release Shift
                    ctypes.windll.user32.keybd_event(0x10, 0, 2, 0)
            else:
                msg = 'Symbol {} is not supported.'.format(symbol)
                logging.warning(msg)

        ctypes.windll.user32.SetForegroundWindow(self.wid2)
        for symbol in cmd:
            time.sleep(0.001) # BUG doesn't work otherwise
            sendkey(symbol)
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
                self.opened_windows[hwnd] = (str(pid)[8:-1], buff.value)
                # logging.debug('0x{} {:>6s} {}'\
                #     .format(hex(hwnd)[2:].zfill(8), str(pid)[8:-1], buff.value))
            return True

        ctypes.windll.user32.EnumWindows(enum_proc, 0)
        return self.opened_windows

    # Align CGX and browser windows
    # CAE window is already aligned
    # NOTE CGX 'wpos' and 'wsize' works worse than 'ctypes'
    def align(self):

        # Align master
        self.mw.setGeometry(0, 0, math.floor(self.w/3), self.h)
        # ctypes.windll.user32.SetProcessDPIAware() # account for scaling
        # if self.wid1 is not None:
        #     ok = forceFocus(self.wid1)
        #     logging.debug(str(ok))

        #     ctypes.windll.user32.ShowWindow(self.wid1, 9) # restore window
        #     time.sleep(0.3)
        #     ctypes.windll.user32.MoveWindow(self.wid1, 0, 0,
        #         math.floor(self.w/3), self.h, True)
        # else:
        #     logging.error('Master WID is None.')

        # Align slave
        if self.wid2 is not None:
            ctypes.windll.user32.ShowWindow(self.wid2, 9) # restore window
            time.sleep(0.3)
            ctypes.windll.user32.MoveWindow(self.wid2,
                math.ceil(self.w/3), 0,
                math.floor(self.w*2/3), self.h, True)
        else:
            logging.error('Slave WID is None.')


# Print current windows list
def test1():
    if os.name == 'posix':
        wc = WindowConnectionLinux(None)
    elif os.name == 'nt':
        wc = WindowConnectionWindows(None)
    else:
        logging.error('Unsupported OS.')
        raise SystemExit

    wc.log_opened_windows()

# Open web-browser, get its wid and align
def test2():
    app = QtWidgets.QApplication(sys.argv) # create application
    p = path.Path()
    s = settings.Settings(p)

    url = 'file://' + p.doc + '/NODE.html'
    webbrowser.open(url, 1)
    browsers = webbrowser._tryorder
    if browsers is None:
        print('ERROR! No web browsers!')
        return

    mw = window.MasterWindow(p, s)
    for slave_title in reversed(browsers):
        slave_title = slave_title.replace('-browser', '')
        mw.create_connection(2, slave_title)
        wid2 = mw.connections[2].wid2
        msg = '{} {}'.format(slave_title, wid2)
        print(msg)
        if wid2 is not None and s.align_windows:
            mw.connections[2].align()
            break

    app.exec() # execute application

# Run test
if __name__ == '__main__':
    start = time.perf_counter() # start time
    logging.basicConfig(level=logging.NOTSET, format='%(message)s')
    clean.screen()

    # test1()
    test2()

    clean.cache()
    print('\nTotal {:.1e} seconds.\n'\
        .format(time.perf_counter()-start)) # spent time
