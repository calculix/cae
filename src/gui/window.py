#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, July 2020
Distributed under GNU General Public License v3.0

Main window class. Here also we keep links to another
app windows like CGX and web browser. Those links (WIDs)
are needed to align application windows.

For documentation on used functions refer to:
http://python-xlib.sourceforge.net/doc/html/python-xlib_21.html
https://github.com/python-xlib/python-xlib
https://github.com/asweigart/pyautogui """

# TODO Allow app to run even without wids

# Standard modules
import os
import sys
import time
import logging
import subprocess
import inspect
import math
import webbrowser
if 'nt' in os.name:
    import ctypes
    from ctypes import wintypes

# External modules
from PyQt5 import QtWidgets, uic
if 'posix' in os.name:
    try:
        import Xlib
        from Xlib import display, protocol, X, XK
        from Xlib.ext.xtest import fake_input
    except:
        msg = 'Please, install Xlib with command:\n' \
            + 'pip3 install Xlib'
        sys.exit(msg)

# My modules
import gui


# Main window
class Window(QtWidgets.QMainWindow):

    # Create main window
    """
    p - Path
    s - Settings
    """
    def __init__(self, p, s):
        self.p = p
        self.s = s
        QtWidgets.QMainWindow.__init__(self) # create main window
        uic.loadUi(p.cae_xml, self) # load form
        self.size = QtWidgets.QDesktopWidget().availableGeometry()

        # Handler to show logs in the CAE's textEdit
        gui.log.add_text_handler(self.textEdit)

        # Window ID - to pass keycodes to
        self.wid1 = None # cae window
        self.wid2 = None # cgx window
        self.wid3 = None # keyword dialog
        self.wid4 = None # help (web browser)
        self.cgx_process = None # running CGX process to send commands to
        self.keyboardMapping = None

        # INP | FRD - corresponds to opened file
        self.mode = None

    # Close opened CGX (if any) and open a new one
    # Get window ID, align windows and post to CGX
    def run_cgx(self, cmd):
        if not os.path.isfile(self.s.path_cgx):
            msg = 'Wrong path to CGX:\n{}\n' \
                + 'Configure it in File->Settings.'
            logging.error(msg.format(self.s.path_cgx))
            return

        gui.cgx.kill(self)
        self.cgx_process = subprocess.Popen(cmd.split(),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT)
        logging.debug('CGX PID={}'.format(self.cgx_process.pid))
        self.wid2 = self.get_wid('CalculiX GraphiX') # could be None
        
        # Start stdout reading and logging thread
        csr = gui.log.CgxStdoutReader(
            self.cgx_process.stdout, 'read_cgx_stdout_', self)
        csr.start()

        if self.s.align_windows:
            self.align()

        # Read config to align model to iso view
        file_name = os.path.join(self.p.config, 'iso.fbd')
        if os.path.isfile(file_name):
            self.post('read ' + file_name)
        else:
            logging.error('No config file iso.fbd')

        # Read config to register additional colors
        # Those colors are needed to paint sets and surfaces
        file_name = os.path.join(self.p.config, 'colors.fbd')
        if os.path.isfile(file_name):
            self.post('read ' + file_name)
        else:
            logging.error('No config file colors.fbd')

        # Caller fuction name: cgx_inp | cgx_frd
        self.mode = inspect.stack()[1].function

    # Open links from the Help menu
    def help(self, url):
        logging.info('Going to\n' + url)
        if not webbrowser.open(url, new=2):
            logging.warning('Can\'t open url: ' + url)


# Common wrapper for Window.get_wid()
def wid_wrapper(w):
    def wrap(method):
        def fcn(w, title):
            start = time.perf_counter() # start time
            wid = None
            while wid is None:
                time.sleep(0.1) # wait for window to start
                wid = method(w, title)
                if time.perf_counter() - start > 5:
                    msg = 'Can\'t get \'{}\' window.'.format(title)
                    logging.error(msg)
                    w.log_window_list()
                    sys.exit(msg)
            msg = '{} WID=0x{}'.format(title, hex(wid)[2:].zfill(8))
            logging.debug(msg)
            return wid
        return fcn
    return wrap


# Common wrapper for Window.post()
def post_wrapper(w):
    def wrap(method):
        def fcn(w, cmd):
            if w.cgx_process is not None \
                and w.cgx_process.poll() is None\
                and w.wid2 is not None:

                # Drop previously pressed keys
                method(w, '\n')

                # Post command
                method(w, cmd + '\n')

                # Flush CGX buffer to get continuous output
                method(w, ' ')

                return
        return fcn
    return wrap


class Linux_window(Window):

    def initialize(self):
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

        self.wid1 = self.get_wid('CalculiX CAE')

    # Get window by title and return its ID
    @wid_wrapper(Window)
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
            # TODO It also could be, for example, a browser window
            # with opened web page about CalculiX CAE
            if title.lower() in wname.lower():
                wid = _id_
                break
        return wid

    def log_window_list(self):
        ids = self.root.get_full_property(
            self.d.intern_atom('_NET_CLIENT_LIST'),
            X.AnyPropertyType).value
        msg = 'Window list:'
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
            msg += '\n0x{} {: 6d} {}'\
                .format(hex(_id_)[2:].zfill(8), pid, wname)
        logging.debug(msg)

    # Activate window, send message to CGX, deactivate
    @post_wrapper(Window)
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
                logging.warning('WARNING! Symbol {} is not supported.'.format(symbol))

        win = self.d.create_resource_object('window', self.wid2)
        for symbol in cmd:
            sendkey(win, symbol)
        self.d.sync()

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
        self.d.sync()

    # Align CGX and browser windows
    # CAE window is already aligned in __init__()
    def align(self):
        for wid in [self.wid1, self.wid3]: # align CAE and keyword dialog
            if wid is not None:
                win = self.d.create_resource_object('window', wid)
                win.set_input_focus(Xlib.X.RevertToNone, Xlib.X.CurrentTime)
                self.send_hotkey('Super_L', 'Down') # restore window
                win.configure(x=0, y=0,
                    width=math.floor(self.size.width()/3),
                    height=self.size.height())
        # if self.wid2 is not None:
        #     self.post('wpos {} {}'\
        #         .format(math.ceil(self.size.width()/3), 0))
        #     self.post('wsize {} {}'\
        #         .format(math.floor(self.size.width()*2/3), self.size.height()))
        for wid in [self.wid2, self.wid4]: # align CGX and webbrowser
            if wid is not None:
                win = self.d.create_resource_object('window', wid)
                win.set_input_focus(Xlib.X.RevertToNone, Xlib.X.CurrentTime)
                self.send_hotkey('Super_L', 'Down') # restore window
                win.configure(x=math.ceil(self.size.width()/3), y=0,
                    width=math.floor(self.size.width()*2/3),
                    height=self.size.height())
        self.d.sync()

    """
        # TODO w.get_attributes('colormap'):
        # self.d.create_resource_object('colormap', ...)

        def change_colormap(self):
            self.colormap = self.screen.default_colormap
            col = self.colormap.alloc_color(0, 0, 0)
            print(col)
            self.colormap.store_colors([col])
            for c in self.colormap.query_colors([0, 1, 2, 3, 4, 5, 6, 7]):
                print(c)
            DefaultVisual
            print(self.screen.default_colormap)
            c = win.create_colormap() # Xlib.xobject.colormap.Colormap
            win.change_attributes(colormap=c)
            win.configure(colormap=c)
            win.set_wm_colormap_windows()
            for cmap in win.list_installed_colormaps():
                print(cmap)
            for cmap in win.list_installed_colormaps():
                print(cmap)
    """


class Windows_window(Window):

    def initialize(self):
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

        self.wid1 = self.get_wid('CalculiX CAE')

    # If window found, its ID is returned
    @wid_wrapper(Window)
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
                if title in buff.value:
                    # TODO It also could be, for example, a browser window
                    # with opened web page about CalculiX CAE
                    self.wid = hwnd
            return True

        ctypes.windll.user32.EnumWindows(enum_proc, 0)
        return self.wid

    def log_window_list(self):
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
                logging.debug('0x{} {:>6s} {}'\
                    .format(hex(hwnd)[2:].zfill(8), str(pid)[8:-1], buff.value))
            return True

        ctypes.windll.user32.EnumWindows(enum_proc, 0)

    # Activate window, send message to CGX, deactivate
    @post_wrapper(Window)
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
                logging.warning('WARNING! Symbol {} is not supported.'.format(symbol))

        ctypes.windll.user32.SetForegroundWindow(self.wid2)
        for symbol in cmd:
            time.sleep(0.001) # BUG doesn't work otherwise
            sendkey(symbol)
        ctypes.windll.user32.SetForegroundWindow(self.wid1)

    # Align CGX and browser windows
    # CAE window is already aligned in __init__()
    # NOTE CGX 'wpos' and 'wsize' works worse than 'ctypes'
    def align(self):
        ctypes.windll.user32.SetProcessDPIAware() # account for scaling
        for wid in [self.wid1, self.wid3]:
            if wid is not None:
                ctypes.windll.user32.ShowWindow(wid, 9) # restore window
                ctypes.windll.user32.MoveWindow(wid, 0, 0,
                    math.floor(self.size.width()/3), self.size.height(), True)
        for wid in [self.wid2, self.wid4]:
            if wid is not None:
                ctypes.windll.user32.ShowWindow(wid, 9) # restore window
                ctypes.windll.user32.MoveWindow(wid,
                    math.ceil(self.size.width()/3), 0,
                    math.floor(self.size.width()*2/3),
                    self.size.height(), True)
