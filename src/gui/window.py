#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, May 2020
Distributed under GNU General Public License v3.0

Main window class. Get windows names and align them.
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
import threading
import inspect
import math
import webbrowser
if 'nt' in os.name:
    import ctypes

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
        self.s = s
        QtWidgets.QMainWindow.__init__(self) # create main window
        uic.loadUi(p.cae_xml, self) # load form
        if self.s.align_windows:
            size = QtWidgets.QDesktopWidget().screenGeometry(-1)
            self.setGeometry(0, 0, size.width()/3, size.height())

        # Handler to show logs in the CAE's textEdit
        gui.log.add_text_handler(self.textEdit)

        # When logger is ready - check if settings read correctly
        if hasattr(s, 'error_path'):
            logging.error('Error path in settings file: ' +\
                s.error_path + '. Loading default values.')

        # Window ID - to pass keycodes to
        self.wid1 = None # cae window
        self.wid2 = None # cgx window
        self.wid3 = None # help (web browser)
        self.process = None # running process to send commands to
        self.keyboardMapping = None
        self.last_command = None

        # INP | FRD - corresponds to opened file
        self.mode = None

        # self.toolBar.setParent(None) # hide toolbar

    # Close opened CGX (if any)
    # open a new one and get window ID
    def run_cgx(self, cmd):

        def start_cgx(cmd):    
            self.process = subprocess.Popen(cmd.split(),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT)

            self.wid2 = self.get_wid('CalculiX GraphiX') # could be None
            if self.wid2 is None:
                msg = 'ERROR! Can\'t get {} window ID.'\
                    .format('CalculiX GraphiX')
                logging.error(msg)
                sys.exit(msg)
            gui.log.read_output(self.process.stdout)
            if self.s.align_windows:
                self.align()
            self.flush_cgx_cache()

        if not os.path.isfile(self.s.path_cgx):
            msg = 'Wrong path to CGX:\n{}\n' \
                + 'Configure it in File->Settings.'
            logging.error(msg.format(self.s.path_cgx))
            return

        gui.cgx.kill()
        start_cgx(cmd)
        self.register_cgx_colors()

        # Caller fuction name: cgx_inp | cgx_frd
        self.mode = inspect.stack()[1].function

    # The only way to get last output from CGX
    def flush_cgx_cache(self):
        self.post(' ')

    # Register new colors to use in CGX
    # https://colourco.de/
    def register_cgx_colors(self):

        # Convert HTML color notation to RGB
        def hex_to_rgb(value):
            value = value.lstrip('#')
            rgb = []
            for i in range(0, 6, 2):
                amp = int(value[i:i+2], 16) / 255
                rgb.append(round(amp, 2))
            return rgb

        blue = [
            '#19324A',
            '#254C71',
            '#306699',
            '#3B80C1',
            '#6099CF']
        pink = [
            '#331132',
            '#5A1D58',
            '#81297E',
            '#A934A5',
            '#C848C4']
        for i in range(5):
            self.post('col blue{} {} {} {}'\
                .format(i, *hex_to_rgb(blue[i])))
            self.post('col pink{} {} {} {}'\
                .format(i, *hex_to_rgb(pink[i])))

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
                time.sleep(0.1)
                wid = method(w, title)
                if time.perf_counter() - start > 1:
                    logging.error('Can\'t get \'{}\' window.'\
                        .format(title))
                    return None
            msg = '{} WID={}'.format(title, wid)
            logging.debug(msg)
            return wid
        return fcn
    return wrap


# Common wrapper for Window.post()
def post_wrapper(w):
    def wrap(method):
        def fcn(w, cmd):
            if w.process is not None \
                and w.process.poll() is None\
                and w.wid2 is not None:
                if len(cmd):
                    if cmd != ' ':
                        w.last_command = cmd
                    method(w, cmd)
                    method(w, ' ')
                    return
                else:
                    logging.warning('Empty command.')
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

    # Get window by title and return its ID
    @wid_wrapper(Window)
    def get_wid(self, title):
        ids = self.root.get_full_property(
            self.d.intern_atom('_NET_CLIENT_LIST'),
            X.AnyPropertyType).value
        for wid in ids:
            w = self.d.create_resource_object('window', wid)
            try:
                wname = w.get_full_property(
                    self.d.intern_atom('_NET_WM_NAME'), X.AnyPropertyType).value.decode()
                pid = w.get_full_property(
                    self.d.intern_atom('_NET_WM_PID'), X.AnyPropertyType).value[0]
            except:
                wname = w.get_wm_name()
                pid = 0
            if title.lower() in wname.lower():
                logging.debug('0x{} {: 6d} {}'.format(hex(wid)[2:].zfill(8), pid, wname))
                return wid
        return None

    """
        def get_wid_old(self, title):
            def get_window_by_title(title, window=None):
                if window is None:
                    window = self.root
                try:
                    for w in window.query_tree().children:
                        wm_name = w.get_wm_name()
                        # print(title, wm_name, w.get_wm_class())
                        if wm_name is not None \
                            and len(wm_name) \
                            and title in wm_name:
                                return w
                        w = get_window_by_title(title, w)
                        if w is not None:
                            return w
                except:
                    pass
                return None

            w = get_window_by_title(title)
            try: return int(w.id)
            except: return None
    """

    # Activate window, send message to CGX, deactivate
    @post_wrapper(Window)
    def post(self, cmd):

        # Create X event
        def event(w, e, keycode, case):
            return e(
                time=Xlib.X.CurrentTime,
                root=self.root,
                window=w,
                same_screen=0, child=Xlib.X.NONE,
                root_x=0, root_y=0, event_x=0, event_y=0,
                state=case, # 0=lowercase, 1=uppercase
                detail=keycode)

        # Key press + release
        def sendkey(w, symbol):
            if symbol in self.keyboardMapping:
                c = self.keyboardMapping[symbol][0]
                keysym = XK.string_to_keysym(c)
                code = self.d.keysym_to_keycode(keysym)
                case = self.keyboardMapping[symbol][1]

                e = event(w, Xlib.protocol.event.KeyPress, code, case)
                w.send_event(e, propagate=True)
                e = event(w, Xlib.protocol.event.KeyRelease, code, case)
                w.send_event(e, propagate=True)

            else:
                logging.warning('WARNING! Symbol {} is not supported.'.format(symbol))

        w2 = self.d.create_resource_object('window', self.wid2)
        for symbol in cmd + '\n':
            sendkey(w2, symbol)

        def get_window_by_id(wid, window=None):
            if window is None:
                window = self.root
            try:
                for w in window.query_tree().children:
                    if w.id == wid:
                        return w
                    w = get_window_by_id(wid, w)
                    if w is not None:
                        return w
            except Xlib.error.BadWindow:
                logging.error('Bad window')
            return None

        w1 = get_window_by_id(self.wid1)
        if w1 is not None:
            self.d.set_input_focus(w1, Xlib.X.RevertToNone, Xlib.X.CurrentTime)
            self.d.sync()

    # Key press + release
    def send_key_combination(self, key1, key2):
        c1 = self.keyboardMapping[key1][0]
        keysym1 = XK.string_to_keysym(c1)
        code1 = self.d.keysym_to_keycode(keysym1)

        c2 = self.keyboardMapping[key2][0]
        keysym2 = XK.string_to_keysym(c2)
        code2 = self.d.keysym_to_keycode(keysym2)

        fake_input(self.d, X.KeyPress, code1)
        fake_input(self.d, X.KeyPress, code2)
        fake_input(self.d, X.KeyRelease, code2)
        fake_input(self.d, X.KeyRelease, code1)
        time.sleep(0.3)
        self.d.sync()

    # Align CGX and browser windows
    # NOTE CGX 'wpos' and 'wsize' works better than 'w2configure'
    def align(self):
        width = self.screen.width_in_pixels
        height = self.screen.height_in_pixels
        if self.wid2 is not None:
            self.post('wpos {} {}'.format(math.ceil(width/3), 0))
            self.post('wsize {} {}'.format(math.floor(width*2/3), height))
        if self.wid3 is not None:
            w3 = self.d.create_resource_object('window', self.wid3)
            w3.set_input_focus(Xlib.X.RevertToNone, Xlib.X.CurrentTime)
            self.send_key_combination('Super_L', 'Down')
            w3.configure(x=math.ceil(width/3), y=0,
                width=math.floor(width*2/3), height=height)
        time.sleep(0.3)
        self.d.sync()

    # TODO w.get_attributes('colormap'):
    """
        def change_colormap(self):
            self.colormap = self.screen.default_colormap
            col = self.colormap.alloc_color(0, 0, 0)
            print(col)
            self.colormap.store_colors([col])
            for c in self.colormap.query_colors([0, 1, 2, 3, 4, 5, 6, 7]):
                print(c)
            DefaultVisual
            print(self.screen.default_colormap)
            c = w1.create_colormap() # Xlib.xobject.colormap.Colormap
            w1.change_attributes(colormap=c)
            w1.configure(colormap=c)
            w1.set_wm_colormap_windows()
            for cmap in w1.list_installed_colormaps():
                print(cmap)
            for cmap in w2.list_installed_colormaps():
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

    # If window found, its ID is returned
    @wid_wrapper(Window)
    def get_wid(self, title):
        hwnd = ctypes.windll.user32.FindWindowW(None, title)
        if hwnd: return int(hwnd)
        else: return None

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
        for symbol in cmd + '\n':
            time.sleep(0.001) # BUG doesn't work otherwise
            sendkey(symbol)
        ctypes.windll.user32.SetForegroundWindow(self.wid1)

    # Align CGX and browser windows
    # NOTE CGX 'wpos' and 'wsize' works worse than 'ctypes'
    def align(self):
        width = ctypes.windll.user32.GetSystemMetrics(0)
        height = ctypes.windll.user32.GetSystemMetrics(1) - 56
        if self.wid2 is not None:
            ctypes.windll.user32.MoveWindow(self.wid2,
                math.ceil(width/3), 0,
                math.floor(width*2/3), height, True)
        time.sleep(0.3)
