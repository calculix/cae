#!/usr/bin/env python3
# -*- coding: utf-8 -*-


""" © Ihor Mirzov, May 2020
Distributed under GNU General Public License v3.0

Main window class. """


# Standard modules
import os
import sys
import time
import logging
import subprocess
import threading
import inspect
import math
if 'nt' in os.name:
    import ctypes

# External modules
from PyQt5 import QtWidgets, uic, QtCore, QtGui
if 'posix' in os.name:
    try:
        import Xlib
        from Xlib import display, protocol, X, XK
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

        # Handler to show logs in the CAE's textEdit
        gui.log.add_text_handler(self.textEdit)

        # When logger is ready - check if settings read correctly
        if hasattr(s, 'error_path'):
            logging.error('Error path in settings file: ' +\
                s.error_path + '. Loading default values.')

        # Window ID - to pass keycodes to
        self.wid1 = None # cae window
        self.wid2 = None # cgx window
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

    # TODO No immediate reaction on manual writing into CGX window
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
    def help(self, link):
        url = QtCore.QUrl(link)
        logging.info('Going to ' + link)
        if not QtGui.QDesktopServices.openUrl(url):
            logging.warning('Can\'t open url: ' + link)


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
                    logging.error('ERROR! Can\'t get {} window.'\
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
        self.d = display.Display()
        self.screen = self.d.screen()
        self.root = self.screen.root

        # 0:lowercase, 1:shifted
        self.keyboardMapping = {
            '\t':(0x17, 0),
            '\n':(0x24, 0),
            ' ': (0x41, 0),
            '!': (0xa, 1),
            '#': (0xc, 1),
            '%': (0xe, 1),
            '$': (0xd, 1),
            '&': (0x10, 1),
            '\"':(0x30, 1),
            '\'':(0x30, 0),
            '(': (0xbb, 1),
            ')': (0xbc, 1),
            '*': (0x11, 1),
            '=': (0x15, 0),
            '+': (0x15, 1),
            ',': (0x3b, 0),
            '-': (0x14, 0),
            '.': (0x3c, 0),
            '/': (0x3d, 0),
            ':': (0x2f, 1),
            ';': (0x2f, 0),
            '<': (0x5e, 0),
            '>': (0x3c, 1),
            '?': (0x3d, 1),
            '@': (0xb, 1),
            '[': (0x22, 0),
            ']': (0x23, 0),
            '^': (0xf, 1),
            '_': (0x14, 1),
            '`': (0x31, 0),
            '{': (0x22, 1),
            '\\':(0x33, 0),
            '|': (0x33, 1),
            '}': (0x23, 1),
            '~': (0x31, 1),
            '0': (0x13, 0),
            '1': (0xa, 0),
            '2': (0xb, 0),
            '3': (0xc, 0),
            '4': (0xd, 0),
            '5': (0xe, 0),
            '6': (0xf, 0),
            '7': (0x10, 0),
            '8': (0x11, 0),
            '9': (0x12, 0),
            'a': (0x26, 0),
            'b': (0x38, 0),
            'c': (0x36, 0),
            'd': (0x28, 0),
            'e': (0x1a, 0),
            'f': (0x29, 0),
            'g': (0x2a, 0),
            'h': (0x2b, 0),
            'i': (0x1f, 0),
            'j': (0x2c, 0),
            'k': (0x2d, 0),
            'l': (0x2e, 0),
            'm': (0x3a, 0),
            'n': (0x39, 0),
            'o': (0x20, 0),
            'p': (0x21, 0),
            'q': (0x18, 0),
            'r': (0x1b, 0),
            's': (0x27, 0),
            't': (0x1c, 0),
            'u': (0x1e, 0),
            'v': (0x37, 0),
            'w': (0x19, 0),
            'x': (0x35, 0),
            'y': (0x1d, 0),
            'z': (0x34, 0),
            'A': (0x26, 1),
            'B': (0x38, 1),
            'C': (0x36, 1),
            'D': (0x28, 1),
            'E': (0x1a, 1),
            'F': (0x29, 1),
            'G': (0x2a, 1),
            'H': (0x2b, 1),
            'I': (0x1f, 1),
            'J': (0x2c, 1),
            'K': (0x2d, 1),
            'L': (0x2e, 1),
            'M': (0x3a, 1),
            'N': (0x39, 1),
            'O': (0x20, 1),
            'P': (0x21, 1),
            'Q': (0x18, 1),
            'R': (0x1b, 1),
            'S': (0x27, 1),
            'T': (0x1c, 1),
            'U': (0x1e, 1),
            'V': (0x37, 1),
            'W': (0x19, 1),
            'X': (0x35, 1),
            'Y': (0x1d, 1),
            'Z': (0x34, 1)}

    # Get window by title and return its ID
    @wid_wrapper(Window)
    def get_wid(self, title):
        def get_window_by_title(title, window=None):
            if window is None:
                window = self.root
            try:
                for w in window.query_tree().children:
                    if w.get_wm_name() == title:
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

    # Activate window, send message to CGX, deactivate
    @post_wrapper(Window)
    def post(self, cmd):

        # Create X event
        def event(w, e, keycode, case):
            return e(
                time=X.CurrentTime,
                root=self.root,
                window=w,
                same_screen=0, child=X.NONE,
                root_x=0, root_y=0, event_x=0, event_y=0,
                state=case, # 0=lowercase, 1=uppercase
                detail=keycode)

        # Key press + release
        def sendkey(w, symbol):
            if symbol in self.keyboardMapping:
                code = self.keyboardMapping[symbol][0]
                case = self.keyboardMapping[symbol][1]
                e1 = event(w, protocol.event.KeyPress, code, case)
                w.send_event(e1, propagate=True)
                e2 = event(w, protocol.event.KeyRelease, code, case)
                w.send_event(e2, propagate=True)
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
            self.d.set_input_focus(w1, X.RevertToNone, X.CurrentTime)
            self.d.sync()
        else:
            print('ERROR! Window is None.')

    # Align CAE and CGX windows
    # NOTE CGX 'wpos' and 'wsize' works better than 'w2configure'
    def align(self):
        width = self.screen.width_in_pixels
        height = self.screen.height_in_pixels
        w1 = self.d.create_resource_object('window', self.wid1)
        w1.configure(x=0, y=0,
            width=math.ceil(width/3), height=height)
        # w2 = self.d.create_resource_object('window', self.wid2)
        # w2.configure(x=math.ceil(width/3), y=0,
        #     width=math.floor(width*2/3), height=height)
        self.d.sync()
        self.post('wpos {} {}'.format(math.ceil(width/3), 0))
        self.post('wsize {} {}'.format(math.floor(width*2/3), height))
        time.sleep(0.3)

    # TODO colormaps for CGX
    def change_colormap(self):
        pass
        # self.colormap = self.screen.default_colormap
        # col = self.colormap.alloc_color(0, 0, 0)
        # print(col)
        # self.colormap.store_colors([col])
        # for c in self.colormap.query_colors([0, 1, 2, 3, 4, 5, 6, 7]):
        #     print(c)
        # DefaultVisual
        # print(self.screen.default_colormap)
        # c = w1.create_colormap() # Xlib.xobject.colormap.Colormap
        # w1.change_attributes(colormap=c)
        # w1.configure(colormap=c)
        # w1.set_wm_colormap_windows()
        # for cmap in w1.list_installed_colormaps():
        #     print(cmap)
        # for cmap in w2.list_installed_colormaps():
        #     print(cmap)


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
            time.sleep(0.001) # BUG doesn't work without print
            sendkey(symbol)
        ctypes.windll.user32.SetForegroundWindow(self.wid1)

    # Align CAE and CGX windows
    # NOTE CGX 'wpos' and 'wsize' works worse than 'ctypes'
    def align(self):
        width = ctypes.windll.user32.GetSystemMetrics(0)
        height = ctypes.windll.user32.GetSystemMetrics(1) - 56
        ctypes.windll.user32.MoveWindow(self.wid1,
            0, 0,
            math.ceil(width/3), height, True)
        ctypes.windll.user32.MoveWindow(self.wid2,
            math.ceil(width/3), 0,
            math.floor(width*2/3), height, True)
        time.sleep(0.3)
