#!/usr/bin/env python3
# -*- coding: utf-8 -*-


""" © Ihor Mirzov, April 2020
Distributed under GNU General Public License v3.0

Methods for CGX window. """

import os
try:
    import psutil
except:
    msg = 'Please, install psutil with command:\n'\
        + 'pip3 install psutil'
    sys.exit(msg)

def kill(parent=None):
    if parent is None:
        parent = psutil.Process()
    for ch in parent.children():
        kill(parent=ch)
        ch.terminate()
        if not ch.pid in psutil.pids():
            print('Killed PID', ch.pid)

def paint_elsets(w, elsets):
    colors = 'rgbymntk'
    i = 0
    for i in range(len(elsets)):
        if elsets[i].upper() == 'ALL':
            elsets.pop(i)
            break
    if len(elsets) > 1:
        for elset in elsets:
            w.post('plus e {} {}'.format(elset, colors[i]))
            i = (i + 1) % len(colors)
