#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, June 2020
Distributed under GNU General Public License v3.0

Methods for CGX window. """

# Standard modules
import os
import time
import logging
import traceback

# Kill all CGX processes
def kill(w):
    w.wid2 = None
    if w.cgx_process is not None:
        count = 0
        while w.cgx_process.poll() is None:
            try:
                w.cgx_process.kill()
            except:
                logging.error(traceback.format_exc())
            time.sleep(0.1)
            count += 1
            if count >= 10:
                break
        if w.cgx_process.poll() is None:
            msg = 'Can not kill CGX, PID={}.'.format(w.cgx_process.pid)
            logging.warning(msg)
        else:
            logging.info('Killed CGX, PID={}.'.format(w.cgx_process.pid))
            w.cgx_process = None

# def paint_elsets_old(w, elsets):
#     colors = 'rgbymntk'
#     i = 0
#     for i in range(len(elsets)):
#         if elsets[i].upper() == 'ALL':
#             elsets.pop(i)
#             break
#     if len(elsets) > 1:
#         for elset in elsets:
#             w.post('plus e {} {}'.format(elset, colors[i]))
#             i = (i + 1) % len(colors)

# Paint element sets in CGX
def paint_elsets(w, m):
    w.post('plot e all')
    w.post('minus e all')
    elsets = [e.name for e in m.Mesh.elsets.values()]
    i = 0
    for elset in elsets:
        if elset.upper() == 'ALL':
            continue
        w.post('plus e {} blue{}'.format(elset, i))
        i = (i + 1) % 5

# Paint surfaces in CGX
def paint_surfaces(w, m):
    w.post('plot e all')
    w.post('minus e all')
    surfaces = [s.name for s in m.Mesh.surfaces.values()]
    i = 0
    for surf in surfaces:
        if surf.upper() == 'ALL':
            continue
        w.post('plus f {} pink{}'.format(surf, i))
        i = (i + 1) % 5
