#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, July 2020
Distributed under GNU General Public License v3.0

Methods for GraphiX window. """

# Standard modules
import os
import time
import logging
import traceback

# TODO Write a class and move w.run_cgx() here

# Kill all CGX processes
def kill(w):
    w.wid2 = None
    if w.cgx_process is not None:
        count = 0
        while w.cgx_process.poll() is None:
            try:
                if os.name == 'nt':
                    os.system('TASKKILL /F /PID {} /T'.format(w.cgx_process.pid))
                else:
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
            logging.debug('Killed CGX, PID={}.'.format(w.cgx_process.pid))
            w.cgx_process = None

"""
    # Paint element sets in default CGX colors
    def paint_elsets_old(w, elsets):
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
"""

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

# Open INP model in GraphiX
def open_inp(w, inp_file, has_nodes=0):
    if os.path.isfile(inp_file):
        kill(w) # close old CGX
        if not has_nodes:
            logging.warning('Empty mesh, CGX will not start!')
            return
        w.run_cgx('-c ' + inp_file)
    else:
        logging.error('File not found:\n' + inp_file)

# Open FRD results in GraphiX
def open_frd(w, frd_file):
    if os.path.isfile(frd_file):
        w.run_cgx('-o ' + frd_file)
    else:
        logging.error('File not found:\n' \
            + frd_file \
            + '\nSubmit analysis first.')
