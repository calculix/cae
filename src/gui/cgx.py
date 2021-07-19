#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

The module represents CGX menu.
It contains functions for CalculiX GraphiX window.

NOTE paint_elsets_old() is not used.
"""

# Standard modules
import os
import sys
import logging

# My modules
sys_path = os.path.abspath(__file__)
sys_path = os.path.dirname(sys_path)
sys_path = os.path.join(sys_path, '..')
sys_path = os.path.normpath(sys_path)
sys_path = os.path.realpath(sys_path)
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
import path
from gui.window import factory
from model import m


def read_fbd_file(basename):
    file_name = os.path.join(path.p.config, basename)
    if os.path.isfile(file_name):
        factory.connection.post('read ' + file_name)
    else:
        logging.error('No config file ' + basename)


def paint_elsets_old(elsets):
    """Paint element sets in default CGX colors."""
    colors = 'rgbymntk'
    i = 0
    for i in range(len(elsets)):
        if elsets[i].upper() == 'ALL':
            elsets.pop(i)
            break
    if len(elsets) > 1:
        for elset in elsets:
            factory.connection.post('plus e {} {}'.format(elset, colors[i]))
            i = (i + 1) % len(colors)


"""Menu CGX."""


def paint_elsets():
    """Paint element sets in CGX when INP is opened."""
    if not (path.p.path_cgx + ' -c ') in factory.sw.cmd:
        msg = 'Please, open INP model to paint elsets.'
        logging.warning(msg)
        return

    factory.connection.post('plot e all')
    factory.connection.post('minus e all')
    elsets = [e.name for e in m.Mesh.elsets.values()]
    i = 0
    for elset in elsets:
        if elset.upper() == 'ALL':
            continue
        factory.connection.post('plus e {} blue{}'.format(elset, i))
        i = (i + 1) % 5


def paint_surfaces():
    """Paint surfaces in CGX when INP is opened."""
    if not (path.p.path_cgx + ' -c ') in factory.sw.cmd:
        msg = 'Please, open INP model to paint surfaces.'
        logging.warning(msg)
        return

    factory.connection.post('plot e all')
    factory.connection.post('minus e all')
    surfaces = [s.name for s in m.Mesh.surfaces.values()]
    i = 0
    for surf in surfaces:
        if surf.upper() == 'ALL':
            continue
        factory.connection.post('plus f {} pink{}'.format(surf, i))
        i = (i + 1) % 5


def open_inp(inp_file, has_nodes=0):
    """Open INP model in GraphiX."""
    if not os.path.isfile(path.p.path_cgx):
        logging.error('CGX not found in ' + path.p.path_cgx)
        raise SystemExit # the best way to exit

    if os.path.isfile(inp_file):
        factory.kill_slave() # close old CGX
        if not has_nodes:
            logging.warning('Empty mesh, CGX will not start!')
            return
        cmd = path.p.path_cgx + ' -c ' + inp_file
        factory.run_slave(cmd)
        read_fbd_file('cgx_start.fbd')
        read_fbd_file('cgx_iso.fbd')
        read_fbd_file('cgx_colors.fbd')
    else:
        logging.error('File not found:\n' + inp_file)


def open_frd(frd_file):
    """Open FRD results in GraphiX."""
    if not os.path.isfile(path.p.path_cgx):
        logging.error('CGX not found in ' + path.p.path_cgx)
        raise SystemExit # the best way to exit

    if os.path.isfile(frd_file):
        cmd = path.p.path_cgx + ' -o ' + frd_file
        factory.run_slave(cmd)
        read_fbd_file('cgx_start.fbd')
        read_fbd_file('cgx_iso.fbd')
    else:
        logging.error('File not found:\n' \
            + frd_file \
            + '\nSubmit analysis first.')


def cmap(colormap):
    """Set custom colormap when FRD is opened."""
    if not (path.p.path_cgx + ' -o ') in factory.sw.cmd:
        msg = 'Please, open FRD model to set colormap.'
        logging.warning(msg)
        return
    factory.connection.post('cmap ' + colormap)
