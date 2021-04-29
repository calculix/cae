#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Methods for GraphiX window. """

# Standard modules
import os
import logging

"""
    # Paint element sets in default CGX colors
    def paint_elsets_old(f, elsets):
        colors = 'rgbymntk'
        i = 0
        for i in range(len(elsets)):
            if elsets[i].upper() == 'ALL':
                elsets.pop(i)
                break
        if len(elsets) > 1:
            for elset in elsets:
                f.connections[1].post('plus e {} {}'.format(elset, colors[i]))
                i = (i + 1) % len(colors)
"""

# Paint element sets in CGX
def paint_elsets(f, m):
    f.connections[1].post('plot e all')
    f.connections[1].post('minus e all')
    elsets = [e.name for e in m.Mesh.elsets.values()]
    i = 0
    for elset in elsets:
        if elset.upper() == 'ALL':
            continue
        f.connections[1].post('plus e {} blue{}'.format(elset, i))
        i = (i + 1) % 5

# Paint surfaces in CGX
def paint_surfaces(f, m):
    f.connections[1].post('plot e all')
    f.connections[1].post('minus e all')
    surfaces = [s.name for s in m.Mesh.surfaces.values()]
    i = 0
    for surf in surfaces:
        if surf.upper() == 'ALL':
            continue
        f.connections[1].post('plus f {} pink{}'.format(surf, i))
        i = (i + 1) % 5

# Open INP model in GraphiX
def open_inp(f, inp_file, has_nodes=0):
    if os.path.isfile(inp_file):
        f.kill_slave() # close old CGX
        if not has_nodes:
            logging.warning('Empty mesh, CGX will not start!')
            return
        cmd = f.p.path_cgx + ' -c ' + inp_file
        slave_title = 'CalculiX GraphiX'
        f.run_slave(cmd, slave_title)
        read_fbd_file(f, 'cgx_start.fbd')
        read_fbd_file(f, 'cgx_iso.fbd')
        read_fbd_file(f, 'cgx_colors.fbd')
    else:
        logging.error('File not found:\n' + inp_file)

# Open FRD results in GraphiX
def open_frd(f, frd_file):
    if os.path.isfile(frd_file):
        cmd = f.p.path_cgx + ' -o ' + frd_file
        slave_title = 'CalculiX GraphiX'
        f.run_slave(cmd, slave_title)
        read_fbd_file(f, 'cgx_start.fbd')
        read_fbd_file(f, 'cgx_iso.fbd')
    else:
        logging.error('File not found:\n' \
            + frd_file \
            + '\nSubmit analysis first.')

def read_fbd_file(f, basename):
    file_name = os.path.join(f.p.config, basename)
    if os.path.isfile(file_name):
        f.connections[1].post('read ' + file_name)
    else:
        logging.error('No config file ' + basename)
