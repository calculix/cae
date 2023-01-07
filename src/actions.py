#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2023
Distributed under GNU General Public License v3.0

Main window actions - all processed signals.
"""

# Standard modules
import os

# External modules
from PyQt5 import QtWidgets

# My modules
from path import p
from settings import s
from gui import cgx
from gui.window import wf

if wf.mw is None:
    raise SystemExit

# File actions
from importer import i
wf.mw.action_file_import.triggered.connect(lambda: i.import_file(None))
wf.mw.action_file_settings.triggered.connect(s.open)
wf.mw.action_file_exit.triggered.connect(QtWidgets.qApp.quit)

# Job actions
from gui.job import j
from model.kom import KWT
wf.mw.action_job_write_input.triggered.connect(
    lambda: j.write_input(KWT.get_inp_code_as_lines()))
wf.mw.action_job_write_input.triggered.connect(
    lambda: wf.mw.setWindowTitle('CalculiX Advanced Environment - ' + j.name))
wf.mw.action_job_edit_inp.triggered.connect(lambda: j.open_inp())
wf.mw.action_job_subroutine.triggered.connect(lambda: j.open_subroutine())
wf.mw.action_job_rebuild_ccx.triggered.connect(lambda: j.rebuild_ccx())
wf.mw.action_job_submit.triggered.connect(lambda: j.submit())
wf.mw.action_job_view_log.triggered.connect(lambda: j.view_log())
wf.mw.action_job_export_vtu.triggered.connect(j.export_vtu)
wf.mw.action_job_paraview.triggered.connect(lambda: j.open_paraview())

# CGX actions
from model import m
wf.mw.action_cgx_paint_elsets.triggered.connect(
    lambda: cgx.paint_elsets())
wf.mw.action_cgx_paint_surfaces.triggered.connect(
    lambda: cgx.paint_surfaces())
wf.mw.action_cgx_inp.triggered.connect(lambda: cgx.open_inp(j.inp, len(m.Mesh.nodes)))
wf.mw.action_cgx_frd.triggered.connect(lambda: cgx.open_frd(j.frd))
wf.mw.action_cgx_cmap_classic.triggered.connect(lambda: cgx.cmap('classic'))
wf.mw.action_cgx_cmap_inferno.triggered.connect(lambda: cgx.cmap('inferno'))
wf.mw.action_cgx_cmap_turbo.triggered.connect(lambda: cgx.cmap('turbo'))
wf.mw.action_cgx_cmap_viridis.triggered.connect(lambda: cgx.cmap('viridis'))

# Help actions
wf.mw.action_help_readme.triggered.connect(
    lambda: wf.mw.help('https://github.com/calculix/cae'))
wf.mw.action_help_examples.triggered.connect(
    lambda: wf.mw.help('https://github.com/calculix/examples'))
wf.mw.action_help_issues.triggered.connect(
    lambda: wf.mw.help('https://github.com/calculix/cae/issues'))

# treeView actions
# if hasattr(wf.mw, 'treeView'):
from gui.tree import t
wf.mw.keyPressEvent = t.keyPressEvent
wf.mw.treeView.doubleClicked.connect(t.doubleClicked)
wf.mw.treeView.clicked.connect(t.clicked)
wf.mw.treeView.customContextMenuRequested.connect(t.rightClicked)
wf.mw.treeView.expanded.connect(t.expanded_or_collapsed)
wf.mw.treeView.collapsed.connect(t.expanded_or_collapsed)

# ToolBar actions
wf.mw.action_view_minus_x.triggered.connect(lambda: wf.connection.post('rot -x'))
wf.mw.action_view_minus_y.triggered.connect(lambda: wf.connection.post('rot -y'))
wf.mw.action_view_minus_z.triggered.connect(lambda: wf.connection.post('rot -z'))
wf.mw.action_view_plus_x.triggered.connect(lambda: wf.connection.post('rot x'))
wf.mw.action_view_plus_y.triggered.connect(lambda: wf.connection.post('rot y'))
wf.mw.action_view_plus_z.triggered.connect(lambda: wf.connection.post('rot z'))
wf.mw.action_view_frame.triggered.connect(lambda: wf.connection.post('frame'))

# Workaround for iso view
# NOTE Three rotation posts to CGX window doesn't work in Windows
# So one may use .fbd commands
def action_view_iso():
    file_name = os.path.join(p.config, 'cgx_iso.fbd')
    cgx.read_fbd_file(file_name)

wf.mw.action_view_iso.triggered.connect(action_view_iso)
wf.mw.action_view_line.triggered.connect(lambda: wf.connection.post('view elem off'))
wf.mw.action_view_line.triggered.connect(lambda: wf.connection.post('view line'))
wf.mw.action_view_fill.triggered.connect(lambda: wf.connection.post('view elem off'))
wf.mw.action_view_fill.triggered.connect(lambda: wf.connection.post('view fill'))
wf.mw.action_view_elem.triggered.connect(lambda: wf.connection.post('view fill'))
wf.mw.action_view_elem.triggered.connect(lambda: wf.connection.post('view elem'))
