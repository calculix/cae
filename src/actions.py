#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Main window actions - all processed signals.
"""

# Standard modules
import os
import logging

# External modules
from PyQt5 import QtWidgets

# My modules
import path
import settings
import gui
from gui.window import factory


# File actions
from importer import i
factory.mw.action_file_import.triggered.connect(lambda: i.import_file(None))
factory.mw.action_file_settings.triggered.connect(settings.s.open)
factory.mw.action_file_exit.triggered.connect(QtWidgets.qApp.quit)

# Job actions
from gui.job import j
from model.kom import KOM
factory.mw.action_job_write_input.triggered.connect(
    lambda: j.write_input(KOM.get_inp_code_as_lines()))
factory.mw.action_job_write_input.triggered.connect(
    lambda: factory.mw.setWindowTitle('CalculiX Advanced Environment - ' + j.name))
factory.mw.action_job_edit_inp.triggered.connect(lambda: j.open_inp())
factory.mw.action_job_subroutine.triggered.connect(lambda: j.open_subroutine())
factory.mw.action_job_rebuild_ccx.triggered.connect(lambda: j.rebuild_ccx())
factory.mw.action_job_submit.triggered.connect(lambda: j.submit())
factory.mw.action_job_view_log.triggered.connect(lambda: j.view_log())
factory.mw.action_job_export_vtu.triggered.connect(j.export_vtu)
factory.mw.action_job_paraview.triggered.connect(lambda: j.open_paraview())

# CGX actions
from model import m
factory.mw.action_cgx_paint_elsets.triggered.connect(
    lambda: gui.cgx.paint_elsets())
factory.mw.action_cgx_paint_surfaces.triggered.connect(
    lambda: gui.cgx.paint_surfaces())
factory.mw.action_cgx_inp.triggered.connect(lambda: gui.cgx.open_inp(j.inp, len(m.Mesh.nodes)))
factory.mw.action_cgx_frd.triggered.connect(lambda: gui.cgx.open_frd(j.frd))
factory.mw.action_cgx_cmap_classic.triggered.connect(lambda: gui.cgx.cmap('classic'))
factory.mw.action_cgx_cmap_inferno.triggered.connect(lambda: gui.cgx.cmap('inferno'))
factory.mw.action_cgx_cmap_turbo.triggered.connect(lambda: gui.cgx.cmap('turbo'))
factory.mw.action_cgx_cmap_viridis.triggered.connect(lambda: gui.cgx.cmap('viridis'))

# Help actions
factory.mw.action_help_readme.triggered.connect(
    lambda: factory.mw.help('https://github.com/calculix/cae'))
factory.mw.action_help_examples.triggered.connect(
    lambda: factory.mw.help('https://github.com/calculix/examples'))
factory.mw.action_help_issues.triggered.connect(
    lambda: factory.mw.help('https://github.com/calculix/cae/issues'))

# treeView actions
if hasattr(factory.mw, 'treeView'):
    from gui.tree import t
    factory.mw.keyPressEvent = t.keyPressEvent
    factory.mw.treeView.doubleClicked.connect(t.doubleClicked)
    factory.mw.treeView.clicked.connect(t.clicked)
    factory.mw.treeView.customContextMenuRequested.connect(t.rightClicked)
    factory.mw.treeView.expanded.connect(t.expanded_or_collapsed)
    factory.mw.treeView.collapsed.connect(t.expanded_or_collapsed)

# ToolBar actions
factory.mw.action_view_minus_x.triggered.connect(lambda: factory.connection.post('rot -x'))
factory.mw.action_view_minus_y.triggered.connect(lambda: factory.connection.post('rot -y'))
factory.mw.action_view_minus_z.triggered.connect(lambda: factory.connection.post('rot -z'))
factory.mw.action_view_plus_x.triggered.connect(lambda: factory.connection.post('rot x'))
factory.mw.action_view_plus_y.triggered.connect(lambda: factory.connection.post('rot y'))
factory.mw.action_view_plus_z.triggered.connect(lambda: factory.connection.post('rot z'))
factory.mw.action_view_frame.triggered.connect(lambda: factory.connection.post('frame'))

# Workaround for iso view
# Three rotation posts to CGX window doesn't work in Windows
# So one may use .fbd commands
def action_view_iso():
    file_name = os.path.join(path.p.config, 'cgx_iso.fbd')
    if not os.path.isfile(file_name):
        logging.error('No config file cgx_iso.fbd')
        return
    factory.connection.post('read ' + file_name)

factory.mw.action_view_iso.triggered.connect(action_view_iso)
factory.mw.action_view_line.triggered.connect(lambda: factory.connection.post('view elem off'))
factory.mw.action_view_line.triggered.connect(lambda: factory.connection.post('view line'))
factory.mw.action_view_fill.triggered.connect(lambda: factory.connection.post('view elem off'))
factory.mw.action_view_fill.triggered.connect(lambda: factory.connection.post('view fill'))
factory.mw.action_view_elem.triggered.connect(lambda: factory.connection.post('view fill'))
factory.mw.action_view_elem.triggered.connect(lambda: factory.connection.post('view elem'))
