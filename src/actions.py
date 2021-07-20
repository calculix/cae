#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Main window actions - all processed signals.

w - Master window
t - Tree
j - Job
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
from model import m
from model.kom import KOM


def actions(t, j, i):
    w = factory.mw
    w.keyPressEvent = t.keyPressEvent

    # File actions
    w.action_file_import.triggered.connect(lambda: i.import_file(None))
    w.action_file_settings.triggered.connect(settings.s.open)
    w.action_file_exit.triggered.connect(QtWidgets.qApp.quit)

    # Job actions
    w.action_job_write_input.triggered.connect(
        lambda: j.write_input(KOM.get_inp_code_as_lines()))
    w.action_job_write_input.triggered.connect(
        lambda: w.setWindowTitle('CalculiX Advanced Environment - ' + j.name))
    w.action_job_edit_inp.triggered.connect(lambda: j.open_inp())
    w.action_job_subroutine.triggered.connect(lambda: j.open_subroutine())
    w.action_job_rebuild_ccx.triggered.connect(lambda: j.rebuild_ccx())
    w.action_job_submit.triggered.connect(lambda: j.submit())
    w.action_job_view_log.triggered.connect(lambda: j.view_log())
    w.action_job_export_vtu.triggered.connect(j.export_vtu)
    w.action_job_paraview.triggered.connect(lambda: j.open_paraview())

    # CGX actions
    w.action_cgx_paint_elsets.triggered.connect(
        lambda: gui.cgx.paint_elsets())
    w.action_cgx_paint_surfaces.triggered.connect(
        lambda: gui.cgx.paint_surfaces())
    w.action_cgx_inp.triggered.connect(lambda: gui.cgx.open_inp(j.inp, len(m.Mesh.nodes)))
    w.action_cgx_frd.triggered.connect(lambda: gui.cgx.open_frd(j.frd))
    w.action_cgx_cmap_classic.triggered.connect(lambda: gui.cgx.cmap('classic'))
    w.action_cgx_cmap_inferno.triggered.connect(lambda: gui.cgx.cmap('inferno'))
    w.action_cgx_cmap_turbo.triggered.connect(lambda: gui.cgx.cmap('turbo'))
    w.action_cgx_cmap_viridis.triggered.connect(lambda: gui.cgx.cmap('viridis'))

    # Help actions
    w.action_help_readme.triggered.connect(
        lambda: w.help('https://github.com/calculix/cae'))
    w.action_help_examples.triggered.connect(
        lambda: w.help('https://github.com/calculix/examples'))
    w.action_help_issues.triggered.connect(
        lambda: w.help('https://github.com/calculix/cae/issues'))

    # treeView actions
    w.treeView.doubleClicked.connect(t.doubleClicked)
    w.treeView.clicked.connect(t.clicked)
    w.treeView.customContextMenuRequested.connect(t.rightClicked)
    w.treeView.expanded.connect(t.expanded_or_collapsed)
    w.treeView.collapsed.connect(t.expanded_or_collapsed)

    # ToolBar actions
    w.action_view_minus_x.triggered.connect(lambda: factory.connection.post('rot -x'))
    w.action_view_minus_y.triggered.connect(lambda: factory.connection.post('rot -y'))
    w.action_view_minus_z.triggered.connect(lambda: factory.connection.post('rot -z'))
    w.action_view_plus_x.triggered.connect(lambda: factory.connection.post('rot x'))
    w.action_view_plus_y.triggered.connect(lambda: factory.connection.post('rot y'))
    w.action_view_plus_z.triggered.connect(lambda: factory.connection.post('rot z'))
    w.action_view_frame.triggered.connect(lambda: factory.connection.post('frame'))

    # Workaround for iso view
    # Three rotation posts to CGX window doesn't work in Windows
    # So one may use .fbd commands
    def action_view_iso():
        file_name = os.path.join(path.p.config, 'cgx_iso.fbd')
        if not os.path.isfile(file_name):
            logging.error('No config file cgx_iso.fbd')
            return
        factory.connection.post('read ' + file_name)
    w.action_view_iso.triggered.connect(action_view_iso)
    w.action_view_line.triggered.connect(lambda: factory.connection.post('view elem off'))
    w.action_view_line.triggered.connect(lambda: factory.connection.post('view line'))
    w.action_view_fill.triggered.connect(lambda: factory.connection.post('view elem off'))
    w.action_view_fill.triggered.connect(lambda: factory.connection.post('view fill'))
    w.action_view_elem.triggered.connect(lambda: factory.connection.post('view fill'))
    w.action_view_elem.triggered.connect(lambda: factory.connection.post('view elem'))
