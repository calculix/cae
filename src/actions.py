#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Main window actions - all processed signals. """

import os
import logging

from PyQt5 import QtWidgets

import gui
import importer

"""
s - Settings
w - Window
m - Model
t - Tree
j - Job
"""
def actions(p, s, w, m, t, j, i):
    w.keyPressEvent = t.keyPressEvent

    # File actions
    w.action_file_import.triggered.connect(lambda: i.import_file(None))
    w.action_file_settings.triggered.connect(s.open)
    w.action_file_exit.triggered.connect(QtWidgets.qApp.quit)

    # Job actions
    w.action_job_write_input.triggered.connect(
        lambda: j.write_input(m.KOM.get_inp_code_as_lines()))
    w.action_job_write_input.triggered.connect(
        lambda: w.setWindowTitle('CalculiX Advanced Environment - ' + j.name))
    w.action_job_edit_inp.triggered.connect(lambda: j.edit_inp())
    w.action_job_subroutine.triggered.connect(lambda: j.open_subroutine())
    w.action_job_rebuild_ccx.triggered.connect(lambda: j.rebuild_ccx())
    w.action_job_submit.triggered.connect(lambda: j.submit())
    w.action_job_view_log.triggered.connect(lambda: j.view_log())
    w.action_job_export_vtu.triggered.connect(j.export_vtu)
    w.action_job_paraview.triggered.connect(lambda: j.open_paraview())

    # CGX actions
    w.action_cgx_paint_elsets.triggered.connect(
        lambda: gui.cgx.paint_elsets(w, m))
    w.action_cgx_paint_surfaces.triggered.connect(
        lambda: gui.cgx.paint_surfaces(w, m))
    w.action_cgx_inp.triggered.connect(lambda: gui.cgx.open_inp(w, j.inp, len(m.Mesh.nodes)))
    w.action_cgx_frd.triggered.connect(lambda: gui.cgx.open_frd(w, j.frd))
    w.action_cgx_cmap_classic.triggered.connect(lambda: w.connections[1].post('cmap classic'))
    w.action_cgx_cmap_inferno.triggered.connect(lambda: w.connections[1].post('cmap inferno'))
    w.action_cgx_cmap_turbo.triggered.connect(lambda: w.connections[1].post('cmap turbo'))
    w.action_cgx_cmap_viridis.triggered.connect(lambda: w.connections[1].post('cmap viridis'))

    # Help actions
    w.action_help_readme.triggered.connect(
        lambda: w.help('file://{}/README.pdf'.format(p.app_home_dir)))
    w.action_help_examples.triggered.connect(
        lambda: w.help('https://github.com/calculix/examples'))
    w.action_help_issues.triggered.connect(
        lambda: w.help('https://github.com/calculix/cae/issues'))

    # treeView actions
    w.treeView.doubleClicked.connect(t.doubleClicked)
    w.treeView.clicked.connect(t.clicked)
    w.treeView.customContextMenuRequested.connect(t.rightClicked)
    w.treeView.expanded.connect(t.treeViewExpanded)
    w.treeView.collapsed.connect(t.treeViewCollapsed)

    # ToolBar actions
    w.action_view_minus_x.triggered.connect(lambda: w.connections[1].post('rot -x'))
    w.action_view_minus_y.triggered.connect(lambda: w.connections[1].post('rot -y'))
    w.action_view_minus_z.triggered.connect(lambda: w.connections[1].post('rot -z'))
    w.action_view_plus_x.triggered.connect(lambda: w.connections[1].post('rot x'))
    w.action_view_plus_y.triggered.connect(lambda: w.connections[1].post('rot y'))
    w.action_view_plus_z.triggered.connect(lambda: w.connections[1].post('rot z'))
    w.action_view_frame.triggered.connect(lambda: w.connections[1].post('frame'))

    # Workaround for iso view
    # Three rotation posts to CGX window doesn't work in Windows
    # So one may use .fbd commands
    def action_view_iso():
        file_name = os.path.join(p.config, 'iso.fbd')
        if not os.path.isfile(file_name):
            logging.error('No config file iso.fbd')
            return
        w.connections[1].post('read ' + file_name)
    w.action_view_iso.triggered.connect(action_view_iso)
    w.action_view_line.triggered.connect(lambda: w.connections[1].post('view elem off'))
    w.action_view_line.triggered.connect(lambda: w.connections[1].post('view line'))
    w.action_view_fill.triggered.connect(lambda: w.connections[1].post('view elem off'))
    w.action_view_fill.triggered.connect(lambda: w.connections[1].post('view fill'))
    w.action_view_elem.triggered.connect(lambda: w.connections[1].post('view fill'))
    w.action_view_elem.triggered.connect(lambda: w.connections[1].post('view elem'))
