#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, December 2019
    Distributed under GNU General Public License v3.0

    Main window actions - all processed signals.
"""


from PyQt5 import QtWidgets
import cae


"""
s - Settings
w - Window
m - Model
t - Tree
j - Job
"""
def actions(s, w, m, t, j):
    w.keyPressEvent = t.keyPressEvent

    # File actions
    w.action_file_import.triggered.connect(lambda: cae.import_file(s, w, m, t, j))
    w.action_file_settings.triggered.connect(s.open)
    w.action_file_exit.triggered.connect(QtWidgets.qApp.quit)

    # Job actions
    w.action_job_write_input.triggered.connect(
        lambda: j.write_input(m.KOM.get_INP_code_as_lines()))
    w.action_job_write_input.triggered.connect(
        lambda: w.setWindowTitle('CalculiX CAE - ' + j.name))
    w.action_job_edit_inp.triggered.connect(lambda: j.edit_inp(s))
    w.action_job_subroutine.triggered.connect(lambda: j.open_subroutine(s))
    w.action_job_rebuild_ccx.triggered.connect(lambda: j.rebuild_ccx(s))
    w.action_job_submit.triggered.connect(lambda: j.submit(s))
    w.action_job_view_log.triggered.connect(lambda: j.view_log(s))
    w.action_job_cgx_inp.triggered.connect(lambda: j.cgx_inp(s, w))
    w.action_job_cgx_frd.triggered.connect(lambda: j.cgx_frd(s, w))
    w.action_job_export_vtu.triggered.connect(j.export_vtu)
    w.action_job_paraview.triggered.connect(lambda: j.open_paraview(s))

    # Help actions
    w.action_help_readme.triggered.connect(
        lambda: w.help('https://github.com/imirzov/ccx_cae#calculix-cae'))
    w.action_help_issues.triggered.connect(
        lambda: w.help('https://github.com/imirzov/ccx_cae/issues'))

    # treeView actions
    w.treeView.doubleClicked.connect(t.doubleClicked)
    w.treeView.clicked.connect(t.clicked)
    w.treeView.customContextMenuRequested.connect(t.rightClicked)
    w.treeView.expanded.connect(t.treeViewExpanded)
    w.treeView.collapsed.connect(t.treeViewCollapsed)
