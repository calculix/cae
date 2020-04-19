#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, December 2019
    Distributed under GNU General Public License v3.0

    Main window actions - all processed signals.
"""



from PyQt5 import QtWidgets
import ie


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
    w.action_file_import.triggered.connect(lambda: ie.importFile(s, w, m, t, j))
    w.action_file_settings.triggered.connect(s.open)
    w.action_file_exit.triggered.connect(QtWidgets.qApp.quit)

    # Job actions
    w.action_job_write_input.triggered.connect(
        lambda: ie.writeInput(j, m.KOM.get_INP_code_as_lines()))
    w.action_job_write_input.triggered.connect(
        lambda: w.setWindowTitle('CalculiX CAE - ' + j.name))
    w.action_job_edit_inp.triggered.connect(lambda: j.editINP(s))
    w.action_job_open_subroutine.triggered.connect(lambda: j.openSubroutine(s))
    w.action_job_rebuild_ccx.triggered.connect(lambda: j.rebuildCCX(s))
    w.action_job_submit.triggered.connect(lambda: j.submit(s))
    w.action_job_view_log.triggered.connect(lambda: j.viewLog(s))
    w.action_job_open_cgx.triggered.connect(lambda: j.openCGX(s, w))
    w.action_job_export_vtu.triggered.connect(j.exportVTU)
    w.action_job_open_paraview.triggered.connect(lambda: j.openParaView(s))

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
