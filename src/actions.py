#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, December 2019
    Distributed under GNU General Public License v3.0

    Main window actions.
"""



from PyQt5 import QtWidgets
from ie import importFile, writeInput



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
    w.action_file_import.triggered.connect(lambda: importFile(s, w, m, t, j))
    w.action_file_settings.triggered.connect(s.open)
    w.action_file_exit.triggered.connect(QtWidgets.qApp.quit)

    # Job actions
    w.action_job_write_input.triggered.connect(lambda: writeInput(m, j))
    w.action_job_edit_inp.triggered.connect(lambda: j.editINP(s))
    w.action_job_open_subroutine.triggered.connect(lambda: j.openSubroutine(s))
    w.action_job_rebuild_ccx.triggered.connect(lambda: j.rebuildCCX(s))
    w.action_job_submit.triggered.connect(lambda: j.submit(s))
    w.action_job_view_log.triggered.connect(lambda: j.viewLog(s))
    w.action_job_open_cgx.triggered.connect(lambda: j.openCGX(s))
    w.action_job_export_vtu.triggered.connect(j.exportVTU)
    w.action_job_open_paraview.triggered.connect(lambda: j.openParaView(s))

    # Help actions
    w.action_help_readme.triggered.connect(lambda:
            w.help('https://github.com/imirzov/ccx_cae#calculix-cae'))
    w.action_help_issues.triggered.connect(lambda:
            w.help('https://github.com/imirzov/ccx_cae/issues'))

    # VTK actions
    if s.show_vtk:
        # w.actionSelectionNodes.triggered.connect(w.VTK.actionSelectionNodes)
        # w.actionSelectionElements.triggered.connect(w.VTK.actionSelectionElements)
        # w.actionSelectionClear.triggered.connect(w.VTK.actionSelectionClear)
        w.actionViewParallel.triggered.connect(w.VTK.actionViewParallel)
        w.actionViewPerspective.triggered.connect(w.VTK.actionViewPerspective)
        w.actionViewFront.triggered.connect(w.VTK.actionViewFront)
        w.actionViewBack.triggered.connect(w.VTK.actionViewBack)
        w.actionViewTop.triggered.connect(w.VTK.actionViewTop)
        w.actionViewBottom.triggered.connect(w.VTK.actionViewBottom)
        w.actionViewLeft.triggered.connect(w.VTK.actionViewLeft)
        w.actionViewRight.triggered.connect(w.VTK.actionViewRight)
        w.actionViewIso.triggered.connect(w.VTK.actionViewIso)
        w.actionViewFit.triggered.connect(w.VTK.actionViewFit)
        w.actionViewWireframe.triggered.connect(w.VTK.actionViewWireframe)
        w.actionViewSurface.triggered.connect(w.VTK.actionViewSurface)
        w.actionViewSurfaceWithEdges.triggered.connect(w.VTK.actionViewSurfaceWithEdges)

    # treeView actions
    w.treeView.doubleClicked.connect(t.doubleClicked)
    w.treeView.clicked.connect(t.clicked)
    w.treeView.customContextMenuRequested.connect(t.rightClicked)
    w.treeView.expanded.connect(t.treeViewExpanded)
    w.treeView.collapsed.connect(t.treeViewCollapsed)
