# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, December 2019
    Distributed under GNU General Public License v3.0

    Main window actions.
"""

from ie import importFile, writeInput
from PyQt5 import QtWidgets


class Actions:


    def __init__(self, s, mw, m, t):
        mw.treeView.keyPressEvent = mw.keyPressEvent

        # File actions
        mw.action_file_import.triggered.connect(lambda: importFile(s, mw, m, t))
        mw.action_file_settings.triggered.connect(s.open)
        mw.action_file_exit.triggered.connect(QtWidgets.qApp.quit)

        # Job actions
        mw.action_job_write_input.triggered.connect(lambda: writeInput(m))
        mw.action_job_edit_inp.triggered.connect(lambda: m.job.editINP(s))
        mw.action_job_open_subroutine.triggered.connect(lambda: m.job.openSubroutine(s))
        mw.action_job_rebuild_ccx.triggered.connect(lambda: m.job.rebuildCCX(s))
        mw.action_job_submit.triggered.connect(lambda: m.job.submit(s))
        mw.action_job_view_log.triggered.connect(lambda: m.job.viewLog(s))
        mw.action_job_open_cgx.triggered.connect(lambda: m.job.openCGX(s))
        mw.action_job_export_vtu.triggered.connect(m.job.exportVTU)
        mw.action_job_open_paraview.triggered.connect(lambda: m.job.openParaView(s))

        # Help actions
        mw.action_help_readme.triggered.connect(lambda:
                mw.help('https://github.com/imirzov/ccx_cae#calculix-cae'))
        mw.action_help_yahoo.triggered.connect(lambda:
                mw.help('https://groups.yahoo.com/neo/groups/CALCULIX/conversations/topics/15616'))
        mw.action_help_issues.triggered.connect(lambda:
                mw.help('https://github.com/imirzov/ccx_cae/issues'))

        # VTK actions
        if s.show_vtk:
            # mw.actionSelectionNodes.triggered.connect(mw.VTK.actionSelectionNodes)
            # mw.actionSelectionElements.triggered.connect(mw.VTK.actionSelectionElements)
            # mw.actionSelectionClear.triggered.connect(mw.VTK.actionSelectionClear)
            mw.actionViewParallel.triggered.connect(mw.VTK.actionViewParallel)
            mw.actionViewPerspective.triggered.connect(mw.VTK.actionViewPerspective)
            mw.actionViewFront.triggered.connect(mw.VTK.actionViewFront)
            mw.actionViewBack.triggered.connect(mw.VTK.actionViewBack)
            mw.actionViewTop.triggered.connect(mw.VTK.actionViewTop)
            mw.actionViewBottom.triggered.connect(mw.VTK.actionViewBottom)
            mw.actionViewLeft.triggered.connect(mw.VTK.actionViewLeft)
            mw.actionViewRight.triggered.connect(mw.VTK.actionViewRight)
            mw.actionViewIso.triggered.connect(mw.VTK.actionViewIso)
            mw.actionViewFit.triggered.connect(mw.VTK.actionViewFit)
            mw.actionViewWireframe.triggered.connect(mw.VTK.actionViewWireframe)
            mw.actionViewSurface.triggered.connect(mw.VTK.actionViewSurface)
            mw.actionViewSurfaceWithEdges.triggered.connect(mw.VTK.actionViewSurfaceWithEdges)

        # treeView actions
        mw.treeView.doubleClicked.connect(t.doubleClicked)
        mw.treeView.clicked.connect(t.clicked)
        mw.treeView.customContextMenuRequested.connect(t.rightClicked)
        mw.treeView.expanded.connect(t.treeViewExpanded)
        mw.treeView.collapsed.connect(t.treeViewCollapsed)
