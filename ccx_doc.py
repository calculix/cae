# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    Dialog window to show PDF help during keyword edit/create.
"""


import os
from PyQt5 import QtWebEngineWidgets, QtCore


class Window(QtWebEngineWidgets.QWebEngineView):

    def __init__(self, item):
        super(Window, self).__init__()

        script = os.path.realpath(__file__)
        folder = os.path.dirname(script)
        PDFJS = 'file://'+folder+'/doc/pdfjs-2.0.943-dist/web/viewer.html'
        PDF = 'file://'+folder+'/doc/ccx_2.15.pdf'
        self.load(QtCore.QUrl.fromUserInput('%s?file=%s' % (PDFJS, PDF)))
