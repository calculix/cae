#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2023
Distributed under GNU General Public License v3.0

Model unites KOM object and parsers.
The FEM model definition will be here.
See scheme on architecture.odp.
"""


class Model:

    def __init__(self):
        # Variable names below should exactly represent KOM group names
        self.Mesh = None
        self.Interactions = None
        self.Constraints = None


m = Model()
