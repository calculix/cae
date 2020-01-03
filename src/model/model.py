#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, December 2019
    Distributed under GNU General Public License v3.0

    Model unites KOM object, Job object and parsers.
    See scheme on ccx_cae.odp
"""


import logging, os
from model.kom import KOM


class Model:


    def __init__(self):
        self.KOM = KOM() # empty KOM w/o implementations

        # Variable names below should exactly represent KOM group names
        self.Mesh = None
        self.Interactions = None
        self.Constraints = None
