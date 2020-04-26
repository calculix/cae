#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, December 2019
Distributed under GNU General Public License v3.0

Model unites KOM object and parsers.
See scheme on architecture.odp. """

from . import parsers
from . import job
from . import kom

import logging, os

class Model:

    def __init__(self):

        # Empty KOM w/o implementations
        self.KOM = kom.KOM()

        # Variable names below should exactly represent KOM group names
        self.Mesh = None
        self.Interactions = None
        self.Constraints = None
