#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Model unites KOM object and parsers.
See scheme on architecture.odp. """

import os
import logging

from . import parsers
from . import job
from . import kom

class Model:

    def __init__(self):

        # Empty KOM w/o implementations
        self.KOM = None

        # Variable names below should exactly represent KOM group names
        self.Mesh = None
        self.Interactions = None
        self.Constraints = None
