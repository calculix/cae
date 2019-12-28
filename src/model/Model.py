# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, December 2019
    Distributed under GNU General Public License v3.0

    Model unites KOM object, Job object and parsers.
    See scheme on ccx_cae.odp
"""


import logging, os
from model.KOM import KOM
from model.Job import Job


class Model:


    def __init__(self, settings, path_to_inp):
        self.mesh = None # mesh from .inp-file - will be parsed in ie.py
        self.KOM = KOM() # empty KOM w/o implementations
        # TODO try to reorder to omit double job calling/renaming
        self.job = Job(settings, path_to_inp) # create job object
