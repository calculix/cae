# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, October 2019
    Distributed under GNU General Public License v3.0

    Main window class.
"""


import logging, os
from kom import KOM
from job import Job


# Main window
class Model:


    # Create main window
    def __init__(self, settings, path_to_inp):
        self.mesh = None # mesh from .inp-file - will be parsed in cae_ie.py
        # self.IE = IE(self) # import/export of .inp-file
        self.KOM = KOM() # empty KOM w/o implementations
        # TODO try to reorder to omit double job calling/renaming
        self.job = Job(settings, path_to_inp) # create job object
