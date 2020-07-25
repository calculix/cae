#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, July 2020
Distributed under GNU General Public License v3.0

Converts UNV file from Salome to CalculiX INP mesh:
python3 ./src/unv2ccx.py ./examples/116.unv

Reads UNV_file, creates an internal FEM object,
then writes the INP_file. """

import os
import sys
import argparse
import logging

sys_path = os.path.dirname(__file__)
sys.path.append(sys_path)
import clean
import UNVParser
import INPWriter


class Converter:

    def __init__(self, unv_file_name):
        self.unv_file_name = os.path.normpath(unv_file_name)
        self.inp_file_name = self.unv_file_name[:-4]+'.inp'

    def run(self):

        # Parse UNV file
        relpath = os.path.relpath(self.unv_file_name,
            start=os.path.dirname(__file__))
        logging.info('Parsing ' + relpath)
        fem = UNVParser.UNVParser(self.unv_file_name).parse()

        # Write INP file
        relpath = os.path.relpath(self.inp_file_name,
            start=os.path.dirname(__file__))
        logging.info('Writing ' + relpath)
        INPWriter.write(fem, self.inp_file_name)


if __name__ == '__main__':
    clean.screen()

    # Configure logging
    logging.basicConfig(level=logging.INFO,
        format='%(levelname)s: %(message)s')

    # Command line parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str,
        help='UNV file name with extension')
    args = parser.parse_args()

    # Create converter and run it
    unv2ccx = Converter(args.filename)
    unv2ccx.run()

    clean.cache()
