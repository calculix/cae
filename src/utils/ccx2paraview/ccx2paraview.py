#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, February 2020
Distributed under GNU General Public License v3.0

Converts CalculiX .frd resutls file to ASCII .vtk or XML .vtu format:
python3 ./src/ccx2paraview.py ./examples/other/Ihor_Mirzov_baffle_2D.frd vtk
python3 ./src/ccx2paraview.py ./examples/other/Ihor_Mirzov_baffle_2D.frd vtu

TODO It would be a killer feature if Paraview could
visualize gauss point results from the dat file...
https://public.kitware.com/pipermail/paraview/2013-January/027121.html

TODO Parse DAT files - there are lots of results

TODO XDMF format:
https://github.com/calculix/ccx2paraview/issues/6

"""

import os
import sys
import logging
import argparse

sys_path = os.path.dirname(__file__)
sys.path.append(sys_path)
import FRDParser
import VTKWriter
import VTUWriter
import PVDWriter
import clean


class Converter:

    def __init__(self, file_name, fmt):
        self.file_name = file_name
        self.fmt = fmt

    def run(self):

        # Parse FRD-file
        relpath = os.path.relpath(self.file_name,
            start=os.path.dirname(__file__))
        logging.info('Parsing ' + relpath)
        p = FRDParser.Parse(self.file_name)

        # If file contains mesh data
        if p.node_block and p.elem_block:
            times = sorted(set([b.value for b in p.result_blocks]))
            l = len(times)
            if l:
                msg = '{} time increment{}'.format(l, 's'*min(1, l-1))
                logging.info(msg)

                """ If model has many time steps - many output files
                will be created. Each output file's name should contain
                increment number padded with zero """
                counter = 1
                times_names = {} # {increment time: file name, ...}
                for t in sorted(times):
                    if l > 1:
                        ext = '.{:0{width}}.{}'\
                            .format(counter, self.fmt, width=len(str(l)))
                        file_name = self.file_name.replace('.frd', ext)
                    else:
                        ext = '.{}'.format(self.fmt)
                        file_name = self.file_name.replace('.frd', ext)
                    times_names[t] = file_name
                    counter += 1

                # For each time increment generate separate .vt* file
                # Output file name will be the same as input
                for t, file_name in times_names.items():
                    relpath = os.path.relpath(file_name,
                        start=os.path.dirname(__file__))
                    logging.info('Writing ' + relpath)
                    if self.fmt == 'vtk':
                        VTKWriter.writeVTK(p, file_name, t)
                    if self.fmt == 'vtu':
                        VTUWriter.writeVTU(p, file_name, t)

                # Write ParaView Data (PVD) for series of VTU files
                if l > 1 and self.fmt == 'vtu':
                    PVDWriter.writePVD(self.file_name\
                        .replace('.frd', '.pvd'), times_names)

            else:
                logging.warning('No time increments!')
                file_name = self.file_name[:-3] + self.fmt
                if self.fmt == 'vtk':
                    VTKWriter.writeVTK(p, file_name, None)
                if self.fmt == 'vtu':
                    VTUWriter.writeVTU(p, file_name, None)
        else:
            logging.warning('File is empty!')


if __name__ == '__main__':
    clean.screen()

    # Configure logging
    logging.basicConfig(level=logging.INFO,
        format='%(levelname)s: %(message)s')

    # Command line parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str,
        help='FRD file name with extension')
    parser.add_argument('format', type=str,
        help='output format: vtu or vtk')
    args = parser.parse_args()

    # Create converter and run it
    if args.format in ['vtk', 'vtu']:
        ccx2paraview = Converter(args.filename, args.format)
        ccx2paraview.run()
    else:
        msg = 'ERROR! Wrong format \"{}\". '.format(args.format) \
            + 'Choose one of: vtk, vtu.'
        print(msg)

    clean.cache()
