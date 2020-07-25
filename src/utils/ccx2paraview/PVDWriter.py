#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, October 2019
Distributed under GNU General Public License v3.0

Writes ParaView Data (PVD) for series of VTU files.
https://www.paraview.org/Wiki/ParaView/Data_formats#PVD_File_Format """

import os

# times and names lists should be the same length
def writePVD(file_name, times_names):
    with open(file_name, 'w') as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
        f.write('\t<Collection>\n')

        # Now Paraview doesn't support timestep formatting
        for t, file_name in times_names.items():
            f.write('\t\t<DataSet file="{}" timestep="{}"/>\n'\
                .format(os.path.basename(file_name), t))

        f.write('\t</Collection>\n')
        f.write('</VTKFile>')
