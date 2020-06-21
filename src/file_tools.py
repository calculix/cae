#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, May 2020
Distributed under GNU General Public License v3.0

Methods to read files. """

import os
import logging

# Recurcively read all the lines of the file and its includes
def read_lines(INP_file):
    INP_file = os.path.abspath(INP_file)
    if not os.path.isfile(INP_file):
        msg_text = 'File not found: ' + INP_file
        logging.error(msg_text)
        return []

    lines = []
    with open(INP_file, 'r') as f:
        for line in f.readlines():
            line = line.strip()

            # Skip comments and empty lines
            if line.startswith('**') or len(line)==0:
                continue

            lines.append(line)

            # Append lines from include file
            if line.upper().startswith('*INCLUDE'):
                inc_file = line.split('=')[1].strip()
                inc_file = os.path.normpath(
                    os.path.join(os.path.dirname(INP_file), inc_file))
                lines.extend(read_lines(inc_file))

    return lines

# Test run
if __name__=='__main__':
    d = os.path.dirname(__file__)
    d = os.path.join(d, '..', 'examples')
    INP_file = os.path.join(d, 'Ihor_Mirzov_baffle_2D.inp')
    INP_file = os.path.normpath(INP_file)
    print(os.path.dirname(INP_file))
    read_lines(INP_file)
