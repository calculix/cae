#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, 2019-2020
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
    with open(INP_file, 'r', errors='ignore') as f:
        for line in f.readlines():
            line = line.strip()
            lines.append(line)

            # Append lines from include file
            if line.upper().startswith('*INCLUDE'):
                inc_file = line.split('=')[1].strip()
                inc_file = os.path.normpath(
                    os.path.join(os.path.dirname(INP_file), inc_file))
                lines.extend(read_lines(inc_file))

    return lines

# Run test
if __name__ == '__main__':
    d = os.path.dirname(__file__)
    d = os.path.join(d, '..', 'examples')
    INP_file = os.path.join(d, 'Ihor_Mirzov_baffle_2D.inp')
    INP_file = os.path.normpath(INP_file)
    print(os.path.dirname(INP_file))
    read_lines(INP_file)
