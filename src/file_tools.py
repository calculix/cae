#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, January 2020
    Distributed under GNU General Public License v3.0

    Methods to read files.
"""


import os, logging


# Recurcively read all the lines of the file and its includes
def readLines(INP_file, include=False):
    lines = []
    INP_file = os.path.abspath(INP_file) # full path
    if os.path.isfile(INP_file):
        with open(INP_file, 'rb') as f:
            line = readByteLine(f)
            while line != None:

                # Skip comments and empty lines
                if (not line.startswith('**')) and len(line):
                    lines.append(line)

                    # Append lines from include file
                    if include and line.upper().startswith('*INCLUDE'):
                        inc_file = line.split('=')[1].strip()
                        inc_file = os.path.join(os.path.dirname(INP_file),
                                        os.path.basename(inc_file)) # file name with path
                        lines.extend(readLines(inc_file))

                line = readByteLine(f)
    else:
        msg_text = 'File not found: ' + INP_file
        logging.error(msg_text)

    return lines


# Read byte line and decode: return None after EOF
def readByteLine(f):

    # Check EOF
    byte = f.read(1)
    if not byte:
        return None

    # Convert first byte
    try:
        line = byte.decode()
    except UnicodeDecodeError:
        line = ' ' # replace endecoded symbols with space

    # Continue reading until EOF or new line
    while byte != b'\n':
        byte = f.read(1)
        if not byte:
            return line.strip() # EOF
        try:
            line += byte.decode()
        except UnicodeDecodeError:
            line += ' ' # replace endecoded symbols with space

    return line.strip()
