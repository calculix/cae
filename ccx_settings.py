# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    Application's global settings
"""


import os, sys, logging, re


class Settings:


    def __init__(self):
        self.file_name = os.path.join(os.path.dirname(sys.argv[0]), 'ccx_settings.env') # full path
        f = open(self.file_name).read()
        self.lines = f.split('\n')
        exec(f)

    def save(self):
        with open(self.file_name, 'w') as f:
            new_line = False
            for line in self.lines:
                match = re.search('^self\.(\S+)\s*=\s*(.+)', line)
                if match:
                    param = match.group(1)
                    value = match.group(2)
                    if value.startswith('\'') and value.endswith('\''):
                        value = '\'' + getattr(self, param) + '\''
                    else:
                        value = getattr(self, param)
                    line = 'self.{} = {}'.format(param, value)
                if new_line:
                    f.write('\n')
                f.write(line)
                new_line = True
        logging.info('Settings saved.')
