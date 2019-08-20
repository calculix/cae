# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    Application's global settings
"""


import logging


class Settings:


    def __init__(self):
        self.file_name = 'ccx_settings.env'
        f = open(self.file_name).read()
        self.lines = f.split('\n')
        exec(f)

    def save(self):
        with open(self.file_name, 'w') as f:
            counter = 0
            for line in self.lines:
                if line.startswith('self.'):
                    param, value = line[5:].split('=')
                    param = param.strip()
                    value = value.strip()
                    if value.startswith('\'') and value.endswith('\''):
                        value = '\'' + getattr(self, param) + '\''
                    else:
                        value = getattr(self, param)
                    line = 'self.{} = {}'.format(param, value)
                if counter:
                    f.write('\n')
                f.write(line)
                counter += 1
        logging.info('Settings saved.')
