# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, UJV Rez, June 2019.
    Distributed under GNU General Public License, version 2.

    CalculiX CAE
    Data object model based on CalculiX keywords hierarchy
"""


class ccx_dom:

    dom = {0:'Model'}
    groups = []
    keywords = []

    # Read CalculiX keywords hierarchy
    def __init__(self):
        try:
            # Read the whole file
            with open('ccx_dom.txt', 'r') as f:
                lines = f.readlines()

            # Analyze keywords hierarchy
            for line in lines:

                # Skip comments and empty lines
                if line.strip() == '': continue 
                if line.strip().startswith('**'): continue

                level = 1
                line = line.rstrip() # cut '\n'
                while line.startswith('    '): # exactly 4 spaces
                    level += 1
                    line = line[4:]
                self.dom[level] = line
                # print(level, self.dom[level-1], line)

                if line.startswith('*'):
                    # keyword
                    keywords.append(keyword(line))
                    # print('keyword', line)
                    pass
                elif line.startswith('-'):
                    # argument
                    # print('argument', line)
                    pass
                else:
                    # group
                    # print('group', line)
                    pass
        except:
            print('\nError with keywords hierarchy!')


# Group of keywords
class group:
    name = ''
    keywords = []
    groups = []
    def __init__(line):
        self.name = line.strip()


# *AMPLITUDE, *BOUNDARY, *STEP etc.
class keyword:
    name = None
    arguments = [] # list of keyword's arguments
    parent = None
    groups = []
    keywords = []

    def __init__(line):
        pass


# Keyword's argument
class argument:
    name = ''
    values = [] # list of strings
    required = False


if __name__ == '__main__':
    ccx_dom()