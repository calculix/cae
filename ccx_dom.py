# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, UJV Rez, June 2019.
    Distributed under GNU General Public License, version 2.

    CalculiX CAE
    Data object model based on CalculiX keywords hierarchy.
    Keywords with all parameters are read from ccx_dom.txt.
"""


class ccx_dom:

    # dom = {0:'Model'} # collections.OrderedDict


    # Read CalculiX keywords hierarchy
    def __init__(self):

        self.items = [] # list of groups and keywords

        # Read the whole file
        with open('ccx_dom.txt', 'r') as f:
            lines = f.readlines()

        # Analyze keywords hierarchy
        level_of_parent_item = -1
        parent_item = None
        for line in lines:

            # Skip comments and empty lines
            if line.strip() == '': continue 
            if line.strip().startswith('**'): continue

            # Define level
            level = 0 # level of current item = amount of tabs
            line = line.rstrip() # cut '\n'
            while line.startswith('\t'):
                level += 1
                line = line[1:] # cut 1 tab

            # self.dom[level] = line
            # print(level, self.dom[level-1], line)

            if line.startswith('*'):
                # Keyword
                item = keyword(line)
                self.items.append(item)

                # Update parent if level decreased
                if level < level_of_parent_item:
                    parent_item = item

            elif line.startswith('-'):
                # Argument
                pass
            else:
                # Group
                item = group(line)
                self.items.append(item)

                # Update parent if level decreased
                if level < level_of_parent_item:
                    parent_item = item
            
            # Update parent for the next itration
            level_of_parent_item = level


# Group of keywords, like 'Properties', 'Constraints', etc.
class group:
    item_type = 'group'
    def __init__(line):
        self.name = line.strip()
        self.items = [] # list of groups and keywords


# *AMPLITUDE, *BOUNDARY, *STEP etc.
class keyword:
    item_type = 'keyword'
    def __init__(line):
        self.name = line.strip()
        self.arguments = [] # list of keyword's arguments
        self.items = [] # list of groups and keywords


# Keyword's argument
class argument:
    def __init__(line, keyword):
        line = line[1:] # cut hyphen
        if ':' in line:
            pass
        else:
            pass

        self.name = line.strip()
        self.values = [] # list of strings
        self.required = False


if __name__ == '__main__':
    ccx_dom()