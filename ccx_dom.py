# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, UJV Rez, June 2019.
    Distributed under GNU General Public License, version 2.

    Data object model based on CalculiX keywords hierarchy.
    Keywords with all arguments are read from ccx_dom.txt.
"""


# CalculiX keywords hierarchy - data object model
class DOM:

    # Read CalculiX keywords hierarchy
    def __init__(self):

        # Top level element in tree
        parent_items = {}

        # Analyze keywords hierarchy
        with open('ccx_dom.txt', 'r') as f:
            for line in f.readlines(): # read the whole file and iterate over line

                # Skip comments and empty lines
                if line.strip() == '': continue 
                if line.strip().startswith('**'): continue

                # Define padding level
                level = 0 # level of current item = amount of tabs
                line = line.rstrip() # cut '\n'
                while line.startswith('\t'):
                    level += 1
                    line = line[1:] # cut 1 tab

                # Define item for current padding level
                if line.startswith('*'):
                    item = keyword(line, level)
                elif line.startswith('-'):
                    item = argument(line, level)
                else:
                    item = group(line, level)

                # Update parent item
                parent_items[level] = item
                if level > 0:
                    parent_items[level-1].addItem(item)

        self.root = parent_items[0] # group 'Model'


# Group of keywords, like 'Properties', 'Constraints', etc.
class group:

    item_type = 'group' # needed to distinguish from 'keyword'

    def __init__(self, line, level):
        self.name = line.strip()
        self.items = [] # list of groups and keywords
        self.level = level # needed for padding in printAll()

    def addItem(self, item):
        self.items.append(item)
    
    def printAll(self):
        string = '\t' * self.level
        string += self.name
        # print(string)
        for item in self.items:
            string += item.printAll()
        return string


# *AMPLITUDE, *BOUNDARY, *STEP etc.
class keyword:

    item_type = 'keyword' # needed to distinguish from 'group'
    # TODO activate keyword after dialog edition

    def __init__(self, line, level):
        self.name = line.strip()
        self.items = [] # list of groups, keywords and arguments
        self.level = level # needed for padding in printAll()
    
    def addItem(self, item):
        self.items.append(item)

    def printAll(self):
        string = self.name + '\n'
        for item in self.items:
            if item.item_type != 'argument':
                continue
            string += item.printAll()
        return string


# Keyword's argument
class argument:

    """
        name        - argument's name
        required    - required or optional
        values      - list of possible values
    """
    items = []
    item_type = 'argument' # needed to distinguish from 'group' and 'keyword'

    def __init__(self, line, level):
        
        self.level = level # needed for padding in printAll()

        line = line[1:] # cut hyphen

        # Required or optional
        self.required = False
        if 'required' in line:
            self.required = True

        self.values = [] # list of strings
        left_part = line
        if ':' in line:
            # arguments : values
            left_part, right_part = line.split(':')

            # Define argument's values
            if '|' in right_part:
                # if few values are present
                self.values = [v for v in right_part.split('|')]
            else:
                # one value only
                self.values = [right_part]

        # Define argument's name
        if ',' in left_part:
            # TODO improve to parse '|'
            
            # if 'required' is present
            self.name = left_part.split(',')[0]
        else:
            # argument name only
            self.name = left_part


    def printAll(self):
        string = '-' + self.name

        if self.required:
            string += ',required'

        amount_of_values = len(self.values)
        if amount_of_values:
            string += ':'
            for i in range(amount_of_values-1):
                string += self.values[i] + '|'
            string += self.values[amount_of_values-1]

        return string + '\n'
