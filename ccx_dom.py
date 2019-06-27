# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, June 2019.
    Distributed under GNU General Public License, version 2.

    Data object model based on CalculiX keywords hierarchy.
    Keywords with all arguments are parsed from ccx_dom.txt.
    Parser supposts mutually exclusive parameters for keywords.
"""


import ccx_log, copy

# TODO enumerations for item_type


# CalculiX keywords hierarchy - data object model
class DOM:

    """
        logger          -   simply logger from ccx_log
        parent_items    -   last parent for given padding level
        root            -   group 'ROOT' from ccx_dom.inp
    """

    # Read CalculiX keywords hierarchy
    def __init__(self, textEdit):

        # Configure logging
        self.logger = ccx_log.logger(textEdit)

        try:
            # Last parent for given padding level
            parent_items = {}

            # Analyze keywords hierarchy
            with open('ccx_dom.inp', 'r') as f:
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
                        item = keyword(line)
                    elif line.endswith('__group'):
                        item = group(line)

                    # Update parent item
                    parent_items[level] = item
                    if level > 0:
                        parent_items[level-1].addItem(item) # for example, add keyword to group
                        item.setParent(parent_items[level-1]) # set item's parent

            self.root = parent_items[0] # group 'Model'
            self.logger.info('CalculiX object model generated.')
        except:
            self.logger.error('Can\'t generate keywords hierarchy!')


# Needed for inheritance by further classes
class item:

    item_type = ''      # item's type: group/keyword/argument/implementation
    name = ''           # name of item, string
    items = []          # list of children
    parent = None       # item's parent item
    expanded = True     # flag for treeView


    # Recursive function to count keyword implementations in item's descendants
    def countImplementations(self):
        if self.item_type == 'argument':
            return 0

        counter = 0
        for i in self.items:
            if i.item_type == 'implementation':
                counter += 1
            else:
                counter += i.countImplementations()
        return counter


    # Get list of keywords implementations
    def getImplementations(self):
        imps = []
        for item in self.items:
            if item.item_type == 'implementation':
                imps.append(item)
        return imps


    # Add child item to 'items' list
    def addItem(self, item):
        self.items.append(item)


    # Get list of items - non implementations
    def getItems(self):
        imps = []
        for item in self.items:
            if item.item_type != 'implementation':
                imps.append(item)
        return imps


    # Set item's parent
    def setParent(self, item):
        self.parent = item


    # Get item's parent
    def getParent(self):
        return self.parent

    """
        def getPath(self, path=None, parent=None):
            if not path:
                path = [self.name]

            if not parent:
                parent = self.getParent()

            if parent:
                path.append(parent.name)
                print(str(path))
                self.getPath(path, parent)

            return path
    """

    # Print all the branch DOM elements starting from current item
    def writeAll(self, f, level=0):

        # How many implementations has group or keyword
        has = ''
        counter = self.countImplementations()
        if counter:
            has = ' ' + str(counter)

        # Name of parent for current item
        parent = ''
        if self.parent:
            parent = ' parent=' + self.parent.name

        # path = ' '# + str(self.getPath())

        # String to write to file for debug purposes
        string = '\t'*level + self.name + has
        f.write(string + '\n')

        # Organize recursion
        for item in self.items:
            if item.item_type != 'argument':
                item.writeAll(f, level+1)


# Group of keywords, like 'Properties', 'Constraints', etc.
class group(item):

    def __init__(self, line):
        self.item_type = 'group'
        line = line.strip().replace('__group', '')
        self.name = line
        self.items = [] # list of groups and keywords


# *AMPLITUDE, *BOUNDARY, *STEP etc.
class keyword(item):

    def __init__(self, line):
        self.item_type = 'keyword'

        # '__' added after '=' for better code highlighting in vscode
        if line.endswith('__'):
            line = line[:-2]

        self.items = [] # list of arguments

        # Start all arguments from the next line?
        self.from_new_line = False # marks that argument should be printed from the new line
        if '__newline' in line:
            self.from_new_line = True # yes
            line = line.replace('__newline', '')

        # Define name
        if not ',' in line:
            self.name = line
        else:
            lines = line.split(',')
            self.name = lines[0]
            for arg in lines[1:]:
                self.items.append(argument(arg))

        # self.path = [] # TODO path to current keyword in the DOM hierarchy


# Keyword's argument
class argument(item):

    def __init__(self, line):

        self.item_type = 'argument' # needed to distinguish from 'group' and 'keyword'

        # Required or optional
        self.required = False
        if '__required' in line:
            self.required = True

        self.items = [] # list of strings
        left_part = line
        if '=' in line:
            # arguments : values
            left_part, right_part = line.split('=')
            left_part = left_part + '='

            # Define argument's values
            # TODO add '+' splitter for arguments that always go together
            if '|' in right_part:
                # if few values are present
                self.items = [v for v in right_part.split('|')]
            elif len(right_part):
                # one value only
                self.items = [right_part]

        # If '__required' is present - remove it
        if '__required' in left_part:
            left_part = left_part.replace('__required', '')

        # Define argument's name
        self.name = left_part


# Keyword implementation - a piece of INP-code for CalculiX input file
class implementation(item):

    def __init__(self, keyword, INP_code, name=None):
        self.item_type = 'implementation'
        self.items = copy.deepcopy(keyword.getItems())

        # Name of current implementation (of *AMPLITUDE, *STEP, *MATERIAL etc.)
        if name:
            self.name = name # it will be used in edit Dialog
            index = int(name.split('-')[1]) - 1
        else:
            index = len(keyword.getImplementations())
            self.name = keyword.name[1:] + '-' + str(index + 1)

        self.INP_code = INP_code # INP-code for current implementation - list of strings
        keyword.items.insert(index, self) # append implementation to keyword's items
