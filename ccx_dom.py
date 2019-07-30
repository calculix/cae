# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, July 2019.
    Distributed under GNU General Public License, version 2.

    Data object model based on CalculiX keywords hierarchy.
    Keywords with all arguments are parsed from ccx_dom.txt.
    Parser supposts mutually exclusive parameters for keywords.
"""


import copy, re, ccx_cae_log
from enum import Enum


# CalculiX keywords hierarchy - data object model
class DOM:

    """
        parent_items    -   last parent for given padding level
        root            -   group 'ROOT' from ccx_dom.inp
        pathes          -   all possible keywords nesting variants
    """


    # Read CalculiX keywords hierarchy
    def __init__(self):
        self.msg_list = [] # list of messages for logger

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

            # All possible keywords nesting variants - needed for parsing INP_doc
            self.pathes = []
            self.buildPathes(self.root)
            self.pathes.sort(key=self.keyword_counter, reverse=True) # maximum nesting first

            msg_text = 'CalculiX object model generated.'
            msg = ccx_cae_log.msg(ccx_cae_log.msgType.INFO, msg_text)
            self.msg_list.append(msg)
        except:
            msg_text = 'Can\'t generate keywords hierarchy!'
            msg = ccx_cae_log.msg(ccx_cae_log.msgType.ERROR, msg_text)
            self.msg_list.append(msg)


    # Recursively builds all possible pathes to nested keywords in DOM
    def buildPathes(self, parent, path=None):
        if not path:
            path = []
        for item in parent.items:
            if (item.item_type != item_type.ARGUMENT):
                self.buildPathes(item, path + [item])
        if len(path):
            if path not in self.pathes:
                self.pathes.append(path)


    # Get nesting path for each of the parsed INP_doc keyword
    def getPath(self, keyword_chain):

        # Modify keyword_chain to remove duplicated words in the end
        if len(keyword_chain) > 1 and \
            keyword_chain[-1] == keyword_chain[-2]:
            del keyword_chain[-2]
        # print(keyword_chain)

        # Now compare keyword_chain with all self.pathes
        for path in self.pathes:

            # Compare last words
            if path[-1].name.lower() != keyword_chain[-1]:
                continue

            matches = 0
            minimum_j = 0
            for i in range(1, len(path)+1):
                for j in range(1, len(keyword_chain)+1):
                    if path[-i].name.lower() == keyword_chain[-j]:
                        matches += 1
                        minimum_j = len(keyword_chain)-j
                        continue

            # If we found all words from path in keyword_chain = if needed path is found
            if matches >= self.keyword_counter(path):
                # print([item.name for item in path], '\n')
                del keyword_chain[:minimum_j]
                return path


    # Count keywords in path
    def keyword_counter(self, path):
        keyword_counter = 0
        for item in path:
            if item.name.startswith('*'):
                keyword_counter += 1
        return keyword_counter


# Enums for 'item_type' variable
class item_type(Enum):
    GROUP = 0
    KEYWORD = 1
    ARGUMENT = 2
    IMPLEMENTATION = 3


# Needed for inheritance by further classes
class item:

    item_type = ''      # item's type: group/keyword/argument/implementation
    name = ''           # name of item, string
    items = []          # list of children
    parent = None       # item's parent item
    expanded = True     # flag for treeView


    # Recursive function to count keyword implementations in item's descendants
    def countImplementations(self):
        if self.item_type == item_type.ARGUMENT:
            return 0

        counter = 0
        for i in self.items:
            if i.item_type == item_type.IMPLEMENTATION:
                counter += 1
            else:
                counter += i.countImplementations()
        return counter


    # Get list of keywords implementations
    def getImplementations(self):
        imps = []
        for item in self.items:
            if item.item_type == item_type.IMPLEMENTATION:
                imps.append(item)
        return imps


    # Add child item to 'items' list
    def addItem(self, item):
        self.items.append(item)


    # Get list of items - non implementations
    def getItems(self):
        imps = []
        for item in self.items:
            if item.item_type == item_type.GROUP or \
                item.item_type == item_type.KEYWORD:
                imps.append(item)
        return imps


    # Set item's parent
    def setParent(self, item):
        self.parent = item


    # Get item's parent
    def getParent(self):
        return self.parent


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

        # String to write to file for debug purposes
        string = '\t'*level + self.name + has
        f.write(string + '\n')

        # Organize recursion
        for item in self.items:
            if item.item_type != item_type.ARGUMENT:
                item.writeAll(f, level+1)


    # Search item by name among children
    def getItemByName(self, name):
        for item in self.items:
            if item.name == name:
                return item


# Group of keywords, like 'Properties', 'Constraints', etc.
class group(item):

    def __init__(self, line):
        self.item_type = item_type.GROUP
        line = line.strip().replace('__group', '')
        self.name = line
        self.items = [] # list of groups and keywords


# *AMPLITUDE, *BOUNDARY, *STEP etc.
class keyword(item):

    def __init__(self, line):
        self.item_type = item_type.KEYWORD

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


# Keyword's argument
class argument(item):

    def __init__(self, line):

        self.item_type = item_type.ARGUMENT

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
        self.item_type = item_type.IMPLEMENTATION
        self.items = copy.deepcopy(keyword.getItems())
        self.parent = keyword

        # Name of current implementation (of *AMPLITUDE, *STEP, *MATERIAL etc.)
        index = len(keyword.getImplementations())
        if name:
            self.name = name # it will be used in edit Dialog
            # index = int(name.split('-')[1]) - 1
        else:
            try:
                match = re.search('(NAME|ELSET|NSET)\s*=\s*(\w*)', INP_code[0].upper())
                self.name = match.group(2)
            except:
                self.name = keyword.name[1:] + '-' + str(index + 1)

        self.INP_code = INP_code # INP-code for current implementation - list of strings
        keyword.items.insert(index, self) # append implementation to keyword's items
