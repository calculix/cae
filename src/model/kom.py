#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, June 2020
Distributed under GNU General Public License v3.0

CalculiX Keyword Object Model (hierarchy).
Parser supposts mutually exclusive parameters for keywords.
Keywords with all arguments are parsed from config/kom.xml.
On INP import KOM is enriched with keyword implementations. """

# Standard modules
import re
import os
import sys
import time
import logging
import copy
from enum import Enum
import xml.etree.ElementTree as ET
import traceback

# My modules
try:
    # Normal run
    import path
except:
    # Test run
    sys_path = os.path.join(os.path.dirname(__file__), '..')
    os.sys.path.append(sys_path)
    import path
    import clean
    import settings


# Keyword Object Model
class KOM:


    # Read CalculiX keywords hierarchy
    def __init__(self, p, s):
        self.s = s

        # List of all existing keywords
        self.keywords = []

        # All possible keywords nesting variants - needed for parsing INP_doc
        self.paths = []

        try:
            self.root = group(s) # group 'Model' from kom.xml

            # Parse keywords hierarchy and build KOM
            t = ET.parse(p.kom_xml)
            self.build(t.getroot(), self.root)

            self.build_paths(self.root)
            self.paths.sort(key=self.keyword_counter, reverse=True) # maximum nesting first
            # for path in self.paths:
            #     logging.debug(str([item.name for item in path]))

            logging.info('Keywords object model generated.')
        except:
            logging.error('Can\'t generate keywords object model!')
            logging.error(traceback.format_exc())


    # Recursively build Keyword Object Model
    def build(self, xml_branch, parent):
        for xml_child in xml_branch:
            logging.log(5, xml_child.tag + ', ' + str(xml_child.attrib))
            """
            xml_child.tag, xml_child.attrib:

            keyword, {'name': '*RESTART'}
            keyword, {'name': '*TIME POINTS'}
            group, {'name': 'Interactions'}
            group, {'name': 'Constraints'}
            """

            # Create item: group, keyword or argument
            klass = globals()[xml_child.tag]
            item = klass(self.s)
            parent.items.append(item)

            # Set fields for object
            setattr(item, 'parent', parent) # set parent
            for k,v in xml_child.attrib.items():
                setattr(item, k, v)

            # Fill self.keywords
            if xml_child.tag=='keyword':
                self.keywords.append(item)

            # Argument values for QComboBox
            if xml_child.tag=='argument' and xml_child.text:
                values = xml_child.text.strip()
                if len(values):
                    values = values.split('|')
                    item.items.extend(values)

            self.build(xml_child, item)


    # Recursively builds all possible paths to nested keywords in KOM
    def build_paths(self, parent, path=None):
        if not path:
            path = [] # list of groups, keywords and implementations
        for item in parent.items:
            if item.item_type != item_type.ARGUMENT:
                self.build_paths(item, path + [item])
        if len(path):
            if path not in self.paths:
                self.paths.append(path)


    # Get nesting path for each of the parsed keyword
    def get_path(self, keyword_chain):

        # Duplicated words in the end of keyword_chain play no role
        if len(keyword_chain) > 1 and \
            keyword_chain[-1] == keyword_chain[-2]:
            del keyword_chain[-2]
        msg = 'keyword_chain: ' + ', '.join(keyword_chain)

        # Now compare keyword_chain with all self.paths
        for path in self.paths:

            # Compare last words
            if path[-1].name.upper() != keyword_chain[-1].upper():
                continue

            matches = 0
            minimum_j = 0
            for i in range(1, len(path)+1):
                for j in range(1, len(keyword_chain)+1):
                    if path[-i].name.upper() == keyword_chain[-j].upper():
                        matches += 1
                        minimum_j = len(keyword_chain)-j
                        continue

            # If we found all words from path in keyword_chain = if needed path is found
            if matches >= self.keyword_counter(path):
                # msg += str([item.name for item in path]) + '\n'
                del keyword_chain[:minimum_j]
                return path, msg

        return None, msg


    # Count keywords in path
    def keyword_counter(self, path):
        keyword_counter = 0
        for item in path:
            if item.name.startswith('*'):
                keyword_counter += 1
        return keyword_counter


    # Get keyword from self.keywords by its name
    def get_keyword_by_name(self, name, parent=None):
        for kw in self.keywords:
            if kw.name == name:
                return kw


    # Recursively get whole model's INP_code as list of strings (lines)
    # Parent is KOM item, level defines code folding/padding
    def get_inp_code_as_lines(self, parent=None, level=0):
        lines = []
        if not parent:
            parent = self.root
        if parent.item_type == item_type.IMPLEMENTATION:
            level += 1

        # For each group/keyword from KOM
        for item in parent.items:
            if item.item_type == item_type.ARGUMENT:
                continue
            if item.item_type == item_type.IMPLEMENTATION:
                # INP_code is stripped
                lines.append(' '*4*level + item.INP_code[0] + '\n')
                for line in item.INP_code[1:]:
                    lines.append(' '*4*(level+1) + line + '\n')

            # Continue call iterator until dig to implementation
            lines.extend(self.get_inp_code_as_lines(item, level))

        return lines


# Enums for 'item_type' variable
class item_type(Enum):
    GROUP = 0
    KEYWORD = 1
    ARGUMENT = 2
    IMPLEMENTATION = 3


# Needed for inheritance by further classes
class item:
    item_type = ''          # item's type: group/keyword/argument/implementation
    name = ''               # name of item, string
    items = []              # list of children
    parent = None           # item's parent item
    active = False

    # Define if item is active
    def isActive(self):

        # Top groups
        if self.parent and self.parent.name == 'Model':
            self.active = True
            return True

        # Children of active groups
        elif self.parent and self.parent.active and self.parent.item_type in \
            [self.item_type.GROUP, self.item_type.IMPLEMENTATION]:
            self.active = True
            return True

        # Implementations are always active
        elif self.item_type == self.item_type.IMPLEMENTATION:
            self.active = True
            return True

        else:
            self.active = False
            return False

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
    def copyItems(self):
        items = []
        for item in self.items:
            item.active = True
            if item.item_type in [item_type.GROUP, item_type.KEYWORD]:
                item = copy.copy(item)
                item.items = item.copyItems()
                items.append(item)
            if item.item_type == item_type.ARGUMENT:
                item = copy.copy(item)
                items.append(item)
        return items

    # Search item by name among children
    def getItemByName(self, name):
        for item in self.items:
            if item.name == name:
                return item

    # Returns first active preceding/parent keyword (not group)
    def getParentKeywordName(self):
        if (self.parent.item_type == self.item_type.KEYWORD or \
            self.parent.name == 'Model') and self.parent.active:
                return self.parent.name
        else:
            return self.parent.getParentKeywordName()

    # # Returns preceding/parent group (not keyword)
    # def getParentGroupName(self):
    #     if (self.parent.item_type == self.item_type.GROUP or \
    #         self.parent.name == 'Model'):
    #             return self.parent.name
    #     else:
    #         return self.parent.getParentGroupName()


# Group of keywords, like 'Properties', 'Constraints', etc.
class group(item):

    def __init__(self, s):
        self.item_type = item_type.GROUP
        self.items = [] # list of groups and keywords
        self.name = 'Model' # default name (root group)
        self.expanded = s.expanded


# *AMPLITUDE, *BOUNDARY, *STEP etc.
class keyword(item):

    def __init__(self, s):
        self.item_type = item_type.KEYWORD
        self.items = [] # list of arguments
        self.name = ''
        self.from_new_line = False # start all arguments from the next line?
        self.expanded = s.expanded


# Keyword's argument
class argument(item):

    def __init__(self, s):
        self.item_type = item_type.ARGUMENT
        self.items = [] # list of strings
        self.name = ''
        self.form = '' # QCheckBox, QLineEdit, QComboBox
        self.required = False


# Keyword implementation - a piece of INP-code for CalculiX input file
class implementation(item):

    def __init__(self, s, keyword, INP_code, name=None):
        self.item_type = item_type.IMPLEMENTATION
        self.items = keyword.copyItems() # newer use deepcopy!
        for item in self.items:
            item.parent = self
        self.parent = keyword
        self.active = True
        self.expanded = s.expanded

        # Name of current implementation (of *AMPLITUDE, *STEP, *MATERIAL etc.)
        index = len(self.parent.getImplementations())
        if name:
            self.name = name # it will be used in edit Dialog
        else:
            lead_line = INP_code[0]
            match = re.search('(NAME|ELSET|NSET)\s*=\s*([\w\!\#\%\$\&\"\'\(\)\*\=\+\-\.\/\:\;\<\>\?\@\[\]\^\_\`\{\\\|\}\~]*)', lead_line.upper())
            if match:
                self.name = lead_line[match.start(2):match.end(2)]
            else:
                self.name = self.parent.name[1:] + '-' + str(index + 1)

        self.INP_code = INP_code # INP-code for current implementation - list of strings
        self.parent.items.insert(index, self) # append implementation to keyword's items
        if name:
            logging.info('{} \"{}\" updated.'.format(keyword.name, self.name))
        else:
            logging.info('{} \"{}\" created.'.format(keyword.name, self.name))


# Run test
if __name__ == '__main__':
    # from pycallgraph import PyCallGraph
    # from pycallgraph import Config
    # from pycallgraph import GlobbingFilter
    # from pycallgraph.output import GraphvizOutput

    clean.screen()
    logging.basicConfig(level=0, format='%(message)s')
    start = time.perf_counter() # start time
    p = path.Path()
    s = settings.Settings(p)

    # modules = [m[:-3]+'*' for m in os.listdir(p.src) if m.endswith('.py')] + ['Window*']
    # config = Config()
    # config.trace_filter = GlobbingFilter(
    #     include=modules, exclude=['logging*', '*FileFinder'])
    # graphviz = GraphvizOutput(output_file=__file__[:-3]+'.png')
    # with PyCallGraph(output=graphviz, config=config):

    KOM(p, s)
    print('\nTotal {:.1e} seconds'\
        .format(time.perf_counter()-start)) # spent time
