#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" © Ihor Mirzov, July 2020
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


# Keyword Object Model
class KOM:


    # Read CalculiX keywords hierarchy
    def __init__(self, p=None, s=None, kom_xml=None):
        self.s = s # could be None

        # List of all existing keywords
        self.keywords = []
        self.keyword_names = ()

        # All possible keywords nesting variants - needed for parsing inp_doc
        self.paths = []

        # Group 'Model' from kom.xml
        self.root = Group(self.s)

        try:
            if p is not None:
                kom_xml = p.kom_xml
            if kom_xml is None:
                raise Exception

            # Parse keywords hierarchy and build KOM
            t = ET.parse(kom_xml)
            self.build(t.getroot(), self.root)

            self.build_paths(self.root)
            self.paths.sort(key=self.keyword_counter, reverse=True) # maximum nesting first
            # for path in self.paths:
            #     logging.debug(str([item.name for item in path]))

            logging.info('Keywords object model generated.')
        except:
            msg = 'Can\'t generate keywords object model!\n' \
                + traceback.format_exc()
            logging.error(msg)


    # Recursively build Keyword Object Model
    def build(self, xml_branch, parent):
        for xml_child in xml_branch:
            logging.log(5, xml_child.tag + ', ' + str(xml_child.attrib))
            """
            xml_child.tag, xml_child.attrib:

            Keyword, {'name': '*RESTART'}
            Keyword, {'name': '*TIME POINTS'}
            Group, {'name': 'Interactions'}
            Group, {'name': 'Constraints'}
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
            if xml_child.tag=='Keyword':
                self.keywords.append(item)
                self.keyword_names += (item.name, )

            # Argument values for QComboBox
            if xml_child.tag=='Argument' and xml_child.text:
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
            if item.itype != ItemType.ARGUMENT:
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


    # Recursively get whole model's inp_code as list of strings (lines)
    # Parent is KOM item, level defines code folding/padding
    def get_inp_code_as_lines(self, parent=None, level=0):
        lines = []
        if not parent:
            parent = self.root
        if parent.itype == ItemType.IMPLEMENTATION:
            level += 1

        # For each group/keyword from KOM
        for item in parent.items:
            if item.itype == ItemType.ARGUMENT:
                continue
            if item.itype == ItemType.IMPLEMENTATION:

                # inp_code is stripped
                for line in item.inp_code:
                    if line.upper().startswith(self.keyword_names):
                        padding = ' '*4*level
                    else:
                        padding = ' '*4*(level+1)
                    if line.startswith('**'):
                        padding = ''
                    lines.append(padding + line + '\n')

            # Continue call iterator until dig to implementation
            lines.extend(self.get_inp_code_as_lines(item, level))

        return lines


# Enums for 'itype' variable
class ItemType(Enum):
    GROUP = 0
    KEYWORD = 1
    ARGUMENT = 2
    IMPLEMENTATION = 3


# Needed for inheritance by further classes
class Item:
    itype = ''          # item's type: group/keyword/argument/implementation
    name = ''               # name of item, string
    items = []              # list of children
    parent = None           # item's parent item
    active = False

    # Define if item is active
    def is_active(self):

        # Top groups
        if self.parent and self.parent.name == 'Model':
            self.active = True
            return True

        # Children of active groups
        elif self.parent and self.parent.active and self.parent.itype in \
            [ItemType.GROUP, ItemType.IMPLEMENTATION]:
            self.active = True
            return True

        # Implementations are always active
        elif self.itype == ItemType.IMPLEMENTATION:
            self.active = True
            return True

        else:
            self.active = False
            return False

    # Recursive function to count keyword implementations in item's descendants
    def count_implementations(self):
        if self.itype == ItemType.ARGUMENT:
            return 0

        counter = 0
        for i in self.items:
            if i.itype == ItemType.IMPLEMENTATION:
                counter += 1
            else:
                counter += i.count_implementations()
        return counter

    # Get list of keywords implementations
    def get_implementations(self):
        imps = []
        for item in self.items:
            if item.itype == ItemType.IMPLEMENTATION:
                imps.append(item)
        return imps

    # Get list of items - non implementations
    def copy_items(self):
        items = []
        for item in self.items:
            item.active = True
            if item.itype in [ItemType.GROUP, ItemType.KEYWORD]:
                item = copy.copy(item)
                item.items = item.copy_items()
                items.append(item)
            if item.itype == ItemType.ARGUMENT:
                item = copy.copy(item)
                items.append(item)
        return items

    # Search item by name among children
    def get_item_by_name(self, name):
        for item in self.items:
            if item.name == name:
                return item

    # Returns first active preceding/parent keyword (not group)
    def get_parent_keyword_name(self):
        if (self.parent.itype == ItemType.KEYWORD or \
            self.parent.name == 'Model') and self.parent.active:
                return self.parent.name
        else:
            return self.parent.get_parent_keyword_name()

    # # Returns preceding/parent group (not keyword)
    # def get_parent_group_name(self):
    #     if (self.parent.itype == ItemType.GROUP or \
    #         self.parent.name == 'Model'):
    #             return self.parent.name
    #     else:
    #         return self.parent.get_parent_group_name()


# Group of keywords, like 'Properties', 'Constraints', etc.
class Group(Item):

    def __init__(self, s):
        self.itype = ItemType.GROUP
        self.items = [] # list of groups and keywords
        self.name = 'Model' # default name (root group)
        if s is not None:
            self.expanded = s.expanded
        else:
            self.expanded = True


# *AMPLITUDE, *BOUNDARY, *STEP etc.
class Keyword(Item):

    def __init__(self, s):
        self.itype = ItemType.KEYWORD
        self.items = [] # list of arguments
        self.name = ''
        self.from_new_line = False # start all arguments from the next line?
        if s is not None:
            self.expanded = s.expanded
        else:
            self.expanded = True


# Keyword's argument
class Argument(Item):

    def __init__(self, s):
        self.itype = ItemType.ARGUMENT
        self.items = [] # list of strings
        self.name = ''
        self.form = '' # QCheckBox, QLineEdit, QComboBox
        self.required = False


# Keyword implementation - a piece of INP-code for CalculiX input file
class Implementation(Item):

    def __init__(self, s, keyword, inp_code, name=None):
        self.itype = ItemType.IMPLEMENTATION
        self.items = keyword.copy_items() # newer use deepcopy!
        for item in self.items:
            item.parent = self
        self.parent = keyword
        self.active = True
        if s is not None:
            self.expanded = s.expanded
        else:
            self.expanded = True

        # Name of current implementation (of *AMPLITUDE, *STEP, *MATERIAL etc.)
        index = len(self.parent.get_implementations())
        if name:
            self.name = name # it will be used in edit Dialog
        else:
            i = 0
            while inp_code[i].startswith('**'):
                i += 1
            lead_line = inp_code[i]

            match = re.search('(NAME|ELSET|NSET)\s*=\s*([\w\!\#\%\$\&\"\'\(\)\*\=\+\-\.\/\:\;\<\>\?\@\[\]\^\_\`\{\\\|\}\~]*)', lead_line.upper())
            if match:
                self.name = lead_line[match.start(2):match.end(2)]
            else:
                self.name = self.parent.name[1:] + '-' + str(index + 1)

        self.inp_code = inp_code # INP-code for current implementation - list of strings
        self.parent.items.insert(index, self) # append implementation to keyword's items
        if name:
            logging.info('{} \"{}\" updated.'.format(keyword.name, self.name))
        else:
            logging.info('{} {} created.'.format(keyword.name, self.name))


# Run test
if __name__ == '__main__':
    # from pycallgraph import PyCallGraph
    # from pycallgraph import Config
    # from pycallgraph import GlobbingFilter
    # from pycallgraph.output import GraphvizOutput

    # Clean screen
    os.system('cls' if os.name=='nt' else 'clear')

    logging.basicConfig(level=0, format='%(message)s')
    start = time.perf_counter() # start time

    # modules = [m[:-3]+'*' for m in os.listdir(p.src) if m.endswith('.py')] + ['Window*']
    # config = Config()
    # config.trace_filter = GlobbingFilter(
    #     include=modules, exclude=['logging*', '*FileFinder'])
    # graphviz = GraphvizOutput(output_file=__file__[:-3]+'.png')
    # with PyCallGraph(output=graphviz, config=config):

    os.chdir(os.path.dirname(__file__))
    kom_xml = '../../config/kom.xml'

    # KOM(None, None)
    k = KOM(None, None, kom_xml)
    print('\nTotal {:.1e} seconds'\
        .format(time.perf_counter()-start)) # spent time

    print(k.keyword_names)
