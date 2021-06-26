#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

CalculiX Keyword Object Model (hierarchy).
Parser supposts mutually exclusive parameters for keywords.
Keywords with all arguments are parsed from config/kom.xml.
On INP import KOM is enriched with keyword implementations.
"""

# Standard modules
import re
import os
import sys
import logging
import copy
from enum import Enum
import xml.etree.ElementTree as ET
import traceback

# My modules
sys_path = os.path.abspath(__file__)
sys_path = os.path.dirname(sys_path)
sys_path = os.path.join(sys_path, '..')
sys_path = os.path.normpath(sys_path)
sys_path = os.path.realpath(sys_path)
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
import path
import settings
import tests


# Keyword Object Model
class KOM:

    # Read CalculiX keywords hierarchy
    def __init__(self, kom_xml=path.p.kom_xml):

        # List of all existing keywords
        self.keywords = []
        self.keyword_names = ()

        # Group 'Model' from kom.xml
        self.root = Group()

        # Parse keywords hierarchy and build KOM
        try:
            t = ET.parse(kom_xml)
            self.build(t.getroot(), self.root)
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
            item = klass()
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

    # Get keyword from self.keywords by its name
    def get_keyword_by_name(self, name, parent=None):
        for kw in self.keywords:
            if kw.name == name:
                return kw

    # Recursively get whole model inp_code as list of strings (lines)
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

    def get_top_keyword_by_name(self, parent, kw_name):
        kw = None
        for item in parent.items:
            if item.itype == ItemType.ARGUMENT:
                continue
            if item.itype == ItemType.KEYWORD \
                and item.name == kw_name.upper():
                return item
            else:
                kw = self.get_top_keyword_by_name(item, kw_name)
                if kw is not None:
                    return kw

    # Test parent-child relations in all tree items
    # Print item paths top-downwards and bottom-upwards
    def test(self, parent=None, path=None):
        if parent is None:
            parent = self.root
            path = self.root.name
        for item in parent.items:
            if item.itype == ItemType.ARGUMENT:
                continue
            path_downwards = path + ' -> ' + item.name
            path_upwards = ' -> '.join(item.get_path())
            if path_downwards != path_upwards:
                all_the_same = False
                msg = 'KOM is built with mistakes!\n{}\n{}'\
                    .format(path_downwards, path_upwards)
                logging.error(msg)
                return
            self.test(item, path_downwards)


# Enums for 'itype' variable
class ItemType(Enum):
    GROUP = 0
    KEYWORD = 1
    ARGUMENT = 2
    IMPLEMENTATION = 3


# Needed for inheritance by further classes
class Item:

    def __init__(self):
        self.itype = ''      # item type: group/keyword/argument/implementation
        self.name = ''       # name of item, string
        self.items = []      # list of children
        self.parent = None   # item parent item
        self.active = False

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

    # Recursively copy tree branch to implementation
    def copy_items_to(self, another_item):
        another_item.items = []
        for item in self.items:
            if item.itype in [ItemType.GROUP, ItemType.KEYWORD]:
                copied_item = copy.copy(item) # newer use deepcopy!
                copied_item.parent = another_item
                copied_item.active = True
                another_item.items.append(copied_item)
                item.copy_items_to(copied_item)
            if item.itype == ItemType.ARGUMENT:
                copied_item = copy.copy(item)
                another_item.items.append(copied_item)

    # Recursive function to count keyword implementations in item descendants
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

    # Search item by name among children
    def get_child_by_name(self, name):
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

    # Get item path bottom-upwards
    def get_path(self):
        parents = []
        while self.parent is not None: # root has None parent
            parents.insert(0, self)
            self = self.parent
        parents.insert(0, self) # insert root
        return [p.name for p in parents]


# Group of keywords, like 'Properties', 'Constraints', etc.
class Group(Item):

    def __init__(self):
        super().__init__()
        self.itype = ItemType.GROUP
        self.items = [] # list of groups and keywords
        self.name = 'Model' # default name (root group)
        self.expanded = settings.s.expanded


# *AMPLITUDE, *BOUNDARY, *STEP etc.
class Keyword(Item):

    def __init__(self):
        super().__init__()
        self.itype = ItemType.KEYWORD
        self.items = [] # list of arguments
        self.name = ''
        self.from_new_line = False # start all arguments from the next line?
        self.expanded = settings.s.expanded


# Keyword argument
class Argument(Item):

    def __init__(self):
        super().__init__()
        self.itype = ItemType.ARGUMENT
        self.items = [] # list of strings
        self.name = ''
        self.form = '' # QCheckBox, QLineEdit, QComboBox
        self.required = False


# Keyword implementation - a piece of INP-code for CalculiX input file
class Implementation(Item):

    def __init__(self, keyword, inp_code, name=None):
        super().__init__()

        self.itype = ItemType.IMPLEMENTATION
        keyword.copy_items_to(self)
        self.parent = keyword
        self.active = True
        self.expanded = settings.s.expanded

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
        self.parent.items.insert(index, self) # append implementation to keyword items
        if name:
            logging.info('{} \"{}\" updated.'.format(keyword.name, self.name))
        else:
            logging.info('{} {} created.'.format(keyword.name, self.name))
            # logging.debug(' > '.join(self.get_path()))


# Test all CalculiX keywords 
@tests.test_wrapper()
def test():
    # from pycallgraph import PyCallGraph
    # from pycallgraph import Config
    # from pycallgraph import GlobbingFilter
    # from pycallgraph.output import GraphvizOutput

    # modules = [m[:-3]+'*' for m in os.listdir(p.src) if m.endswith('.py')] + ['Window*']
    # config = Config()
    # config.trace_filter = GlobbingFilter(
    #     include=modules, exclude=['logging*', '*FileFinder'])
    # graphviz = GraphvizOutput(output_file=__file__[:-3]+'.png')
    # with PyCallGraph(output=graphviz, config=config):

    os.chdir(os.path.dirname(__file__))
    k = KOM()

    # Print all CalculiX keywords
    for kw in k.keyword_names:
        print(kw)

    # Test KOM
    # k.test()


# Run test
if __name__ == '__main__':
    test()
