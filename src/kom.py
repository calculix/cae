# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, August 2019
    Distributed under GNU General Public License v3.0

    CalculiX Keyword Object Model.
    Keywords with all arguments are parsed from kom.xml.
    Parser supposts mutually exclusive parameters for keywords.
"""


import re, logging, copy, os, sys, time
from enum import Enum
import xml.etree.ElementTree as ET
from settings import Settings


# Keyword Object Model
class KOM:


    # Read CalculiX keywords hierarchy
    def __init__(self):

        try:
            self.root = group() # group 'Model' from kom.xml

            # List of all existing keywords
            self.keywords = []

            # Analyze keywords hierarchy
            kom_xml = os.path.join(os.path.dirname(sys.argv[0]), 'kom.xml')
            tree = ET.parse(kom_xml)
            self.buildKOM(tree.getroot(), self.root)

            # All possible keywords nesting variants - needed for parsing INP_doc
            self.pathes = []
            self.buildPathes(self.root)
            self.pathes.sort(key=self.keyword_counter, reverse=True) # maximum nesting first
            # for path in self.pathes:
            #     logging.debug(str([item.name for item in path]))

            # # Regenerate all HTML help pages
            # for item in self.keywords:
            #     import dialog
            #     dialog.saveHTML(item)

            logging.info('CalculiX object model generated.')
        except:
            logging.error('Can\'t generate keywords hierarchy!')


    # Recursively build Keyword Object Model
    def buildKOM(self, xml_branch, parent):
        for xml_child in xml_branch:
            # logging.debug(xml_child.tag + ', ' + str(xml_child.attrib))
            """
                xml_child.tag, xml_child.attrib:

                keyword, {'name': '*RESTART'}
                keyword, {'name': '*TIME POINTS'}
                group, {'name': 'Interactions'}
                group, {'name': 'Constraints'}
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
            if xml_child.tag=='keyword':
                self.keywords.append(item)

            # Argument values for QComboBox
            if xml_child.tag=='argument' and xml_child.text:
                values = xml_child.text.strip()
                if len(values):
                    values = values.split('|')
                    item.items.extend(values)
                    # print(parent.name, ":", item.name, values)

            self.buildKOM(xml_child, item)


    # Recursively builds all possible pathes to nested keywords in KOM
    def buildPathes(self, parent, path=None):
        if not path:
            path = [] # list of groups, keywords and implementations
        for item in parent.items:
            if item.item_type != item_type.ARGUMENT:
                self.buildPathes(item, path + [item])
        if len(path):
            if path not in self.pathes:
                self.pathes.append(path)


    # Get nesting path for each of the parsed keyword
    def getPath(self, keyword_chain):

        # Duplicated words in the end of keyword_chain play no role
        if len(keyword_chain) > 1 and \
            keyword_chain[-1] == keyword_chain[-2]:
            del keyword_chain[-2]
        logging.debug('')
        logging.debug('keyword_chain: ' + str(keyword_chain))

        # Now compare keyword_chain with all self.pathes
        for path in self.pathes:

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
                # logging.debug(str([item.name for item in path]) + '\n')
                del keyword_chain[:minimum_j]
                return path


    # Count keywords in path
    def keyword_counter(self, path):
        keyword_counter = 0
        for item in path:
            if item.name.startswith('*'):
                keyword_counter += 1
        return keyword_counter


    # Get keyword from self.keywords by its name
    def getKeywordByName(self, name, parent=None):
        for kw in self.keywords:
            if kw.name == name:
                return kw


# Enums for 'item_type' variable
class item_type(Enum):
    GROUP = 0
    KEYWORD = 1
    ARGUMENT = 2
    IMPLEMENTATION = 3


# Needed for inheritance by further classes
class item:

    # Read application's global settings
    settings = Settings()

    item_type = ''          # item's type: group/keyword/argument/implementation
    name = ''               # name of item, string
    items = []              # list of children
    parent = None           # item's parent item
    expanded = settings.expanded
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


    # Print all the branch KOM elements starting from current item
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

    def __init__(self):
        self.item_type = item_type.GROUP
        self.items = [] # list of groups and keywords
        self.name = 'Model' # default name (root group)


# *AMPLITUDE, *BOUNDARY, *STEP etc.
class keyword(item):

    def __init__(self):
        self.item_type = item_type.KEYWORD
        self.items = [] # list of arguments
        self.name = ''
        self.from_new_line = False # start all arguments from the next line?


# Keyword's argument
class argument(item):

    def __init__(self):
        self.item_type = item_type.ARGUMENT
        self.items = [] # list of strings
        self.name = ''
        self.form = '' # QCheckBox, QLineEdit, QComboBox
        self.required = False


# Keyword implementation - a piece of INP-code for CalculiX input file
class implementation(item):

    def __init__(self, keyword, INP_code, name=None):
        self.item_type = item_type.IMPLEMENTATION
        self.items = keyword.copyItems() # newer use deepcopy!
        for item in self.items:
            item.parent = self
        self.parent = keyword
        self.active = True

        # Name of current implementation (of *AMPLITUDE, *STEP, *MATERIAL etc.)
        index = len(self.parent.getImplementations())
        if name:
            self.name = name # it will be used in edit Dialog
        else:
            lead_line = INP_code[0]
            match = re.search('(NAME|ELSET|NSET)\s*=\s*([\w\-]*)', lead_line.upper())
            if match:
                self.name = lead_line[match.start(2):match.end(2)]
            else:
                self.name = self.parent.name[1:] + '-' + str(index + 1)

        self.INP_code = INP_code # INP-code for current implementation - list of strings
        self.parent.items.insert(index, self) # append implementation to keyword's items
        logging.info(self.name + ' created.')


# Test module
if __name__ == '__main__':
    start = time.perf_counter() # start time
    logging.getLogger().setLevel(logging.DEBUG)
    KOM()
    print('\nTotal {:.1e} seconds'.format(time.perf_counter()-start)) # end time
