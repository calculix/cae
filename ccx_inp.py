# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, UJV Rez, June 2019.
    Distributed under GNU General Public License, version 2.

    INP parser. Imports INP model into CAE.
"""


import ccx_log, ccx_dom


class Parse:

    """
        DOM         -   data object model, generated once per session during start.

        lines    -   is a multiline text = list of strings, could be just contents of the whole
                        INP file or piece of INP code, passed from some keyword's implementation.
    """

    def __init__(self, DOM, textEdit, lines):

        # Configure logging
        self.logger = ccx_log.logger(textEdit)

        # Make lines available for all methods of the class + clear lines
        self.lines = []
        for line in lines:
            line = line.strip() # clear line
            if not(line.startswith('**')): # remove comments
                self.lines.append(line)

        self.parse(DOM.root)


    # Parse each keyword in DOM
    def parse(self, parent):
        for item in parent.items: # for each keyword from DOM
            if (item.item_type == 'group') or (item.item_type == 'keyword'):
                if item.item_type == 'keyword':
                    self.search(item)
                self.parse(item)


    # Search each keyword in lines (INP file or piece of INP code)
    def search(self, item): # item is a keyword object

        for i in range(len(self.lines)):
            # TODO distinguish 'NODE' and 'NODE PRINT'
            if self.lines[i].startswith(item.name): # if starts with keyword

                INP_code = [self.lines[i]] # must be list of strings

                while i+1<len(self.lines) and not self.lines[i+1].startswith('*'): # there will be no comments
                    INP_code.append(self.lines[i+1])
                    i += 1

                # TODO parse item's (keyword's) arguments

                # Create implementation object
                impl = ccx_dom.implementation(item, INP_code)
                # impl.show()
