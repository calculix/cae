# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, June 2019.
    Distributed under GNU General Public License, version 2.

    INP parser. Imports INP model into CAE.
"""


import ccx_log, ccx_dom


class Parse:

    """
        DOM     -   data object model, generated once per session during start.

        INP_doc -   multiline text = list of strings, could be contents of the whole INP
                    file or piece of INP code, passed from some keyword's implementation.
    """

    def __init__(self, DOM, textEdit, INP_doc):

        # Configure logging
        self.logger = ccx_log.logger(textEdit)

        # Make INP_doc available for all methods of the class + clear each line
        self.INP_doc = INP_doc
        # self.INP_doc = []
        # for line in INP_doc:
        #     line = line.strip() # clear line
        #     if not(line.startswith('**')): # remove comments
        #         self.INP_doc.append(line)

        # Start recursive DOM iterator
        self.iterator(DOM.root)


    # Recursively iterate over DOM items, call parser for each keyword
    def iterator(self, parent):
        for item in parent.items: # for each group/keyword from DOM
            if (item.item_type == 'group'):
                self.iterator(item) # continue call iterator until dig to keyword
            if item.item_type == 'keyword':
                self.parser(item) # parse keyword
                self.iterator(item) # and dig to child keywords


    # Search CalculiX'es keywords in the INP_doc lines
    def parser(self, item): # item is a keyword object

        for i in range(len(self.INP_doc)):
            line = self.INP_doc[i].rstrip() # cut '\n'

            # Skip comments
            if line.lstrip().startswith('**'):
                continue

            # Distinguish 'NODE' and 'NODE PRINT'
            if ',' in line:
                keyword_name = line.split(',')[0]
            else:
                keyword_name = line

            # If we found a keyword in the INP_doc
            if keyword_name.strip() == item.name:
                INP_code = [line] # must be list of strings

                while i+1<len(self.INP_doc)\
                    and not self.INP_doc[i+1].lstrip().startswith('*'): # there will be no comments

                    INP_code.append(self.INP_doc[i+1].rstrip()) # cut '\n'
                    i += 1

                # Create implementation object
                impl = ccx_dom.implementation(item, INP_code)
