#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2023
Distributed under GNU General Public License v3.0

Test if all keywords has HTML file in the documentation
and if there are no redundant HTML files.
"""

# Standard modules
import os
import re
import sys
import glob
import unittest

# My modules
sys_path = os.path.abspath(__file__)
sys_path = os.path.dirname(sys_path)
sys_path = os.path.join(sys_path, '..')
sys_path = os.path.normpath(sys_path)
sys_path = os.path.realpath(sys_path)
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
from path import p
from model.kom import KWL


class TestDoc(unittest.TestCase):

    def test_doc(self):
        """Amount of HTML files has to coincide with the amount of keywords."""
        keywords = set()
        for item in KWL.keywords:
            keyword_name = item.name[1:] # cut star
            html_page_name = re.sub(r'[ -]', '_', keyword_name)
            keywords.add(html_page_name)

        files = glob.glob(p.doc + '/*.html')
        files = [os.path.basename(f)[:-5] for f in files]

        for k in keywords:
            if not k in files:
                print('No documentation file for', k)
        for f in files:
            if not f in keywords:
                print('Redundant doc file', f)
        
        self.assertEqual(len(files), len(keywords))


if __name__ == '__main__':
    unittest.main()
