#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2023
Distributed under GNU General Public License v3.0

Test for src/doc.py
"""

# Standard modules
import os
import sys
import re
import unittest

# My modules
sys_path = os.path.abspath(__file__)
sys_path = os.path.dirname(sys_path)
sys_path = os.path.join(sys_path, '..')
sys_path = os.path.normpath(sys_path)
sys_path = os.path.realpath(sys_path)
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
from model.kom import KWL
from path import p


class TestUtilsDoc(unittest.TestCase):

    def test_utils_doc(self):
        """Check if HTML pages are generated for all keywords"""
        keywords = [re.sub(r'[ -]', '_', kw.name[1:]) for kw in KWL.keywords]
        keywords = sorted(set(keywords))
        # print(keywords)
        # print('\nTotal {} keywords'.format(len(keywords)))
        pages = [fn for fn in os.listdir(p.doc) if fn.endswith('.html')]
        pages = sorted(pages)
        # print('Total {} HTML pages'.format(len(pages)))
        # print(pages)
        if len(keywords) > len(pages):
            for page in pages:
                if page[:-5] in keywords:
                    keywords.remove(page[:-5])
            print('\nThose keywords have no HTML pages:')
            print(keywords)
        self.assertTrue(len(keywords) == len(pages))


if __name__ == '__main__':
    unittest.main()
