#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Test for src/model/kom.py
Ensure that all app's folders and files are present.
"""

# Standard modules
import os
import sys
import unittest
import logging

# My modules
sys_path = os.path.abspath(__file__)
sys_path = os.path.dirname(sys_path)
sys_path = os.path.join(sys_path, '..')
sys_path = os.path.normpath(sys_path)
sys_path = os.path.realpath(sys_path)
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
from model.kom import KOM, ItemType


class Test(unittest.TestCase):

    def setUp(self):
        print()


class TestKOM(Test):

    def test_kom1(self):
        """Print all CalculiX keywords"""
        for kw in KOM.keyword_names:
            print(kw)
        self.assertTrue(len(KOM.keyword_names))

    def test_kom2(self):
        """Test parent-child relations in all tree items"""

        def test(parent=None, path=None):
            """Print item paths top-downwards and bottom-upwards"""
            if parent is None:
                parent = KOM.root
                path = KOM.root.name
            for item in parent.items:
                if item.itype == ItemType.ARGUMENT:
                    continue
                path_downwards = path + ' -> ' + item.name
                path_upwards = ' -> '.join(item.get_path())
                self.assertEqual(path_downwards, path_upwards)
                if path_downwards != path_upwards:
                    msg = 'Keyword object model is built with mistakes!\n{}\n{}'\
                        .format(path_downwards, path_upwards)
                    logging.error(msg)
                    return
                test(item, path_downwards)
        test()


if __name__ == '__main__':
    unittest.main()
