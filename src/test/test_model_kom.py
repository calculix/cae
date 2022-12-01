#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2023
Distributed under GNU General Public License v3.0

Test for src/model/kom.py
"""

# Standard modules
import os
import sys
import unittest

# My modules
sys_path = os.path.abspath(__file__)
sys_path = os.path.dirname(sys_path)
sys_path = os.path.join(sys_path, '..')
sys_path = os.path.normpath(sys_path)
sys_path = os.path.realpath(sys_path)
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
from model.kom import KWL, ItemType


class TestKOM(unittest.TestCase):

    def test_kom1(self):
        """Test if CalculiX keywords are generated"""
        self.assertTrue(len(KWL.keywords))

    def test_kom2(self, parent=None, path=None):
        """Test parent-child relations in all tree items"""

        if parent is None:
            parent = KWL.root
            path = KWL.root.name
        for item in parent.items:
            if item.itype == ItemType.ARGUMENT:
                continue
            path_downwards = path + ' -> ' + item.name
            path_upwards = ' -> '.join(item.get_path())
            self.assertEqual(path_downwards, path_upwards)
            self.test_kom2(item, path_downwards)


if __name__ == '__main__':
    unittest.main()
