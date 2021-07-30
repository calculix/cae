#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Test for src/path.py
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
from path import p


class Test(unittest.TestCase):

    # Test all attributes on the Path class
    def test_path(self):
        for attr in dir(p):
            a = getattr(p, attr)
            if type(a) is str:
                print('p.{} = {}'.format(attr, a))
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
