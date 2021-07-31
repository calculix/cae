#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2021
Distributed under GNU General Public License v3.0

Test for src/model/parsers/mesh.py
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
from path import p
from model.parsers.mesh import Mesh


class TestKOM(unittest.TestCase):

    def test_mesh(self):
        """Test if parsed mesh has nodes and elements"""
        inp = os.path.join(p.examples, 'default.inp')
        m = Mesh(inp)
        self.assertTrue(len(m.nodes))
        self.assertTrue(len(m.elements))


if __name__ == '__main__':
    unittest.main()
