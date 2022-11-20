#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""© Ihor Mirzov, 2019-2023
Distributed under GNU General Public License v3.0

Test if all keywords have a corresponding icon.
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
from model.kom import KOM


class TestIcons(unittest.TestCase):

    def test_icons(self):
        """Test if all keywords have a corresponding icon."""
        keywords = set(['icon_logo', 'icon_job'])
        for item in KOM.keywords + KOM.groups:
            item_name = re.sub(r'\*', '', item.name) # cut star
            icon_name = 'icon_' + re.sub(r'[ -]', '_', item_name)
            keywords.add(icon_name.lower())

        files = glob.glob(p.img + '/icon_*.png')
        files = [os.path.basename(f)[:-4] for f in files]

        for k in keywords:
            if not k in files:
                print('No icon for', k)
        for f in files:
            if not f in keywords:
                print('Redundant', f)

        self.assertEqual(len(files), len(keywords))


if __name__ == '__main__':
    unittest.main()
