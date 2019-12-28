# -*- coding: utf-8 -*-


"""
    © Ihor Mirzov, October 2019
    Distributed under GNU General Public License v3.0

    Parses nodes of the finite element mesh from the CalculiX .inp-file.
"""


import re, logging


class Nodes:


    # Initialization
    def __init__(self, INP_code):
        self.nodes = {} # all mesh nodes with coordinates

        # Parse nodes
        try:
            self.parse_nodes(INP_code)

            msg_text = '{} nodes'.format(len(self.nodes))
            # msg_text += ': ' + str(list(self.nodes.keys()))
            logging.info(msg_text)
        except:
            logging.error('Can\'t parse nodes')


    # Parse nodes with coordinates - *NODE keyword
    def parse_nodes(self, lines):
        for i in range(len(lines)):

            # Distinguish 'NODE' and 'NODE PRINT'
            if ',' in lines[i]:
                keyword_name = lines[i].split(',')[0]
            else:
                keyword_name = lines[i]

            if keyword_name.upper() == '*NODE':
                nodes = []
                lead_line = lines[i]
                match = re.search('NSET\s*=\s*([\w\-]*)', lead_line.upper())

                while i+1<len(lines) and not lines[i+1].startswith('*'): # read the whole block and return
                    a = lines[i+1].replace(',', ' ').split() # to avoid redundant commas in the end of line
                    num = int(a[0]) # node number
                    coords = [float(coord) for coord in a[1:]] # node coordinates
                    if len(coords) == 2: # in 2D case add Z coord equal to zero
                        coords.append(0)

                    # Create NODE
                    node = NODE(num, coords)

                    # Check duplicates
                    if num in self.nodes:
                        msg_text = 'Duplicated node {}.'.format(num)
                        logging.warning(msg_text)
                        del node
                    else:
                        if len(coords):
                            nodes.append(node)
                            self.nodes[num] = node

                            # Bounding box
                            for j,coord in enumerate(coords):
                                if coord < self.bounds[j*2]:
                                    self.bounds[j*2] = coord # update min coords values
                                if coord > self.bounds[j*2+1]:
                                    self.bounds[j*2+1] = coord # update max coords values
                        else:
                            msg_text = 'Node {} has no coordinates and will be removed.'.format(num)
                            logging.warning(msg_text)
                            del node

                    i += 1

                # If all nodes are named as a set
                if match:
                    name = lead_line[match.start(1):match.end(1)]
                    self.create_or_extend_nset(name, nodes)

                # do not return to parse few *NODE sections


class NODE:
    """
        1: [ 0.0, -1742.5, 0.0],
        2: [74.8, -1663.7, 0.0],
        ...
    """
    def __init__(self, num, coords):
        self.num = num
        self.coords = coords
