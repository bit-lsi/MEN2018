# -*- coding: utf-8 -*-

"""Testing in python."""

import logging
import unittest
import os
import pandas as pd
import networkx as nx

log = logging.getLogger(__name__)
log.setLevel(10)


from src.pyrcr.rcr_3 import search_node, example_1, example_2, example_3, gene_f, protein_a


class TestRCR(unittest.TestCase):
    """Test RCR Algoritihm."""
    
    def test_rcr_3(self):
        """Test rcr."""
        the_path = os.getcwd()
        vals = the_path+'/data/test.csv'
        result1 = search_node(example_2,gene_f, vals)
        result2 = search_node(example_1,protein_a, vals)
        result3 = search_node(example_3,protein_a, vals)
        
        self.assertEqual(result1[0][1], 'true')
        self.assertEqual(result2[0][1], 'true')
        self.assertEqual(result3[0][1], 'true')