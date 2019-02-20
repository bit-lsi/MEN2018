# -*- coding: utf-8 -*-

"""Testing in python."""

import logging
import unittest
import os
import pandas as pd
import networkx as nx

log = logging.getLogger(__name__)
log.setLevel(10)


#from constants import make_graph_2, protein_e, protein_b, gene_f


#from pyrcr.rcr_3 import search_node

class TestRCR(unittest.TestCase):
    """Test RCR Algoritihm."""
    
    def test_rcr_3(self):
        """Test rcr."""
        the_path = os.getcwd()
        graph_1 = make_graph_1()
        #vals = the_path+'/data/test.csv'
        result1 = search_node(example_2,gene_f, 'test.csv')
        result2 = search_node(example_1,protein_a, 'test.csv')
        result3 = search_node(example_3,protein_a, 'test.csv')
        
        self.assertEqual(result1[0][1], 'true')
        self.assertEqual(result2[0][1], 'true')
        self.assertEqual(result3[0][1], 'true')
        print(result1[0][1])
        print(result2[0][1])
        print(result3[0][1])

a = TestRCR()
print(a.test_rcr_3())
