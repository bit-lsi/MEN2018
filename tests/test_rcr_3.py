# -*- coding: utf-8 -*-

"""Testing in python."""

import logging
import unittest
import os
import pandas as pd
import networkx as nx

log = logging.getLogger(__name__)
log.setLevel(10)

#from .constants import make_graph_1, gene_c,protein_a,rna_d,protein_b
#from .constants import make_graph_2, protein_e, protein_b, gene_f

#from pyrcr.rcr_1 import is_correct
#from pyrcr.rcr_3 import search_node

class TestRCR(unittest.TestCase):
    """Test RCR Algoritihm."""
    
    def test_rcr_3(self):
        """Test rcr."""
        
        result1 = search_node(example_2,gene_f, 'test.csv')
        self.assertEqual(result1[0][1], 'true')
        return result1[0][1]

a = TestRCR()
print(a.test_rcr_3())
