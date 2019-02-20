# -*- coding: utf-8 -*-

"""Testing in pythonx."""

import logging
import unittest
import os
import pandas as pd
import networkx as nx

log = logging.getLogger(__name__)
log.setLevel(10)

from .constants import make_graph_1, gene_c,protein_a,rna_d,protein_b

from pyrcr.rcr_1 import is_correct

class TestRCR(unittest.TestCase):
    """Test RCR Algoritihm."""

    def test_rcr_1(self):
        """Test rcr."""
        
        the_path = os.getcwd()
        graph_1 = make_graph_1()
        vals = pd.read_csv(the_path+'/data/test.csv',',')
        vals.columns = ['val','gene']
        vals = {k:v for k,v in zip(vals['gene'].values,vals['val'].values)}
        
        nodes = {}
        
        for i in graph_1.nodes():
            nodes[i] = vals[i.name]
    
        nx.set_node_attributes(graph_1,nodes,name='data')

        result = is_correct(graph_1, rna_d,protein_b)

        self.assertEqual(result, False)
 