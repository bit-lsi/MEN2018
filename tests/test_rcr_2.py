# -*- coding: utf-8 -*-

"""Testing in python."""

import logging
import unittest
import os

from constants import make_graph_3, protein_a, protein_e, protein_b, rna_d, gene_c, protein_g

from pyrcr.rcr_2 import path_validation

log = logging.getLogger(__name__)
log.setLevel(10)

HERE = os.path.realpath('__file__')
TEST_DIRECTORY = os.path.abspath(os.path.join(HERE, os.pardir))
PROJECT_DIRECTORY = os.path.abspath(os.path.join(TEST_DIRECTORY, os.pardir))

TESTFile = os.path.join(PROJECT_DIRECTORY, 'data','test.csv')

class TestSpia(unittest.TestCase):
    
    def test_rcr_2(self):
        graph = make_graph_3()
        assert path_validation(graph,TESTFile,protein_a,gene_c) == False
        assert path_validation(graph,TESTFile,protein_a,protein_b) == True
        assert path_validation(graph,TESTFile,protein_a,protein_a) == "No path between same nodes"
        assert path_validation(graph,TESTFile,protein_a,rna_d) == "path doesn't exist"
        assert path_validation(graph,TESTFile,gene_c,protein_g) == "path cannot be evaluated"
