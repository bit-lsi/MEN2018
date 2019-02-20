# -*- coding: utf-8 -*-

"""Testing in pythonx."""

import logging
import unittest

log = logging.getLogger(__name__)
log.setLevel(10)

from .constants import make_graph_1, gene_c

from pyrcr.rcr_1 import is_correct

class TestRCR(unittest.TestCase):
    """Test RCR Algoritihm."""

    def test_rcr_1(self):
        """Test rcr."""

        graph_1 = make_graph_1()

        from .constants import my_function

        result = my_function(graph_1, gene_c, optional=data)

        self.assertRaises(count_dict[False], 3)
        self.assertEqual(count_dict[True], 2)
