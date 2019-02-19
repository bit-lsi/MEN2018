# -*- coding: utf-8 -*-

"""Testing in python."""

import logging
import unittest

log = logging.getLogger(__name__)
log.setLevel(10)


class TestSpia(unittest.TestCase):
    """Test RCR Algoritihm."""

    def test_rcr_1(self):
        """Test rcr."""
        pass

        self.assertEqual(1, 1)
