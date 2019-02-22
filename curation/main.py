# -*- coding: utf-8 -*-

"""Run this script to generate all summaries.

Run with ``python main.py compile -r``.
"""

import os

from bel_enrichment import BELSheetsRepository

HERE = os.path.abspath(os.path.dirname(__file__))
repo = BELSheetsRepository(directory=HERE, sheet_suffix='.xls')
main = repo.build_cli()

if __name__ == '__main__':
    main()
