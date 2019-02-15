Curation
========
1. Generate a PubMed identifier list:

.. code-block:: sh

   $ cat pubmeds.tsv | cut -f 1 | grep -e "^\d\+$" > pubmeds.txt
   
   
2. Generate curation sheet and pickles:

.. code-block:: sh

   $ bel-enrichment from-pmids --pmids pubmeds.txt --output curation.tsv --pickle-file statements.pkl
