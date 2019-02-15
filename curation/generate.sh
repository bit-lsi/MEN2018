#!/usr/bin/env bash

# Slice out just the identifiers column and skip the missing entries
cat pubmeds.tsv | cut -f 1 | grep -e "^\d\+$" > pubmeds.txt

# Generate INDRA statements and BEL from those PubMed identifiers
bel-enrichment from_pmids --pmids pubmeds.txt --output curation.tsv --pickle-file statements.pkl 2> generation-log.txt
