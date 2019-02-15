#!/usr/bin/env bash

# Generate INDRA statements and BEL from those PubMed identifiers
cat pubmeds.tsv | cut -f 1 | grep -e "^\d\+$" > articles/pubmeds.txt
bel-enrichment from_pmids --pmids articles/pubmeds.txt --output articles/curation.tsv --pickle-file articles/statements.pkl 2> articles/log.txt

mkdir TREM2
bel-enrichment from_agents -a "TREM2" --output TREM2/curation.tsv --pickle-file TREM2/statements.pkl 2> TREM2/log.txt

mkdir MAOB
bel-enrichment from_agents -a "MAOB" --output MAOB/curation.tsv --pickle-file MAOB/statements.pkl 2> MAOB/log.txt

mkdir IRS1
bel-enrichment from_agents -a "IRS1" --output IRS1/curation.tsv --pickle-file IRS1/statements.pkl 2> IRS1/log.txt

mkdir IGF1R
bel-enrichment from_agents -a "IGF1R" --output IGF1R/curation.tsv --pickle-file IGF1R/statements.pkl 2> IGF1R/log.txt
