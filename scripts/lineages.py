"""
Filter opendata fasta file for genomes that were used for Pango
lineage designations.
Locate all differences from reference genome.
Partition into segments based on ARTIC v3 primers.

cd data
wget https://github.com/cov-lineages/pango-designation/raw/master/lineages.csv
"""

import csv
import lzma


