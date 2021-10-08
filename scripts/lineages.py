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
from seq_utils import iter_fasta

# import tile coordinates file
reader = csv.DictReader(open('data/tile-coords.csv'))
for row in reader:



# parse Pango lineage designations file
reader = csv.reader(open("data/lineages.csv"))
lineages = dict([row for row in reader])

# open stream to xz-compressed FASTA
handle = lzma.open("data/sequences.fasta.xz", 'rt')
for h, s in iter_fasta(handle):
    print(h)
    print(s[:100])
    break
