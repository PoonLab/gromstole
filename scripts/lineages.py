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
from seq_utils import iter_fasta, batch_fasta
import minimap2
import sys

# import tile coordinates file
tiles = []
reader = csv.reader(open('data/tile-coords.csv'))
_ = next(reader)
for row in reader:
    idx, lpos, llen, rpos, rlen = row
    tiles.append(tuple([int(lpos)+int(llen), int(rpos)-1]))

# parse Pango lineage designations file
reader = csv.reader(open("data/lineages.csv"))
lineages = dict([row for row in reader])


# open stream to xz-compressed FASTA
handle = lzma.open("data/sequences.fasta.xz", 'rt')
ref_file = "data/NC_045512.fa"

batcher = batch_fasta(filter(lambda x: x[0] in lineages, iter_fasta(handle)), size=100)
for fasta in batcher:
    mm2 = minimap2.minimap2(fasta, ref=ref_file, stream=True, nthread=1, minlen=29000)
    for qname, diffs, missing in [minimap2.encode_diffs(row) for row in mm2]:
        print(qname, diffs, missing)
        break
    break
