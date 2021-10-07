from minimap2 import minimap2
from seq_utils import convert_fasta, revcomp
import csv

# assume we are running from project root
with open('data/NC_045512.fa') as handle:
    ref = convert_fasta(handle)[0][1]

reader = csv.DictReader(open('data/nCoV-2019.tsv'), delimiter='\t')

# convert to FASTA file
coords = {}
for row in reader:
    name = row['name']
    seq = row['seq']
    right = False
    if 'RIGHT' in name:
        right = True
        seq = revcomp(seq)  # reverse-complement 3' primer

    rpos = ref.index(seq)  # exact match works for reference genome

    pid = name.split('_')[1]
    if pid not in coords:
        coords.update({pid: {
            'left': {'rpos': None, 'len': 0},
            'right': {'rpos': None, 'len': 0}}
        })

    coords[pid]['right' if right else 'left']['rpos'] = rpos
    coords[pid]['right' if right else 'left']['len'] = len(seq)

outfile = open('data/tile-coords.csv', 'w')
outfile.write("primer.index,left.rpos,left.len,right.rpox,right.len\n")
for pid, data in coords.items():
    outfile.write("{},{},{},{},{}\n".format(
        pid, data['left']['rpos'], data['left']['len'],
        data['right']['rpos'], data['right']['len']
    ))

outfile.close()
