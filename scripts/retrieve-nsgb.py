"""
Stream genome and metadata from Nextstrain-curated Genbank feeds
Export as xz-compressed JSON with structure:
  {"B.1": {"mutations": {"~|210|T": 54, ...}, "count": 6441}, ...
"""

from urllib import request
import lzma
import gzip
from seq_utils import iter_fasta
import csv
import subprocess
import re
import argparse
from minimap2 import encode_diffs


def parse_metadata(metaurl):
    """ extract PANGO lineage assignments from metadata """
    lineages = {}  # header : PANGO lineage
    handle = gzip.open(request.urlopen(metaurl), 'rt')
    count = 0
    for row in csv.DictReader(handle, delimiter='\t'):
        if row["strain"] == '?':
            continue
        lineages.update({row['strain']: row["pango_lineage"]})
    return lineages


def batcher(handle, size=100):
    """ Break FASTA data into batches to stream to minimap2 """
    stdin = ''
    for i, record in iter_fasta(handle):
        stdin += '>{0}\n{1}\n'.format(*record)
        if i > 0 and i % size == 0:
            yield stdin
            stdin = ''
    if stdin:
        yield stdin


def minimap2(stdin, refpath, path='minimap2', nthread=1, minlen=29000):
    p = subprocess.Popen(
        [path, '-t', str(nthread), '-a', '--eqx', refpath, '-'], encoding='utf8',
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
    )
    stdout, stderr = p.communicate(stdin)
    for line in stdout.split('\n'):
        if line == '' or line.startswith('@'):
            # split on \n leaves empty line; @ prefix header lines
            continue
        qname, flag, rname, rpos, _, cigar, _, _, _, seq, qual = \
            line.strip('\n').split('\t')[:11]

        if rname == '*' or ((int(flag) & 0x800) != 0):
            # did not map, or supplementary alignment
            continue

        if len(seq) < minlen:
            # reject sequence that is too short
            continue

        # validate CIGAR string
        is_valid = re.match(r'^((\d+)([MIDNSHPX=]))*$', cigar)
        if not is_valid:
            raise RuntimeError('Invalid CIGAR string: {!r}.'.format(cigar))

        rpos = int(rpos) - 1  # convert to 0-index
        yield qname, rpos, cigar, seq, qual


def count_mutations(sequrl, lineages, refpath):
    # stream genomes in batches and extract and count differences
    results = {}
    handle = lzma.open(request.urlopen(sequrl), 'rt')
    for batch in batcher(handle):
        for row in minimap2(batch, refpath):
            qname, diffs, missing = encode_diffs(row)

            # retrieve PANGO lineage assignment
            lineage = lineages.get(qname, None)
            if lineage is None:
                continue

            if lineage not in results:
                results.update({lineage: {'mutations': {}, 'count': 0}})

            results[lineage]['count'] += 1
            for diff in diffs:

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequrl", help="URL to Nextstrain sequences.fasta.xz file",
                        default="https://data.nextstrain.org/files/ncov/open/sequences.fasta.xz")
    parser.add_argument("--metaurl", help="URL to Nextstrain metadata.tsv.gz file",
                        default="https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz")
    parser.add_argument("--mm2bin", default="minimap2", help="path to minimap2 binary executable file")
    parser.add_argument("--refpath", default="data/NC_045512.fa", help="path to reference genome FASTA")
    args = parser.parse_args()

    lineages = parse_metadata(args.metaurl)
    results = count_mutations(args.sequrl, lineages, refpath=args.refpath)
