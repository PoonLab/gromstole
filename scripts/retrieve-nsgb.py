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
import json
import sys


def parse_metadata(metaurl, progress=1e5):
    """ extract PANGO lineage assignments from metadata """
    lineages = {}  # header : PANGO lineage
    handle = gzip.open(request.urlopen(metaurl), 'rt')
    for count, row in enumerate(csv.DictReader(handle, delimiter='\t')):
        if progress and count % progress == 0:
            sys.stdout.write('.')
            sys.stdout.flush()
        if row["strain"] == '?':
            continue
        lineages.update({row['strain']: row["pango_lineage"]})
    sys.stdout.write('\n')
    return lineages


def batcher(handle, size=100, limit=None, minlen=29000):
    """ Break FASTA data into batches to stream to minimap2 """
    stdin = ''
    for i, record in enumerate(iter_fasta(handle)):
        if limit is not None and i == limit:
            break
        qname, seq = record
        if len(seq.replace('N', '')) < minlen:
            continue
        stdin += '>{}\n{}\n'.format(qname, seq)
        if i > 0 and i % size == 0:
            yield stdin
            stdin = ''
    if stdin:
        yield stdin


def minimap2(stdin, refpath, path='minimap2', nthread=1):
    """ Wrapper function on minimap2 """
    p = subprocess.Popen(
        [path, '-t', str(nthread), '-a', '--eqx', refpath, '-'], encoding='utf8',
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
    )
    stdout, stderr = p.communicate(stdin)
    for line in stdout.split('\n'):
        if line == '' or line.startswith('@'):
            # split on \n leaves empty line; @ prefix header lines
            continue
        qname, flag, rname, rpos, _, cigar, _, _, _, seq = \
            line.strip('\n').split('\t')[:10]
        qual = 'I' * len(seq)  # highest quality

        if rname == '*' or ((int(flag) & 0x800) != 0):
            # did not map, or supplementary alignment
            continue

        # validate CIGAR string
        is_valid = re.match(r'^((\d+)([MIDNSHPX=]))*$', cigar)
        if not is_valid:
            raise RuntimeError('Invalid CIGAR string: {!r}.'.format(cigar))

        rpos = int(rpos) - 1  # convert to 0-index
        yield qname, rpos, cigar, seq, qual


def count_mutations(sequrl, lineages, refpath, binpath='minimap2', nthreads=1,
                    minlen=29000, limit=None):
    """
    Stream genomes in batches and extract and count differences

    :param sequrl:  str, URL to remote xz-compressed FASTA file
    :param lineages:  dict, mapping sequence labels to PANGO lineage
    :param refpath:  str, path to FASTA file containing reference genome
    :param binpath:  str, path to minimap2 binary executable
    :param nthreads:  int, number of threads to run minimap2
    :param minlen:  int, minimum genome length
    :param limit:  int, limit number of genomes (for debugging)
    :return:  dict, mutation counts by lineage
    """
    results = {}
    handle = lzma.open(request.urlopen(sequrl), 'rt')
    for batch in batcher(handle, limit=limit, minlen=minlen):
        for row in minimap2(batch, refpath, path=binpath, nthread=nthreads):
            qname, diffs, missing = encode_diffs(row)

            # retrieve PANGO lineage assignment
            lineage = lineages.get(qname, None)
            if lineage is None:
                continue

            if lineage not in results:
                results.update({lineage: {'mutations': {}, 'count': 0}})

            results[lineage]['count'] += 1  # denominator
            for diff in diffs:
                key = '|'.join(map(str, diff))
                if key not in results[lineage]['mutations']:
                    results[lineage]['mutations'].update({key: 0})
                results[lineage]['mutations'][key] += 1

    return results


if __name__ == "__main__":
    # command line interface
    parser = argparse.ArgumentParser()
    parser.add_argument("outfile", type=argparse.FileType('w'), help="Path to write JSON output.")
    parser.add_argument("--sequrl", help="URL to Nextstrain sequences.fasta.xz file",
                        default="https://data.nextstrain.org/files/ncov/open/sequences.fasta.xz")
    parser.add_argument("--metaurl", help="URL to Nextstrain metadata.tsv.gz file; to use a local file, "
                                          "use file:// prefix instead of https://",
                        default="https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz")
    parser.add_argument("--mm2bin", default="minimap2", help="path to minimap2 binary executable file")
    parser.add_argument("-t", "--nthreads", type=int, default=1, help="number of threads for minimap2")
    parser.add_argument("--minlen", default=29000, help="minimum genome length")
    parser.add_argument("--refpath", default="data/NC_045512.fa", help="path to reference genome FASTA")
    parser.add_argument("--limit", type=int, default=None, help="limit number of genomes (for debugging)")
    args = parser.parse_args()

    print("\n"
          "*******************************************************\n"
          "Note, please avoid running this script multiple times a\n"
          "day unless using a local copy of metadata.tsv.gz       \n"
          "(override default --metaurl with file://) and setting  \n"
          "--limit to a small number, to avoid overloading the    \n"
          "nextstrain.org server!                                 \n"
          "*******************************************************\n\n")

    print("Parsing metadata")
    lineages = parse_metadata(args.metaurl)
    print("Streaming genomes and extracting features")
    results = count_mutations(args.sequrl, lineages, refpath=args.refpath, binpath=args.mm2bin,
                              nthreads=args.nthreads, minlen=args.minlen, limit=args.limit)
    print("Writing result to output file {}".format(args.outfile.name))
    json.dump(results, args.outfile, indent=2)
