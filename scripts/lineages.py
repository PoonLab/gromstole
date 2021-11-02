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
import subprocess
import re
import json
import time

def encode_diffs_classic(row, reflen=29903, alphabet='ACGT'):
    """
    Serialize differences of query sequences to reference
    genome, which comprise nucleotide substitutions, in-frame
    indels, and locations of missing data.
    NOTE: runs of 'N's are also represented by 'X' tokens in the CIGAR
    string.
    :param row:  tuple from minimap2(), i.e., (qname, rpos, cigar, seq)
    :param reflen:  length of reference genome
    """
    qname, rpos, cigar, seq = row  # unpack tuple
    diffs = []
    missing = []
    if rpos > 0:
        # incomplete on left
        missing.append(tuple([0, rpos]))
    left = 0  # index for query

    tokens = re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE)
    for length, operator in tokens:
        length = int(length)
        substr = seq[left:(left + length)]
        if operator == 'X':
            # each nucleotide is a separate diff
            if 'N' in substr:
                # for now, assume the whole substring is bs
                missing.append(tuple([rpos, rpos+length]))
            else:
                # assume adjacent mismatches are independent substitutions
                for i, nt in enumerate(substr):
                    if nt in alphabet:
                        diffs.append(tuple(['~', rpos + i, nt]))
                    else:
                        # skip ambiguous base calls, like "R"
                        missing.append(tuple([rpos+i, rpos+i+1]))

            left += length
            rpos += length
        elif operator == 'S':
            # discard soft clip
            left += length
        elif operator == 'I':
            # insertion relative to reference
            diffs.append(tuple(['+', rpos, substr]))
            left += length
        elif operator == 'D':
            # deletion relative to reference
            diffs.append(tuple(['-', rpos, length]))
            rpos += length
        elif operator == '=':
            # exact match
            left += length
            rpos += length
        elif operator == 'H':
            # hard clip, do nothing
            pass
        else:
            print("ERROR: unexpected operator {}".format(operator))
            sys.exit()

    # update missing if sequence incomplete on the right
    if rpos < reflen:
        missing.append(tuple([rpos, reflen]))

    return qname, diffs, missing

def minimap2_fasta(infile, ref, stream=False, path='minimap2', nthread=3, minlen=29000):
    """
    Wrapper function for minimap2.

    :param infile:  file object or StringIO
    :param ref:  str, path to FASTA with reference sequence(s)
    :param stream:  bool, if True then stream data from <infile> object via stdin
    :param path:  str, path to binary executable
    :param nthread:  int, number of threads for parallel execution of minimap2
    :param minlen:  int, filter genomes below minimum length; to accept all, set to 0.

    :yield:  query sequence name, reference index, CIGAR and original
             sequence
    """
    if stream:
        # input from StringIO in memory
        p = subprocess.Popen(
            [path, '-t', str(nthread), '-a', '--eqx', ref, '-'], encoding='utf8',
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
        )
        output, outerr = p.communicate(infile)
        output = output.split('\n')
    else:
        # input read from file
        p = subprocess.Popen(
            [path, '-t', str(nthread), '-a', '--eqx', ref, infile.name],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
        )
        output = map(lambda x: x.decode('utf-8'), p.stdout)

    for line in output:
        if line == '' or line.startswith('@'):
            # split on \n leaves empty line; @ prefix header lines
            continue

        qname, flag, rname, rpos, _, cigar, _, _, _, seq = \
            line.strip('\n').split('\t')[:10]

        if rname == '*' or ((int(flag) & 0x800) != 0):
            # did not map, or supplementary alignment
            continue

        if len(seq) < minlen:
            # sequence too short
            continue

        # validate CIGAR string
        is_valid = re.match(r'^((\d+)([MIDNSHPX=]))*$', cigar)
        if not is_valid:
            raise RuntimeError('Invalid CIGAR string: {!r}.'.format(cigar))

        rpos = int(rpos) - 1  # convert to 0-index
        yield qname, rpos, cigar, seq

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
counter = 0
res = []
tic0 = time.perf_counter()
for fasta in batcher:
    tic = time.perf_counter()
    mm2 = minimap2_fasta(fasta, ref=ref_file, stream=True, nthread=4, minlen=28000)
    counter += 1
    for qname, diffs, missing in [encode_diffs_classic(row) for row in mm2]:
        diffs = ["".join(map(str, x)) for x in diffs]
        res.append({"qname": qname, "lineage": lineages[qname], "diffs": diffs, "missing": missing})
    if not counter % 200:
        toc = time.perf_counter()
        print(f"{toc - tic:0.3f} for loop {counter}, {toc-tic0:0.3f} total.")
    #if counter >= 8:
    #    break

toc0 = time.perf_counter()
print(f"total time: {toc0 - tic0:0.3f}")

serial = json.dumps(res).replace('},', '},\n')
#print(serial)
open("data/pangodiffs.json", "w").write(serial)
