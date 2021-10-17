"""
Adapted from https://github.com/PoonLab/covizu/covizu/minimap2.py
"""

import subprocess
import argparse
import re
import sys
import os
import json

from seq_utils import convert_fasta, SC2Locator


def minimap2_paired(fq1, fq2, ref, path='minimap2', nthread=3):
    """
    Wrapper for minimap2 for processing paired-end read data.
    :param fq1:  str, path to FASTQ R1 file (uncompressed or gzipped)
    :param fq2:  str, path to FASTQ R2 file (uncompressed or gzipped)
    :param ref:  str, path to FASTA file with reference sequence
    :param path:  str, path to minimap2 binary executable
    :param nthread:  int, number of threads to run
    :yield:  tuple (qname, rpos, cigar, seq)
    """
    p = subprocess.Popen(
        [path, '-t', str(nthread), '-ax', 'sr', '--eqx', '--secondary=no', ref, fq1, fq2],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
    )
    output = map(lambda x: x.decode('utf-8'), p.stdout)
    for line in output:
        if line == '' or line.startswith('@'):
            # split on \n leaves empty line; @ prefix header lines
            continue

        qname, flag, rname, rpos, _, cigar, _, _, _, seq, qual = \
            line.strip('\n').split('\t')[:11]

        if rname == '*' or ((int(flag) & 0x800) != 0):
            # did not map, or supplementary alignment
            continue

        rpos = int(rpos) - 1  # convert to 0-index
        yield qname, rpos, cigar, seq, qual


def matchmaker(samfile):
    """
    An iterator that returns pairs of reads sharing a common qname from a SAM file.
    Note that unpaired reads will be left in the cached_rows dictionary and
    discarded.
    FIXME: it would probably be simpler to just pair up rows after headers
    @param samfile: open file handle to a SAM file
    @return: yields a tuple for each read pair with fields split by tab chars:
        ([qname, flag, rname, ...], [qname, flag, rname, ...])
    """
    cached_rows = {}
    for row in samfile:
        qname = row[0]
        old_row = cached_rows.pop(qname, None)
        if old_row is None:
            cached_rows[qname] = row
        else:
            # current row should be the second read of the pair
            yield old_row, row


def encode_diffs(row, reflen=29903, alphabet='ACGT', minq=10):
    """
    Serialize differences of query sequences to reference
    genome, which comprise nucleotide substitutions, in-frame
    indels, and locations of missing data.
    NOTE: runs of 'N's are also represented by 'X' tokens in the CIGAR
    string.
    :param row:  tuple, from minimap2(), i.e., (qname, rpos, cigar, seq)
    :param reflen:  int, length of reference genome
    :param alphabet:  str, accepted character states (nucleotides)
    :param minq:  int, minimum base quality
    """
    qname, rpos, cigar, seq, qual = row  # unpack tuple
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
        subqual = qual[left:(left + length)]

        if operator == 'X':
            # base mismatches
            if 'N' in substr:
                # for now, assume the whole substring is bs
                missing.append(tuple([rpos, rpos+length]))
            else:
                # assume adjacent mismatches are independent substitutions
                for i, nt in enumerate(substr):
                    q = subqual[i]
                    if nt in alphabet and q >= minq:
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


def check_miss(miss):
    """ Restore missing lower or upper bounds for missing ranges """
    if len(miss) < 2:
        if miss[0][0] == 0:
            # missing right bound
            miss.append(tuple([29903, 29903]))
        elif miss[0][1] == 29903:
            # missing left bound
            miss.insert(0, tuple([0, 0]))
        else:
            print("Unexpected case in check_miss(): {}".format(miss))
            sys.exit()
    return miss


def merge_diffs(diff1, diff2, miss1, miss2):
    """
    Merge the differences from reference for paired reads (including coverage).

    If diff is on the overlap of coverage, the diff must appear in both.
    If diff exists where only one read has coverage, always accept.
    If diff is only in one read and both have coverage, reject.
    """
    miss1 = check_miss(miss1)
    miss2 = check_miss(miss2)

    lo1 = miss1[0][1]  # miss1 is [(0, lo1), (hi1, 29903)]
    hi1 = miss1[1][0]
    lo2 = miss2[0][1]
    hi2 = miss2[1][0]

    if miss1 == miss2: 
        # same coverage: only return duplicates
        same_diffs = [value for value in diff1 if value in diff2]
        return same_diffs, [(lo1, hi1)]
    elif lo1 > hi2 or lo2 > hi1:
        # non-overlapping, return all mutations
        all_diffs = list(set(diff1 + diff2))
        return all_diffs, [(lo1, hi1), (lo2, hi2)]
    else: 
        # imperfect overlap
        cov_union = (min([lo1, lo2]), max([hi1, hi2]))
        cov_intersect = (max([lo1, lo2]), min([hi1, hi2]))

        all_diffs = list(set(diff1 + diff2))
        return_diffs = [value for value in diff1 if value in diff2]
        diff_diffs = [value for value in all_diffs if value not in return_diffs]

        diff_loc = [value for value in map(lambda x: x[1], diff_diffs)]
        for i in range(len(diff_diffs)):
            if cov_union[0] <= diff_loc[i] <= cov_union[1] and not cov_intersect[0] <= diff_loc[i] <= cov_intersect[1]:
                return_diffs.append(diff_diffs[i])
        return return_diffs, cov_union


def parse_mm2(mm2, locator, report=10000, stop=1e6, maxpos=29903):
    """
    Iterate over minimap2 output
    :param mm2:  generator, yields tuples of SAM output from minimap2
    :param locator:  SC2locator object, maps mutations to nicer notation
    :param report:  int, console reporting frequency
    :param stop:  int, limit to number of paired reads that mapped to sc2 to report
    :return:
    """
    count = 0
    total_coverage = dict([(pos, 0) for pos in range(maxpos)])
    res = []

    for r1, r2 in matchmaker(mm2):
        _, diff1, miss1 = encode_diffs(r1)
        _, diff2, miss2 = encode_diffs(r2)

        diffs, coverage = merge_diffs(diff1, diff2, miss1, miss2)

        # update coverage stats
        if len(coverage) == 2 and type(coverage[0]) is not tuple:
            coverage = [coverage]  # FIXME: this is a patch!
        for left, right in coverage:
            for pos in range(left, right):
                total_coverage[pos] += 1

        # get indel or AA sub notations
        mut = [locator.parse_mutation(d) for d in diffs]
        res.append({"diff": diffs, "mutations": mut})

        count += 1
        if count % report == 0:
            sys.stderr.write('{}\n'.format(count))
            sys.stderr.flush()
        if count > stop:
            break

    return res, total_coverage


def get_frequencies(res, coverage):
    counts = {}
    for row in res:
        for i, diff in enumerate(row['diff']):
            typ, pos, mut = diff
            label = row['mutations'][i]
            denom = coverage[pos]
            key = ''.join(map(str, [typ, pos, mut]))
            if pos not in counts:
                counts.update({pos: {}})
            if key not in counts[pos]:
                counts[pos].update({key: {'label': label, 'count': 0}})
            counts[pos][key]['count'] += 1/denom
    return counts


def parse_args():
    parser = argparse.ArgumentParser("Wrapper script for minimap2")

    parser.add_argument('fq1', type=str, help="<input> path to FASTQ R1 file.  May be gzip'ed.")
    parser.add_argument('fq2', type=str, help="<input> path to FASTQ R2 file.  May be gzip'ed.")
    parser.add_argument('-o', '--outfile',
                        type=argparse.FileType('w'),
                        required=False,
                        help="<output, optional> path to write output, "
                             "defaults to stdout")
    parser.add_argument('--path', type=str, default='minimap2', help="<option> path to minimap2 executable")

    parser.add_argument('-a', '--align', action='store_true',
                        help="<option> output aligned sequences as FASTA")
    parser.add_argument('-t', '--thread', type=int, default=3, 
                        help="<option> number of threads")

    parser.add_argument('--filter', action='store_true',
                        help="<option> filter problematic sites")

    ref_file = os.path.join("data", "NC_045512.fa")
    vcf_file = os.path.join("data", "problematic_sites_sarsCov2.vcf")

    parser.add_argument('--ref', type=str, default=ref_file,
                        help="<input> path to target FASTA (reference)")
    parser.add_argument("--vcf", type=str, default=vcf_file,
                        help="Path to VCF file of problematic sites in SARS-COV-2 genome. "
                             "Source: https://github.com/W-L/ProblematicSites_SARS-CoV2")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    if args.outfile is None:
        args.outfile = sys.stdout

    mm2 = minimap2_paired(args.fq1, args.fq2, ref=args.ref, nthread=args.thread, path=args.path)
    locator = SC2Locator(ref_file='data/NC_045512.fa')
    res, coverage = parse_mm2(mm2, locator, stop=1e3)
    counts = get_frequencies(res, coverage)

    # serial = json.dumps(res).replace('},', '},\n')
    # args.outfile.write(serial)

    args.outfile.write("mutation,frequency,coverage\n")
    for pos, mutations in counts.items():
        denom = coverage[pos]
        for mut, freq in mutations.items():
            args.outfile.write('{},{},{},{}\n'.format(mut, freq['label'], freq['count'], denom))
