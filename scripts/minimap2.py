"""
Adapted from https://github.com/PoonLab/covizu/covizu/minimap2.py
"""

import subprocess
import argparse
import re
import sys
import os
import json

from seq_utils import convert_fasta, filter_problematic_sites


def apply_cigar(seq, rpos, cigar):
    """
    Use CIGAR to pad sequence with gaps as required to
    align to reference.  Adapted from http://github.com/cfe-lab/MiCall
    """
    is_valid = re.match(r'^((\d+)([MIDNSHPX=]))*$', cigar)

    if not is_valid:
        raise RuntimeError('Invalid CIGAR string: {!r}.'.format(cigar))
    tokens = re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE)
    aligned = '-'*rpos
    left = 0
    for length, operation in tokens:
        length = int(length)
        if operation in 'M=X':
            aligned += seq[left:(left+length)]
            left += length
        elif operation == 'D':
            aligned += '-'*length
        elif operation in 'SI':
            left += length  # soft clip

    return aligned


def minimap2(infile, ref, stream=False, path='minimap2', nthread=3, minlen=29000):
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
            # reject sequence that is too short
            continue

        # validate CIGAR string
        is_valid = re.match(r'^((\d+)([MIDNSHPX=]))*$', cigar)
        if not is_valid:
            raise RuntimeError('Invalid CIGAR string: {!r}.'.format(cigar))

        rpos = int(rpos) - 1  # convert to 0-index
        yield qname, rpos, cigar, seq


def minimap2_paired(fq1, fq2, ref, path='minimap2', nthread=3):
    """

    :param fq1:  str, path to FASTQ R1 file (uncompressed or gzipped)
    :param fq2:  str, path to FASTQ R2 file (uncompressed or gzipped)
    :param ref:  str, path to FASTA file with reference sequence
    :param path:  str, path to minimap2 binary executable
    :param nthread:  int, number of threads to run
    :return:
    """
    p = subprocess.Popen(
        [path, '-t', str(nthread), '-ax', 'sr', '--eqx', ref, fq1, fq2],
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

        rpos = int(rpos) - 1  # convert to 0-index
        yield qname, rpos, cigar, seq


def matchmaker(samfile):
    """
    An iterator that returns pairs of reads sharing a common qname from a SAM file.
    Note that unpaired reads will be left in the cached_rows dictionary and
    discarded.
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


# return aligned sequence?
def output_fasta(iter, outfile, reflen=0):
    """
    Stream output from minimap2 into FASTA file
    of aligned sequences.  CIGAR parsing code adapted from
    http://github.com/cfe-lab/MiCall

    :param iter:  generator from minimap2()
    :param outfile:  open file stream in write mode
    :param reflen:  int, length of reference genome to pad sequences;
                    defaults to no padding.
    """
    for qname, rpos, cigar, seq in iter:
        tokens = re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE)
        aligned = '-' * rpos
        left = 0
        for length, operation in tokens:
            length = int(length)
            if operation in 'M=X':
                aligned += seq[left:(left + length)]
                left += length
            elif operation == 'D':
                aligned += '-' * length
            elif operation in 'SI':
                left += length  # soft clip

        # pad on right
        aligned += '-'*(reflen-len(aligned))
        outfile.write('>{}\n{}\n'.format(qname, aligned))


def stream_fasta(iter, reflen=0):
    """
    Stream output from minimap2 into list of tuples [(header, seq), ... , (header, seq)]
    of aligned sequences.  CIGAR parsing code adapted from
    http://github.com/cfe-lab/MiCall

    :param iter:  generator from minimap2()
    :param reflen:  int, length of reference genome to pad sequences;
                    defaults to no padding.
    :yield:  tuple, header and aligned sequence
    """
    for qname, rpos, cigar, seq in iter:
        tokens = re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE)
        aligned = '-' * rpos
        left = 0
        for length, operation in tokens:
            length = int(length)
            if operation in 'M=X':
                aligned += seq[left:(left + length)]
                left += length
            elif operation == 'D':
                aligned += '-' * length
            elif operation in 'SI':
                left += length  # soft clip

        # pad on right
        aligned += '-'*(reflen-len(aligned))
        yield qname, aligned


def encode_diffs(qname, rpos, cigar, seq, reflen=29903, alphabet='ACGT'):
    """
    Serialize differences of query sequences to reference
    genome, which comprise nucleotide substitutions, in-frame
    indels, and locations of missing data.
    NOTE: runs of 'N's are also represented by 'X' tokens in the CIGAR
    string.
    :param iter:  generator from minimap2()
    :param reflen:  length of reference genome
    """
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


def extract_features(batcher, ref_file, binpath='minimap2', nthread=3, minlen=29000):
    """
    Stream output from JSON.xz file via load_gisaid() into minimap2
    via subprocess.

    :param batcher:  generator, returned by batch_fasta()
    :param ref_file:  str, path to reference genome (FASTA format)
    :param binpath:  str, path to minimap2 binary executable
    :param nthread:  int, number of threads to run minimap2
    :param minlen:  int, minimum genome length

    :yield:  dict, record augmented with genetic differences and missing sites;
    """
    with open(ref_file) as handle:
        reflen = len(convert_fasta(handle)[0][1])

    for fasta, batch in batcher:
        mm2 = minimap2(fasta, ref_file, stream=True, path=binpath, nthread=nthread,
                       minlen=minlen)
        result = list(encode_diffs(mm2, reflen=reflen))
        for row, record in zip(result, batch):
            # reconcile minimap2 output with GISAID record
            qname, diffs, missing = row
            record.update({'diffs': diffs, 'missing': missing})
            yield record


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
    for pairs in matchmaker(mm2):
        for qname, rpos, cigar, seq in pairs:
            diffs = encode_diffs(qname, rpos, cigar, seq)
            print(diffs)
        break

    sys.exit()

    # get length of reference
    refseq = convert_fasta(open(args.ref))[0][1]
    reflen = len(refseq)
    """
        if args.align:
        if args.filter:
            vcf = load_vcf(args.vcf)
            encoded = encode_diffs(mm2, reflen=reflen)
            for qname, diffs, missing in filter_problematic_sites(encoded, mask=vcf):
                seq = apply_features(diffs, missing, refseq=refseq)
                args.outfile.write(">{}\n{}\n".format(qname, seq))
        else:
            output_fasta(mm2, reflen=reflen, outfile=args.outfile)
    """

    # serialize feature vectors as JSON
    res = []
    for row in mm2:
        qname, diffs, missing = encode_diffs(row)
        res.append({'name': qname, 'diffs': diffs, 'missing': missing})
    serial = json.dumps(res).replace('},', '},\n')
    args.outfile.write(serial)

